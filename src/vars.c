#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include "memory.h"
#include "vars.h"

static void
iproc_vars_free (iproc_vars *vars)
{
    if (vars) {
        iproc_actors_unref(vars->receivers);
        iproc_actors_unref(vars->senders);
        if (vars->free_user_data)
            vars->free_user_data(vars->user_data);
        iproc_free(vars);
    }
}

iproc_vars *     
iproc_vars_new (iproc_actors *senders,
                iproc_actors *receivers,
                int64_t       ndynamic,
                void         *user_data,
                void        (*get_sender_vars) (iproc_vars_ctx *ctx),
                void        (*free_user_data)  (void *user_data))
{
    assert(senders);
    assert(receivers);
    assert(ndynamic >= 0);

    iproc_vars *vars = iproc_malloc(sizeof(*vars));

    if (!vars)
        return NULL;

    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    int64_t nstatic = p * q;

    vars->senders = iproc_actors_ref(senders);
    vars->receivers = iproc_actors_ref(receivers);
    vars->nstatic = nstatic;
    vars->ndynamic = ndynamic;
    vars->user_data = user_data;
    vars->get_sender_vars = get_sender_vars;
    vars->free_user_data = free_user_data;
    iproc_refcount_init(&vars->refcount);
    return vars;
}

iproc_vars *
iproc_vars_ref (iproc_vars *vars)
{
    if (vars) {
        iproc_refcount_get(&vars->refcount);
    }
    return vars;
}

static void
iproc_vars_release (iproc_refcount *refcount)
{
    iproc_vars *vars = container_of(refcount, iproc_vars, refcount);
    iproc_vars_free(vars);
}

void
iproc_vars_unref (iproc_vars *vars)
{
    if (!vars)
        return;

    iproc_refcount_put(&vars->refcount, iproc_vars_release);
}

int64_t
iproc_vars_dim (iproc_vars *vars)
{
    assert(vars);
    int64_t nstatic = iproc_vars_nstatic(vars);
    int64_t ndynamic = iproc_vars_ndynamic(vars);
    int64_t dim = nstatic + ndynamic;
    return dim;
}

int64_t
iproc_vars_nstatic (iproc_vars *vars)
{
    assert(vars);
    return vars->nstatic;
}

int64_t
iproc_vars_istatic (iproc_vars *vars,
                    int64_t     i)
{
    assert(vars);
    assert(i >= 0);
    assert(i < iproc_vars_nstatic(vars));

    int64_t ndynamic = iproc_vars_ndynamic(vars);
    int64_t istatic = ndynamic + i;
    return istatic;
}

int64_t
iproc_vars_ndynamic (iproc_vars *vars)
{
    assert(vars);
    return vars->ndynamic;
}

int64_t
iproc_vars_idynamic (iproc_vars *vars,
                     int64_t     i)
{
    assert(vars);
    assert(i >= 0);
    assert(i < iproc_vars_ndynamic(vars));

    int64_t idynamic = i;
    return idynamic;
}

int64_t
iproc_vars_nsender (iproc_vars *vars)
{
    assert(vars);
    iproc_actors *senders = iproc_vars_senders(vars);
    return iproc_actors_size(senders);
}

int64_t
iproc_vars_nreceiver (iproc_vars *vars)
{
    assert(vars);
    iproc_actors *receivers = iproc_vars_receivers(vars);
    return iproc_actors_size(receivers);
}

iproc_actors *
iproc_vars_senders (iproc_vars *vars)
{
    assert(vars);
    iproc_actors *senders = vars->senders;
    return senders;
}

iproc_actors *
iproc_vars_receivers (iproc_vars *vars)
{
    assert(vars);
    iproc_actors *receivers = vars->receivers;
    return receivers;
}


void
iproc_vars_sender0_mul (double        alpha,
                        iproc_trans   trans,
                        iproc_vars   *vars,
                        int64_t       isend,
                        iproc_vector *x,
                        double        beta,
                        iproc_vector *y)
{
    assert(vars);
    assert(isend >= 0);
    assert(isend < iproc_vars_nsender(vars));
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_vars_dim(vars));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_vars_nreceiver(vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_vars_nreceiver(vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_vars_dim(vars));

    /* y := beta y */
    if (beta == 0.0) {
        iproc_vector_set_all(y, 0.0);
    } else if (beta != 1.0) {
        iproc_vector_scale(y, beta);
    }

    if (iproc_vars_nstatic(vars) == 0)
        return;

    iproc_actors *senders = iproc_vars_senders(vars);
    iproc_actors *receivers = iproc_vars_receivers(vars);
    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    int64_t ix_begin = iproc_vars_istatic(vars, 0);
    int64_t nstatic = iproc_vars_nstatic(vars);
    iproc_vector *s = iproc_actors_traits(senders, isend);
    iproc_vector *z = iproc_vector_new(q);

    if (trans == IPROC_TRANS_NOTRANS) {
        iproc_vector_view xsub = iproc_vector_subvector(x, ix_begin, nstatic);

        /* z := alpha t(x) s */
        iproc_matrix_view xmat = iproc_matrix_view_vector(&xsub.vector, p, q);
        iproc_matrix_mul(alpha, IPROC_TRANS_TRANS, &xmat.matrix, s, 0.0, z);

        /* y := y + R z */
        iproc_actors_mul(1.0, IPROC_TRANS_NOTRANS, receivers, z, 1.0, y);
    } else {
        /* z := alpha t(R) x */
        iproc_actors_mul(alpha, IPROC_TRANS_TRANS, receivers, x, 0.0, z);

        /* y := y + s \otimes z */
        iproc_vector_view ysub = iproc_vector_subvector(y, ix_begin, nstatic);
        iproc_matrix_view ymat = iproc_matrix_view_vector(&ysub.vector, p, q);
        iproc_matrix_view smat = iproc_matrix_view_vector(s, p, 1);
        iproc_matrix_view zmat = iproc_matrix_view_vector(z, 1, q);
        iproc_matrix_matmul(1.0, IPROC_TRANS_NOTRANS, &smat.matrix, &zmat.matrix,
                            1.0, &ymat.matrix);
    }

    iproc_vector_unref(z);
}


void
iproc_vars_sender0_muls (double          alpha,
                         iproc_trans     trans,
                         iproc_vars     *vars,
                         int64_t         isend,
                         iproc_svector  *x,
                         double          beta,
                         iproc_vector   *y)
{
    assert(vars);
    assert(isend >= 0);
    assert(isend < iproc_vars_nsender(vars));
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_vars_dim(vars));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_vars_nreceiver(vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_vars_nreceiver(vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_vars_dim(vars));

    /* y := beta y */
    if (beta == 0.0) {
        iproc_vector_set_all(y, 0.0);
    } else if (beta != 1.0) {
        iproc_vector_scale(y, beta);
    }

    if (iproc_vars_nstatic(vars) == 0)
        return;

    iproc_actors *senders = iproc_vars_senders(vars);
    iproc_actors *receivers = iproc_vars_receivers(vars);
    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    int64_t ix_begin = iproc_vars_istatic(vars, 0);
    int64_t nstatic = iproc_vars_nstatic(vars);
    int64_t ix_end = ix_begin + nstatic;
    iproc_vector *s = iproc_actors_traits(senders, isend);
    iproc_vector *z = iproc_vector_new(q);

    if (trans == IPROC_TRANS_NOTRANS) {
        /* z := alpha t(x) s 
         *
         * z[j] = alpha * { \sum_i (x[i,j] * s[i]) }
         */
        iproc_vector_set_all(z, 0.0);
        int64_t inz, nnz = iproc_svector_nnz(x);
        for (inz = 0; inz < nnz; inz++) {
            int64_t ix = iproc_svector_nz(x, inz);

            if (ix < ix_begin)
                continue;

            if (ix >= ix_end)
                break;
            
            imaxdiv_t ij = imaxdiv(ix - ix_begin, p);
            int64_t i = ij.rem;  /* ix % p */
            int64_t j = ij.quot; /* ix / p */
            double x_ij = iproc_svector_nz_val(x, inz);
            double s_i = iproc_vector_get(s, i);
            iproc_vector_inc(z, j, x_ij * s_i);
        }
        iproc_vector_scale(z, alpha);
                             
        /* y := y + R z */
        iproc_actors_mul(1.0, IPROC_TRANS_NOTRANS, receivers, z, 1.0, y);
    } else {
        /* z := alpha t(R) x */
        iproc_actors_muls(alpha, IPROC_TRANS_TRANS, receivers, x, 0.0, z);

        /* y := y + s \otimes z */
        iproc_vector_view ysub = iproc_vector_subvector(y, ix_begin, nstatic);
        iproc_matrix_view ymat = iproc_matrix_view_vector(&ysub.vector, p, q);
        iproc_matrix_view smat = iproc_matrix_view_vector(s, p, 1);
        iproc_matrix_view zmat = iproc_matrix_view_vector(z, 1, q);
        iproc_matrix_matmul(1.0, IPROC_TRANS_NOTRANS, &smat.matrix, &zmat.matrix,
                            1.0, &ymat.matrix);
    }

    iproc_vector_unref(z);
}

