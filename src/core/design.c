#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include "memory.h"
#include "design.h"


static void
iproc_design_free0 (iproc_design *design)
{
    if (design->free_user_data) {
        design->free_user_data(design->user_data);
    }
}

static void
iproc_design_clear_svectors (iproc_array *svectors)
{
    if (svectors) {
        int64_t i, n = iproc_array_size(svectors);
        for (i = 0; i < n; i++) {
            iproc_svector *x = iproc_array_index(svectors,
                                                 iproc_svector *,
                                                 i);
            iproc_svector_unref(x);
        }
        iproc_array_set_size(svectors, 0);
    }
}


static void
iproc_design_free_svectors (iproc_array *svectors)
{
    if (svectors) {
        iproc_design_clear_svectors(svectors);
        iproc_array_unref(svectors);
    }
}

static void
iproc_design_ctx_free_dealloc (iproc_design_ctx *ctx)
{
    if (ctx) {
        assert(iproc_array_size(ctx->dxs) == 0);
        iproc_array_unref(ctx->dxs);
        iproc_free(ctx);
    }
}

static void
iproc_design_free_ctxs (iproc_array *ctxs)
{
    if (ctxs) {
        int64_t i, n = iproc_array_size(ctxs);
        for (i = 0; i < n; i++) {
            iproc_design_ctx *ctx = iproc_array_index(ctxs,
                                                      iproc_design_ctx *,
                                                      i);
            iproc_design_ctx_free_dealloc(ctx);
        }
        iproc_array_unref(ctxs);
    }
}

static void
iproc_design_free_vars (iproc_array *vars)
{
    if (vars) {
        int64_t i, n = iproc_array_size(vars);
        for (i = 0; i < n; i++) {
            iproc_design_var *var = iproc_array_index(vars,
                                                      iproc_design_var *,
                                                      i);
            if (var->free)
                var->free(var);
        }
        iproc_array_unref(vars);
    }
}

static void
iproc_design_free (iproc_design *design)
{

    if (design) {
        iproc_design_free0(design);
        
        
        iproc_design_free_svectors(design->svectors);        
        iproc_design_free_ctxs(design->ctxs);
        iproc_design_free_vars(design->vars);
        iproc_actors_unref(design->receivers);
        iproc_actors_unref(design->senders);
        iproc_free(design);
    }
}

iproc_design *
iproc_design_new (iproc_actors *senders,
                  iproc_actors *receivers)
{
    assert(senders);
    assert(receivers);

    iproc_design *design = iproc_malloc(sizeof(*design));
    
    if (!design)
        return NULL;

    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    int64_t dim0 = p * q;

    design->senders = iproc_actors_ref(senders);
    design->receivers = iproc_actors_ref(receivers);
    design->vars = iproc_array_new(sizeof(iproc_design_var *));
    design->dim0 = dim0;
    design->dim = dim0;
    design->ctxs = iproc_array_new(sizeof(iproc_design_ctx *));
    design->svectors = iproc_array_new(sizeof(iproc_svector *));
    iproc_refcount_init(&design->refcount);

    if (!(design->vars && design->ctxs && design->svectors)) {
        iproc_design_free(design);
        design = NULL;
    }

    return design;
}

iproc_design *
iproc_design_ref (iproc_design *design)
{
    if (design) {
        iproc_refcount_get(&design->refcount);
    }
    return design;
}

static void
iproc_design_release (iproc_refcount *refcount)
{
    iproc_design *design = container_of(refcount, iproc_design, refcount);
    iproc_design_free(design);
}

void
iproc_design_unref (iproc_design *design)
{
    if (!design)
        return;

    iproc_refcount_put(&design->refcount, iproc_design_release);
}

int64_t
iproc_design_dim (iproc_design *design)
{
    assert(design);
    int64_t nstatic = iproc_design_nstatic(design);
    int64_t ndynamic = iproc_design_ndynamic(design);
    int64_t dim = nstatic + ndynamic;
    return dim;
}

void
iproc_design_append (iproc_design     *design,                                           
                     iproc_design_var *var)
{
    assert(design);
    assert(var);
    assert(var->dim >= 0);
    assert(var->get_dxs);

    // the old svectors are invalid since the dimension has changed
    iproc_design_clear_svectors(design->svectors);
    iproc_array_append(design->vars, &var);
    design->dim += var->dim;
}

void
iproc_design_mul0 (double        alpha,
                   iproc_trans   trans,
                   iproc_design *design,
                   int64_t       isend,
                   iproc_vector *x,
                   double        beta,
                   iproc_vector *y)
{
    assert(design);
    assert(isend >= 0);
    assert(isend < iproc_design_nsender(design));
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_design_dim(design));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_design_nreceiver(design));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_design_nreceiver(design));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_design_dim(design));

    /* y := beta y */
    if (beta == 0.0) {
        iproc_vector_set_all(y, 0.0);
    } else if (beta != 1.0) {
        iproc_vector_scale(y, beta);
    }

    if (iproc_design_nstatic(design) == 0)
        return;

    iproc_actors *senders = iproc_design_senders(design);
    iproc_actors *receivers = iproc_design_receivers(design);
    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    int64_t ix_begin = iproc_design_istatic(design, 0);
    int64_t nstatic = iproc_design_nstatic(design);
    iproc_vector *s = iproc_actors_get(senders, isend);
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
iproc_design_muls0 (double          alpha,
                         iproc_trans     trans,
                         iproc_design     *design,
                         int64_t         isend,
                         iproc_svector  *x,
                         double          beta,
                         iproc_vector   *y)
{
    assert(design);
    assert(isend >= 0);
    assert(isend < iproc_design_nsender(design));
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_design_dim(design));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_design_nreceiver(design));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_design_nreceiver(design));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_design_dim(design));

    /* y := beta y */
    if (beta == 0.0) {
        iproc_vector_set_all(y, 0.0);
    } else if (beta != 1.0) {
        iproc_vector_scale(y, beta);
    }

    if (iproc_design_nstatic(design) == 0)
        return;

    iproc_actors *senders = iproc_design_senders(design);
    iproc_actors *receivers = iproc_design_receivers(design);
    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    int64_t ix_begin = iproc_design_istatic(design, 0);
    int64_t nstatic = iproc_design_nstatic(design);
    int64_t ix_end = ix_begin + nstatic;
    iproc_vector *s = iproc_actors_get(senders, isend);
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
            double x_ij = iproc_svector_nz_get(x, inz);
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


////////////////////////////////


iproc_design *     
iproc_design_new0 (iproc_actors *senders,
                   iproc_actors *receivers,
                   int64_t       ndynamic,
                   void         *user_data,
                   void        (*get_sdesign_vars) (iproc_design_ctx *ctx),
                   void        (*free_user_data)  (void *user_data))
{
    assert(ndynamic >= 0);
    
    iproc_design *design = iproc_design_new(senders, receivers);
    
    if (!design)
        return NULL;
    
    design->nstatic = design->dim0;
    design->ndynamic = ndynamic;
    design->dim += ndynamic;
    design->user_data = user_data;
    design->get_sdesign_vars = get_sdesign_vars;
    design->free_user_data = free_user_data;
    
    return design;
}

/////// make these static ?

int64_t
iproc_design_nstatic (iproc_design *design)
{
    assert(design);
    return design->dim0;
}

int64_t
iproc_design_istatic (iproc_design *design,
                      int64_t     i)
{
    assert(design);
    assert(i >= 0);
    assert(i < iproc_design_nstatic(design));
    
    return i;
}

int64_t
iproc_design_ndynamic (iproc_design *design)
{
    assert(design);
    return design->ndynamic;
}

int64_t
iproc_design_idynamic (iproc_design *design,
                       int64_t       i)
{
    assert(design);
    assert(i >= 0);
    assert(i < iproc_design_ndynamic(design));
    
    int64_t idynamic = i + design->dim0;
    return idynamic;
}

int64_t
iproc_design_nsender (iproc_design *design)
{
    assert(design);
    iproc_actors *senders = iproc_design_senders(design);
    return iproc_actors_size(senders);
}

int64_t
iproc_design_nreceiver (iproc_design *design)
{
    assert(design);
    iproc_actors *receivers = iproc_design_receivers(design);
    return iproc_actors_size(receivers);
}

iproc_actors *
iproc_design_senders (iproc_design *design)
{
    assert(design);
    iproc_actors *senders = design->senders;
    return senders;
}

iproc_actors *
iproc_design_receivers (iproc_design *design)
{
    assert(design);
    iproc_actors *receivers = design->receivers;
    return receivers;
}




