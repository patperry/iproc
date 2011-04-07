#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include "memory.h"
#include "design.h"


static void
iproc_design_clear_svectors (struct darray *svectors)
{
    if (svectors) {
        int64_t i, n = darray_size(svectors);
        for (i = 0; i < n; i++) {
            iproc_svector *x = darray_index(svectors,
                                            iproc_svector *,
                                            i);
            iproc_svector_unref(x);
        }
        darray_resize(svectors, 0);
    }
}


static void
iproc_design_free_svectors (struct darray *svectors)
{
    if (svectors) {
        iproc_design_clear_svectors(svectors);
        darray_free(svectors);
    }
}

static void
iproc_design_ctx_free_dealloc (iproc_design_ctx *ctx)
{
    if (ctx) {
        assert(darray_size(ctx->dxs) == 0);
        darray_free(ctx->dxs);
        iproc_free(ctx);
    }
}

static void
iproc_design_free_ctxs (struct darray *ctxs)
{
    if (ctxs) {
        int64_t i, n = darray_size(ctxs);
        for (i = 0; i < n; i++) {
            iproc_design_ctx *ctx = darray_index(ctxs,
                                                 iproc_design_ctx *,
                                                 i);
            iproc_design_ctx_free_dealloc(ctx);
        }
        darray_free(ctxs);
    }
}

static void
iproc_design_free_vars (struct darray *vars)
{
    if (vars) {
        int64_t i, n = darray_size(vars);
        for (i = 0; i < n; i++) {
            iproc_design_var *var = darray_index(vars,
                                                 iproc_design_var *,
                                                 i);
            if (var->free)
                var->free(var);
        }
        darray_free(vars);
    }
}

static void
iproc_design_free (iproc_design *design)
{

    if (design) {
        iproc_design_free_svectors(design->svectors);        
        iproc_design_free_ctxs(design->ctxs);
        iproc_design_free_vars(design->vars);
        iproc_actors_unref(design->receivers);
        iproc_actors_unref(design->senders);
        iproc_free(design);
    }
}

void
iproc_design_var_init  (iproc_design_var *var,
                        int64_t           dim,
                        void (*get_dxs) (iproc_design_var *,
                                         iproc_design_ctx *ctx,
                                         int64_t),
                        void (*free)    (iproc_design_var *))
{
    assert(var);
    assert(dim >= 0);

    var->dim = dim;
    var->get_dxs = get_dxs;
    var->free = free;
}


iproc_design *
iproc_design_new (iproc_actors *senders,
                  iproc_actors *receivers,
                  bool          has_reffects)
{
    assert(senders);
    assert(receivers);

    iproc_design *design = iproc_malloc(sizeof(*design));
    
    if (!design)
        return NULL;

    int64_t nreceivers = iproc_actors_size(receivers);
    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);

    design->senders = iproc_actors_ref(senders);
    design->receivers = iproc_actors_ref(receivers);
    design->has_reffects = has_reffects;
    design->vars = darray_new(iproc_design_var *);
    design->ireffects = 0;
    design->nreffects = has_reffects ? nreceivers : 0;
    design->istatic = design->ireffects + design->nreffects;
    design->nstatic = p * q;
    design->idynamic = design->istatic + design->nstatic;
    design->ndynamic = 0;
    design->dim = design->idynamic + design->ndynamic;
    design->ctxs = darray_new(iproc_design_ctx *);
    design->svectors = darray_new(iproc_svector *);
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
    return design->dim;
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
    darray_push_back(design->vars, &var);
    design->ndynamic += var->dim;
    design->dim += var->dim;
}


static void
iproc_design_mul0_reffects (double        alpha,
                            iproc_trans   trans,
                            iproc_design *design,
                            iproc_vector *x,
                            iproc_vector *y)
{
    if (!design->has_reffects)
        return;
    
    int64_t off = design->ireffects;
    int64_t dim = design->nreffects;
    
    if (trans == IPROC_TRANS_NOTRANS) {
        iproc_vector_view xsub = iproc_vector_subvector(x, off, dim);
        iproc_vector_acc(y, alpha, &xsub.vector);
    } else {
        iproc_vector_view ysub = iproc_vector_subvector(y, off, dim);
        iproc_vector_acc(&ysub.vector, alpha, x);
    }
}

static void
iproc_design_muls0_reffects (double         alpha,
                             iproc_trans    trans,
                             iproc_design  *design,
                             iproc_svector *x,
                             iproc_vector  *y)
{
    if (!design->has_reffects)
        return;

    int64_t off = design->ireffects;
    int64_t dim = design->nreffects;
    int64_t end = off + dim;
    
    if (trans == IPROC_TRANS_NOTRANS) {
        int64_t nnz = iproc_svector_nnz(x);
        int64_t inz = iproc_svector_find_nz(x, off);
        
        if (inz < 0)
            inz = ~inz;
        
        while (inz < nnz) {
            int64_t i = iproc_svector_nz(x, inz);
            if (i >= end)
                break;
            
            double x_i = iproc_svector_nz_get(x, inz);
            iproc_vector_inc(y, i, alpha * x_i);
            
            inz++;
        }
    } else {
        iproc_vector_view ysub = iproc_vector_subvector(y, off, dim);
        iproc_vector_sacc(&ysub.vector, alpha, x);
    }
}


static void
iproc_design_mul0_static (double        alpha,
                          iproc_trans   trans,
                          iproc_design *design,
                          int64_t       isend,
                          iproc_vector *x,
                          iproc_vector *y)
{
    if (design->nstatic == 0)
        return;

    iproc_actors *senders = iproc_design_senders(design);
    iproc_actors *receivers = iproc_design_receivers(design);
    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    int64_t ix_begin = design->istatic;
    int64_t nstatic = design->nstatic;
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


static void
iproc_design_muls0_static (double         alpha,
                           iproc_trans    trans,
                           iproc_design  *design,
                           int64_t        isend,
                           iproc_svector *x,
                           iproc_vector  *y)
{
    if (design->nstatic == 0)
        return;

    iproc_actors *senders = iproc_design_senders(design);
    iproc_actors *receivers = iproc_design_receivers(design);
    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    int64_t ix_begin = design->istatic;
    int64_t nstatic = design->nstatic;
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
    
    iproc_design_mul0_reffects(alpha, trans, design, x, y);
    iproc_design_mul0_static(alpha, trans, design, isend, x, y);
}


void
iproc_design_muls0 (double         alpha,
                    iproc_trans    trans,
                    iproc_design  *design,
                    int64_t        isend,
                    iproc_svector *x,
                    double         beta,
                    iproc_vector  *y)
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
    
    iproc_design_muls0_reffects(alpha, trans, design, x, y);
    iproc_design_muls0_static(alpha, trans, design, isend, x, y);
}


/////// make these static ?

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




