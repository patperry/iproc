#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include "memory.h"
#include "frame.h"

static void
iproc_frame_free (iproc_frame *frame)
{
    if (frame) {
        iproc_actors_unref(frame->receivers);
        iproc_actors_unref(frame->senders);
        if (frame->free_user_data)
            frame->free_user_data(frame->user_data);
        iproc_free(frame);
    }
}

iproc_frame *     
iproc_frame_new (iproc_actors *senders,
                iproc_actors *receivers,
                int64_t       ndynamic,
                void         *user_data,
                void        (*get_sender_frame) (iproc_frame_ctx *ctx),
                void        (*free_user_data)  (void *user_data))
{
    assert(senders);
    assert(receivers);
    assert(ndynamic >= 0);

    iproc_frame *frame = iproc_malloc(sizeof(*frame));

    if (!frame)
        return NULL;

    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    int64_t nstatic = p * q;

    frame->senders = iproc_actors_ref(senders);
    frame->receivers = iproc_actors_ref(receivers);
    frame->nstatic = nstatic;
    frame->ndynamic = ndynamic;
    frame->user_data = user_data;
    frame->get_sender_frame = get_sender_frame;
    frame->free_user_data = free_user_data;
    iproc_refcount_init(&frame->refcount);
    return frame;
}

iproc_frame *
iproc_frame_ref (iproc_frame *frame)
{
    if (frame) {
        iproc_refcount_get(&frame->refcount);
    }
    return frame;
}

static void
iproc_frame_release (iproc_refcount *refcount)
{
    iproc_frame *frame = container_of(refcount, iproc_frame, refcount);
    iproc_frame_free(frame);
}

void
iproc_frame_unref (iproc_frame *frame)
{
    if (!frame)
        return;

    iproc_refcount_put(&frame->refcount, iproc_frame_release);
}

int64_t
iproc_frame_dim (iproc_frame *frame)
{
    assert(frame);
    int64_t nstatic = iproc_frame_nstatic(frame);
    int64_t ndynamic = iproc_frame_ndynamic(frame);
    int64_t dim = nstatic + ndynamic;
    return dim;
}

int64_t
iproc_frame_nstatic (iproc_frame *frame)
{
    assert(frame);
    return frame->nstatic;
}

int64_t
iproc_frame_istatic (iproc_frame *frame,
                    int64_t     i)
{
    assert(frame);
    assert(i >= 0);
    assert(i < iproc_frame_nstatic(frame));

    int64_t ndynamic = iproc_frame_ndynamic(frame);
    int64_t istatic = ndynamic + i;
    return istatic;
}

int64_t
iproc_frame_ndynamic (iproc_frame *frame)
{
    assert(frame);
    return frame->ndynamic;
}

int64_t
iproc_frame_idynamic (iproc_frame *frame,
                     int64_t     i)
{
    assert(frame);
    assert(i >= 0);
    assert(i < iproc_frame_ndynamic(frame));

    int64_t idynamic = i;
    return idynamic;
}

int64_t
iproc_frame_nsender (iproc_frame *frame)
{
    assert(frame);
    iproc_actors *senders = iproc_frame_senders(frame);
    return iproc_actors_size(senders);
}

int64_t
iproc_frame_nreceiver (iproc_frame *frame)
{
    assert(frame);
    iproc_actors *receivers = iproc_frame_receivers(frame);
    return iproc_actors_size(receivers);
}

iproc_actors *
iproc_frame_senders (iproc_frame *frame)
{
    assert(frame);
    iproc_actors *senders = frame->senders;
    return senders;
}

iproc_actors *
iproc_frame_receivers (iproc_frame *frame)
{
    assert(frame);
    iproc_actors *receivers = frame->receivers;
    return receivers;
}


void
iproc_frame_sender0_mul (double        alpha,
                        iproc_trans   trans,
                        iproc_frame   *frame,
                        int64_t       isend,
                        iproc_vector *x,
                        double        beta,
                        iproc_vector *y)
{
    assert(frame);
    assert(isend >= 0);
    assert(isend < iproc_frame_nsender(frame));
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_frame_dim(frame));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_frame_nreceiver(frame));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_frame_nreceiver(frame));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_frame_dim(frame));

    /* y := beta y */
    if (beta == 0.0) {
        iproc_vector_set_all(y, 0.0);
    } else if (beta != 1.0) {
        iproc_vector_scale(y, beta);
    }

    if (iproc_frame_nstatic(frame) == 0)
        return;

    iproc_actors *senders = iproc_frame_senders(frame);
    iproc_actors *receivers = iproc_frame_receivers(frame);
    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    int64_t ix_begin = iproc_frame_istatic(frame, 0);
    int64_t nstatic = iproc_frame_nstatic(frame);
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
iproc_frame_sender0_muls (double          alpha,
                         iproc_trans     trans,
                         iproc_frame     *frame,
                         int64_t         isend,
                         iproc_svector  *x,
                         double          beta,
                         iproc_vector   *y)
{
    assert(frame);
    assert(isend >= 0);
    assert(isend < iproc_frame_nsender(frame));
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_frame_dim(frame));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_frame_nreceiver(frame));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_frame_nreceiver(frame));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_frame_dim(frame));

    /* y := beta y */
    if (beta == 0.0) {
        iproc_vector_set_all(y, 0.0);
    } else if (beta != 1.0) {
        iproc_vector_scale(y, beta);
    }

    if (iproc_frame_nstatic(frame) == 0)
        return;

    iproc_actors *senders = iproc_frame_senders(frame);
    iproc_actors *receivers = iproc_frame_receivers(frame);
    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    int64_t ix_begin = iproc_frame_istatic(frame, 0);
    int64_t nstatic = iproc_frame_nstatic(frame);
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

