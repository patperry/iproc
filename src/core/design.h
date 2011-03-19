#ifndef _IPROC_DESIGN_H
#define _IPROC_DESIGN_H

#include <stdbool.h>
#include <stdint.h>
#include "actors.h"
#include "history.h"
#include "matrix.h"
#include "refcount.h"
#include "svector.h"
#include "vector.h"


/* I. Introduction
 * ---------------
 *
 * A "design" holds the covariates for the interaction process.  It can
 * be thought of as a J-by-p matrix, indexed by time and sender.
 * 
 *         X[t,i] = [ x[t,i,1], ..., x[t,i,J] ]^T;
 * 
 * that is, the j-th row of X[t,i] is x[t,i,j].
 *
 * We represent the design as
 *
 *         X[t,i] = X[0,i] + dX[t,i].
 *
 *
 * II. Static part
 * ---------------
 * The matrix X[0,i] is referred to as the static design matrix.  This matrix
 * has a special structure which enables efficent multiplication.  Let
 * s[i] be the vector of sender i's traits and let
 * R = [ r[1], ..., r[J] ]^T be the matrix of the receivers' traits.  Then,
 *
 *         X[0,i] = [ kronecker(R, s[i]^T)  0 ]
 *
 * Suppose that s is p-dimensional and R is J-by-q.  If Y is a p-by-q
 * matrix, then
 *
 *         kronecker(R, s[i]^T) * vec(Y) = R * (Y^T * s[i]).
 *
 * If y is J-dimensional, then
 *
 *         (kronecker(R, s[i]^T))^T * y = vec(s[i]^T * (R^T * y)^T).
 *
 *
 * III. Dynamic part
 * -----------------
 * The dynamic part dX[t,i] is computed by client-supplied functions.  These
 * functions in general will depend on the history of the process.  The only
 * restriction is that dX[0,i] = 0.  Write
 *
 *         dX[t,i] = [  0  dX{1}[t,i]  dX{2}[t,i]  ...  dX{K}[t,i] ]
 *
 * The presumption is that each dX{k}[t,i] is sparse, in the sence that
 * (a) most of the rows of dX{k}[t,i] are zero, and (b) the nonzero rows
 * of dX{k}[t,i] are sparse vectors.  For each dX{k}[t,i], the client
 * must provide the following:
 *
 *         A. `dim` : the number of columns in dX{k}[t,i];
 *
 *         B. `get_dx(t, i, dX, offset, udata)` : a function that
 *            copies the rows of dX{k}[t,i] into dX[:,offset:(offset+dim)];
 *
 *         C. `free_udata` : a funciton to free the extra client data;
 *
 *         D. `udata` : the extra client data.
 *
 * The matrix dX{k}[t,i] is represented as an array of { j, dx[t,i,j] }
 * pairs, where dx[t,i,j] is a sparse vector.
 */

typedef struct _iproc_design     iproc_design;
typedef struct _iproc_design_var iproc_design_var;
typedef struct _iproc_design_ctx iproc_design_ctx;

typedef struct _iproc_design_dx  iproc_design_dx; // private


struct _iproc_design {
    iproc_actors  *senders;
    iproc_actors  *receivers;
    iproc_array   *vars;
    int64_t        dim0;
    int64_t        dim;
    
    iproc_array   *ctxs;
    iproc_array   *svectors;
    iproc_refcount refcount;    
};

/* dX[t,i] */
struct _iproc_design_ctx {
    iproc_design  *design;
    iproc_history *history;
    int64_t        isend;
    iproc_array   *dxs;
    iproc_refcount refcount;
};

/* dX{k}[t,i] */
struct _iproc_design_var {
    int64_t dim;
    void (*get_dxs) (iproc_design_var *var, iproc_design_ctx *ctx, int64_t offset);
    void (*free)    (iproc_design_var *var);
};


/* dx[t,i,j] */
struct _iproc_design_dx {
    int64_t        jrecv; // NOTE: layout is important here
    iproc_svector *dx;
};

void               iproc_design_var_init  (iproc_design_var *var,
                                           int64_t           dim,
                                           void (*get_dxs) (iproc_design_var *,
                                                            iproc_design_ctx *ctx,
                                                            int64_t),
                                           void (*free)    (iproc_design_var *));


iproc_design *     iproc_design_new       (iproc_actors     *senders,
                                           iproc_actors     *receivers);
iproc_design *     iproc_design_ref       (iproc_design     *design);
void               iproc_design_unref     (iproc_design     *design);

int64_t            iproc_design_dim       (iproc_design     *design);

// NOTE: a call to `append` invalidates existing design_ctxs
void               iproc_design_append    (iproc_design     *design,
                                           iproc_design_var *var);


void               iproc_design_mul0      (double            alpha,
                                           iproc_trans       trans,
                                           iproc_design     *design,
                                           int64_t           isend,
                                           iproc_vector     *x,
                                           double            beta,
                                           iproc_vector     *y);
void               iproc_design_muls0     (double            alpha,
                                           iproc_trans       trans,
                                           iproc_design     *design,
                                           int64_t           isend,
                                           iproc_svector    *x,
                                           double            beta,
                                           iproc_vector     *y);


iproc_design_ctx * iproc_design_ctx_new   (iproc_design     *design,
                                           int64_t           isend,
                                           iproc_history    *h);
iproc_design_ctx * iproc_design_ctx_ref   (iproc_design_ctx *ctx);
void               iproc_design_ctx_unref (iproc_design_ctx *ctx);


void               iproc_design_ctx_mul   (double            alpha,
                                           iproc_trans       trans,
                                           iproc_design_ctx *ctx,
                                           iproc_vector     *x,
                                           double            beta,
                                           iproc_vector     *y);
void               iproc_design_ctx_muls  (double            alpha,
                                           iproc_trans       trans,
                                           iproc_design_ctx *ctx,
                                           iproc_svector    *x,
                                           double            beta,
                                           iproc_vector     *y);


iproc_svector *    iproc_design_ctx_dx    (iproc_design_ctx *ctx,
                                           int64_t           jrecv,
                                           bool              null_ok);

/* number of nonzero rows in dX */
int64_t            iproc_design_ctx_nnz   (iproc_design_ctx *ctx);

/* inz-th nonzero row of dX */
iproc_svector *    iproc_design_ctx_nz    (iproc_design_ctx *ctx,
                                           int64_t           inz,
                                           int64_t          *jrecv);


void               iproc_design_ctx_dmul  (double            alpha,
                                           iproc_trans       trans,
                                           iproc_design_ctx *ctx,
                                           iproc_vector     *x,
                                           double            beta,
                                           iproc_svector    *y);
void               iproc_design_ctx_dmuls (double            alpha,
                                           iproc_trans       trans,
                                           iproc_design_ctx *ctx,
                                           iproc_svector    *x,
                                           double            beta,
                                           iproc_svector    *y);


int64_t          iproc_design_nsender     (iproc_design     *design);
int64_t          iproc_design_nreceiver   (iproc_design     *design);
iproc_actors *   iproc_design_senders     (iproc_design     *design);
iproc_actors *   iproc_design_receivers   (iproc_design     *design);

int64_t          iproc_design_nstatic     (iproc_design     *design);
int64_t          iproc_design_istatic     (iproc_design     *design,
                                           int64_t           i);
int64_t          iproc_design_ndynamic    (iproc_design     *design);
int64_t          iproc_design_idynamic    (iproc_design     *design,
                                           int64_t           i);



#endif /* _IPROC_DESIGN_H */
