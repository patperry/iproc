#ifndef _VARS_H
#define _VARS_H

#include "design.h"
#include "frame.h"

enum var_class {
	VAR_RECV_VAR,
	VAR_SEND_VAR
};

struct design_var;
struct frame_var;

struct var_type {
	enum var_class var_class;

	void (*init) (struct design_var * dv, const struct design * d);
	void (*deinit) (struct design_var * dv);

	void (*frame_init) (struct frame_var * fv, struct frame * f);
	void (*frame_deinit) (struct frame_var * fv);

	struct frame_callbacks callbacks;
};

struct design_var {
	const struct var_type *type;
	ssize_t dim;
	ssize_t index;
	void *udata;
};

struct frame_var {
	struct design_var *design;
	void *udata;
};

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i] = 1{ j -> i in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct var_type *SEND_VAR_IRECV;

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i] = #{ j -> i in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct var_type *SEND_VAR_NRECV;

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i,j] = 1{ j -> i in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct var_type *RECV_VAR_IRECV;

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i,j] = #{ j -> i in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct var_type *RECV_VAR_NRECV;

extern const struct var_type *RECV_VAR_NSEND;


#endif /* _VARS_H */
