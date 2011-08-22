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

	void (*init) (struct design_var * dv, const struct design * d,
		      void *params);
	void (*deinit) (struct design_var * dv);

	void (*frame_init) (struct frame_var * fv, struct frame * f);
	void (*frame_deinit) (struct frame_var * fv);

	struct frame_callbacks callbacks;
};

struct design_var {
	const struct var_type *type;
	char **names;
	ssize_t dim;
	ssize_t dyn_index;
	void *udata;
};

struct frame_var {
	struct design_var *design;
	void *udata;
};

/* Indicator 1{ j -> i in (-Infty, t) } */
extern const struct var_type *RECV_VAR_IRECV;

/* Indicator 1{ i -> j in (-Infty, t) } */
extern const struct var_type *RECV_VAR_ISEND;

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i,j] = #{ j -> i in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct var_type *RECV_VAR_NRECV;

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i,j] = #{ i -> j in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct var_type *RECV_VAR_NSEND;

extern const struct var_type *RECV_VAR_IRECV2;
extern const struct var_type *RECV_VAR_NRECV2;

extern const struct var_type *RECV_VAR_ISEND2;
extern const struct var_type *RECV_VAR_NSEND2;

extern const struct var_type *RECV_VAR_ISIB;
extern const struct var_type *RECV_VAR_NSIB;

extern const struct var_type *RECV_VAR_ICOSIB;
extern const struct var_type *RECV_VAR_NCOSIB;


char **var_names_alloc(char *name, size_t len, ssize_t n);
char **var_names_alloc2(char *name, size_t len, ssize_t m, ssize_t n);
void var_names_free(char **names);



#endif /* _VARS_H */
