#ifndef _VARS_H
#define _VARS_H

/* forward declarations */
struct design;
struct design_var;
struct frame;			
struct frame_event;
struct frame_var;

enum var_class {
	VAR_RECV_VAR,
	VAR_SEND_VAR
};

struct var_type {
	enum var_class var_class;
	uint8_t event_mask;
	
	void (*init) (struct design_var *dv, const struct design *d);
	void (*deinit) (struct design_var *dv);
	
	void (*frame_init) (struct frame_var *fv, struct frame *f);
	void (*frame_deinit) (struct frame_var *fv);
	void (*frame_clear) (struct frame_var *fv);
	
	void (*handle_event) (struct frame_var * fv,
			      const struct frame_event * e, struct frame *f);
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


#endif /* _VARS_H */

