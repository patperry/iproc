#ifndef VAR_H
#define VAR_H

#include <stdarg.h>
#include <stddef.h>
#include "deltaset.h"
#include "uintset.h"
#include "history.h"


struct var2 {
	struct var_meta meta;
	struct design2 *design;
	size_t index;
};

struct kvar2 {
	struct var2 var;
	double *xi, *xj;
	size_t dimi, dimj;
};

struct tvar2 {
	struct var2 var;
	const struct tvar2_type *type;
	void *udata;
};

struct tvar2_type {
	void (*init) (struct tvar2 *tv, const char *name, struct history *h, va_list ap);
	void (*deinit) (struct tvar2 * tv, struct history *h);
};


extern const struct tvar_type *VAR_IRECVTOT;
extern const struct tvar_type *VAR_ISENDTOT;

extern const struct tvar_type *VAR_NRECVTOT;
extern const struct tvar_type *VAR_NSENDTOT;


/* Indicator 1{ j -> i in (-Infty, t) } */
extern const struct tvar2_type *VAR2_IRECV;

/* Indicator 1{ i -> j in (-Infty, t) } */
extern const struct tvar2_type *VAR2_ISEND;

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i,j] = #{ j -> i in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct tvar2_type *VAR2_NRECV;

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i,j] = #{ i -> j in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct tvar2_type *VAR2_NSEND;

extern const struct tvar2_type *VAR2_IRECV2;
extern const struct tvar2_type *VAR2_NRECV2;

extern const struct tvar2_type *VAR2_ISEND2;
extern const struct tvar2_type *VAR2_NSEND2;

extern const struct tvar2_type *VAR2_ISIB;
extern const struct tvar2_type *VAR2_NSIB;

extern const struct tvar2_type *VAR2_ICOSIB;
extern const struct tvar2_type *VAR2_NCOSIB;


#endif /* VAR_H */
