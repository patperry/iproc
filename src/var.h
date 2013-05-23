#ifndef VAR_H
#define VAR_H

#include <stdarg.h>
#include <stddef.h>
#include "uintset.h"
#include "history.h"


#define VAR_RANK_MAX 8


enum var_type {
	VAR_TYPE_TRAIT,
	VAR_TYPE_TVAR
};

struct var_meta {
	const char *name;
	enum var_type type;
	size_t dims[VAR_RANK_MAX];
	size_t rank;
	size_t size;

	struct uintset changed;
	int cleared;
};



void var_meta_init(struct var_meta *meta, const char *name, enum var_type type,
		   const size_t *dims, size_t rank);
void var_meta_deinit(struct var_meta *meta);

#define VAR_NAME_OPEN0	"["
#define VAR_NAME_CLOSE0	"]"
#define VAR_NAME_SEP0	","
#define	VAR_NAME_IFMT0	"%zd"
#define VAR_NAME_IONE0	1

#define VAR_NAME_FMT0	{ VAR_NAME_OPEN0, VAR_NAME_CLOSE0, VAR_NAME_SEP0, \
			  VAR_NAME_IFMT0, VAR_NAME_IONE0 }

struct var_name_fmt {
	const char *open;
	const char *close;
	const char *sep;
	const char *ifmt;
	int ione;
};

char *alloc_var_name(const struct var_name_fmt *fmt,
		     const struct var_meta *meta, size_t i);
int sprint_var_name(char *str, const struct var_name_fmt *fmt,
		    const struct var_meta *meta, size_t i);
int snprint_var_name(char *str, size_t size, const struct var_name_fmt *fmt,
		     const struct var_meta *meta, size_t i);

struct var {
	struct var_meta meta;
	struct design *design;
	size_t index;
};

struct tvar {
	struct var var;
	const struct tvar_type *type;
	void *udata;
};

void var_change(struct var *v, size_t i);
void var_clear(struct var *v);
void var_delta_clear(struct var *v);


struct tvar_type {
	void (*init) (struct tvar *tv, const char *name, struct history *h, va_list ap);
	void (*deinit) (struct tvar * tv, struct history *h);
	void (*update) (struct tvar *tv);
};



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
