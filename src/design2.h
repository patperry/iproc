#ifndef DESIGN2_H
#define DESIGN2_H

#include "blas.h"
#include "sblas.h"
#include <stdarg.h>


struct design2 {
	struct frame *frame;
	size_t count1, count2;

	struct dmatrix traits;
	struct var2 **trait_vars;
	size_t ntrait, ntrait_max;

	size_t tvar_dim;
	struct tvar2 **tvars;
	size_t ntvar, ntvar_max;
	struct vpattern pat_buf;
	
	size_t *ir;
	size_t *jc;
	double *dx;	// transpose of dX[t]
	size_t nnz, nnz_max;	
	
	struct design2_observer *observers;
	size_t nobs, nobs_max;
};

struct var2 {
	const struct design2 *design;
	enum var_type type;
	const char *name;
	size_t dim;
	size_t index;
};

struct tvar2 {
	struct var2 var;
	const struct tvar2_type *type;
	void *udata;
};

struct tvar2_type {
	void (*init) (struct tvar2 *tv, const struct design2 *d, va_list ap);
	void (*deinit) (struct tvar2 * tv, const struct design2 *d);
};

struct design2_callbacks {
	void (*update) (void *udata, struct design2 *d, size_t i, size_t j,
			    const double *delta, const struct vpattern *pat);
	void (*update_var) (void *udata, struct design2 *d, const struct var2 *v, size_t i, size_t j,
			    const double *delta, const struct vpattern *pat);
};

struct design2_observer {
	void *udata;
	struct design2_callbacks callbacks;
};

struct coefs2 {
	double *all;
	double *traits;
	double *tvars;
	size_t dim;
};


void design2_init(struct design2 *d, struct frame *f, size_t count1, size_t count2);
void design2_deinit(struct design2 *d);

/* observers */
void design2_add_observer(struct design2 *d, void *udata,
			const struct design2_callbacks *callbacks);
void design2_remove_observer(struct design2 *d, void *udata);

/* properties */
static inline struct frame *design2_frame(const struct design2 *d);
static inline size_t design2_count1(const struct design2 *d);
static inline size_t design2_count2(const struct design2 *d);
const struct var2 *design2_var(const struct design2 *d, const char *name);

/* traits */
static inline size_t design2_trait_dim(const struct design2 *d);
static inline const struct dmatrix *design2_traits(const struct design2 *d);
const char *design2_trait_name(const struct design2 *d, size_t k);
const struct var2 *design2_add_trait(struct design2 *d, const char *name, const double *x);
void design2_add_traits(struct design2 *d, size_t ntrait, const char * const *names, const struct dmatrix *x);

//void design2_add_kron(struct design2 *d, const struct var *i,
//			const struct var *j);

void design2_traits_mul(double alpha, const struct design2 *d, size_t i,
			const double *x, double beta, double *y);
void design2_traits_tmul(double alpha, const struct design2 *d, size_t i, const double *x, double beta, double *y);
void design2_traits_axpy(double alpha, const struct design2 *d, size_t i, size_t j, double *y);


/* tvars */
static inline size_t design2_tvar_dim(const struct design2 *d);
static inline const double *design2_tvars(const struct design2 *d, size_t i, size_t j);
const char *design2_tvar_name(const struct design2 *d, size_t k);
const struct var2 *design2_add_tvar(struct design2 *d, const char *name, const struct tvar2_type *type, ...);

static inline void design2_tvars_get(const struct design2 *d, size_t i, const double **dxp, const size_t **jp,
				     size_t *nzp);

void design2_tvars_mul(double alpha, const struct design2 *d, size_t i,
		      const double *x, double beta, double *y);
void design2_tvars_tmul(double alpha, const struct design2 *d, size_t i, const double *x, double beta, double *y);
void design2_tvars_axpy(double alpha, const struct design2 *d, size_t i, size_t j, double *y);


/* coefs2 */
void coefs2_init(struct coefs2 *c, const struct design2 *d);
void coefs2_deinit(struct coefs2 *c);

void design2_mul(double alpha, const struct design2 *d, size_t i,
		 const struct coefs2 *c, double beta, double *y);
void design2_tmul(double alpha, const struct design2 *d, size_t i,
		  const double *x, double beta, struct coefs2 *c);
void design2_axpy(double alpha, const struct design2 *d, size_t i, size_t j,
		  struct coefs2 *c);


/* internal functions (for use by tvar callbacks) */
void design2_update(struct design2 *d, const struct var2 *v, size_t i, size_t j,
		    const double *delta, const struct vpattern *pat);


/* inline function definitions */
struct frame *design2_frame(const struct design2 *d)
{
	return d->frame;
}


size_t design2_count1(const struct design2 *d)
{
	return d->count1;
}

size_t design2_count2(const struct design2 *d)
{
	return d->count2;
}

size_t design2_trait_dim(const struct design2 *d)
{
	return d->ntrait;
}


const struct dmatrix *design2_traits(const struct design2 *d)
{
	return &d->traits;
}


size_t design2_tvar_dim(const struct design2 *d)
{
	return d->tvar_dim;
}


const double *design2_tvars(const struct design2 *d, size_t i, size_t j)
{
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));
	
	size_t ir0 = d->ir[i];
	size_t ir1 = d->ir[i+1];
	size_t nz = ir1 - ir0;
	size_t *indx = d->jc + ir0;
	struct vpattern pat = vpattern_make(indx, nz);
	
	ptrdiff_t ix = vpattern_find(&pat, j);

	if (ix < 0)
		return NULL;

	return d->dx + (ir0 + ix) * d->tvar_dim;
}


void design2_tvars_get(const struct design2 *d, size_t i, const double **dxp, const size_t **jp, size_t *nzp)
{
	assert(i < design2_count1(d));

	size_t ir0 = d->ir[i];
	size_t ir1 = d->ir[i+1];
	
	*nzp = ir1 - ir0;
	*jp = d->jc + ir0;
	*dxp = d->dx + ir0 * d->tvar_dim;
}


#endif /* DESIGN2_H */
