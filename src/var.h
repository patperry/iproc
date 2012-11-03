#ifndef VAR_H
#define VAR_H

#include <stdarg.h>
#include <stddef.h>
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
};



void var_meta_init(struct var_meta *meta, const char *name, enum var_type type,
		   const size_t *dims, size_t rank);
void var_meta_deinit(struct var_meta *meta);



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

struct tvar_type {
	void (*init) (struct tvar *tv, const char *name, struct history *h, va_list ap);
	void (*deinit) (struct tvar * tv, struct history *h);
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

#endif /* VAR_H */
