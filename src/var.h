#ifndef VAR_H
#define VAR_H

#include <stdarg.h>


#define VAR_RANK_MAX 8


enum var_type {
	VAR_TYPE_TRAIT,
	VAR_TYPE_TVAR
};

struct var {
	struct design *design;
	enum var_type type;
	const char *name;
	size_t dims[VAR_RANK_MAX];
	size_t rank;
	size_t size;
	size_t index;
};

struct tvar {
	struct var var;
	const struct tvar_type *type;
	void *udata;
};

struct tvar_type {
	void (*init) (struct tvar *tv, struct design *d, va_list ap);
	void (*deinit) (struct tvar * tv, struct design *d);
};


#endif /* VAR_H */
