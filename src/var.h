#ifndef VAR_H
#define VAR_H

#include <stddef.h>


#define VAR_NAME_OPEN0	"["
#define VAR_NAME_CLOSE0	"]"
#define VAR_NAME_SEP0	","
#define	VAR_NAME_IFMT0	"%zd"
#define VAR_NAME_IONE0	1

#define VAR_NAME_FMT0	{ VAR_NAME_OPEN0, VAR_NAME_CLOSE0, VAR_NAME_SEP0, \
			  VAR_NAME_IFMT0, VAR_NAME_IONE0 }

#define VAR_RANK_MAX 8


enum var_type {
	VAR_TYPE_TRAIT,
	VAR_TYPE_TVAR
};


struct var_name_fmt {
	const char *open;
	const char *close;
	const char *sep;
	const char *ifmt;
	int ione;
};


struct var_meta {
	enum var_type type;
	const char *name;
	size_t dims[VAR_RANK_MAX];
	size_t rank;
	size_t size;
};

void var_meta_init(struct var_meta *meta, enum var_type type, const char *name,
		   const size_t *dims, size_t rank);
void var_meta_deinit(struct var_meta *meta);


char *alloc_var_name(const struct var_name_fmt *fmt,
		     const struct var_meta *meta, size_t i);
int sprint_var_name(char *str, const struct var_name_fmt *fmt,
		    const struct var_meta *meta, size_t i);
int snprint_var_name(char *str, size_t size, const struct var_name_fmt *fmt,
		     const struct var_meta *meta, size_t i);


#endif /* VAR_H */
