#include "port.h"
#include "coreutil.h"
#include "xalloc.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "var.h"


static size_t compute_size(const size_t *dims, size_t rank);
static void get_indices(const size_t *dims, size_t rank, size_t i, size_t *ix);


void var_meta_init(struct var_meta *meta, enum var_type type, const char *name,
		   const size_t *dims, size_t rank)
{
	meta->type = type;
	meta->name = name ? xstrdup(name) : NULL;
	meta->rank = rank;
	memcpy(meta->dims, dims, meta->rank * sizeof(size_t));
	meta->size = compute_size(meta->dims, meta->rank);
}


void var_meta_deinit(struct var_meta *meta)
{
	free((void*)(meta->name));
}


char *alloc_var_name(const struct var_name_fmt *fmt,
		     const struct var_meta *meta, size_t i)
{
	size_t len = snprint_var_name(NULL, 0, fmt, meta, i);
	char *str = xmalloc(len + 1);
	sprint_var_name(str, fmt, meta, i);
	return str;
}


int sprint_var_name(char *str, const struct var_name_fmt *fmt,
		    const struct var_meta *meta, size_t i)
{
	assert(i < meta->size);
	return snprint_var_name(str, SIZE_MAX, fmt, meta, i);
}

#define PRINT(f, a) \
	do { \
		size_t l = (size_t)snprintf(str, size, (f), (a)); \
		len += l; \
		if (l >= size) { \
			str = NULL; \
			size = 0; \
		} else { \
			str += l; \
			size -= l; \
		} \
	} while (0)

int snprint_var_name(char *str, size_t size, const struct var_name_fmt *fmt,
		     const struct var_meta *meta, size_t i)
{
	assert(i < meta->size);

	struct var_name_fmt fmt0 = VAR_NAME_FMT0;
	if (!fmt)
		fmt = &fmt0;

	size_t j, rank = meta->rank;
	size_t indices[VAR_RANK_MAX];
	size_t len = 0;

	get_indices(meta->dims, meta->rank, i, indices);

	if (fmt->ione) {
		for (j = 0; j < rank; j++) {
			indices[j]++;
		}
	}

	PRINT("%s", meta->name ? meta->name : "(NULL)");

	if (rank > 0) {
		PRINT("%s", fmt->open);
		PRINT(fmt->ifmt, indices[0]);

		for (j = 1; j < rank; j++) {
			PRINT("%s", fmt->sep);
			PRINT(fmt->ifmt, indices[j]);
		}

		PRINT("%s", fmt->close);
	}

	return (int)len;
}

#undef PRINT


size_t compute_size(const size_t *dims, size_t rank)
{
	size_t size = 1;
	size_t i;

	for (i = 0; i < rank; i++) {
		size *= dims[i];
	}

	return size;
}

void get_indices(const size_t *dims, size_t rank, size_t i, size_t *indices)
{
	while (rank > 0) {
		size_t td = compute_size(dims + 1, rank - 1);
		*indices++ = i / td;

		dims++;
		rank--;
		i %= td;
	}
}
