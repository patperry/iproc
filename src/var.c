#include "port.h"
#include "xalloc.h"
#include <stdlib.h>
#include <string.h>
#include "var.h"

static size_t compute_size(const size_t *dims, size_t rank);


void var_meta_init(struct var_meta *meta, const char *name, enum var_type type,
		   const size_t *dims, size_t rank)
{
	meta->name = xstrdup(name);
	meta->type = type;
	meta->rank = rank;
	memcpy(meta->dims, dims, meta->rank * sizeof(size_t));
	meta->size = compute_size(meta->dims, meta->rank);
}


void var_meta_deinit(struct var_meta *meta)
{
	free((void*)(meta->name));
}


size_t compute_size(const size_t *dims, size_t rank)
{
	size_t size = 1;
	size_t i;

	for (i = 0; i < rank; i++) {
		size *= dims[i];
	}

	return size;
}

