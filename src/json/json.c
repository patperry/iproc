#include "port.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "json.h"

#define YG(gen) \
	do { \
		if ((err = gen) != yajl_gen_status_ok) \
			return err; \
	} while (0)
#define YSTR(s) ((const unsigned char *)s)

#define VECTOR_DIM  "dim"
#define VECTOR_DATA "data"

#define SVECTOR_DIM   "dim"
#define SVECTOR_COUNT "count"
#define SVECTOR_PATTERN "pattern"
#define SVECTOR_DATA "data"

#define MATRIX_NROW "nrow"
#define MATRIX_NCOL "ncol"
#define MATRIX_DATA "data"

yajl_gen_status yajl_gen_ieee754(yajl_gen hand, double val)
{
	switch (fpclassify(val)) {
	case FP_INFINITE:
		return yajl_gen_double(hand, val > 0 ? DBL_MAX : -DBL_MAX);
	case FP_NAN:
		return yajl_gen_null(hand);
	default:
		return yajl_gen_double(hand, val);
	}
}

yajl_gen_status yajl_gen_vector(yajl_gen hand, const struct vector * x)
{
	assert(x);

	yajl_gen_status err = yajl_gen_status_ok;
	ssize_t i, n = vector_dim(x);

	YG(yajl_gen_array_open(hand));
	for (i = 0; i < n; i++) {
		double val = vector_item(x, i);
		YG(yajl_gen_ieee754(hand, val));
	}
	YG(yajl_gen_array_close(hand));

	return err;
}

yajl_gen_status yajl_gen_svector(yajl_gen hand, const struct svector *x)
{
	assert(x);
	yajl_gen_status err = yajl_gen_status_ok;
	
	size_t inz, nnz = svector_count(x);
	const size_t *index = svector_index_ptr(x);
	const double *data = svector_data_ptr(x);
	
	YG(yajl_gen_map_open(hand));
	
	YG(yajl_gen_string(hand, YSTR(SVECTOR_DIM), strlen(SVECTOR_DIM)));
	YG(yajl_gen_integer(hand, svector_dim(x)));
	
	YG(yajl_gen_string(hand, YSTR(SVECTOR_COUNT), strlen(SVECTOR_COUNT)));
	YG(yajl_gen_integer(hand, nnz));

	YG(yajl_gen_string(hand, YSTR(SVECTOR_PATTERN), strlen(SVECTOR_PATTERN)));
	YG(yajl_gen_array_open(hand));
	for (inz = 0; inz < nnz; inz++) {
		YG(yajl_gen_integer(hand, index[inz] + 1));
	}
	YG(yajl_gen_array_close(hand));
	
	YG(yajl_gen_string(hand, YSTR(SVECTOR_DATA), strlen(SVECTOR_DATA)));
	YG(yajl_gen_array_open(hand));
	for (inz = 0; inz < nnz; inz++) {
		YG(yajl_gen_ieee754(hand, data[inz]));
	}
	YG(yajl_gen_array_close(hand));
	
	YG(yajl_gen_map_close(hand));
	
	return err;
}

yajl_gen_status yajl_gen_matrix(yajl_gen hand, const struct matrix * a)
{
	assert(a);

	yajl_gen_status err = yajl_gen_status_ok;
	ssize_t i, m = matrix_nrow(a);
	ssize_t j, n = matrix_ncol(a);

	YG(yajl_gen_map_open(hand));

	YG(yajl_gen_string(hand, YSTR(MATRIX_NROW), strlen(MATRIX_NROW)));
	YG(yajl_gen_integer(hand, m));

	YG(yajl_gen_string(hand, YSTR(MATRIX_NCOL), strlen(MATRIX_NCOL)));
	YG(yajl_gen_integer(hand, n));

	YG(yajl_gen_string(hand, YSTR(MATRIX_DATA), strlen(MATRIX_DATA)));
	YG(yajl_gen_array_open(hand));
	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			double val = matrix_item(a, i, j);
			YG(yajl_gen_ieee754(hand, val));
		}
	}
	YG(yajl_gen_array_close(hand));

	YG(yajl_gen_map_close(hand));

	return err;

}
