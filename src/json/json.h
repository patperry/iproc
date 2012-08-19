#ifndef _JSON_H
#define _JSON_H

#include "blas.h"
#include <yajl/yajl_gen.h>

yajl_gen_status yajl_gen_ieee754(yajl_gen hand, double val);
yajl_gen_status yajl_gen_vector(yajl_gen hand, size_t n, const double *x);
yajl_gen_status yajl_gen_matrix(yajl_gen hand, size_t m, size_t n, const double *a);

#endif /* _JSON_H */
