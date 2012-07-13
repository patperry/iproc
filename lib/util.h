#ifndef UTIL_H
#define UTIL_H

#include <stddef.h>

void free2(void **ptrs, size_t len);
char **xstrdup2(const char *const *strs, size_t len);

void vector_shift(size_t n, double alpha, double *x);
void vector_exp(size_t n, double *x);
double vector_logsumexp(size_t n, double *x);
double vector_max(size_t n, double *x);


#endif /* UTIL_H */
