#ifndef UTIL_H
#define UTIL_H

#include <stddef.h>

void free2(void **ptrs, size_t len);
char **xstrdup2(const char *const *strs, size_t len);



#endif /* UTIL_H */
