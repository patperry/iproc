#ifndef F77_H
#define F77_H

#include <stddef.h>

#ifndef F77_FUNC
# define F77_FUNC(f) f ## _
#endif

#ifndef f77int
# ifdef HAVE_BLAS64
#  define f77int long long int
# else
#  define f77int int
# endif
#endif

void f77int_pack(ptrdiff_t *base, size_t len);
void f77int_unpack(ptrdiff_t *base, size_t len);

#endif /* F77_H */
