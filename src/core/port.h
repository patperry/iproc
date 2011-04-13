#ifndef _PORT_H
#define _PORT_H

/* Portabiliy types and macros
 *
 * This file should be included in all ".c" files.
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#if defined(__APPLE__)
# include "macconfig.h"
#endif

#ifdef HAVE_INTTYPES_H
# include <inttypes.h>
#endif
#ifdef HAVE_STDBOOL_H
# include <stdbool.h>
#endif
#ifdef HAVE_STDINT_H
# include <stdint.h>
#endif
#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>		/* ssize_t */
#endif

#include <limits.h>

/* http://www.gnu.org/software/autoconf/manual/autoconf.html#index-HAVE_005fSTDBOOL_005fH-624 */
#ifdef HAVE_STDBOOL_H
#				/* bool defined in stdbool.h */
#else
# ifndef HAVE__BOOL
#  ifdef __cplusplus
typedef bool _Bool;
#  else
#   define _Bool signed char
#  endif
# endif
# define bool _Bool
# define false 0
# define true 1
# define __bool_true_false_are_defined 1
#endif

#ifdef HAVE_INTPTR_T
/* intptr_t defined in stdint.h or inttypes.h */
#elif SIZEOF_VOID_P <= SIZEOF_INT
typedef int intptr_t;
#elif SIZEOF_VOID_P <= SIZEOF_LONG
typedef long intptr_t;
#elif SIZEOF_VOID_P <= SIZEOF_LONG_LONG
typedef long long intptr_t;
#else
# error "Need a typedef for intptr_t in config.h"
#endif /* HAVE_UINTPTR_T */

#ifdef HAVE_INTTYPES_H
/* imaxabs already defined */
#elif  HAVE_IMAXABS
/* imaxabs already defined */
#else
intmax_t imaxabs(intmax_t i);
#endif

/* http://www.mail-archive.com/bug-gnulib@gnu.org/msg02492.html */
#ifndef SSIZE_MAX
# define SSIZE_MAX  (~ (-1L << (SIZEOF_SIZE_T * CHAR_BIT - 1)))
#endif
#ifndef SSIZE_MIN
# define SSIZE_MIN  (-SSIZE_MAX-1)
#endif

#ifndef SSIZE_FMT
# if SIZEOF_SIZE_T == SIZEOF_LONG
#  define SSIZE_FMT "ld"
# elif SIZEOF_SIZE_T == SIZEOF_LONG_LONG
#  define SSIZE_FMT "lld"
# elif SIZEOF_SIZE_T == SIZEOF_INT
#  define SSIZE_FMT "d"
# else
#  error "Cannot determine format string for ssize_t; \
define SSIZE_FMT in config.h"
# endif
#endif

#endif /* _PORT_H */
