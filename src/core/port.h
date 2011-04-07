#ifndef _PORT_H
#define _PORT_H

/* Portabiliy types and macros */


#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

/* http://www.gnu.org/software/autoconf/manual/autoconf.html#index-HAVE_005fSTDBOOL_005fH-624 */
#ifdef HAVE_STDBOOL_H
# include <stdbool.h>
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



#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h> /* ssize_t */
#endif

/* http://www.mail-archive.com/bug-gnulib@gnu.org/msg02492.html */
#include <limits.h>
#ifndef SSIZE_MAX
# define SSIZE_MAX  (~ (-1L << (SIZEOF_SIZE_T * CHAR_BIT - 1)))
#endif


#endif /* _PORT_H */
