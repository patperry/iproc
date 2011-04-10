#ifndef _MACCONFIG_H
#define _MACCONFIG_H
     /*
      * Modified version of "pyconfig.h" from Python-3.2.
      *
      * This file moves some of the autoconf magic to compile-time
      * when building on MacOSX. This is needed for building 4-way
      * universal binaries and for 64-bit universal binaries because
      * the values redefined below aren't configure-time constant but
      * only compile-time constant in these scenarios.
      */

#if defined(__APPLE__)


# undef HAVE_INTTYPES_H
# undef HAVE_STDBOOL_H
# undef HAVE_STDINT_H
# undef HAVE_SYS_TYPES_H

# undef HAVE__BOOL
# undef HAVE_INTPTR_T

# undef SIZEOF_INT
# undef SIZEOF_LONG
# undef SIZEOF_LONG_LONG
# undef SIZEOF_SIZE_T
# undef SIZEOF_VOID_P
# undef SIZEOF_INTPTR_T

# undef WORDS_BIGENDIAN
# undef DOUBLE_IS_ARM_MIXED_ENDIAN_IEEE754
# undef DOUBLE_IS_BIG_ENDIAN_IEEE754
# undef DOUBLE_IS_LITTLE_ENDIAN_IEEE754



# define HAVE_INTTYPES_H   1
# define HAVE_STDBOOL_H    1
# define HAVE_STDINT_H     1
# define HAVE_SYS_TYPES_H  1


# define HAVE__BOOL    1
# define HAVE_INTPTR_T 1

# ifdef __LP64__
#  define SIZEOF_INT              4
#  define SIZEOF_LONG             8
#  define SIZEOF_LONG_LONG        8
#  define SIZEOF_SIZE_T           8
#  define SIZEOF_VOID_P           8
#  define SIZEOF_INTPTR_T         8
# else
#  define SIZEOF_INT              4
#  define SIZEOF_LONG             4
#  define SIZEOF_LONG_LONG        8
#  define SIZEOF_SIZE_T           4
#  define SIZEOF_VOID_P           4
#  define SIZEOF_INTPTR_T         4
# endif


# ifdef __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
#  define DOUBLE_IS_BIG_ENDIAN_IEEE754
# else
#  define DOUBLE_IS_LITTLE_ENDIAN_IEEE754
# endif


#endif /* defined(__APPLE__) */

#endif /* _MACCONFIG_H */
