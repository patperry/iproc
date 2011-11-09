#ifndef PORT_H
#define PORT_H

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#if defined(__APPLE__)
# undef WORDS_BIGENDIAN
# ifdef __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#endif

#endif /* PORT_H */
