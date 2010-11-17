ACLOCAL_AMFLAGS    = -I m4 --install

lib_LTLIBRARIES    = libla.la
include_HEADERS    = libla/index.h \
                     libla/vector.h


libla_la_SOURCES   = libla/blas-private.h \
                     libla/index.h \
                     libla/vector.h \
                     libla/vector.c

libla_la_CPPFLAGS  = $(GLIB_CFLAGS) \
                     $(AM_CPPFLAGS)

libla_la_LIBADD    = $(BLAS_LIBS) \
                     $(GLIB_LIBS) \
                     $(LIBS) \
                     $(FLIBS)
