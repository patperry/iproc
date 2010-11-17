ACLOCAL_AMFLAGS        = -I m4 --install

lib_LTLIBRARIES        = libla.la
include_HEADERS        = libla/ieee754.h \
                         libla/index.h \
                         libla/vector.h

libla_la_SOURCES       = libla/blas-private.h \
                         libla/ieee754.h \
                         libla/ieee754.c \
                         libla/index.h \
                         libla/vector.h \
                         libla/vector.c

libla_la_CPPFLAGS      = $(GLIB_CFLAGS) \
                         $(AM_CPPFLAGS)

libla_la_LIBADD        = $(BLAS_LIBS) \
                         $(GLIB_LIBS) \
                         $(LIBS) \
                         $(FLIBS)

TESTS                  = tests/lacheck
check_PROGRAMS         = tests/lacheck

tests_lacheck_SOURCES  = tests/ieee754.c \
                         tests/main.c \
                         tests/vector.c

tests_lacheck_CPPFLAGS = $(CHECK_CFLAGS) \
                         $(GLIB_CFLAGS) \
                         $(AM_CPPFLAGS)

tests_lacheck_LDADD    = $(top_builddir)/libla.la \
                         $(CHECK_LIBS) \
                         $(GLIB_LIBS)
