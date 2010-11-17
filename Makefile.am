ACLOCAL_AMFLAGS       = -I m4 --install
                      
lib_LTLIBRARIES       = libla.la
include_HEADERS       = libla/ieee754.h \
                        libla/index.h \
                        libla/vector.h
                      
libla_la_SOURCES      = libla/blas-private.h \
                        libla/ieee754.h \
                        libla/ieee754.c \
                        libla/index.h \
                        libla/vector.h \
                        libla/vector.c
                      
libla_la_CPPFLAGS     = $(GLIB_CFLAGS) \
                        $(AM_CPPFLAGS)
                      
libla_la_LIBADD       = $(BLAS_LIBS) \
                        $(GLIB_LIBS) \
                        $(LIBS) \
                        $(FLIBS)

TESTS                 = tests/vector
check_PROGRAMS        = tests/vector

tests_vector_SOURCES  = tests/vector.c

tests_vector_CPPFLAGS = $(CHECK_CFLAGS) \
                        $(AM_CPPFLAGS)

tests_vector_LDADD    = $(top_builddir)/libla.la \
                        $(CHECK_LIBS)

