ACLOCAL_AMFLAGS = -I m4 --install

include_HEADERS = \
	libla/index.h \
	libla/vector.h

noinst_HEADERS = \
	libla/blas-private.h

lib_LTLIBRARIES  = libla.la
libla_la_SOURCES = \
	libla/index.h \
	libla/vector.h \
	libla/vector.c

libla_la_CPPFLAGS = $(GLIB_CFLAGS) $(AM_CPPFLAGS)
libla_la_LIBADD = $(BLAS_LIBS) $(GLIB_LIBS) $(LIBS) $(FLIBS)
