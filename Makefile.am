ACLOCAL_AMFLAGS      = -I m4 --install

lib_LTLIBRARIES      = libiproc.la
include_HEADERS      = iproc/ieee754.h \
                       iproc/vector.h

libiproc_la_SOURCES  = iproc/blas-private.h \
                       iproc/ieee754.h \
                       iproc/ieee754.c \
                       iproc/vector.h \
                       iproc/vector.c

libiproc_la_CPPFLAGS = $(GLIB_CFLAGS) \
                       $(AM_CPPFLAGS)

libiproc_la_LIBADD   = $(BLAS_LIBS) \
                       $(GLIB_LIBS) \
                       $(LIBS) \
                       $(FLIBS)

TESTS                = check-iproc
check_PROGRAMS       = check-iproc

check_iproc_SOURCES  = tests/checkutils.h \
                       tests/ieee754.c \
                       tests/main.c \
                       tests/vector.c

check_iproc_CPPFLAGS = -Itests \
                       $(CHECK_CFLAGS) \
                       $(GLIB_CFLAGS) \
                       $(AM_CPPFLAGS)

check_iproc_LDADD    = $(top_builddir)/libiproc.la \
                       $(CHECK_LIBS) \
                       $(GLIB_LIBS)
