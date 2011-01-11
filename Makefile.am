ACLOCAL_AMFLAGS = -I m4 --install

lib_LTLIBRARIES = \
	libiproc.la

include_HEADERS = \
	iproc/actors.h \
	iproc/ieee754.h \
	iproc/loglik.h \
	iproc/matrix.h \
	iproc/model.h \
	iproc/vars.h \
	iproc/vector.h

libiproc_la_SOURCES = \
	iproc/actors.h \
	iproc/actors.c \
	iproc/array.h \
	iproc/array.c \
	iproc/blas-private.h \
	iproc/events.h \
	iproc/events.c \
	iproc/history.h \
	iproc/history.c \
	iproc/ieee754.h \
	iproc/ieee754.c \
	iproc/loglik.h \
	iproc/loglik.c \
	iproc/matrix.h \
	iproc/matrix.c \
	iproc/model.h \
	iproc/model.c \
	iproc/refcount.h \
	iproc/refcount.c \
	iproc/sloglik.h \
	iproc/sloglik.c \
	iproc/svector.h \
	iproc/svector.c \
	iproc/memory.h \
	iproc/utils.h \
	iproc/vars.h \
	iproc/vars.c \
	iproc/vector.h \
	iproc/vector.c

libiproc_la_LIBADD = \
	$(BLAS_LIBS) \
	$(LIBS) \
	$(FLIBS)

TESTS = \
	check-iproc

check_PROGRAMS = \
	check-iproc

check_iproc_SOURCES = \
	tests/actors.c \
	tests/array.c \
	tests/checkutils.h \
	tests/events.c \
	tests/history.c \
	tests/ieee754.c \
	tests/main.c \
	tests/svector.c \
	tests/vector.c

check_iproc_CPPFLAGS = \
	-Itests \
	$(CHECK_CFLAGS) \
	$(AM_CPPFLAGS)

check_iproc_LDADD = \
	$(top_builddir)/libiproc.la \
	$(CHECK_LIBS)
