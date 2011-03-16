ACLOCAL_AMFLAGS = -I m4 --install

lib_LTLIBRARIES = \
	libiproc.la

libiproc_la_SOURCES =

libiproc_la_LIBADD = \
	libcore.la

noinst_LTLIBRARIES = \
	libcore.la \
	libf2c.la \
	libminpack.la

libcore_la_SOURCES = \
	src/core/actors.c \
	src/core/actors.h \
	src/core/array.c \
	src/core/array.h \
	src/core/blas-private.h \
	src/core/design-ctx.c \
	src/core/design.c \
	src/core/design.h \
	src/core/events.c \
	src/core/events.h \
	src/core/fit.c \
	src/core/fit.h \
	src/core/hash.c \
	src/core/hash.h \
	src/core/history.c \
	src/core/history.h \
	src/core/ieee754.c \
	src/core/ieee754.h \
	src/core/loglik.c \
	src/core/loglik.h \
	src/core/logsumexp.c \
	src/core/logsumexp.h \
	src/core/matrix.c \
	src/core/matrix.h \
	src/core/memory.h \
	src/core/message-iter.c \
	src/core/messages.c \
	src/core/messages.h \
	src/core/model-ctx.c \
	src/core/model.c \
	src/core/model.h \
	src/core/refcount.c \
	src/core/refcount.h \
	src/core/sloglik.c \
	src/core/sloglik.h \
	src/core/svector.c \
	src/core/svector.h \
	src/core/utils.h \
	src/core/vector.c \
	src/core/vector.h \
	src/core/vrecip.c \
	src/core/vrecip.h

libcore_la_LIBADD = \
	libminpack.la \
	$(LAPACK_LIBS) \
	$(BLAS_LIBS) \
	$(LIBS) \
	$(FLIBS)

libf2c_la_SOURCES = \
	src/libf2c/f2c.h \
	src/libf2c/s_cmp.c \
	src/libf2c/s_copy.c

libminpack_la_SOURCES = \
	src/minpack/dcsrch.c \
	src/minpack/dcstep.c

libminpack_la_CPPFLAGS = \
	-I./src/libf2c \
	$(AM_CPPFLAGS)

libminpack_la_LIBADD = \
	$(BLAS_LIBS) \
	$(LIBS) \
	$(FLIBS)

dist_noinst_DATA = \
	@RPACKAGE@

rpkg_extra_sources = \
	$(libcore_la_SOURCES) \
	$(libf2c_la_SOURCES) \
	$(libminpack_la_SOURCES)

@RPACKAGE@: $(rpkg_extra_sources)
	cp $(rpkg_extra_sources)  interfaces/R/src
	R CMD build interfaces/R

TESTS_ENVIRONMENT = R CMD check --no-codoc --no-examples --no-manual --no-multiarch
TESTS = @RPACKAGE@

clean-local:
	rm -rf iproc.Rcheck
	rm -f @RPACKAGE@
	rm -f $(addprefix interfaces/R/src/,$(notdir $(rpkg_extra_sources)))
