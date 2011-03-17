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

rpkg_sources = \
	interfaces/R/DESCRIPTION \
	interfaces/R/LICENSE \
	interfaces/R/NAMESPACE \
	interfaces/R/data/enron.rda \
	interfaces/R/inst/extdata/enron/employees.tsv \
	interfaces/R/inst/extdata/enron/messages.tsv \
	interfaces/R/inst/extdata/enron/recipients.tsv \
	interfaces/R/inst/unitTests/common.loglik.R \
	interfaces/R/inst/unitTests/common.loglik.all.R \
	interfaces/R/inst/unitTests/common.loglik.groups.R \
	interfaces/R/inst/unitTests/common.loglik.recip.R \
	interfaces/R/inst/unitTests/common.model.R \
	interfaces/R/inst/unitTests/common.model.all.R \
	interfaces/R/inst/unitTests/common.model.groups.R \
	interfaces/R/inst/unitTests/common.model.recip.R \
	interfaces/R/inst/unitTests/runit.actors.R \
	interfaces/R/inst/unitTests/runit.cursor.R \
	interfaces/R/inst/unitTests/runit.design.all.R \
	interfaces/R/inst/unitTests/runit.design.groups.R \
	interfaces/R/inst/unitTests/runit.design.recip.R \
	interfaces/R/inst/unitTests/runit.enron.R \
	interfaces/R/inst/unitTests/runit.loglik.all.loops.R \
	interfaces/R/inst/unitTests/runit.loglik.all.noloops.R \
	interfaces/R/inst/unitTests/runit.loglik.groups.loops.R \
	interfaces/R/inst/unitTests/runit.loglik.groups.noloops.R \
	interfaces/R/inst/unitTests/runit.loglik.recip.loops.R \
	interfaces/R/inst/unitTests/runit.loglik.recip.noloops.R \
	interfaces/R/inst/unitTests/runit.messages.R \
	interfaces/R/inst/unitTests/runit.model.all.loops.R \
	interfaces/R/inst/unitTests/runit.model.all.noloops.R \
	interfaces/R/inst/unitTests/runit.model.groups.loops.R \
	interfaces/R/inst/unitTests/runit.model.groups.noloops.R \
	interfaces/R/inst/unitTests/runit.model.recip.loops.R \
	interfaces/R/inst/unitTests/runit.model.recip.noloops.R \
	interfaces/R/R/actors.R \
	interfaces/R/R/cursor.R \
	interfaces/R/R/design.R \
	interfaces/R/R/enron.R \
	interfaces/R/R/fit.R \
	interfaces/R/R/loglik.R \
	interfaces/R/R/messages.R \
	interfaces/R/R/methods.R \
	interfaces/R/R/model.R \
	interfaces/R/R/runit.R \
	interfaces/R/R/utils.R \
	interfaces/R/src/Makevars \
	interfaces/R/src/Makevars.win \
	interfaces/R/src/config.h \
	interfaces/R/src/r-actors.c \
	interfaces/R/src/r-actors.h \
	interfaces/R/src/r-cursor.c \
	interfaces/R/src/r-cursor.h \
	interfaces/R/src/r-design-udata.c \
	interfaces/R/src/r-design.c \
	interfaces/R/src/r-design.h \
	interfaces/R/src/r-fit.c \
	interfaces/R/src/r-fit.h \
	interfaces/R/src/r-init.c \
	interfaces/R/src/r-loglik.c \
	interfaces/R/src/r-loglik.h \
	interfaces/R/src/r-messages.c \
	interfaces/R/src/r-messages.h \
	interfaces/R/src/r-model.c \
	interfaces/R/src/r-model.h \
	interfaces/R/src/r-utils.c \
	interfaces/R/src/r-utils.h \
	interfaces/R/tests/doRUnit.R

rpkg_extra_sources = \
	$(libcore_la_SOURCES) \
	$(libf2c_la_SOURCES) \
	$(libminpack_la_SOURCES)

@RPACKAGE@: $(rpkg_sources) $(rpkg_extra_sources)
	cp $(rpkg_extra_sources)  interfaces/R/src
	R CMD build interfaces/R

TESTS_ENVIRONMENT = R CMD check --no-codoc --no-examples --no-manual --no-multiarch
TESTS = @RPACKAGE@

clean-local:
	rm -rf iproc.Rcheck
	rm -f @RPACKAGE@
	rm -f $(addprefix interfaces/R/src/,$(notdir $(rpkg_extra_sources)))
