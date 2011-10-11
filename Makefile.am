ACLOCAL_AMFLAGS = -I m4 --install

AM_CPPFLAGS = \
	-Ilib/core/src \
	-Ilib/linalg/lib \
	-Ilib/linalg/src \
	-Ilib/include \
	-Ilib \
	-Isrc/core \
	-Isrc

lib_LTLIBRARIES = \
	libiproc.la

noinst_LTLIBRARIES = \
	libcore.la \
	libdSFMT.la \
	libf2c.la \
	libjson.la \
	libminpack.la \
	libyajl.la

libiproc_la_SOURCES = \
	src/design.c \
	src/design.h \
	src/frame.c \
	src/frame.h \
	src/messages.c \
	src/messages.h \
	src/recv_boot.c \
	src/recv_boot.h \
	src/recv_fit.c \
	src/recv_fit.h \
	src/recv_loglik.c \
	src/recv_loglik.h \
	src/recv_model.c \
	src/recv_model.h \
	src/recv_resid.c \
	src/recv_resid.h \
	src/rv-icosib.c \
	src/rv-irecv.c \
	src/rv-irecv2.c \
	src/rv-isend.c \
	src/rv-isend2.c \
	src/rv-isib.c \
	src/rv-ncosib.c \
	src/rv-nrecv.c \
	src/rv-nrecv2.c \
	src/rv-nsend.c \
	src/rv-nsend2.c \
	src/rv-nsib.c \
	src/strata.c \
	src/strata.h \
	src/vars.c \
	src/vars.h

libiproc_la_CPPFLAGS = \
	-I./src/core \
	-I./lib \
	$(AM_CPPFLAGS)

libiproc_la_LIBADD = \
	libcore.la \
	libdSFMT.la

libcore_la_SOURCES = \
	lib/core/src/coreutil.h \
	lib/core/src/hash.h \
	lib/core/src/hashset.c \
	lib/core/src/hashset.h \
	lib/core/src/pqueue.c \
	lib/core/src/pqueue.h \
	lib/core/src/timsort-impl.h \
	lib/core/src/timsort.c \
	lib/core/src/timsort.h \
	lib/core/src/xalloc-die.c \
	lib/core/src/xalloc.c \
	lib/core/src/xalloc.h \
	lib/linalg/lib/f77.c \
	lib/linalg/lib/f77.h \
	lib/linalg/lib/f77blas.h \
	lib/linalg/lib/f77lapack.h \
	lib/linalg/lib/xerbla.c \
	lib/linalg/src/blas.c \
	lib/linalg/src/blas.h \
	lib/linalg/src/lapack.c \
	lib/linalg/src/lapack.h \
	lib/linalg/src/sblas.c \
	lib/linalg/src/sblas.h \
	src/core/array.c \
	src/core/array.h \
	src/core/bfgs.c \
	src/core/bfgs.h \
	src/core/compare.h \
	src/core/delegate.h \
	src/core/intmap.c \
	src/core/intmap.h \
	src/core/intset.c \
	src/core/intset.h \
	src/core/ieee754.c \
	src/core/ieee754.h \
	src/core/linalg.c \
	src/core/linalg.h \
	src/core/linesearch.c \
	src/core/linesearch.h \
	src/core/logsumexp.c \
	src/core/logsumexp.h \
	src/core/macconfig.h \
	src/core/matrix.c \
	src/core/matrix.h \
	src/core/port.h \
	src/core/refcount.c \
	src/core/refcount.h \
	src/core/slice.c \
	src/core/slice.h \
	src/core/svector.c \
	src/core/svector.h \
	src/core/util.h \
	src/core/util.c \
	src/core/vector.c \
	src/core/vector.h

libcore_la_LIBADD = \
	libminpack.la \
	$(LAPACK_LIBS) \
	$(BLAS_LIBS) \
	$(LIBS) \
	$(FLIBS) \
	$(LTLIBOBJS)

libjson_la_SOURCES = \
	src/json/json.c \
	src/json/json.h

libjson_la_LIBADD = \
	libcore.la \
	libyajl.la \
	$(LIBS) \
	$(FLIBS) \
	$(LTLIBOBJS)

libjson_la_CPPFLAGS = \
	$(libcore_la_CPPFLAGS) \
	$(libyajl_la_CPPFLAGS) \
	$(AM_CPPFLAGS)

libdSFMT_la_SOURCES = \
	lib/dSFMT/dSFMT-params.h \
	lib/dSFMT/dSFMT-params11213.h \
	lib/dSFMT/dSFMT-params1279.h \
	lib/dSFMT/dSFMT-params132049.h \
	lib/dSFMT/dSFMT-params19937.h \
	lib/dSFMT/dSFMT-params216091.h \
	lib/dSFMT/dSFMT-params2203.h \
	lib/dSFMT/dSFMT-params4253.h \
	lib/dSFMT/dSFMT-params44497.h \
	lib/dSFMT/dSFMT-params521.h \
	lib/dSFMT/dSFMT-params86243.h \
	lib/dSFMT/dSFMT.c \
	lib/dSFMT/dSFMT.h

libdSFMT_la_CPPFLAGS = \
	-Wmissing-prototypes -Wall \
	-O3 -finline-functions -fomit-frame-pointer -DNDEBUG \
	-fno-strict-aliasing --param max-inline-insns-single=1800 \
	-std=c99 \
	-DDSFMT_MEXP=19937 \
	$(AM_CPPFLAGS)

libdSFMT_la_LIBADD = \
	$(LIBS) \
	$(LTLIBOBJS)

libf2c_la_SOURCES = \
	lib/libf2c/f2c.h \
	lib/libf2c/s_cmp.c \
	lib/libf2c/s_copy.c

libminpack_la_SOURCES = \
	lib/minpack/dcsrch.c \
	lib/minpack/dcstep.c

libminpack_la_CPPFLAGS = \
	-I./lib/libf2c \
	$(AM_CPPFLAGS)

libminpack_la_LIBADD = \
	libf2c.la \
	$(BLAS_LIBS) \
	$(LIBS) \
	$(FLIBS)

libyajl_la_SOURCES = \
	lib/yajl/src/api/yajl_common.h \
	lib/yajl/src/api/yajl_gen.h \
	lib/yajl/src/api/yajl_parse.h \
	lib/yajl/src/api/yajl_tree.h \
	lib/yajl/src/yajl.c \
	lib/yajl/src/yajl_alloc.c \
	lib/yajl/src/yajl_alloc.h \
	lib/yajl/src/yajl_buf.c \
	lib/yajl/src/yajl_buf.h \
	lib/yajl/src/yajl_bytestack.h \
	lib/yajl/src/yajl_encode.c \
	lib/yajl/src/yajl_encode.h \
	lib/yajl/src/yajl_gen.c \
	lib/yajl/src/yajl_lex.c \
	lib/yajl/src/yajl_lex.h \
	lib/yajl/src/yajl_parser.c \
	lib/yajl/src/yajl_parser.h \
	lib/yajl/src/yajl_tree.c

libyajl_la_CPPFLAGS = \
	-I./lib/yajl/src \
	-I./lib/yajl/src/api \
	$(AM_CPPFLAGS)


## dist_noinst_DATA = \
##	@RPACKAGE@

## rpkg_sources = \
##	interfaces/R/DESCRIPTION \
##	interfaces/R/LICENSE \
##	interfaces/R/NAMESPACE \
##	interfaces/R/data/enron.rda \
##	interfaces/R/inst/extdata/enron/employees.tsv \
##	interfaces/R/inst/extdata/enron/messages.tsv \
##	interfaces/R/inst/extdata/enron/recipients.tsv \
##	interfaces/R/inst/unitTests/runit.actors.R \
##	interfaces/R/R/actors.R \
##	interfaces/R/R/design.R \
##	interfaces/R/R/enron.R \
##	interfaces/R/R/fit.R \
##	interfaces/R/R/iproc.R \
##	interfaces/R/R/loglik.R \
##	interfaces/R/R/messages.R \
##	interfaces/R/R/methods.R \
##	interfaces/R/R/model.R \
##	interfaces/R/R/runit.R \
##	interfaces/R/R/utils.R \
##	interfaces/R/src/Makevars \
##	interfaces/R/src/Makevars.win \
##	interfaces/R/src/config.h \
##	interfaces/R/src/r-actors.c \
##	interfaces/R/src/r-actors.h \
##	interfaces/R/src/r-design.c \
##	interfaces/R/src/r-design.h \
##	interfaces/R/src/r-fit.c \
##	interfaces/R/src/r-fit.h \
##	interfaces/R/src/r-init.c \
##	interfaces/R/src/r-messages.c \
##	interfaces/R/src/r-messages.h \
##	interfaces/R/src/r-utils.c \
##	interfaces/R/src/r-utils.h \
##	interfaces/R/tests/doRUnit.R

##rpkg_extra_sources = \
##	$(libiproc_la_SOURCES) \
##	$(libcore_la_SOURCES) \
##	$(libf2c_la_SOURCES) \
##	$(libminpack_la_SOURCES)

##@RPACKAGE@: $(rpkg_sources) $(rpkg_extra_sources)
##	$(AM_V_GEN)cp $(rpkg_extra_sources)  interfaces/R/src
##	R CMD build interfaces/R

check_LTLIBRARIES = \
	tests/libcmockery.la \
	tests/libenron.la

tests_libcmockery_la_SOURCES = \
	tests/cmockery/cmockery.c \
	tests/cmockery/cmockery.h

tests_libcmockery_la_CPPFLAGS = \
	-I./tests/cmockery \
	$(AM_CPPFLAGS)

tests_libenron_la_SOURCES = \
	tests/enron/employees.c \
	tests/enron/enron.h \
	tests/enron/messages.c

tests_libenron_la_CPPFLAGS = \
	-I./lib \
	-I./src \
	-I./src/core \
	$(AM_CPPFLAGS)

tests_libenron_la_LIBADD = \
	libiproc.la \
	libjson.la \
	libyajl.la


cmockery_cppflags = -Itests/cmockery $(AM_CPPFLAGS)
cmockery_ldadd    = tests/libcmockery.la

check_PROGRAMS = # empty

check_PROGRAMS += tests/array_test
tests_array_test_CPPFLAGS = -Isrc/core $(cmockery_cppflags)
tests_array_test_LDADD    = libcore.la $(cmockery_ldadd)

check_PROGRAMS += tests/intmap_test
tests_intmap_test_CPPFLAGS = -Isrc/core $(cmockery_cppflags)
tests_intmap_test_LDADD    = libcore.la $(cmockery_ldadd)

check_PROGRAMS += tests/design_test
tests_design_test_CPPFLAGS = -Isrc -Isrc/core -Itests/enron $(cmockery_cppflags)
tests_design_test_LDADD    = libiproc.la tests/libenron.la $(cmockery_ldadd)

check_PROGRAMS += tests/frame_test
tests_frame_test_CPPFLAGS = -Isrc -Isrc/core -Itests/enron $(cmockery_cppflags)
tests_frame_test_LDADD    = libiproc.la tests/libenron.la $(cmockery_ldadd)

check_PROGRAMS += tests/model_test
tests_model_test_CPPFLAGS = -Isrc -Isrc/core -Itests/enron $(cmockery_cppflags)
tests_model_test_LDADD    = libiproc.la tests/libenron.la $(cmockery_ldadd)

check_PROGRAMS += tests/recv_loglik_test
tests_recv_loglik_test_CPPFLAGS = -Isrc -Isrc/core -Itests/enron $(cmockery_cppflags)
tests_recv_loglik_test_LDADD    = libiproc.la tests/libenron.la $(cmockery_ldadd)

check_PROGRAMS += tests/enronfit
tests_enronfit_SOURCES  = tests/enron/enronfit.c
tests_enronfit_CPPFLAGS = -Isrc/core -Isrc/json -Itests/enron $(AM_CPPFLAGS)
tests_enronfit_LDADD    = libiproc.la tests/libenron.la

TESTS = \
	$(check_PROGRAMS)

##check-local:
##	R CMD check --no-codoc --no-examples --no-manual --no-multiarch @RPACKAGE@

##clean-local:
##	rm -rf iproc.Rcheck
##	rm -f @RPACKAGE@
##	rm -f $(addprefix interfaces/R/src/,$(notdir $(rpkg_extra_sources)))
