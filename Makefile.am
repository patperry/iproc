ACLOCAL_AMFLAGS = -I m4 --install

AM_CPPFLAGS = \
	-Ilib/core/src \
	-Ilib/linalg/lib \
	-Ilib/linalg/src \
	-Ilib/include \
	-Ilib \
	-Isrc

LIBADD = \
	$(LIBOBJS) \
	$(ALLOCA)

LDADD = \
	$(LAPACK_LIBS) \
	$(BLAS_LIBS) \
	$(LIBADD)


noinst_LIBRARIES = \
	libcore.a \
	libdSFMT.a \
	libiproc.a \
	libyajl.a


libiproc_a_SOURCES = \
	src/design.c \
	src/design_alg.c \
	src/design_prod.c \
	src/design.h \
	src/design2.c \
	src/design2_alg.c \
	src/design2_prod.c \
	src/design2.h \
	src/history.c \
	src/history.h \
	src/recv_fit.c \
	src/recv_fit.h \
	src/recv_loglik.c \
	src/recv_loglik.h \
	src/recv_model.c \
	src/recv_model.h \
	src/recv_resid.c \
	src/recv_resid.h \
	src/send_model.c \
	src/send_model.h \
	src/v-irecvtot.c \
	src/v-isendtot.c \
	src/v-itot-source.h \
	src/v-nrecvtot.c \
	src/v-nsendtot.c \
	src/v-ntot-source.h \
	src/v2-idyad-source.h \
	src/v2-irecv.c \
	src/v2-isend.c \
	src/v2-ndyad-source.h \
	src/v2-nrecv.c \
	src/v2-nsend.c \
	src/v2-ntriad-source.h \
	src/v2-ncosib.c \
	src/v2-nsend2.c \
	src/v2-nsib.c \
	src/v2-nrecv2.c \
	src/var.c \
	src/var.h


libiproc_a_LIBADD = \
	libcore.a \
	libdSFMT.a \
	$(LIBADD)


libcore_a_SOURCES = \
	lib/port.h \
	lib/catdist.c \
	lib/catdist.h \
	lib/catdist1.c \
	lib/catdist1.h \
	lib/constr.c \
	lib/constr.h \
	lib/deltaset.c \
	lib/deltaset.h \
	lib/linesearch.h \
	lib/logsumexp.c \
	lib/logsumexp.h \
	lib/mlogit.c \
	lib/mlogit.h \
	lib/mlogitaug.c \
	lib/mlogitaug.h \
	lib/newton.c \
	lib/newton.h \
	lib/strata.c \
	lib/strata.h \
	lib/util.c \
	lib/util.h \
	lib/version.c \
	lib/version.h \
	lib/core/src/coreutil.c \
	lib/core/src/coreutil.h \
	lib/core/src/hash.h \
	lib/core/src/hashset.c \
	lib/core/src/hashset.h \
	lib/core/src/intset.c \
	lib/core/src/intset.h \
	lib/core/src/ieee754.c \
	lib/core/src/ieee754.h \
	lib/core/src/pqueue.c \
	lib/core/src/pqueue.h \
	lib/core/src/timsort-impl.h \
	lib/core/src/timsort.c \
	lib/core/src/timsort.h \
	lib/core/src/uintset.c \
	lib/core/src/uintset.h \
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
	lib/linalg/src/matrixutil.c \
	lib/linalg/src/matrixutil.h \
	lib/linalg/src/sblas.c \
	lib/linalg/src/sblas.h

libcore_a_LIBADD = \
	$(LIBADD)


libdSFMT_a_SOURCES = \
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

libdSFMT_a_CFLAGS = \
	-Wmissing-prototypes -Wall \
	-O3 -finline-functions -fomit-frame-pointer \
	-fno-strict-aliasing --param max-inline-insns-single=1800 \
	-std=c99 \
	$(AM_CFLAGS)

libdSFMT_a_CPPFLAGS = \
	-DNDEBUG \
	-DDSFMT_MEXP=19937 \
	$(AM_CPPFLAGS)


libyajl_a_SOURCES = \
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

libyajl_a_CFLAGS = \
	-Wno-sign-compare \
	$(AM_CFLAGS)

libyajl_a_CPPFLAGS = \
	-I./lib/yajl/src \
	-I./lib/yajl/src/api \
	$(AM_CPPFLAGS)


dist_noinst_DATA = \
	@RPACKAGE@

rpkg_sources = \
	R/DESCRIPTION \
	R/LICENSE \
	R/NAMESPACE \
	R/data/enron.rda \
	R/demo/00Index \
	R/demo/enron.R \
	R/inst/tests/test-Mesg.R \
	R/R/Mesg.R \
	R/R/fit.R \
	R/R/recv.frame.R \
	R/R/recv.matrix.R \
	R/R/recv.model.R \
	R/R/summary.recv.model.R \
	R/R/utils.R \
	R/src/Makevars \
	R/src/Makevars.win \
	R/src/r-recv_model.c \
	R/src/r-xalloc-die.c \
	R/tests/test-all.R

rpkg_extra_sources = \
	$(libiproc_a_SOURCES) \
	$(libcore_a_SOURCES)

@RPACKAGE@: $(rpkg_sources) $(rpkg_extra_sources)
	$(AM_V_GEN)cp $(rpkg_extra_sources)  R/src
	rm -f R/src/xerbla.c R/src/xalloc-die.c
	R CMD build R

check_LIBRARIES = \
	tests/libcmockery.a \
	tests/libenron.a \
	tests/libfixtures.a

tests_libcmockery_a_SOURCES = \
	tests/cmockery/cmockery.c \
	tests/cmockery/cmockery.h \
	tests/testutil.c \
	tests/testutil.h

tests_libcmockery_a_LIBADD = \
	libcore.a

tests_libcmockery_a_CPPFLAGS = \
	-Itests/cmockery \
	$(AM_CPPFLAGS)

tests_libenron_a_SOURCES = \
	tests/enron/actors.c \
	tests/enron/actors.h \
	tests/enron/messages.c \
	tests/enron/messages.h \
	tests/enron/paths.h

tests_libenron_a_CPPFLAGS = \
	-Itests \
	$(AM_CPPFLAGS)

tests_libfixtures_a_SOURCES = \
	tests/fixtures/design.c \
	tests/fixtures/design.h \
	tests/fixtures/design2.c \
	tests/fixtures/design2.h \
	tests/fixtures/history.c \
	tests/fixtures/history.h \
	tests/fixtures/recv_model.c \
	tests/fixtures/recv_model.h \
	tests/fixtures/send_model.c \
	tests/fixtures/send_model.h

tests_libfixtures_a_LIBADD = \
	tests/libenron.a \
	libcore.a \
	libiproc.a

tests_libfixtures_a_CPPFLAGS = \
	-Itests \
	$(AM_CPPFLAGS)


check_LDADD = \
	$(noinst_LIBRARIES) \
	$(check_LIBRARIES) \
	$(LDADD)

check_CPPFLAGS = \
	-Isrc \
	-Itests \
	-Itests/cmockery \
	$(AM_CPPFLAGS)

check_PROGRAMS = # empty

check_PROGRAMS += tests/catdist_test
tests_catdist_test_CPPFLAGS = $(check_CPPFLAGS)
tests_catdist_test_LDADD    = $(check_LDADD)

check_PROGRAMS += tests/catdist1_test
tests_catdist1_test_CPPFLAGS = $(check_CPPFLAGS)
tests_catdist1_test_LDADD    = $(check_LDADD)

check_PROGRAMS += tests/mlogit_test
tests_mlogit_test_CPPFLAGS = $(check_CPPFLAGS)
tests_mlogit_test_LDADD    = $(check_LDADD)

check_PROGRAMS += tests/history_test
tests_history_test_CPPFLAGS = $(check_CPPFLAGS)
tests_history_test_LDADD    = $(check_LDADD)

check_PROGRAMS += tests/design_test
tests_design_test_CPPFLAGS = $(check_CPPFLAGS)
tests_design_test_LDADD    = $(check_LDADD)

check_PROGRAMS += tests/send_model_test
tests_send_model_test_CPPFLAGS = $(check_CPPFLAGS)
tests_send_model_test_LDADD    = $(check_LDADD)

check_PROGRAMS += tests/design2_test
tests_design2_test_CPPFLAGS = $(check_CPPFLAGS)
tests_design2_test_LDADD    = $(check_LDADD)

check_PROGRAMS += tests/recv_model_test
tests_recv_model_test_CPPFLAGS = $(check_CPPFLAGS)
tests_recv_model_test_LDADD    = $(check_LDADD)

check_PROGRAMS += tests/recv_loglik_test
tests_recv_loglik_test_CPPFLAGS = $(check_CPPFLAGS)
tests_recv_loglik_test_LDADD    = $(check_LDADD)

check_PROGRAMS += tests/enronfit
tests_enronfit_SOURCES  = tests/enronfit.c
tests_enronfit_CPPFLAGS = $(check_CPPFLAGS) $(HDF5_CPPFLAGS)
tests_enronfit_LDADD    = $(check_LDADD) $(HDF5_LDFLAGS)

TESTS = \
	$(check_PROGRAMS)

##check-local:
##	R CMD check --no-codoc --no-examples --no-manual --no-multiarch @RPACKAGE@

clean-local:
	rm -rf iproc.Rcheck
	rm -f @RPACKAGE@
	rm -f $(addprefix R/src/,$(notdir $(rpkg_extra_sources)))
	rm -f R/src/*.o R/src/*.so
