ACLOCAL_AMFLAGS = -I m4 --install

AM_CPPFLAGS = \
	-Ilib/core/src \
	-Ilib/linalg/lib \
	-Ilib/linalg/src \
	-Ilib/include \
	-Ilib \
	-Isrc/core \
	-Isrc/json \
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
	libf2c.a \
	libiproc.a \
	libjson.a \
	libminpack.a \
	libyajl.a


libiproc_a_SOURCES = \
	src/design.c \
	src/design.h \
	src/design2.c \
	src/design2.h \
	src/dv-irecv.c \
	src/dv-isend.c \
	src/dv-nrecv.c \
	src/dv-nsend.c \
	src/frame.c \
	src/frame.h \
	src/history.c \
	src/history.h \
	src/messages.c \
	src/messages.h \
	src/recv_model.c \
	src/recv_model.h \
	src/strata.c \
	src/strata.h \
	src/vars.c \
	src/vars.h


libiproc_a_LIBADD = \
	libcore.a \
	libdSFMT.a
	$(LIBADD)


libcore_a_SOURCES = \
	lib/port.h \
	lib/catdist.c \
	lib/catdist.h \
	lib/util.c \
	lib/util.h \
	lib/core/src/coreutil.c \
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
	lib/linalg/src/matrixutil.c \
	lib/linalg/src/matrixutil.h \
	lib/linalg/src/sblas.c \
	lib/linalg/src/sblas.h \
	src/core/ieee754.c \
	src/core/ieee754.h \
	src/core/linesearch.c \
	src/core/linesearch.h \
	src/core/logsumexp.c \
	src/core/logsumexp.h

libcore_a_LIBADD = \
	libminpack.a \
	$(LIBADD)


libjson_a_SOURCES = \
	src/json/json.c \
	src/json/json.h

libjson_a_LIBADD = \
	libcore.a \
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

libdSFMT_a_CPPFLAGS = \
	-Wmissing-prototypes -Wall \
	-O3 -finline-functions -fomit-frame-pointer -DNDEBUG \
	-fno-strict-aliasing --param max-inline-insns-single=1800 \
	-std=c99 \
	-DDSFMT_MEXP=19937 \
	$(AM_CPPFLAGS)


libf2c_a_SOURCES = \
	lib/libf2c/f2c.h \
	lib/libf2c/s_cmp.c \
	lib/libf2c/s_copy.c

libminpack_a_SOURCES = \
	lib/minpack/dcsrch.c \
	lib/minpack/dcstep.c

libminpack_a_CPPFLAGS = \
	-I./lib/libf2c \
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

libyajl_a_CPPFLAGS = \
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
##	$(libiproc_a_SOURCES) \
##	$(libcore_a_SOURCES) \
##	$(libf2c_a_SOURCES) \
##	$(libminpack_a_SOURCES)

##@RPACKAGE@: $(rpkg_sources) $(rpkg_extra_sources)
##	$(AM_V_GEN)cp $(rpkg_extra_sources)  interfaces/R/src
##	R CMD build interfaces/R

check_LIBRARIES = \
	tests/libcmockery.a \
	tests/libenron.a

tests_libcmockery_a_SOURCES = \
	tests/cmockery/cmockery.c \
	tests/cmockery/cmockery.h \
	tests/testutil.c \
	tests/testutil.h

tests_libcmockery_a_LIBADD = \
	libcore.a

tests_libcmockery_a_CPPFLAGS = \
	-I./tests/cmockery \
	$(AM_CPPFLAGS)

tests_libenron_a_SOURCES = \
	tests/enron/employees.c \
	tests/enron/enron.h \
	tests/enron/messages.c

tests_libenron_a_LIBADD = \
	libiproc.a


check_LDADD = \
	$(noinst_LIBRARIES) \
	$(check_LIBRARIES) \
	$(LDADD)

check_CPPFLAGS = \
	-Itests/enron \
	-Itests/cmockery \
	$(AM_CPPFLAGS)

check_PROGRAMS = # empty

## check_PROGRAMS += tests/design_test
## tests_design_test_CPPFLAGS = $(check_CPPFLAGS)
## tests_design_test_LDADD    = $(check_LDADD)

## check_PROGRAMS += tests/frame_test
## tests_frame_test_CPPFLAGS = $(check_CPPFLAGS)
## tests_frame_test_LDADD    = $(check_LDADD)

check_PROGRAMS += tests/catdist_test
tests_catdist_test_CPPFLAGS = $(check_CPPFLAGS)
tests_catdist_test_LDADD    = $(check_LDADD)


## check_PROGRAMS += tests/mlogit_test
## tests_mlogit_test_CPPFLAGS = $(check_CPPFLAGS)
## tests_mlogit_test_LDADD    = $(check_LDADD)

## check_PROGRAMS += tests/mlogit1_test
## tests_mlogit1_test_CPPFLAGS = $(check_CPPFLAGS)
## tests_mlogit1_test_LDADD    = $(check_LDADD)


## check_PROGRAMS += tests/model_test
## tests_model_test_CPPFLAGS = $(check_CPPFLAGS)
## tests_model_test_LDADD    = $(check_LDADD)

#check_PROGRAMS += tests/recv_loglik_test
#tests_recv_loglik_test_CPPFLAGS = $(check_CPPFLAGS)
#tests_recv_loglik_test_LDADD    = $(check_LDADD)

#check_PROGRAMS += tests/enronfit
#tests_enronfit_SOURCES  = tests/enron/enronfit.c
#tests_enronfit_CPPFLAGS = $(check_CPPFLAGS)
#tests_enronfit_LDADD    = $(check_LDADD)

TESTS = \
	$(check_PROGRAMS)

##check-local:
##	R CMD check --no-codoc --no-examples --no-manual --no-multiarch @RPACKAGE@

##clean-local:
##	rm -rf iproc.Rcheck
##	rm -f @RPACKAGE@
##	rm -f $(addprefix interfaces/R/src/,$(notdir $(rpkg_extra_sources)))
