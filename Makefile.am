ACLOCAL_AMFLAGS = -I m4 --install

lib_LTLIBRARIES = \
	libiproc.la

noinst_LTLIBRARIES = \
	libcore.la \
	libf2c.la \
	libjson.la \
	libminpack.la \
	libyajl.la

libiproc_la_SOURCES = \
	src/actors.c \
	src/actors.h \
	src/design.c \
	src/design.h \
	src/event-trace.c \
	src/event-trace.h \
	src/frame.c \
	src/frame.h \
	src/history.c \
	src/history.h \
	src/messages.c \
	src/messages.h \
	src/recv_fit.c \
	src/recv_fit.h \
	src/recv_loglik.c \
	src/recv_loglik.h \
	src/recv_model.c \
	src/recv_model.h \
	src/rv-irecv.c \
	src/rv-nrecv.c \
	src/strata.c \
	src/strata.h \
	src/sv-irecv.c \
	src/sv-nrecv.c \
	src/vars.h

libiproc_la_CPPFLAGS = \
	-I./src/core \
	$(AM_CPPFLAGS)

libiproc_la_LIBADD = \
	libcore.la

libcore_la_SOURCES = \
	lib/timsort/timsort.h \
	lib/timsort/timsort-impl.h \
	lib/timsort/timsort.c \
	src/core/array.c \
	src/core/array.h \
	src/core/bfgs.c \
	src/core/bfgs.h \
	src/core/blas-private.h \
	src/core/compare.h \
	src/core/delegate.h \
	src/core/hash.c \
	src/core/hash.h \
	src/core/hashset.c \
	src/core/hashset.h \
	src/core/intmap.c \
	src/core/intmap.h \
	src/core/intset.c \
	src/core/intset.h \
	src/core/ieee754.c \
	src/core/ieee754.h \
	src/core/lapack-private.h \
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
	src/core/pqueue.c \
	src/core/pqueue.h \
	src/core/refcount.c \
	src/core/refcount.h \
	src/core/slice.c \
	src/core/slice.h \
	src/core/sparsetable.c \
	src/core/sparsetable.h \
	src/core/svector.c \
	src/core/svector.h \
	src/core/util.h \
	src/core/util.c \
	src/core/vector.c \
	src/core/vector.h

libcore_la_CPPFLAGS = \
	-I./src/core \
	-I./lib/timsort \
	$(AM_CPPFLAGS)

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
	lib/yajl/yajl_common.h \
	lib/yajl/yajl_gen.h \
	lib/yajl/yajl_parse.h \
	lib/yajl/yajl_tree.h \
	lib/yajl/yajl_version.h \
	lib/yajl/api/yajl_common.h \
	lib/yajl/api/yajl_gen.h \
	lib/yajl/api/yajl_parse.h \
	lib/yajl/api/yajl_tree.h \
	lib/yajl/yajl.c \
	lib/yajl/yajl_alloc.c \
	lib/yajl/yajl_alloc.h \
	lib/yajl/yajl_buf.c \
	lib/yajl/yajl_buf.h \
	lib/yajl/yajl_bytestack.h \
	lib/yajl/yajl_encode.c \
	lib/yajl/yajl_encode.h \
	lib/yajl/yajl_gen.c \
	lib/yajl/yajl_lex.c \
	lib/yajl/yajl_lex.h \
	lib/yajl/yajl_parser.c \
	lib/yajl/yajl_parser.h \
	lib/yajl/yajl_tree.c \
	lib/yajl/yajl_version.c

libyajl_la_CPPFLAGS = \
	-I./lib/ \
	-I./lib/yajl \
	$(AM_CPPFLAGS)


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
	interfaces/R/inst/unitTests/runit.actors.R \
	interfaces/R/R/actors.R \
	interfaces/R/R/design.R \
	interfaces/R/R/enron.R \
	interfaces/R/R/fit.R \
	interfaces/R/R/iproc.R \
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
	interfaces/R/src/r-design.c \
	interfaces/R/src/r-design.h \
	interfaces/R/src/r-fit.c \
	interfaces/R/src/r-fit.h \
	interfaces/R/src/r-init.c \
	interfaces/R/src/r-messages.c \
	interfaces/R/src/r-messages.h \
	interfaces/R/src/r-utils.c \
	interfaces/R/src/r-utils.h \
	interfaces/R/tests/doRUnit.R

rpkg_extra_sources = \
	$(libiproc_la_SOURCES) \
	$(libcore_la_SOURCES) \
	$(libf2c_la_SOURCES) \
	$(libminpack_la_SOURCES)

@RPACKAGE@: $(rpkg_sources) $(rpkg_extra_sources)
	$(AM_V_GEN)cp $(rpkg_extra_sources)  interfaces/R/src
	R CMD build interfaces/R

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

check_PROGRAMS += tests/hashset_test
tests_hashset_test_CPPFLAGS = -Isrc/core $(cmockery_cppflags)
tests_hashset_test_LDADD    = libcore.la $(cmockery_ldadd)

check_PROGRAMS += tests/pqueue_test
tests_pqueue_test_CPPFLAGS = -Isrc/core $(cmockery_cppflags)
tests_pqueue_test_LDADD    = libcore.la $(cmockery_ldadd)

check_PROGRAMS += tests/intmap_test
tests_intmap_test_CPPFLAGS = -Isrc/core $(cmockery_cppflags)
tests_intmap_test_LDADD    = libcore.la $(cmockery_ldadd)

check_PROGRAMS += tests/actors_test
tests_actors_test_CPPFLAGS = -Isrc -Isrc/core -Itests/enron $(cmockery_cppflags)
tests_actors_test_LDADD    = libiproc.la tests/libenron.la $(cmockery_ldadd)

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
tests_enronfit_CPPFLAGS = -Isrc -Isrc/core -Isrc/json -Itests/enron -Ilib
tests_enronfit_LDADD    = libiproc.la tests/libenron.la

TESTS = \
	$(check_PROGRAMS)

check-local:
	R CMD check --no-codoc --no-examples --no-manual --no-multiarch @RPACKAGE@

clean-local:
	rm -rf iproc.Rcheck
	rm -f @RPACKAGE@
	rm -f $(addprefix interfaces/R/src/,$(notdir $(rpkg_extra_sources)))
