#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>

#include <design.h>


#include "coreutil.h"
#include "testutil.h"
#include "cmockery.h"
#include "lapack.h"
#include "matrixutil.h"
#include "sblas.h"
#include "xalloc.h"

#include "fixtures/history.h"
#include "fixtures/design.h"


struct fixture {
	struct history_fixture h;
	struct design_fixture  d;
};

struct fixture ctx;

#define NRECV	ctx.h.nrecv
#define NSEND	ctx.h.nsend



struct test {
	struct history h;
	struct design d;
};

struct test test;

#define DESIGN	&test.d
#define HISTORY	&test.h



static double vector_dist(size_t n, double *x, double *y)
{
	size_t i;
	double dist = 0;

	for (i = 0; i < n; i++) {
		double d = x[i] - y[i];
		dist += d * d;
	}

	return sqrt(dist);
}

static void fixture_setup_enron()
{
	print_message("Enron employees\n");
	print_message("---------------\n");
	history_fixture_setup_enron(&ctx.h);
	design_fixture_setup_enron(&ctx.d);
	design_fixture_add_traits_enron(&ctx.d);
}

static void fixture_teardown()
{
	design_fixture_teardown(&ctx.d);
	history_fixture_teardown(&ctx.h);
	print_message("\n\n");
}


static void setup()
{
	history_test_setup(&test.h, &ctx.h);
	design_test_setup(&test.d, &test.h, &ctx.d);
}

static void teardown()
{
	design_test_teardown(&test.d);
	history_test_teardown(&test.h);
}


static void matrix_assign_reffects(double *x, size_t tdx,
				   size_t nrecv,
				   int has_reffects)
{
	if (!has_reffects)
		return;

	size_t i;

	matrix_dzero(nrecv, nrecv, x, tdx);
	for (i = 0; i < nrecv; i++) {
		x[i * tdx + i] = 1.0;
	}
}

static void matrix_init_design0(double **px, const struct design *d)
{
	const double *traits = design_trait_matrix(d);
	//int has_reffects = design_has_effects(d);
	int has_reffects = 0;
	
	size_t nrecv = design_count(d);
	size_t pr = design_trait_dim(d);
	size_t ireff = 0;
	size_t nreff = has_reffects ? nrecv : 0;
	size_t istat = ireff + nreff;
	size_t nstat = pr;
	size_t dim =  istat + nstat;
	double *x = xmalloc(nrecv * dim * sizeof(double));
	size_t tdx = MAX(1, dim);
	
	*px = x;
	double *xreff = x + ireff;
	matrix_assign_reffects(xreff, tdx, nrecv, has_reffects);
	double *xstat = x + istat;
	lapack_dlacpy(LA_COPY_ALL, nstat, nrecv, traits, MAX(1, nstat), xstat, tdx);
}

static void test_mul0()
{

	double *a;
	size_t tda;
	
	double *x, *y, *y1;	
	size_t i, n, p;
	
	n = design_count(DESIGN);
	p = design_trait_dim(DESIGN);
	
	x = xmalloc(p * sizeof(double));
	y = xmalloc(n * sizeof(double));
	y1 = xmalloc(n * sizeof(double));
	
	matrix_init_design0(&a, DESIGN);
	tda = MAX(1, p);
	
	for (i = 0; i < p; i++) {
		x[i] =  (7 * i) % 5 - 2.0;
	}
		
	design_traits_mul(1.0, DESIGN, x, 0.0, y);
	blas_dgemv(BLAS_TRANS, p, n, 1.0, a, tda, x, 1, 0.0, y1, 1);
	assert_true(vector_dist(n, y, y1) == 0.0);		
		
	design_traits_mul(1.0, DESIGN, x, 1.0, y);
	blas_dgemv(BLAS_TRANS, p, n, 1.0, a, tda, x, 1, 1.0, y1, 1);
	assert_true(vector_dist(n, y, y1) == 0.0);
		
	design_traits_mul(1.0, DESIGN, x, -1.0, y);
	blas_dgemv(BLAS_TRANS, p, n, 1.0, a, tda, x, 1, -1.0, y1, 1);	
	assert_true(vector_dist(n, y, y1) == 0.0);		
		
	design_traits_mul(2.0, DESIGN, x, 2.0, y);
	blas_dgemv(BLAS_TRANS, p, n, 2.0, a, tda, x, 1, 2.0, y1, 1);	
	assert_true(vector_dist(n, y, y1) == 0.0);		
		
	free(a);
	free(y1);
	free(y);
	free(x);
}


static void test_tmul0()
{
	
	double *a;
	size_t tda;
	double *x, *y, *y1;	
	size_t i, n, p;
	
	n = design_count(DESIGN);
	p = design_trait_dim(DESIGN);
	
	x = xmalloc(n * sizeof(double));
	y = xmalloc(p * sizeof(double));
	y1 = xmalloc(p * sizeof(double));
	matrix_init_design0(&a, DESIGN);
	tda = MAX(1, p);
		
	for (i = 0; i < n; i++) {
		x[i] = (3 * i + 1) % 5 - 2.0;
	}
		
	design_traits_tmul(1.0, DESIGN, x, 0.0, y);
	blas_dgemv(BLAS_NOTRANS, p, n, 1.0, a, tda, x, 1, 0.0, y1, 1);
	assert_true(vector_dist(p, y, y1) == 0.0);		
		
	design_traits_tmul(1.0, DESIGN, x, 1.0, y);
	blas_dgemv(BLAS_NOTRANS, p, n, 1.0, a, tda, x, 1, 1.0, y1, 1);
	assert_true(vector_dist(p, y, y1) == 0.0);
		
	design_traits_tmul(1.0, DESIGN, x, -1.0, y);
	blas_dgemv(BLAS_NOTRANS, p, n, 1.0, a, tda, x, 1, -1.0, y1, 1);	
	assert_true(vector_dist(p, y, y1) == 0.0);
		
	design_traits_tmul(2.0, DESIGN, x, 2.0, y);
	blas_dgemv(BLAS_NOTRANS, p, n, 2.0, a, tda, x, 1, 2.0, y1, 1);
	assert_true(vector_dist(p, y, y1) == 0.0);
		
	free(a);
	free(y1);
	free(y);
	free(x);
}

/*
static void test_tmuls0()
{
	struct dmatrix matrix;
	double *x, *y, *y1;
	struct vpattern pat;
	size_t i, n, p;

	n = design_count(design);
	p = design_dim(design);
	x = xmalloc(n * sizeof(x[0]));
	pat.indx = xmalloc(n * sizeof(pat.indx[0]));

	y = xmalloc(p * sizeof(double));
	y1 = xmalloc(p * sizeof(double));
	matrix_init_design0(&matrix, design);

	pat.nz = 0;
	for (i = 0; i < n; i++) {
		if (i % 3 == 0) {
			x[pat.nz] = (3 * i + 1) % 5 - 2.0;
			pat.indx[pat.nz] = i;
			pat.nz++;
		}
	}

	design_muls0(1.0, BLAS_TRANS, design, x, &pat, 0.0, y);
	sblas_dgemvi(BLAS_TRANS, n, p, 1.0, &matrix, x, &pat, 0.0, y1);
	assert_true(vector_dist(p, y, y1) == 0.0);

	design_muls0(1.0, BLAS_TRANS, design, x, &pat, 1.0, y);
	sblas_dgemvi(BLAS_TRANS, n, p, 1.0, &matrix, x, &pat, 1.0, y1);
	assert_true(vector_dist(p, y, y1) == 0.0);

	design_muls0(1.0, BLAS_TRANS, design, x, &pat, -1.0, y);
	sblas_dgemvi(BLAS_TRANS, n, p, 1.0, &matrix, x, &pat, -1.0, y1);
	assert_true(vector_dist(p, y, y1) == 0.0);

	design_muls0(2.0, BLAS_TRANS, design, x, &pat, 2.0, y);
	sblas_dgemvi(BLAS_TRANS, n, p, 2.0, &matrix, x, &pat, 2.0, y1);
	assert_true(vector_dist(p, y, y1) == 0.0);

	free(matrix.data);
	free(y1);
	free(y);
	free(pat.indx);
	free(x);
}


static void test_muls0()
{

	struct dmatrix matrix;
	double *x, *y, *y1;
	struct vpattern pat;
	size_t i, n, p;

	n = design_count(design);
	p = design_dim(design);
	x = xmalloc(n * sizeof(x[0]));
	pat.indx = xmalloc(n * sizeof(pat.indx[0]));

	y = xmalloc(n * sizeof(double));
	y1 = xmalloc(n * sizeof(double));
	matrix_init_design0(&matrix, design);

	pat.nz = 0;
	for (i = 0; i < p; i++) {
		if (i % 2 == 0) {
			x[pat.nz] = (7 * i) % 5 - 2.0;
			pat.indx[pat.nz] = i;
			pat.nz++;
		}
	}

	design_muls0(1.0, BLAS_NOTRANS, design, x, &pat, 0.0, y);
	sblas_dgemvi(BLAS_NOTRANS, n, p, 1.0, &matrix, x, &pat, 0.0, y1);
	assert_true(vector_dist(n, y, y1) == 0.0);

	design_muls0(1.0, BLAS_NOTRANS, design, x, &pat, 1.0, y);
	sblas_dgemvi(BLAS_NOTRANS, n, p, 1.0, &matrix, x, &pat, 1.0, y1);
	assert_true(vector_dist(n, y, y1) == 0.0);

	design_muls0(1.0, BLAS_NOTRANS, design, x, &pat, -1.0, y);
	sblas_dgemvi(BLAS_NOTRANS, n, p, 1.0, &matrix, x, &pat, -1.0, y1);
	assert_true(vector_dist(n, y, y1) == 0.0);

	design_muls0(2.0, BLAS_NOTRANS, design, x, &pat, 2.0, y);
	sblas_dgemvi(BLAS_NOTRANS, n, p, 2.0, &matrix, x, &pat, 2.0, y1);
	assert_true(vector_dist(n, y, y1) == 0.0);

	free(matrix.data);
	free(y1);
	free(y);
	free(pat.indx);
	free(x);
}
*/

static void test_irecvtot()
{
	double window = 2000;
	const struct var *v = design_add_tvar(DESIGN, "IRecvTot", VAR_IRECVTOT, window);

	history_advance(HISTORY, 988116420);

	size_t imsgs, nmsgs = history_count(HISTORY);
	for (imsgs = nmsgs; imsgs > 0; imsgs--) {
		const struct message *msg = history_item(HISTORY, imsgs - 1);

		if (msg->time < history_time(HISTORY) - window)
			break;

		size_t ito, nto = msg->nto;
		for (ito = 0; ito < nto; ito++) {
			size_t to = msg->to[ito];
			const double *x = design_tvar(DESIGN, v, to);
			assert_true(x);
			assert_real_identical(x[0], 1.0);
		}
	}

	assert_true(v);
}


int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, fixture_setup_enron),
		unit_test_setup_teardown(test_mul0, setup, teardown),		
		unit_test_setup_teardown(test_tmul0, setup, teardown),
		//unit_test_setup_teardown(test_muls0, setup, teardown),
		//unit_test_setup_teardown(test_tmuls0, setup, teardown),
		unit_test_setup_teardown(test_irecvtot, setup, teardown),
		unit_test_teardown(enron_suite, fixture_teardown),
	};
	return run_tests(tests);
}


