
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include "cmockery.h"

#include <float.h> // DBL_DIG, DBL_MANT_DIG
#include <math.h> // fabs
#include <stdlib.h> // rand
#include <string.h> // strlen
#include "ieee754.h" // double_identical

#include "testutil.h"


#define STR(x) #x
#define STRSTR(x) STR(x)
#define LargestRealTypePrintfFormat "%." STRSTR(DBL_DIG + 1) "g"

void setup_fixture(const char *name)
{
	size_t i, n = strlen(name);
	print_message("%s\n", name);
	for (i = 0; i < n; i++) {
		print_message("-");
	}
	print_message("\n");
}


void teardown_fixture()
{
	print_message("\n\n");
}

double runif(double min, double max)
{
	double u = (double)rand() / RAND_MAX;
	double x =  min + u * (max - min);
	return x;
}

static int reals_approx_display_error(const double left, const double right)
{
	const int minprec = DBL_MANT_DIG / 2 + 1;
	const int p = double_eqrel(left, right);
	const int rel_approx = p >= minprec;
	const int abs_approx = (fabs(left) <= DBL_EPSILON / 2
				&& fabs(right) <= DBL_EPSILON / 2);
	const int approx = rel_approx || abs_approx;
	
	if (!approx) {
		print_error(LargestRealTypePrintfFormat  " != "
			    LargestRealTypePrintfFormat " (%d bits < %d bits)\n",
			    left, right, p, minprec);
	}
	return approx;
	
}

static int reals_eqrel_display_error(const int precision, const double left, const double right)
{
	const int p = double_eqrel(left, right);
	const int eqrel = p >= precision;
	if (!eqrel) {
		print_error(LargestRealTypePrintfFormat  " != "
			    LargestRealTypePrintfFormat " (%d bits < %d bits)\n",
			    left, right, p, precision);
	}
	return eqrel;
	
}

static int reals_identical_display_error(const double left, const double right)
{
	const int identical = double_identical(left, right);
	if (!identical) {
		print_error(LargestRealTypePrintfFormat  " != "
				LargestRealTypePrintfFormat "\n", left, right);
	}
	return identical;

}


void _assert_real_approx(const double a, const double b,
			 const char * const file, const int line)
{
	if (!reals_approx_display_error(a, b)) {
		_fail(file, line);
	}
}

void _assert_real_eqrel(int precision, const double a, const double b,
			const char * const file, const int line)
{
	if (!reals_eqrel_display_error(precision, a, b)) {
		_fail(file, line);
	}
}

void _assert_real_identical(const double a, const double b, const char * const file, const int line)
{
	if (!reals_identical_display_error(a, b)) {
		_fail(file, line);
	}
}

