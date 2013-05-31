#ifndef TESTUTIL_H
#define TESTUTIL_H

#include "blas.h"


// print a message with the name of the fixture
void setup_fixture(const char *name);

// print a message indicating the end of the fixture
void teardown_fixture();

// generate a uniformly distributed variate in the range
// [min, max] using the system's random number generator (rand)
double runif(double min, double max);


// Assert that the two given real numbers are approximately equal, otherwise fail.
#define assert_real_approx(a, b) \
	_assert_real_approx((double)(a), \
			    (double)(b), \
			    __FILE__, __LINE__)


// Assert that the two given real numbers are equal relative to the
// given precision, otherwise fail.
#define assert_real_eqrel(bits_of_precision, a, b) \
	_assert_real_eqrel((bits_of_precision), \
			   (double)(a), \
			   (double)(b), \
			   __FILE__, __LINE__)

// Assert that the two given real numbers are identical, otherwise fail.
#define assert_real_identical(a, b) \
	_assert_real_identical((double)(a), \
	                       (double)(b), \
	                       __FILE__, __LINE__)

// Assert that the two given real vectors are identical, otherwise fail.
#define assert_vec_identical(a, b, n) \
	_assert_vec_identical((a), (b), (n), \
			      __FILE__, __LINE__)

// Assert that the two given symmetric matrices are approximately equal,
// otherwise fail.
#define assert_sym_approx(a, b, uplo, n) \
	_assert_sym_approx((a), (b), (uplo), (n), \
			   __FILE__, __LINE__)


void _assert_real_approx(const double a, const double b,
			 const char * const file, const int line);

void _assert_real_eqrel(int precision, const double a, const double b,
			const char * const file, const int line);

void _assert_real_identical(const double a, const double b,
			    const char * const file, const int line);

void _assert_vec_identical(const double *a, const double *b, size_t n,
			   const char * const file, const int line);

void _assert_sym_approx(const double *a, const double *b, enum blas_uplo uplo,
			size_t n, const char * const file, const int line);

#endif // TESTUTIL_H
