#ifndef _TESTUTIL_H
#define _TESTUTIL_H


// print a message with the name of the fixture
void setup_fixture(const char *name);

// print a message indicating the end of the fixture
void teardown_fixture();



// Assert that the two given real numbers are equal, otherwise fail.
#define assert_real_identical(a, b) \
	_assert_real_identical((double)(a), \
			   (double)(b), \
			   __FILE__, __LINE__)

void _assert_real_identical(const double a, const double b,
		const char * const file, const int line);


#endif // _TESTUTIL_H
