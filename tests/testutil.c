
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include "cmockery.h"

#include <float.h> // DBL_DIG
#include <string.h> // strlen
#include "ieee754.h" // double_identical

#include "testutil.h"


#define STR(x) #x
#define STRSTR(x) STR(x)
#define LargestRealTypePrintfFormat "%." STRSTR(DBL_DIG) "g"

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


static int reals_identical_display_error(const double left, const double right)
{
	const int identical = double_identical(left, right);
	if (!identical) {
		print_error(LargestRealTypePrintfFormat  " != "
				LargestRealTypePrintfFormat "\n", left, right);
	}
	return identical;

}


void _assert_real_identical(const double a, const double b, const char * const file, const int line)
{
	if (!reals_identical_display_error(a, b)) {
		_fail(file, line);
	}
}
