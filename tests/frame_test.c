#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include "cmockery.h"

#include "enron.h"
#include "design.h"
#include "frame.h"


static struct design design;
static struct frame frame;

int main(int argc, char **argv)
{
	UnitTest tests[1];
	return run_tests(tests);
}
