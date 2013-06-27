#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "xalloc.h"


void xalloc_die(void)
{
	error("memory allocation failed");
}
