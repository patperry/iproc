#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include <stdio.h>
#include "f77.h"

void xerbla(char *srname, f77int *info);
void xerbla_(char *srname, f77int *info);

void xerbla(char *srname, f77int *info)
{
	int i = *info;
	printf("** On entry to %6s, parameter number %2d had an illegal value\n",
		srname, i);
}


void xerbla_(char *srname, f77int *info)
{
	xerbla(srname, info);
}

