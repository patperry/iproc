
#include "port.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "xalloc.h"

char **var_names_alloc(char *name, size_t len, size_t n)
{
	char **res = xcalloc(n + 1, sizeof(res[0]));
	size_t i;

	for (i = 0; i < n; i++) {
		size_t digits = 1 + (i + 1) / 10;
		size_t len1 = len + 2 + digits + 1;
		res[i] = xcalloc(len1, sizeof(res[i][0]));
		size_t nout =
		    snprintf(res[i], len1, "%s(%" SSIZE_FMT ")", name, i + 1);
		assert(nout + 1 == len1);
		(void)nout;
	}
	res[n] = NULL;

	return res;
}

char **var_names_alloc2(char *name, size_t len, size_t m, size_t n)
{
	char **res = xcalloc(m * n + 1, sizeof(res[0]));
	size_t i, j, ix;

	for (j = 0; j < n; j++) {
		size_t jdigits = 1 + (j + 1) / 10;
		for (i = 0; i < m; i++) {
			size_t idigits = 1 + (i + 1) / 10;
			size_t len1 = len + 3 + idigits + jdigits + 1;

			ix = i + j * m;
			res[ix] = xcalloc(len1, sizeof(res[ix][0]));
			size_t nout =
			    snprintf(res[ix], len1,
				     "%s(%" SSIZE_FMT ",%" SSIZE_FMT ")", name,
				     i + 1, j + 1);
			assert(nout + 1 == len1);
			(void)nout;
		}
	}
	res[m * n] = NULL;
	return res;
}

void var_names_free(char **names)
{
	if (names) {
		size_t i = 0;
		while (names[i]) {
			free(names[i]);
			i++;
		}
		free(names);
	}
}
