
#include "util.h"
#include "xalloc.h"
#include <stdlib.h>



void free2(void **ptrs, size_t len)
{
	size_t i;
	
	if (!ptrs)
		return;
	
	for (i = len; i > 0; i--) {
		free(ptrs[i-1]);
	}
	free(ptrs);
}


char **xstrdup2(const char *const *strs, size_t len)
{
	if (!strs)
		return NULL;

	char **res = xmalloc(len * sizeof(res[0]));
	size_t i;

	for (i = 0; i < len; i++) {
		res[i] = xstrdup(strs[i]);
	}
	return res;
}

