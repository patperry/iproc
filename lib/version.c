#include "port.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "coreutil.h"
#include "xalloc.h"
#include "version.h"



void version_init(struct version *v)
{
	v->w = NULL;
	v->nw = 0;
	v->nwmax = 0;
}


void version_deinit(struct version *v)
{
	free(v->w);
}


void version_update(struct version *v)
{
	if (v->counter == SIZE_MAX) {
		size_t i, n = v->nw;
		for (i = 0; i < n; i++) {
			v->w[i]->counter = 0;
		}
		v->counter = 0;
	}

	v->counter++;
}


void version_watch_init(struct version_watch *w, struct version *v)
{
	if (needs_grow(v->nw + 1, &v->nwmax)) {
		v->w = xrealloc(v->w, v->nwmax * sizeof(w));
	}
	v->w[v->nw++] = w;
	w->counter = v->counter;
}


void version_watch_deinit(struct version_watch *w, struct version *v)
{
	size_t i, n = v->nw;
	for (i = n; i > 0; i--) {
		if (v->w[i - 1] == w) {
			memmove(v->w + i - 1, v->w + i, (n - i) * sizeof(w));
			break;
		}
	}
}


