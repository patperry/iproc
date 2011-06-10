#include "port.h"
#include <assert.h>
#include "refcount.h"

#define REFCOUNT_MAX INT_MAX


void refcount_init(struct refcount *refcount)
{
	assert(refcount);
	refcount_set(refcount, 1);
}

void refcount_deinit(struct refcount *refcount)
{
	assert(refcount);
	assert(refcount->count == 1);
}

void refcount_set(struct refcount *refcount, int count)
{
	assert(refcount);
	assert(count > 0);
	assert(count <= REFCOUNT_MAX);
	refcount->count = count;
}

void refcount_get(struct refcount *refcount)
{
	assert(refcount);
	assert(refcount->count < REFCOUNT_MAX);
	refcount->count++;
}

bool refcount_put(struct refcount *refcount,
		  void (*release) (struct refcount * refcount))
{
	assert(refcount);
	assert(refcount->count > 0);

	refcount->count--;
	if (refcount->count == 0) {
		if (release)
			release(refcount);
		return true;
	} else {
		return false;
	}
}
