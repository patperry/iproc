#include "port.h"
#include <assert.h>
#include "refcount.h"

bool refcount_init(struct refcount *refcount)
{
	assert(refcount);
	refcount_set(refcount, 1);
	return true;
}

void refcount_deinit (struct refcount *refcount)
{
	assert(refcount);
}

void refcount_set(struct refcount *refcount, int count)
{
	assert(refcount);
	assert(count > 0);
	assert(count <= REFCOUNT_MAX);
	refcount->count = count;
}

bool refcount_get(struct refcount *refcount)
{
	assert(refcount);

	if (refcount->count == REFCOUNT_MAX)
		return false;

	refcount->count++;
	return true;
}

bool refcount_put(struct refcount * refcount,
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
