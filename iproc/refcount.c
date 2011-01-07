#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <iproc/refcount.h>

void
iproc_refcount_init (iproc_refcount *refcount)
{
    assert(refcount);
    iproc_refcount_set(refcount, 1);
}

void
iproc_refcount_set (iproc_refcount *refcount,
                    int             count)
{
    assert(refcount);
    assert(count > 0);
    refcount->count = count;
}


void
iproc_refcount_get (iproc_refcount *refcount)
{
    assert(refcount);
    (refcount->count)++;
}

int
iproc_refcount_put (iproc_refcount *refcount,
                    void (*release) (iproc_refcount *refcount))
{
    assert(refcount);
    assert(release);

    (refcount->count)--;
    if (refcount->count == 0) {
        release(refcount);
        return 1;
    } else {
        return 0;
    }
}