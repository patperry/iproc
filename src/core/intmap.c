#include "port.h"
#include <assert.h>
#include "intmap.h"



bool _intmap_init (struct intmap *m, size_t elt_size)
{
    assert(m);
    assert(elt_size > 0);

    if (intset_init(&m->keys)) {
        if (_darray_init(&m->vals, elt_size)) {
            return m;
        }
        intset_deinit(&m->keys);
    }
    
    return NULL;
}


bool intmap_init_copy (struct intmap *m, const struct intmap *src)
{
    assert(m);
    assert(src);

    if (_intmap_init(m, intmap_elt_size(src))) {
        if (intmap_assign_copy(m, src)) {
            return m;
        }
        intmap_deinit(m);
    }
    return NULL;
}


void intmap_deinit (struct intmap *m)
{
    assert(m);
    
    darray_deinit(&m->vals);
    intset_deinit(&m->keys);
}


bool intmap_assign_copy (struct intmap *m, const struct intmap *src)
{
    assert(m);
    assert(src);
    
    if (intmap_reserve(m, intmap_size(src))) {
        intset_assign_copy(&m->keys, &src->keys);
        darray_assign_copy(&m->vals, &src->vals);
        return m;
    }
    
    return NULL;
}


void * intmap_copy_values_to (const struct intmap *m, void *dst)
{
    assert(m);
    assert(dst || intmap_size(m) == 0);
    
    return darray_copy_to(&m->vals, dst);
}


intptr_t * intmap_copy_keys_to (const struct intmap *m, intptr_t *dst)
{
    assert(m);
    assert(dst);
    
    return intset_copy_to(&m->keys, dst);
}


void intmap_clear (struct intmap *m)
{
    assert(m);
    intset_clear(&m->keys);
    darray_clear(&m->vals);
}


bool intmap_reserve (struct intmap *m, ssize_t n)
{
    assert(m);
    assert(n >= 0);
    
    if (intset_reserve(&m->keys, n)
        && darray_reserve(&m->vals, n)) {
        return m;
    }
    
    return NULL;
}


bool intmap_contains (const struct intmap *m, intptr_t key)
{
    assert(m);
    return intset_contains(&m->keys, key);
}


ssize_t intmap_index (const struct intmap *m, intptr_t key)
{
    assert(m);
    return intset_index(&m->keys, key);
}


void * intmap_lookup (const struct intmap *m, intptr_t key)
{
    assert(m);
    return intmap_lookup_with(m, key, NULL);
}


void * intmap_lookup_with (const struct intmap *m, intptr_t key, const void *val0)
{
    assert(m);
    struct intmap_pos pos;
    void *val;
    
    if ((val = intmap_find(m, key, &pos))) {
        return val;
    }
    
    return (void *)val0;
}


bool intmap_add (struct intmap *m, intptr_t key, const void *val)
{
    assert(m);
    assert(val);
    return intmap_add_all(m, &key, val, 1);
}


bool intmap_add_all (struct intmap *m, const intptr_t *keys, const void *vals, ssize_t n)
{
    assert(m);
    assert(keys || n == 0);
    assert(vals || n == 0);

    size_t elt_size = intmap_elt_size(m);
    struct intmap_pos pos;
    ssize_t i;
    
    if (!intmap_reserve(m, intmap_size(m) + n))
        return false;
    
    for (i = 0; i < n; i++, vals += elt_size) {
        if (intmap_find(m, keys[i], &pos)) {
            darray_set(&m->vals, pos.key.index, vals);
        } else {
            intset_insert(&m->keys, &pos.key);
            darray_insert(&m->vals, pos.key.index, vals);
        }
    }
    
    return true;
}


void intmap_remove (struct intmap *m, intptr_t key)
{
    assert(m);
    intmap_remove_all(m, &key, 1);
}


void intmap_remove_all (struct intmap *m, const intptr_t *keys, ssize_t n)
{
    assert(m);
    assert(keys || n == 0);
    
    struct intmap_pos pos;
    ssize_t i;
    
    for (i = 0; i < n; i++) {
        if (intmap_find(m, keys[i], &pos)) {
            intset_erase(&m->keys, &pos.key);
            darray_erase(&m->vals, pos.key.index);
        }
    }
}


void * intmap_find (const struct intmap *m, intptr_t key, struct intmap_pos *pos)
{
    assert(m);
    assert(pos);
    
    if (intset_find(&m->keys, key, &pos->key)) {
        return darray_ptr(&m->vals, pos->key.index);
    }
    
    return NULL;
}


void * intmap_insert (struct intmap *m, const struct intmap_pos *pos, const void *val)
{
    assert(m);
    assert(pos);
    
    if (!intmap_reserve(m, intmap_size(m) + 1))
        return NULL;
    
    intset_insert(&m->keys, &pos->key);
    return darray_insert(&m->vals, pos->key.index, val);
}


void intmap_erase (struct intmap *m, const struct intmap_pos *pos)
{
    assert(m);
    assert(pos);
    
    intset_erase(&m->keys, &pos->key);
    darray_erase(&m->vals, pos->key.index);
}

