#include "port.h"
#include <assert.h>
#include "hashset.h"


#define HASH_INDEX(h) ((ssize_t)(SIZEOF_SIZE_T <= 4 ? (h) >> 1 : (h)))


static struct hashset_bucket * hashset_bucket_init (const struct hashset *s,
                                                    struct hashset_bucket *b)
{
    assert(s);
    assert(b);

    
    if (darray_assign_copy(&b->values, &(s->empty_bucket).values)) {
        return b;
    }
    
    return NULL;
}


static void hashset_bucket_deinit (const struct hashset *s,
                                   struct hashset_bucket *b)
                                   
{
    assert(s);    
    assert(b);

    struct darray *values = &b->values;
    destroy_fn destroy = s->destroy;
    
    if (destroy) {
        size_t elt_size = darray_elt_size(values);
        const void *end = darray_end(values);
        void *ptr;
        
        for (ptr = darray_begin(values); ptr != end; ptr += elt_size) {
            destroy(ptr);
        }
    }
    
    darray_deinit(values);
}


struct hashset * _hashset_init (struct hashset *s, hash_fn hash, equals_fn equal,
                                size_t elt_size)
{
    assert(s);
    assert(hash);
    assert(equal);
    assert(elt_size > 0);
    
    return _hashset_init_with_destroy(s, hash, equal, NULL, elt_size);
}


struct hashset * _hashset_init_with_destroy (struct hashset *s, hash_fn hash,
                                             equals_fn equal, destroy_fn destroy,
                                             size_t elt_size)
{
    assert(s);
    assert(hash);
    assert(equal);
    assert(elt_size > 0);

    if (sarray_init(&s->buckets, struct hashset_bucket)) {
        if (_darray_init(&(s->empty_bucket).values, elt_size)) {
            s->hash = hash;
            s->equal = equal;
            s->destroy = destroy;
            return s;
        }
        sarray_deinit(&s->buckets);
    }

    return NULL;
}


void hashset_deinit (struct hashset *s)
{
    assert(s);
    assert(darray_size(&(s->empty_bucket).values) == 0);
    
    hashset_clear(s);
    darray_deinit(&(s->empty_bucket).values);
    sarray_deinit(&s->buckets);
}


void hashset_clear (struct hashset *s)
{
    assert(s);
    
    struct sarray *buckets = &s->buckets;
    const struct hashset_bucket *buckets_end = sarray_end(buckets);
    struct hashset_bucket *b;
    
    for (b = sarray_begin(buckets); b != buckets_end; b++) {
        hashset_bucket_deinit(s, b);
    }
}


ssize_t hashset_size (const struct hashset *s)
{
    assert(s);
    
    ssize_t n = 0;
    struct hashset_bucket *end = sarray_end(&s->buckets);
    struct hashset_bucket *b;
    
    for (b = sarray_begin(&s->buckets); b != end; b++) {
        n += darray_size(&b->values);
    }
    
    return n;
}


const void * hashset_find (const struct hashset *s, const void *key)
{
    assert(s);

    return hashset_find_with(s, key, NULL);
}


const void * hashset_find_with (const struct hashset *s, const void *key,
                                const void *val0)
{
    assert(s);

    ssize_t index = HASH_INDEX(s->hash(key));
    const struct hashset_bucket *b = sarray_find(&s->buckets, index);

    if (b) {
        ssize_t pos = darray_find_index(&b->values, key, s->equal);
        if (pos >= 0) {
            return darray_ptr(&b->values, pos);
        }
    }
    
    return val0;
}


void * hashset_add (struct hashset *s, const void *val)
{
    assert(s);
    assert(val);

    ssize_t index = HASH_INDEX(s->hash(val));
    struct hashset_bucket *b = &sarray_index_with(&s->buckets,
                                                  struct hashset_bucket,
                                                  index,
                                                  &s->empty_bucket);
    ssize_t pos = darray_find_index(&b->values, val, s->equal);
    if (pos >= 0) {
        if (s->destroy) {
            s->destroy(darray_ptr(&b->values, pos));
        }
        darray_set(&b->values, pos, val);
        return (void *)val + hashset_elt_size(s);
    } else {
        return darray_push_back(&b->values, val);
    }
}


void * hashset_add_all (struct hashset *s, const void *ptr, ssize_t n)
{
    assert(s);
    assert(ptr || n == 0);
    assert(n >= 0);

    ssize_t i;
    
    for (i = 0; i < n; i++) {
        ptr = hashset_add(s, ptr);
    }

    return (void *)ptr;
}


void hashset_remove (struct hashset *s, const void *key)
{
    assert(s);
    
    ssize_t index = HASH_INDEX(s->hash(key));
    struct hashset_bucket *b;
    ssize_t pos;
    
    if ((b = sarray_find(&s->buckets, index))
        && ((pos = darray_find_index(&b->values, key, s->equal)) >= 0)) {
        
        if (s->destroy)
            s->destroy(darray_ptr(&b->values, pos));
            
        darray_erase(&b->values, pos);
        
        if (darray_empty(&b->values)) {
            darray_deinit(&b->values);
            
            sarray_remove(&s->buckets,
                          b - 
                          (struct hashset_bucket *)sarray_begin(&s->buckets));
        }
        
    }
    
}


void hashset_remove_all (struct hashset *s, const void *ptr, ssize_t n)
{
    assert(s);
    assert(ptr || n == 0);
    assert(n >= 0);
    

    size_t elt_size = hashset_elt_size(s);
    const void *end = ptr + n * elt_size;

    for (; ptr < end; ptr += elt_size) {
        hashset_remove(s, ptr);
    }
}





