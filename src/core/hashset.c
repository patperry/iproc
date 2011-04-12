#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <sys/queue.h>
#include "hashset.h"

struct hashset_node {
    SLIST_ENTRY(hashset_node) nodes;
};

struct hashset_bucket {
    SLIST_HEAD(nodes_head, hashset_node) collisions;
};


static bool bucket_init (struct hashset *s, struct hashset_bucket *b)
{
    SLIST_INIT(&b->collisions);
    return b;
}

static void bucket_deinit (struct hashset *s, struct hashset_bucket *b)
{
    return;
}

static size_t bucket_size (const struct hashset *s)
{
    return intmap_elt_size(&s->buckets);
}


static size_t node_size (const struct hashset *s)
{
    return bucket_size(s);
}


static struct hashset_node * node_alloc (struct hashset *s)
{
    if (trashstack_empty(&s->nodes)) {
        return malloc(node_size(s));        
    } else {
        return trashstack_pop(&s->nodes);
    }
}


static void node_free (struct hashset *s, struct hashset_node *n)
{
    trashstack_push(&s->nodes, n);
}


static void * node_val (const struct hashset *s, const struct hashset_node *n)
{
    return (void *)n + s->elt_offset;
}


static void * bucket_val (const struct hashset *s, const struct hashset_bucket *b)
{
    return node_val(s, (const struct hashset_node *)b);
}



bool _hashset_init (struct hashset *s, hash_fn hash, equals_fn equals,
                    size_t elt_size, size_t elt_offset)
{
    assert(s);
    assert(hash);
    assert(equals);
    assert(elt_size > 0);
    assert(elt_offset >= sizeof(struct hashset_node));

    size_t bucket_size = elt_offset + elt_size;
    if (_intmap_init(&s->buckets, bucket_size)) {
        s->hash = hash;
        s->equals = equals;
        s->elt_size = elt_size;
        s->elt_offset = elt_offset;
        trashstack_init(&s->nodes);
        return s;
    }
    
    return NULL;
}

/*
bool hashset_init_copy (struct hashset *s, const struct hashset *src)
{
    assert(s);
    assert(src);
    
    if (_hashset_init(s, src->hash, src->equals, src->elt_size, src->elt_offset)) {
        if (hashset_assign_copy(s, src)) {
            return s;
        }
        hashset_deinit(s);
    }
    return NULL;
}
 */


static void hashset_deinit_nodes (struct trashstack *nodes)
{
    struct node *n;
    
    while (!trashstack_empty(nodes)) {
        n = trashstack_pop(nodes);
        free(n);
    }
    
}


void hashset_deinit (struct hashset *s)
{
    hashset_clear(s);
    intmap_deinit(&s->buckets);
    hashset_deinit_nodes(&s->nodes);
}


void hashset_clear (struct hashset *s)
{
    struct hashset_bucket *buckets, *b;
    struct hashset_node *c;
    ssize_t i, n = intmap_size(&s->buckets);
    
    buckets = intmap_vals_begin(&s->buckets);
    
    for (i = 0; i < n; i++) {
        b = &buckets[i];
        while (!SLIST_EMPTY(&b->collisions)) {
            c = SLIST_FIRST(&b->collisions);
            SLIST_REMOVE_HEAD(&b->collisions, nodes);
            node_free(s, c);
        }
    }
    
    intmap_clear(&s->buckets);
}


// assume worst case: each added value is a collision
bool hashset_reserve (struct hashset *s, ssize_t n)
{
    assert(s);
    assert(n >= 0);
    
    if (!intmap_reserve(&s->buckets, n))
        return false;
    
    ssize_t i, n0 = trashstack_size(&s->nodes);
        
    for (i = n0; i < n; i++) {
        struct hashset_node *node = malloc(sizeof(*node));
        if (!node)
            return false;
        trashstack_push(&s->nodes, node);
    }

    return true;
}


void * hashset_add (struct hashset *s, const void *val)
{
    assert(s);
    assert(val);

    struct hashset_pos pos;
    
    if (hashset_find(s, val, &pos)) {
        return hashset_replace(s, &pos, val);
    } else {
        return hashset_insert(s, &pos, val);
    }
}


bool hashset_add_all (struct hashset *s, const void *vals, ssize_t n)
{
    assert(s);
    assert(vals || n == 0);
    
    size_t elt_size = hashset_elt_size(s);
    ssize_t i;
    
    for (i = 0; i < n; i++, vals += elt_size) {
        if (!hashset_add(s, vals))
            goto undo;
    }
    return true;
undo:
    for (vals -= elt_size; i > 0; i--, vals -= elt_size) {
        hashset_remove(s, vals);
    }
    return false;
}


void hashset_remove (struct hashset *s, const void *key)
{
    assert(s);
    assert(key);
    
    struct hashset_pos pos;

    if (hashset_find(s, key, &pos)) {
        hashset_erase(s, &pos);
    }
}


void hashset_remove_all (struct hashset *s, const void *keys, ssize_t n)
{
    assert(s);
    assert(keys || n == 0);
    assert(n >= 0);

    size_t elt_size = hashset_elt_size(s);
    const void *end = keys + n * elt_size;

    for (; keys < end; keys += elt_size) {
        hashset_remove(s, keys);
    }
}


void * hashset_find (const struct hashset *s, const void *key,
                     struct hashset_pos *pos)
{
    assert(s);
    assert(key);
    assert(pos);

    struct hashset_bucket *bucket = NULL;
    struct hashset_node   *node   = NULL;
    void                  *val    = NULL;

    bucket = intmap_find(&s->buckets, hashset_hash(s, key), &pos->key);
    
    // bucket exists
    if (bucket) {
        val = bucket_val(s, bucket);
    
        // key is in bucket
        if (hashset_equals(s, key, val)) {
            pos->node.b = bucket;
            goto out;
        }
    
        SLIST_FOREACH(node, &bucket->collisions, nodes) {
            val = node_val(s, node);

            // key is in node
            if (hashset_equals(s, key, val)) {
                pos->node.n = node;
                goto out;
            }
        }
        pos->node.n= NULL;
    }
    val = NULL;
out:
    pos->bucket = bucket;
    return val;
}


void * hashset_insert (struct hashset *s, struct hashset_pos *pos, const void *val)
{
    assert(s);
    assert(pos);
    assert(!pos->bucket || !pos->node.n);

    struct hashset_bucket *bucket = pos->bucket;
    struct hashset_node   *node;
    void *ptr = NULL;
    
    // create a bucket
    if (!bucket) {
        if ((bucket = intmap_insert(&s->buckets, &pos->key, NULL))) {
            bucket_init(s, bucket);
            ptr = bucket_val(s, bucket);
            pos->bucket = bucket;
            goto success;
        }

    // create a node
    } else {
        if ((node = node_alloc(s))) {
            ptr = node_val(s, node);
            pos->node.n = node;
            goto success;
        }
    }
    return NULL;

success:
    if (val) {
        assert(hashset_hash(s, val) == pos->key.key.value);
        assert(pos->node.b == pos->bucket || hashset_equals(s, val, bucket_val(s, pos->node.b)));
        assert(pos->node.b != pos->bucket || hashset_equals(s, val, node_val(s, pos->node.n)));
        memcpy(ptr, val, hashset_elt_size(s));
    }
    return ptr;
}


void * hashset_replace (struct hashset *s, struct hashset_pos *pos, const void *val)
{
    assert(s);
    assert(pos);
    assert(pos->bucket);
    assert(pos->node.n || pos->node.b == pos->bucket);

    void *dst;
    
    if (!val) {
        hashset_erase(s, pos);
        return NULL;
    } else if (pos->node.b == pos->bucket) {
        dst = bucket_val(s, pos->bucket);
    } else {
        dst = node_val(s, pos->node.n);
    }
    
    memcpy(dst, val, hashset_elt_size(s));
    return dst;
}


void hashset_erase (struct hashset *s, struct hashset_pos *pos)
{
    assert(s);
    assert(pos);
    assert(pos->bucket);
    assert(pos->node.n || pos->node.b == pos->bucket);
    struct hashset_bucket *bucket = pos->bucket;
    
    // remove a node
    if (pos->node.b != pos->bucket) { 
        SLIST_REMOVE(&bucket->collisions, pos->node.n, hashset_node, nodes);
        node_free(s, pos->node.n);
        pos->node.n = NULL;

    // remove a bucket     
    } else if (SLIST_EMPTY(&bucket->collisions)) {
        bucket_deinit(s, bucket);
        intmap_erase(&s->buckets, &pos->key);
        pos->bucket = NULL;

    // move a node to a bucket    
    } else { 
        struct hashset_node *node = SLIST_FIRST(&bucket->collisions);
        SLIST_REMOVE_HEAD(&bucket->collisions, nodes);
        memcpy(bucket_val(s, bucket), node_val(s, node), hashset_elt_size(s));
        node_free(s, node);
    }
}


