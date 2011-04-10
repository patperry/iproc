#ifndef _IPROC_HASH_H
#define _IPROC_HASH_H

#include <stddef.h>

typedef uint32_t (*hash_fn) (const void *key);

static inline void hash_combine  (uint32_t *seedp, uint32_t hash);
static inline void hash_finalize (uint32_t *hashp);

uint32_t memory_hash (const void *ptr, ssize_t n);
uint32_t double_hash (const void *val);



/* based on http://www.azillionmonkeys.com/qed/hash.html */
void hash_combine (uint32_t *seedp, uint32_t hash)
{
    union { uint32_t value; int16_t word[2]; } v;
    uint32_t seed, tmp;
    
    seed = *seedp;
    v.value = hash;
    
    seed += v.word[0];
    tmp   = (v.word[1] << 11) ^ seed;
    hash  = (seed << 16) ^ tmp;
    seed += seed >> 11;
    
    *seedp = seed;
}


/* force "avalanching" of final 127 bits */
void hash_finalize (uint32_t *hashp)
{
    uint32_t hash = *hashp;
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;
    
    *hashp = hash;
}


#endif /* _IPROC_HASH_H */
