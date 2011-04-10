#ifndef _COMPARE_H
#define _COMPARE_H


typedef bool (*equal_fn)   (const void *px, const void *py);
typedef int  (*compare_fn) (const void *px, const void *py);


#define DEFINE_COMPARE_FN(name, t) \
        static inline int name (const void *px, const void *py) \
        {                        \
            t x = *(t *)px;      \
            t y = *(t *)py;      \
                                 \
            if (x < y) {         \
                return -1;       \
            } else if (x > y) {  \
                return +1;       \
            } else {             \
                return 0;        \
            }                    \
        }

#define DEFINE_EQUAL_FN(name, t) \
        static inline bool name (const void *px, const void *py) \
        {                        \
            t x = *(t *)px;      \
            t y = *(t *)py;      \
                                 \
            return (x == y);     \
        }

#define DEFINE_COMPARE_AND_EQUAL_FN(compare_name, equal_name, t) \
        DEFINE_COMPARE_FN(compare_name, t); \
        DEFINE_EQUAL_FN(equal_name, t);


DEFINE_COMPARE_AND_EQUAL_FN(int64_compare,  int64_equal,  int64_t);
DEFINE_COMPARE_AND_EQUAL_FN(uint64_compare, uint64_equal, uint64_t);
DEFINE_COMPARE_AND_EQUAL_FN(size_compare,   size_equal,   size_t);
DEFINE_COMPARE_AND_EQUAL_FN(ssize_compare,  ssize_equal,  ssize_t);


/* compare using uint64_t instead of double to avoid dealing with NaN
 * comparisons; this relies on IEEE doubles being 8 bytes and lexicographically
 * ordered, and uint64_t having the same endianness and alignment as double */
#define double_compare uint64_compare
#define double_equal   uint64_equal


#endif /* _COMPARE_H */
