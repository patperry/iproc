#ifndef _COMPARE_H
#define _COMPARE_H

typedef bool (*equals_fn) (const void *px, const void *py, void *udata);
typedef int (*compare_fn) (const void *px, const void *py, void *udata);

/* Example use:
 * 
 * DEFINE_COMPARE_AND_EQUALS_FN(int_compare, int_equals, int)
 */

#define DEFINE_COMPARE_FN(name, t) \
        static inline int name (const void *px, const void *py, void *udata)	\
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

#define DEFINE_EQUALS_FN(name, t) \
        static inline bool name (const void *px, const void *py, void *udata)	\
        {                        \
            t x = *(t *)px;      \
            t y = *(t *)py;      \
                                 \
            return (x == y);     \
        }

#define DEFINE_COMPARE_AND_EQUALS_FN(compare_name, equal_name, t) \
        DEFINE_COMPARE_FN(compare_name, t) \
        DEFINE_EQUALS_FN(equal_name, t)

#endif /* _COMPARE_H */
