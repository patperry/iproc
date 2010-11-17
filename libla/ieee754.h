
#ifndef _LA_IEEE754_H
#define _LA_IEEE754_H

#include <stdint.h>

int      la_identical (double   x,
                       double   y);
double   la_nextup    (double   x);
double   la_nextdown  (double   x);
double   la_ieeemean  (double   x,
                       double   y);
double   la_mknan     (uint64_t payload);
uint64_t la_getnan    (double   x);
int      la_feqrel    (double   x,
                       double   y);

#endif /* _LA_IEEE754_H */
