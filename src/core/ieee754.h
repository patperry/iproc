
#ifndef _IEEE754_H
#define _IEEE754_H

#define  MAX_NAN_PAYLOAD  UINT64_C(0x0007FFFFFFFFFFFF)

bool double_identical(double x, double y);
double double_nextup(double x);
double double_nextdown(double x);
double double_ieeemean(double x, double y);
double double_mknan(uint64_t payload);
uint64_t double_getnan(double x);
int double_eqrel(double x, double y);

bool double_equals(const void *x, const void *y, void *udata);
int double_compare(const void *x, const void *y, void *udata);
int double_rcompare(const void *x, const void *y, void *udata);

#endif /* _IEEE754_H */
