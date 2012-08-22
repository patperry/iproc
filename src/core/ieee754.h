
#ifndef _IEEE754_H
#define _IEEE754_H

#include <stdint.h>

#define  MAX_NAN_PAYLOAD  UINT64_C(0x0007FFFFFFFFFFFF)

#define SQRT_DBL_EPSILON	1.4901161193847656e-08
#define ROOT3_DBL_EPSILON	6.0554544523933429e-06
#define ROOT4_DBL_EPSILON	1.2207031250000000e-04
#define ROOT5_DBL_EPSILON	7.4009597974140505e-04
#define ROOT6_DBL_EPSILON	2.4607833005759251e-03
#define LOG_DBL_MAX		7.0978271289338397e+02


int double_identical(double x, double y);
double double_nextup(double x);
double double_nextdown(double x);
double double_ieeemean(double x, double y);
double double_mknan(uint64_t payload);
uint64_t double_getnan(double x);
int double_eqrel(double x, double y);

int double_equals(const void *x, const void *y);
int double_compare(const void *x, const void *y);
int double_rcompare(const void *x, const void *y);

#endif /* _IEEE754_H */
