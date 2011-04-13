
#ifndef _IPROC_IEEE754_H
#define _IPROC_IEEE754_H

#include <stdint.h>

#define  IPROC_MAX_NAN_PAYLOAD  UINT64_C(0x0007FFFFFFFFFFFF)

int iproc_identical(double x, double y);
double iproc_nextup(double x);
double iproc_nextdown(double x);
double iproc_ieeemean(double x, double y);
double iproc_mknan(uint64_t payload);
uint64_t iproc_getnan(double x);
int iproc_feqrel(double x, double y);

#endif /* _IPROC_IEEE754_H */
