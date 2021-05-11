#ifndef __MATH_FUNCTIONS_FROM_STAN_H__
#define __MATH_FUNCTIONS_FROM_STAN_H__

#include <Rinternals.h>

void partial_logSumExp(long double*, R_xlen_t, long double,
                       long double*, int, long double*);

// This function helps reduce the floating point rounding error.
static inline void KahanSum(long double* tot, long double x, long double* c)
{
  long double t, y;
  y = x - *c;
  t = *tot + y;
  *c = t - *tot - y;
  *tot = t;
}

#endif
