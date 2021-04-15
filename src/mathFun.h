#ifndef __MATH_FUNCTIONS_FROM_STAN_H__
#define __MATH_FUNCTIONS_FROM_STAN_H__

#include <Rinternals.h>

// Auxiliary functions roughly copied from the stan header files

long double partial_logSumExp(long double* fun, R_xlen_t til, long double maxA);

#endif
