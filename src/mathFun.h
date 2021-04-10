#ifndef __MATH_FUNCTIONS_FROM_STAN_H__
#define __MATH_FUNCTIONS_FROM_STAN_H__

#include <Rinternals.h>
#include <Rmath.h>
#include "math.h"

// Auxiliary functions roughly copied from the stan header files

long double Rf_log1mexp(long double);

static inline long double log1m(long double x)
{return log1p(-x);}

static inline long double log1m_exp(long double x)
{return (x > -0.693147 ? log(-expm1(x)) : log1m(exp(x)));}

static inline long double log_diff_exp(long double x, long double y)
{return (R_FINITE(y) ? x + Rf_log1mexp(y - x) : x);}

static inline long double logit(long double logp)
{return (R_FINITE(logp) ? -log_diff_exp(-logp, 0) : logp);}

long double partial_logSumExp(long double* fun, R_xlen_t til, long double maxA);

#endif
