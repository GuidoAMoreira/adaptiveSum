#ifndef __MATH_FUNCTIONS_FROM_STAN_H__
#define __MATH_FUNCTIONS_FROM_STAN_H__

#include <R.h>
#include <Rinternals.h>
#include "string.h"

// Auxiliary functions roughly copied from the stan header files

static inline double log1m(double x)
{return log1p(-x);}

static inline double log1m_exp(double x)
{return (x>-0.693147?log(-expm1(x)):log1m(exp(x)));}

static inline double log_diff_exp(double x, double y)
{return (R_FINITE(y)?x + log1m_exp(y-x):x);}

static inline double log1p_exp(double x)
{return log1p(exp(x));}

static inline double logit(double logp)
{return (R_FINITE(logp)?-log_diff_exp(-logp,0):logp);}

#endif
