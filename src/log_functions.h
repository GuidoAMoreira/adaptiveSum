#ifndef __LOG_FUNCTIONS_ADAPTIVE_SUM_H__
#define __LOG_FUNCTIONS_ADAPTIVE_SUM_H__
#include <Rinternals.h>

// To add functions implementations, they must be written here with the
// same parameters as the existent ones. Additionally, they must be assigned
// a "number" in the adapt_sum_callPrecomp function in the adapt_sum_compiled.c
// file switch statement. Finally, a condition must be included in the adapt_sum
// wrapper function in R in the adapt_sum.R file.

long double negbin_marginal(R_xlen_t k, double *Theta);

long double noObs(R_xlen_t k, double *Theta);

long double COMP(R_xlen_t k, double *Theta);

long double dR0(R_xlen_t k, double *Theta);

long double powerLawDiff(R_xlen_t k, double *Theta);

long double negbin_sentinel(R_xlen_t k, double *Theta);

long double poisson_sentinel(R_xlen_t k, double *Theta);

long double weird_series_constL(R_xlen_t k, double *Theta);

long double weird_series(R_xlen_t k, double *Theta);

#endif
