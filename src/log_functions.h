#ifndef __LOG_FUNCTIONS_ADAPTIVE_SUM_H__
#define __LOG_FUNCTIONS_ADAPTIVE_SUM_H__
#include <R.h>
#include <Rinternals.h>
#include "mathFun.h"

// To add functions implementations, they must be written here with the
// same parameters as the existent ones. Additionally, they must be assigned
// a "number" in the adapt_sum_callPrecomp function in the adapt_sum_compiled.c
// file switch statement. Finally, a condition must be included in the adapt_sum
// wrapper function in R in the adapt_sum.R file.

double Rf_dnbinom(double x, double size, double prob, int give_log);
double Rf_dbinom(double x, double n, double p, int give_log);
double Rf_dpois(double x, double lambda, int give_log);

double negbin_marginal(R_xlen_t k, double *Theta);

double noObs(R_xlen_t k, double *Theta);

double COMP(R_xlen_t k, double *Theta);

double dR0(R_xlen_t k, double *Theta);

double powerLawDiff(R_xlen_t k, double *Theta);

double negbin_sentinel(R_xlen_t k, double *Theta);

double poisson_sentinel(R_xlen_t k, double *Theta);

double weird_series_constL(R_xlen_t k, double *Theta);

double weird_series(R_xlen_t k, double *Theta);

#endif
