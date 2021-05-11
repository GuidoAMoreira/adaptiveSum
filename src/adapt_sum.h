#ifndef __ADAPTIVE_APPROXIMATION_OF_INFINITE_SUM_H__
#define __ADAPTIVE_APPROXIMATION_OF_INFINITE_SUM_H__

#include <Rinternals.h>

#define l2 0.69314718055994530941723212145817656831

double Rf_logspace_sum(const double*, int);
double Rf_logspace_sub(double, double);
double Rf_log1pexp(double);

// Some macros
static inline long double logz(double loga, double logap1)
{return logap1 + loga - Rf_logspace_sub(loga, logap1);}

static inline long double delta(double logz, double loga, double log1ml)
{
  double ls = loga - log1ml;
  return (logz > ls ? Rf_logspace_sub(logz, ls) : Rf_logspace_sub(ls, logz));
}

static inline double feval(SEXP lF, SEXP rho)
{return REAL(eval(lF, rho))[0];}

static inline SEXP retFun(double res, R_xlen_t mI)
{
  SEXP out = PROTECT(allocVector(VECSXP,2));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(res));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(mI));

  UNPROTECT(1);
  return out;
}

//// Adaptive functions ////
// Calculating the adaptive sum method with known L
SEXP adapt_sum(SEXP logFun, SEXP params, SEXP epsilon, SEXP maxIter_, SEXP logL_,
               SEXP n0_, SEXP rho);

// Calculating the adaptive sum method with known L - pre-compiled code
SEXP adapt_sum_precomp(long double logFun(R_xlen_t k, double *Theta),
                       double *params, double eps,
                       R_xlen_t maxIter, double logL, R_xlen_t n0);

// Wrapper function for the pre-compiled code
SEXP adapt_sum_callPrecomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                           SEXP logL, SEXP n0);

//// Naive functions ////
// Calculating the adaptive sum method with known L
SEXP naive_sum(SEXP logFun, SEXP params, SEXP epsilon, SEXP maxIter_,
               SEXP n0_, SEXP rho);

// Calculating the adaptive sum method with known L - pre-compiled code
SEXP naive_sum_precomp(long double logFun(R_xlen_t k, double *Theta),
                       double *params, double eps,
                       R_xlen_t maxIter, R_xlen_t n0);

// Wrapper function for the pre-compiled code
SEXP naive_sum_callPrecomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                           SEXP n0);

#endif
