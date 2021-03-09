#ifndef __ADAPTIVE_APPROXIMATION_OF_INFINITE_SUM_H__
#define __ADAPTIVE_APPROXIMATION_OF_INFINITE_SUM_H__

#include <R.h>
#include <Rinternals.h>
#include "mathFun.h"

// Some macros
static inline double logz(double loga, double logap1)
{return logap1 - log_diff_exp(0,logap1 - loga);}

static inline double delta(double logz, double loga, double logr)
{
  double ls = loga + logr;
  return (logz>ls?log_diff_exp(logz,ls):log_diff_exp(ls,logz));
}

static inline double feval(SEXP lF, SEXP rho)
{return REAL(eval(lF, rho))[0];}

static inline SEXP retFun(SEXP res, SEXP mI)
{
  SEXP out = PROTECT(allocVector(VECSXP,2));
  SET_VECTOR_ELT(out,0,res);
  SET_VECTOR_ELT(out,1,mI);

  UNPROTECT(1);
  return out;
}

// Calculating the adaptive sum method with known L
SEXP adapt_sum(SEXP logFun, SEXP params, SEXP epsilon, SEXP maxIter_, SEXP logL_,
               SEXP n0_, SEXP rho);

// Calculating the adaptive sum method with known L - pre-compiled code
SEXP adapt_sum_precomp(double logFun(R_xlen_t k, double *Theta),
                       double *params, double eps,
                       R_xlen_t maxIter, double logL, R_xlen_t n0);

// Wrapper function for the pre-compiled code

SEXP adapt_sum_callPrecomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                           SEXP logL, SEXP n0);

#endif
