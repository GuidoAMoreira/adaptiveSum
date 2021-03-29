#include "adapt_sum.h"
#include "log_functions.h"

SEXP adapt_sum_precomp(double logFun(R_xlen_t k, double *Theta),
                       double *params, double eps,
                       R_xlen_t maxIter, double logL, R_xlen_t n0)
{
  // Declaration
  SEXP result = PROTECT(allocVector(REALSXP, maxIter+1)),
    i_SEXP = PROTECT(allocVector(INTSXP,1));
  R_xlen_t n = 0;
  double *logFunVal = REAL(result), lEps = log(eps) + log(2), logap1;

  // Find the maximum
  logFunVal[n] = logFun(n0, params);
  do
  {logFunVal[++n] = logFun(++n0, params);}
  while (!R_FINITE(logFunVal[n]) ||
    (logFunVal[n] >= logFunVal[n - 1] &&
    n <= (maxIter - 1)));

  // If too many iterations
  if (n == maxIter)
  {
    logFunVal[maxIter] = 1;
    INTEGER(i_SEXP)[0] = maxIter;
    UNPROTECT(2);
    return retFun(result, i_SEXP);
  }

  // Calculate the tail
  do
  {logFunVal[++n] = logFun(++n0, params);}
  while ((delta(logz(logFunVal[n - 1], logFunVal[n]),
                  logFunVal[n], logL) >= lEps) & (n <= (maxIter - 1)));

  // If too many iterations
  if (n == maxIter)
  {
    logFunVal[maxIter] = 1;
    INTEGER(i_SEXP)[0] = maxIter;
    UNPROTECT(2);
    return retFun(result, i_SEXP);
  }

  logap1 = logFunVal[n];
  logFunVal[n] = logz(logFunVal[n - 1], logap1) - log(2);
  logFunVal[n + 1] = logap1 - log_diff_exp(0, logL) - log(2);

  INTEGER(i_SEXP)[0] = n;
  UNPROTECT(2);
  return retFun(result, i_SEXP);
}

SEXP adapt_sum_callPrecomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                           SEXP logL, SEXP n0)
{
  unsigned int funSelect = INTEGER(lF)[0];

  switch (funSelect)
  {
  case 1:
    return adapt_sum_precomp(negbin_marginal,REAL(params),REAL(epsilon)[0],
                             INTEGER(maxIter)[0],REAL(logL)[0],INTEGER(n0)[0]);
  case 2:
    return adapt_sum_precomp(noObs,REAL(params),REAL(epsilon)[0],
                             INTEGER(maxIter)[0],REAL(logL)[0],INTEGER(n0)[0]);
  case 3:
    return adapt_sum_precomp(COMP,REAL(params),REAL(epsilon)[0],
                             INTEGER(maxIter)[0],REAL(logL)[0],INTEGER(n0)[0]);
  case 4:
    return adapt_sum_precomp(dR0,REAL(params),REAL(epsilon)[0],
                             INTEGER(maxIter)[0],REAL(logL)[0],INTEGER(n0)[0]);
  case 5:
    return adapt_sum_precomp(powerLawDiff,REAL(params),REAL(epsilon)[0],
                             INTEGER(maxIter)[0],REAL(logL)[0],INTEGER(n0)[0]);
  case 6:
    return adapt_sum_precomp(negbin_sentinel,REAL(params),REAL(epsilon)[0],
                             INTEGER(maxIter)[0],REAL(logL)[0],INTEGER(n0)[0]);
  case 7:
    return adapt_sum_precomp(poisson_sentinel,REAL(params),REAL(epsilon)[0],
                             INTEGER(maxIter)[0],REAL(logL)[0],INTEGER(n0)[0]);
  case 8:
    return adapt_sum_precomp(weird_series_constL,REAL(params),REAL(epsilon)[0],
                             INTEGER(maxIter)[0],REAL(logL)[0],INTEGER(n0)[0]);
  case 9:
    return adapt_sum_precomp(weird_series,REAL(params),REAL(epsilon)[0],
                             INTEGER(maxIter)[0],REAL(logL)[0],INTEGER(n0)[0]);
  default:
    error("No implemented logFunction found.");
  }

}
