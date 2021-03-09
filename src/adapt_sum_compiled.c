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
  double *logFunVal = REAL(result), lEps = log(eps) + log(2);
  double logR = logit(logL);

  // Find the maximum
  logFunVal[n] = logFun(n0, params);
  do
  {logFunVal[++n] = logFun(++n0, params);}
  while (!R_FINITE(logFunVal[n]) ||
    (logFunVal[n] - logFunVal[n - 1] >= 0 &&
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
  while ((delta(logz(logFunVal[n-1], logFunVal[n]),
                  logFunVal[n-1], logR) > lEps) & (n <= (maxIter - 1)));

  // If too many iterations
  if (n == maxIter)
  {
    logFunVal[maxIter] = 1;
    INTEGER(i_SEXP)[0] = maxIter;
    UNPROTECT(2);
    return retFun(result, i_SEXP);
  }

  logFunVal[n] = logz(logFunVal[n-1], logFunVal[n]) - log(2);
  logFunVal[n+1] = logFunVal[n-1] + logR - log(2);

  INTEGER(i_SEXP)[0] = n;
  UNPROTECT(2);
  return retFun(result,i_SEXP);
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
  default:
    error("No implemented logFunction found.");
  }

}
