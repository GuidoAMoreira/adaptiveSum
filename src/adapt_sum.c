#include "adapt_sum.h"

SEXP adapt_sum(SEXP logFun, SEXP params, SEXP epsilon, SEXP maxIter_,
               SEXP logL_, SEXP n0_, SEXP rho)
{
  // Error checking
  if(!isReal(params)) error("'params' must be a vector");
  if(!isEnvironment(rho)) error("'rho' should be an environment");

  // Setting up
  R_xlen_t maxIter = INTEGER(maxIter_)[0], n0 = INTEGER(n0_)[0], n = 0;
  SEXP result = PROTECT(allocVector(REALSXP, maxIter+1)),
    i_SEXP = PROTECT(allocVector(INTSXP,1));
  double *logFunVal = REAL(result);
  double lEps = log(REAL(epsilon)[0]) + log(2), logL = REAL(logL_)[0],
    logR = logit(logL);
  defineVar(install("Theta"), params, rho);

  INTEGER(i_SEXP)[0] = n0;
  defineVar(install("k"), i_SEXP, rho);
  logFunVal[n] = feval(logFun,rho);
  do
  {
    INTEGER(i_SEXP)[0] = ++n0;
    defineVar(install("k"), i_SEXP, rho);
    logFunVal[++n] = feval(logFun,rho);
  } while (!R_FINITE(logFunVal[n]) ||
             (logFunVal[n] - logFunVal[n - 1] >= 0 &&
             n <= (maxIter - 1)));

  // If too many iterations
  if (n == maxIter)
  {
    logFunVal[maxIter] = 1;
    UNPROTECT(2);
    return retFun(result, maxIter_);
  }

  do
  {
    INTEGER(i_SEXP)[0] = ++n0;
    defineVar(install("k"), i_SEXP, rho);
    logFunVal[++n] = feval(logFun,rho);
  } while ((delta(logz(logFunVal[n-1], logFunVal[n]),
                 logFunVal[n-1], logR) > lEps) & (n <= (maxIter - 1)));

  // If too many iterations
  if (n == maxIter)
  {
    logFunVal[maxIter] = 1;
    UNPROTECT(2);
    return retFun(result, maxIter_);
  }

  logFunVal[n] = logz(logFunVal[n-1], logFunVal[n]) - log(2);
  logFunVal[n+1] = logFunVal[n-1] + logR - log(2);

  INTEGER(i_SEXP)[0] = n;
  UNPROTECT(2);
  return retFun(result,i_SEXP);
}
