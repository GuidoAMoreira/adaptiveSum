#include "adapt_sum.h"

SEXP naive_sum(SEXP logFun, SEXP params, SEXP epsilon, SEXP maxIter_,
               SEXP n0_, SEXP rho)
{
  // Error checking
  if(!isReal(params)) error("'params' must be a vector");
  if(!isEnvironment(rho)) error("'rho' should be an environment");

  // Setting up
  R_xlen_t maxIter = INTEGER(maxIter_)[0], n0 = INTEGER(n0_)[0], n = 0;
  SEXP result = PROTECT(allocVector(REALSXP, maxIter+1)),
    i_SEXP = PROTECT(allocVector(INTSXP,1));
  double *logFunVal = REAL(result);
  double lEps = log(REAL(epsilon)[0]);
  defineVar(install("Theta"), params, rho);

  // Finding function max. Only check convergence after max is reached
  INTEGER(i_SEXP)[0] = n0;
  defineVar(install("k"), i_SEXP, rho);
  logFunVal[n] = feval(logFun,rho);
  do
  {
    INTEGER(i_SEXP)[0] = ++n0;
    defineVar(install("k"), i_SEXP, rho);
    logFunVal[++n] = feval(logFun,rho);
  } while (!R_FINITE(logFunVal[n]) ||
    (logFunVal[n] >= logFunVal[n - 1] &&
    n <= (maxIter - 1)));

  // If too many iterations
  if (n == maxIter)
  {
    logFunVal[maxIter] = 1;
    UNPROTECT(2);
    return retFun(result, maxIter_);
  }

  // Now for the convergence checking
  do
  {
    INTEGER(i_SEXP)[0] = ++n0;
    defineVar(install("k"), i_SEXP, rho);
    logFunVal[++n] = feval(logFun,rho);
  } while ((logFunVal[n] >= lEps) & (n < maxIter));

  // If too many iterations
  if (n == maxIter)
  {
    logFunVal[maxIter] = 1;
    UNPROTECT(2);
    return retFun(result, maxIter_);
  }

  INTEGER(i_SEXP)[0] = n + 1;
  UNPROTECT(2);
  return retFun(result, i_SEXP);
}
