#include <Rinternals.h>
#include "adapt_sum.h"

SEXP naive_sum(SEXP logFun, SEXP params, SEXP epsilon, SEXP maxIter_,
               SEXP n0_, SEXP rho)
{
  // Error checking
  if(!isReal(params)) error("'params' must be a vector");
  if(!isEnvironment(rho)) error("'rho' should be an environment");

  // Setting up
  R_xlen_t maxIter = INTEGER(maxIter_)[0], n0 = INTEGER(n0_)[0], n = 0;
  long double logFunVal[maxIter + 1], lEps = log(REAL(epsilon)[0]), maxA,
    total;
  defineVar(install("Theta"), params, rho);

  // Finding function max. Only check convergence after max is reached
  defineVar(install("k"), Rf_ScalarInteger(n0), rho);
  logFunVal[n] = feval(logFun,rho);
  maxA = logFunVal[n];
  do
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun,rho);
  } while (!R_FINITE(logFunVal[n]) ||
    (logFunVal[n] >= logFunVal[n - 1] &&
    n <= (maxIter - 1)));

  // If too many iterations
  if (n == maxIter)
    return retFun(logFunVal[n - 1] +
                  log1p(partial_logSumExp(logFunVal, maxIter - 1,
                                          logFunVal[n - 1])),
                                          maxIter);

  // I know which is the max due to the stop criteria.
  // Assumed local max = global max.
  maxA = logFunVal[n - 1];
  total = partial_logSumExp(logFunVal, n - 2, maxA);
  total += exp(logFunVal[n] - maxA);

  // Now for the convergence checking. Only loop once.
  do
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun,rho);
    total += exp(logFunVal[n] - maxA);
  } while ((logFunVal[n] >= lEps) & (n < maxIter));

  return retFun(maxA + log1p(total), n);
}
