#include <Rinternals.h>
#include "adapt_sum.h"
#include "mathFun.h"

SEXP naive_sum(SEXP logFun, SEXP params, SEXP epsilon, SEXP maxIter_,
               SEXP n0_, SEXP rho)
{
  // Error checking
  if(!isReal(params)) error("'params' must be a vector");
  if(!isEnvironment(rho)) error("'rho' should be an environment");

  // Setting up
  R_xlen_t maxIter = INTEGER(maxIter_)[0], n0 = INTEGER(n0_)[0], n = 0, nMax;
  long double logFunVal[maxIter + 1], lEps = log(REAL(epsilon)[0]), maxA,
    total = 0., totalBack = 0., c = 0., cb = 0.;
  defineVar(install("Theta"), params, rho);

  // Finding function max. Only check convergence after max is reached
  defineVar(install("k"), Rf_ScalarInteger(n0), rho);
  logFunVal[n] = feval(logFun,rho);
  while (!R_FINITE(logFunVal[n]))
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun,rho);
  }

  do
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun,rho);
  }
  while (logFunVal[n] >= logFunVal[n - 1] && n <= (maxIter - 1));

  // If too many iterations
  if (n == maxIter)
  {
    partial_logSumExp(logFunVal, maxIter - 1, logFunVal[n], &c, 0, &total);
    return retFun(logFunVal[n] + log1p(total), maxIter);
  }

  // I know which is the max due to the stop criteria.
  // Assumed local max = global max.
  maxA = logFunVal[n - 1];
  nMax = n;
  if (n > 1)
    partial_logSumExp(logFunVal, n - 2, maxA, &c, 0, &total);

  // Now for the convergence checking.
  do
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun,rho);
  } while ((logFunVal[n] >= lEps) & (n < maxIter));
  partial_logSumExp(&logFunVal[nMax], n - nMax, maxA, &cb, 1, &totalBack);

  return retFun((double)(maxA) + log1pl(total + totalBack), n);
}
