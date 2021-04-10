#include <Rinternals.h>
#include "adapt_sum.h"
#include "mathFun.h"

SEXP adapt_sum(SEXP logFun, SEXP params, SEXP epsilon, SEXP maxIter_,
               SEXP logL_, SEXP n0_, SEXP rho)
{
  // Error checking
  if(!isReal(params)) error("'params' must be a vector");
  if(!isEnvironment(rho)) error("'rho' should be an environment");

  // Setting up
  R_xlen_t maxIter = INTEGER(maxIter_)[0], n0 = INTEGER(n0_)[0], n = 0;
  double log1mL = log_diff_exp(0, REAL(logL_)[0]);
  long double maxA, logFunVal[maxIter + 1],
                             lEps = log(REAL(epsilon)[0]) + log(2), total;
  defineVar(install("Theta"), params, rho);

  defineVar(install("k"), Rf_ScalarInteger(n0), rho);
  logFunVal[n] = feval(logFun,rho);
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
  // Remove the max to add later.
  maxA = logFunVal[n - 1];
  total = partial_logSumExp(logFunVal, n - 2, maxA);

  do
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun,rho);
    total += exp(logFunVal[n - 1] - maxA);
  } while ((delta(logz(logFunVal[n - 1], logFunVal[n]),
                 logFunVal[n], log1mL) >= lEps) & (n < maxIter));

  // Braden bounds
  total += exp(logz(logFunVal[n - 1], logFunVal[n]) - log(2) - maxA) +
    exp(logFunVal[n] - log1mL - log(2) - maxA);

  return retFun(maxA + log1p(total), n);
}
