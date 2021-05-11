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
  R_xlen_t maxIter = INTEGER(maxIter_)[0], n0 = INTEGER(n0_)[0], n = 0, nMax;
  long double log1mL = Rf_logspace_sub(0, REAL(logL_)[0]), maxA,
    logFunVal[maxIter + 1], lEps = log(REAL(epsilon)[0]) + log(2), total = 0.,
    totalBack = 0., c = 0., cb = 0.;
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

  // If too many iterations. Last iter is max.
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

  do
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun,rho);
  } while ( (log1mL ? // if L = 0 there's a simpler convergence check
      (delta(logz(logFunVal[n - 1], logFunVal[n]),
                  logFunVal[n], log1mL) >= lEps) :
               (logFunVal[n - 1] - logFunVal[n] < Rf_log1pexp(logFunVal[n] -
                 lEps))) &
                (n < maxIter));

  // Braden bounds
  KahanSum(&totalBack, expl(logFunVal[n] - log1mL - l2 - maxA), &cb);
  KahanSum(&totalBack, expl(logz(logFunVal[n - 1], logFunVal[n]) - l2 -
    maxA), &cb);
  partial_logSumExp(&logFunVal[nMax], n - nMax - 1, maxA, &cb, 1, &totalBack);

  return retFun((double)(maxA + log1pl(total + totalBack)), n);
}
