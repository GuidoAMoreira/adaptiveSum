#include <Rinternals.h>
#include "mathFun.h"
#include "adapt_sum.h"
#include "log_functions.h"

SEXP adapt_sum_precomp(long double logFun(R_xlen_t k, double *Theta),
                       double *params, double eps,
                       R_xlen_t maxIter, double logL, R_xlen_t n0)
{
  // Declaration
  R_xlen_t n = 0;
  long double maxA, logFunVal[maxIter + 1], lEps = log(eps) + log(2), total,
    log1mL = Rf_logspace_sub(0, logL);

  // Find the maximum
  // Finding function max. Only check convergence after max is reached
  logFunVal[n] = logFun(n0, params);
  while (!R_FINITE(logFunVal[n]))
    logFunVal[++n] = logFun(++n0, params);

  do
    logFunVal[++n] = logFun(++n0, params);
  while (logFunVal[n] >= logFunVal[n - 1] && n <= (maxIter - 1));

  // If too many iterations. Last iter is max.
  if (n == maxIter)
    return retFun(logFunVal[n] +
                  log1p(partial_logSumExp(logFunVal, maxIter - 1,
                                          logFunVal[n])),
                                          maxIter);

  // I know which is the max due to the stop criteria.
  // Assumed local max = global max.
  // 20 added to make calculations with more precision.
  maxA = logFunVal[n - 1] + 20;
  total = partial_logSumExp(logFunVal, n - 1, maxA);

  // Calculate the tail. Only loop once.
  do
  {
    total += exp(logFunVal[n] - maxA);
    logFunVal[++n] = logFun(++n0, params);
  }
  while ((delta(logz(logFunVal[n - 1], logFunVal[n]),
                  logFunVal[n], log1mL) >= lEps) & (n <= (maxIter - 1)));

  // Braden bounds
  total += exp(logz(logFunVal[n - 1], logFunVal[n]) - log(2) - maxA) +
    exp(logFunVal[n] - log1mL - log(2) - maxA);

  return retFun(maxA + log(total), n);
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
