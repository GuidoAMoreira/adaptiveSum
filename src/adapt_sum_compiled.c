#include <Rinternals.h>
#include "mathFun.h"
#include "adapt_sum.h"
#include "log_functions.h"

SEXP adapt_sum_precomp(long double logFun(R_xlen_t k, double *Theta),
                       double *params, double eps,
                       R_xlen_t maxIter, double logL, R_xlen_t n0)
{
  // Declaration
  R_xlen_t n = 0, nMax;
  long double maxA, logFunVal[maxIter + 1], lEps = log(eps) + l2, total = 0.,
    totalBack = 0., log1mL = Rf_logspace_sub(0, logL), c = 0., cb = 0.;

  // Find the maximum
  // Finding function max. Only check convergence after max is reached
  logFunVal[n] = logFun(n0, params);
  while (!R_FINITE(logFunVal[n])) // In case the series starts with inf values.
    logFunVal[++n] = logFun(++n0, params);

  do
    logFunVal[++n] = logFun(++n0, params);
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
    logFunVal[++n] = logFun(++n0, params);
  while  ( (log1mL ? // if L = 0 there's a simpler convergence check
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

SEXP adapt_sum_callPrecomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                           SEXP logL, SEXP n0)
{
  unsigned int funSelect = INTEGER(lF)[0];

  switch (funSelect)
  {
  case 1:
    return adapt_sum_precomp(negbin_marginal, REAL(params), REAL(epsilon)[0],
                             INTEGER(maxIter)[0], REAL(logL)[0], INTEGER(n0)[0]);
  case 2:
    return adapt_sum_precomp(noObs,REAL(params), REAL(epsilon)[0],
                             INTEGER(maxIter)[0], REAL(logL)[0], INTEGER(n0)[0]);
  case 3:
    return adapt_sum_precomp(COMP,REAL(params), REAL(epsilon)[0],
                             INTEGER(maxIter)[0], REAL(logL)[0], INTEGER(n0)[0]);
  case 4:
    return adapt_sum_precomp(dR0,REAL(params), REAL(epsilon)[0],
                             INTEGER(maxIter)[0], REAL(logL)[0], INTEGER(n0)[0]);
  case 5:
    return adapt_sum_precomp(powerLawDiff,REAL(params),REAL(epsilon)[0],
                             INTEGER(maxIter)[0],REAL(logL)[0],INTEGER(n0)[0]);
  case 6:
    return adapt_sum_precomp(negbin_sentinel, REAL(params), REAL(epsilon)[0],
                             INTEGER(maxIter)[0], REAL(logL)[0], INTEGER(n0)[0]);
  case 7:
    return adapt_sum_precomp(poisson_sentinel, REAL(params), REAL(epsilon)[0],
                             INTEGER(maxIter)[0], REAL(logL)[0], INTEGER(n0)[0]);
  case 8:
    return adapt_sum_precomp(weird_series_constL, REAL(params), REAL(epsilon)[0],
                             INTEGER(maxIter)[0], REAL(logL)[0], INTEGER(n0)[0]);
  case 9:
    return adapt_sum_precomp(weird_series, REAL(params), REAL(epsilon)[0],
                             INTEGER(maxIter)[0], REAL(logL)[0], INTEGER(n0)[0]);
  default:
    error("No implemented logFunction found.");
  }

}
