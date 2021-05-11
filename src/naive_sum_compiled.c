#include <Rinternals.h>
#include "adapt_sum.h"
#include "mathFun.h"
#include "log_functions.h"

SEXP naive_sum_precomp(long double logFun(R_xlen_t k, double *Theta),
                       double *params, double eps,
                       R_xlen_t maxIter, R_xlen_t n0)
{
  // Declaration
  R_xlen_t n = 0, nMax;
  long double maxA, lEps = log(eps), logFunVal[maxIter + 1], total = 0.,
    totalBack = 0., c = 0., cb = 0.;

  // Finding function max. Only check convergence after max is reached
  logFunVal[n] = logFun(n0, params);
  while (!R_FINITE(logFunVal[n]))
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

  // Calculate the tail. Only loop once.
  do
    logFunVal[++n] = logFun(++n0, params);
  while ((logFunVal[n] >= lEps) & (n <= (maxIter - 1)));
  partial_logSumExp(&logFunVal[nMax], n - nMax, maxA, &cb, 1, &totalBack);

  return retFun((double)(maxA) + log1pl(total + totalBack), n);
}

SEXP naive_sum_callPrecomp(SEXP lF, SEXP params, SEXP epsilon, SEXP maxIter,
                           SEXP n0)
{
  unsigned int funSelect = INTEGER(lF)[0];

  switch (funSelect)
  {
    case 1:
      return naive_sum_precomp(negbin_marginal,REAL(params),REAL(epsilon)[0],
                               INTEGER(maxIter)[0],INTEGER(n0)[0]);
    case 2:
      return naive_sum_precomp(noObs,REAL(params),REAL(epsilon)[0],
                               INTEGER(maxIter)[0],INTEGER(n0)[0]);
    case 3:
      return naive_sum_precomp(COMP,REAL(params),REAL(epsilon)[0],
                               INTEGER(maxIter)[0],INTEGER(n0)[0]);
    case 4:
      return naive_sum_precomp(dR0,REAL(params),REAL(epsilon)[0],
                                 INTEGER(maxIter)[0],INTEGER(n0)[0]);
    case 5:
      return naive_sum_precomp(powerLawDiff,REAL(params),REAL(epsilon)[0],
                               INTEGER(maxIter)[0],INTEGER(n0)[0]);
    case 6:
      return naive_sum_precomp(negbin_sentinel,REAL(params),REAL(epsilon)[0],
                               INTEGER(maxIter)[0],INTEGER(n0)[0]);
    case 7:
      return naive_sum_precomp(poisson_sentinel,REAL(params),REAL(epsilon)[0],
                               INTEGER(maxIter)[0],INTEGER(n0)[0]);
    case 8:
      return naive_sum_precomp(weird_series_constL,REAL(params),REAL(epsilon)[0],
                               INTEGER(maxIter)[0],INTEGER(n0)[0]);
    case 9:
      return naive_sum_precomp(weird_series,REAL(params),REAL(epsilon)[0],
                               INTEGER(maxIter)[0],INTEGER(n0)[0]);
    default:
      error("No implemented logFunction found.");
  }

}
