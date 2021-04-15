#include <Rinternals.h>
#include "adapt_sum.h"
#include "mathFun.h"
#include "log_functions.h"

SEXP doubling_sum_precomp(long double logFun(R_xlen_t k, double *Theta),
                       double *params, double eps, R_xlen_t N_start, R_xlen_t c,
                       R_xlen_t maxIter, R_xlen_t n0)
{
  // Declaration
  R_xlen_t n = 0, N, N_inc = N_start * c;
  long double maxA, lEps = log(eps), logFunVal[maxIter + 1],
            partial, *checkStart = logFunVal, S;

  logFunVal[n] = logFun(n0, params);
  while (!R_FINITE(logFunVal[n]) && n < (maxIter - 1))
    logFunVal[++n] = logFun(++n0, params);

  // Find the maximum
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
  lEps -= maxA;
  while (n < N_inc && n < maxIter) // Complete the two first checkpoints.
    logFunVal[++n] = logFun(++n0, params);
  N = n == N_inc ? N_start : ((n - N_inc) / N_inc) * N_inc; // Second to last completed checkpoint.
  partial = partial_logSumExp(logFunVal, N, maxA);
  checkStart += N + 1;
  N = n == N_inc ? N_inc - N_start - 1 : N_inc - 1; // Last completed checkpoint.
  S = partial_logSumExp(checkStart, N, maxA);

  // Calculate the tail. Only loop once.
  while (log(S) >= lEps && n < maxIter)
  {
    partial += S;
    for (N = 0, S = 0; N < N_inc; N++)
      S += exp(logFun(++n0, params) - maxA);
    n += N_inc;
  }
  partial += S;

  return retFun(maxA + log(partial), n);
}

SEXP doubling_sum_callPrecomp(SEXP lF, SEXP params, SEXP epsilon, SEXP N_start,
                              SEXP c, SEXP maxIter, SEXP n0)
{
  unsigned int funSelect = INTEGER(lF)[0];

  switch (funSelect)
  {
    case 1:
      return doubling_sum_precomp(negbin_marginal, REAL(params), REAL(epsilon)[0],
                INTEGER(N_start)[0], INTEGER(c)[0], INTEGER(maxIter)[0],
                INTEGER(n0)[0]);
    case 2:
      return doubling_sum_precomp(noObs, REAL(params), REAL(epsilon)[0],
                INTEGER(N_start)[0], INTEGER(c)[0], INTEGER(maxIter)[0],
                INTEGER(n0)[0]);
    case 3:
      return doubling_sum_precomp(COMP, REAL(params), REAL(epsilon)[0],
                INTEGER(N_start)[0], INTEGER(c)[0], INTEGER(maxIter)[0],
                INTEGER(n0)[0]);
    case 4:
      return doubling_sum_precomp(dR0, REAL(params), REAL(epsilon)[0],
                INTEGER(N_start)[0], INTEGER(c)[0], INTEGER(maxIter)[0],
                INTEGER(n0)[0]);
    case 5:
      return doubling_sum_precomp(powerLawDiff, REAL(params), REAL(epsilon)[0],
                INTEGER(N_start)[0], INTEGER(c)[0], INTEGER(maxIter)[0],
                INTEGER(n0)[0]);
    case 6:
      return doubling_sum_precomp(negbin_sentinel, REAL(params), REAL(epsilon)[0],
                INTEGER(N_start)[0], INTEGER(c)[0], INTEGER(maxIter)[0],
                INTEGER(n0)[0]);
    case 7:
      return doubling_sum_precomp(poisson_sentinel, REAL(params), REAL(epsilon)[0],
                INTEGER(N_start)[0], INTEGER(c)[0], INTEGER(maxIter)[0],
                INTEGER(n0)[0]);
    case 8:
      return doubling_sum_precomp(weird_series_constL, REAL(params), REAL(epsilon)[0],
                INTEGER(N_start)[0], INTEGER(c)[0], INTEGER(maxIter)[0],
                INTEGER(n0)[0]);
    case 9:
      return doubling_sum_precomp(weird_series, REAL(params), REAL(epsilon)[0],
                INTEGER(N_start)[0], INTEGER(c)[0], INTEGER(maxIter)[0],
                INTEGER(n0)[0]);
    default:
      error("No implemented logFunction found.");
  }

}
