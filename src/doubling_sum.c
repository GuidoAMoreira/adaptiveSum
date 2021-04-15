#include <Rinternals.h>
#include "adapt_sum.h"
#include "mathFun.h"

SEXP doubling_sum(SEXP logFun, SEXP params, SEXP epsilon, SEXP N_start_,
                  SEXP c_, SEXP maxIter_, SEXP n0_, SEXP rho)
{
  // Error checking
  if(!isReal(params)) error("'params' must be a vector");
  if(!isEnvironment(rho)) error("'rho' should be an environment");

  // Setting up
  R_xlen_t maxIter = INTEGER(maxIter_)[0], n0 =
    INTEGER(n0_)[0], c = INTEGER(c_)[0], N_start =
    INTEGER(N_start_)[0], N, n = 0, N_inc = N_start * c;
  long double logFunVal[maxIter + 1], lEps = log(REAL(epsilon)[0]),
    maxA, partial, *checkStart = logFunVal, S;
  defineVar(install("Theta"), params, rho);

  defineVar(install("k"), Rf_ScalarInteger(n0), rho);
  logFunVal[n] = feval(logFun,rho);
  while (!R_FINITE(logFunVal[n]) && n < (maxIter - 1))
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun,rho);
  }

  // Find the maximum
  do
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun,rho);
  } while (logFunVal[n] >= logFunVal[n - 1] && n <= (maxIter - 1));

  // If too many iterations. Last iter is max.
  if (n == maxIter)
    return retFun(logFunVal[n] +
                  log1p(partial_logSumExp(logFunVal, n - 1,
                                          logFunVal[n])),
                                          maxIter);

  // I know which is the max due to the stop criteria.
  // Assumed local max = global max.
  // 20 added to make calculations with more precision.
  maxA = logFunVal[n - 1] + 20;
  lEps -= maxA;
  while (n < N_inc && n < maxIter) // Complete second bulk.
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun,rho);
  }
  N = n == N_inc ? N_start : ((n - N_inc) / N_inc) * N_inc; // Second to last completed checkpoint.
  partial = partial_logSumExp(logFunVal, N, maxA);
  checkStart += N + 1;
  N = n == N_inc ? N_inc - N_start - 1 : N_inc - 1; // Last completed checkpoint.
  S = partial_logSumExp(checkStart, N, maxA);

  while (log(S) >= lEps && n < maxIter)
  {
    partial += S;
    for (N = 0, S = 0; N < N_inc; N++)
    {
      defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
      S += exp(feval(logFun,rho) - maxA);
    }
    n += N_inc;
  }
  partial += S;

  return retFun(maxA + log(partial), n);
}
