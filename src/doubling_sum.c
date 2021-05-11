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
    maxA, partial = 0., *checkStart = logFunVal, S = 0., cc = 0.,
    ccc = 0., total = 0.;
  defineVar(install("Theta"), params, rho);

  defineVar(install("k"), Rf_ScalarInteger(n0), rho);
  logFunVal[n] = feval(logFun,rho);
  while (!R_FINITE(logFunVal[n]) && n < (maxIter - 1))
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun, rho);
  }

  // Find the maximum
  do
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun, rho);
  }
  while (logFunVal[n] >= logFunVal[n - 1] && n <= (maxIter - 1));

  // If too many iterations. Last iter is max.
  if (n == maxIter)
  {
    partial_logSumExp(logFunVal, maxIter - 1, logFunVal[n], &cc, 0, &total);
    return retFun(logFunVal[n] + log1p(total), maxIter);
  }

  // I know which is the max due to the stop criteria.
  // Assumed local max = global max.
  // 20 added to make calculations with more precision.
  maxA = logFunVal[n - 1];
  lEps -= maxA; // For the convergence checking
  while (n < N_inc && n < maxIter) // Complete the two first checkpoints.
  {
    defineVar(install("k"), Rf_ScalarInteger(++n0), rho);
    logFunVal[++n] = feval(logFun, rho);
  }
  N = n == N_inc ? N_start : ((n - N_inc) / N_inc) * N_inc; // Second to last completed checkpoint.
  partial_logSumExp(logFunVal, N, maxA, &cc, 0, &partial);
  checkStart += N + 1;
  N = n == N_inc ? N_inc - N_start - 1 : N_inc - 1; // Last completed checkpoint.
  cc = 0.;
  partial_logSumExp(checkStart, N, maxA, &cc, 0, &S);

  // Calculate the tail. Only loop once.
  while (log(S) >= lEps && n < maxIter)
  {
    KahanSum(&partial, S, &ccc);
    n0 += N_inc; // Going backwards for numerical stability
    cc = 0.;
    for (N = 0, S = 0.; N < N_inc; N++)
    {
      defineVar(install("k"), Rf_ScalarInteger(n0--), rho);
      KahanSum(&S, exp(feval(logFun, rho) - maxA), &cc);
    }
    n0 += N_inc;
    n += N_inc;
  }
  KahanSum(&partial, S, &ccc);

  return retFun((double)(maxA + logl(partial)), n);
}
