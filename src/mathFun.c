#include <Rinternals.h>
#include "mathFun.h"

// Auxiliary functions roughly copied from the stan header files

void partial_logSumExp(long double* fun, R_xlen_t evals, long double maxA,
                       long double* c, int backwards, long double* res)
{
  if (backwards)
    for (R_xlen_t i = evals; i >= 0; i--)
      KahanSum(res, expl(fun[i] - maxA), c);
  else
    for (R_xlen_t i = 0; i <= evals; i++)
      KahanSum(res, expl(fun[i] - maxA), c);
}
