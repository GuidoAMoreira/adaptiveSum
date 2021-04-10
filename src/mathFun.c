#include <Rinternals.h>

// Auxiliary functions roughly copied from the stan header files

long double partial_logSumExp(long double* fun, R_xlen_t til, long double maxA)
{
  long double total = 0;
  for (R_xlen_t i = 0; i <= til; i++)
    total += exp(fun[i] - maxA);

  return total;
}
