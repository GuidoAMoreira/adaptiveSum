#include "log_functions.h"

// To add functions implementations, they must be written here with the
// same parameters as the existent ones. Additionally, they must be assigned
// a "number" in the adapt_sum_callPrecomp function in the adapt_sum_compiled.c
// file switch statement. Finally, a condition must be included in the adapt_sum
// wrapper function in R in the adapt_sum.R file.

double negbin_marginal(R_xlen_t k, double *Theta)
{
  double x = Theta[3], s = Theta[1];
  return (k<x ? -INFINITY :
            Rf_dnbinom(k, s, s / (s + Theta[0]), 1) +
              Rf_dbinom(x, k, Theta[2], 1));
}

double noObs(R_xlen_t k, double *Theta)
{return k * log1p(-Theta[0]);}

double COMP(R_xlen_t k, double *Theta)
{return k * log(Theta[0]) - Theta[1] * lgamma(k+1);}

double dR0(R_xlen_t k, double *Theta)
{
  double x = Theta[2];
  if (k == 0 || k < x) {return -INFINITY;}
  else
  {
    double wy = Theta[1] * k, wyym1 = wy + k - 1;
    return lgamma(wyym1) - (lgamma(wy) + lgamma(k+1)) + ((k-1) * (log(Theta[0])-
                  log(Theta[1])) - wyym1 * log1p(Theta[0]/Theta[1])) +
                  Rf_dbinom(x, k, Theta[3], 1);
  }
}

double powerLawDiff(R_xlen_t k, double *Theta)
{return (k<Theta[1] ? -INFINITY :
           -Theta[0] * log(k) + log_diff_exp(0,
                           -Theta[2] - Theta[3] * (k - Theta[1])));}

                           double negbin_sentinel(R_xlen_t k, double *Theta)
                           {
                             return Rf_dnbinom(k, Theta[1], Theta[1]/ (Theta[1] + Theta[0]), 1) +
                               k * log1p(-Theta[2]);
                           }

double poisson_sentinel(R_xlen_t k, double *Theta)
{
  return Rf_dpois(k, Theta[0], 1) +
    k * log1p(-Theta[1]);
}

double weird_series_constL(R_xlen_t k, double *Theta)
{return (k==0 ? -INFINITY : -(2*log(k) + k * log(Theta[0])));}

double weird_series(R_xlen_t k, double *Theta)
{return (k==0 ? -INFINITY : lgamma(k+1) - k * log(k));}


