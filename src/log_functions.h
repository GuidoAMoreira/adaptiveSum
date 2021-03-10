#ifndef __LOG_FUNCTIONS_ADAPTIVE_SUM_H__
#define __LOG_FUNCTIONS_ADAPTIVE_SUM_H__
#include <R.h>
#include <Rinternals.h>

double Rf_dnbinom(double x, double size, double prob, int give_log);
double Rf_dbinom(double x, double n, double p, int give_log);

double negbin_marginal(R_xlen_t k, double *Theta)
{
  double x = Theta[3], s = Theta[1];
  return (k<x ? -INFINITY :
            Rf_dnbinom(k, s, s / (s + Theta[0]), 1) +
              Rf_dbinom(x, k, Theta[2], 1));
}

double noObs(R_xlen_t k, double *Theta)
{return k * log(1-Theta[0]);}

double COMP(R_xlen_t k, double *Theta)
{return k * log(Theta[0]) - Theta[1] * lgamma(k+1);}

#endif
