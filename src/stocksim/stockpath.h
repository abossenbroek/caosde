#ifndef _STOCKPATH_H_
#define _STOCKPATH_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

enum { 
   EULER,
   MILSTEIN,
   RK,
   NO_METHOD
};

enum {
   NO_ANTIVAR,
   ANTIVAR
};

void stockPath(gsl_rng *rng_stock, gsl_rng *rng_vol, long int samples, double Dt,
      double sigma_0, double S_0, double xi_0, double mu, double p, double
      alpha, long int N, double r, char num_method, char antivar, double
      *stock_mn, double *vol_mn, double *xi_mn);

#endif

/* vim: set et : tw=80 : spell spelllang=en: */

