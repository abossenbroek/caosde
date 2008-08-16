#ifndef _STOCKPATH_H_
#define _STOCKPATH_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

#include "/Applications/MATLAB_R2008a/extern/include/mex.h"

enum { 
   EULER,
   MILSTEIN,
   RK,
   NO_METHOD
};

void stockPath(gsl_rng *rng_stock, gsl_rng *rng_vol, mwSize samples, double Dt,
      double sigma_0, double S_0, double xi_0, double mu, double p, double
      alpha, mwSize N, double r, char num_method,  double
      *stock_paths, double *vol_paths, double *xi_paths);

#endif

/* vim: set et : tw=80 : spell spelllang=en: */

