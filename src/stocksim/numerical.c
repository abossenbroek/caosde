#include <gsl/gsl_math.h>

#include <math.h>

#include  "numerical.h"


double 
euler_stock(const double stock_t, const double vol_t, const double mu,
      const double Dt, const double phi_stock)
{
   return stock_t + mu * stock_t * Dt + vol_t * stock_t * phi_stock;
}

double 
euler_vol(const double vol_t, const double xi_t, const double p,
      const double Dt, const double phi_vol)
{
   return vol_t - (vol_t - xi_t) * Dt + p * vol_t * phi_vol;
}

double 
milstein_stock(const double stock_t, const double vol_t, const double mu,
      const double Dt, const double phi_stock)
{
   return stock_t + mu * stock_t * Dt + vol_t * stock_t * phi_stock
      + 0.5 * gsl_pow_2(vol_t) * stock_t * (gsl_pow_2(phi_stock) - Dt);
}

double 
milstein_vol(const double vol_t, const double xi_t, const double p,
      const double Dt, const double phi_vol)
{
   return vol_t - (vol_t - xi_t) * Dt + p * vol_t * phi_vol 
      + 0.5 * gsl_pow_2(p) * vol_t * (gsl_pow_2(phi_vol) - Dt);

}

double 
rk_stock(const double stock_t, const double vol_t, const double mu,
      const double Dt, const double phi_stock)
{
   double stock_hat = stock_t + mu * stock_t * Dt + vol_t * stock_t * sqrt(Dt);
   
   return stock_t + mu * stock_t * Dt + vol_t * stock_t * phi_stock
      + 1 / (2 * sqrt(Dt)) * (gsl_pow_2(phi_stock) - Dt) * (vol_t * stock_hat 
            - vol_t * stock_t);
}

double 
rk_vol(const double vol_t, const double xi_t, const double p,
      const double Dt, const double phi_vol)
{
   double vol_hat = vol_t - (vol_t - xi_t) * Dt + p * vol_t * sqrt(Dt);

   return vol_t - (vol_t - xi_t) * Dt + p * vol_t * phi_vol 
      + 1 / (2 * sqrt(Dt)) * (gsl_pow_2(phi_vol) - Dt) * (p * vol_hat 
            - p * vol_t);
}


/* vim: set et : tw=80 : spell spelllang=en: */

