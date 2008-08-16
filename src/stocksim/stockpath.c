/* 
 * File:    stockpaths.c
 *
 * Abstract:
 *    This file permits to compute the paths of a random walk which models
 *    a stock price. The approximation is done using a monte-carlo simulation.
 *    The stock process is a diffusion process where the diffusion term is 
 *    the volatility. The volatility is also a stochastic process.
 *    The processes are defined as:
 *    dS_t = \mu S_t dt + \sigma_t S_t dW^{(1)}_t
 *    and
 *    d\sigma_t = -(\sigma_t - \xi_t)\sigma_t dt + p\sigma_t dW^{(2)}_t
 *    where \xi is:
 *    \xi_t = 1/alpha (\sigma_t - \xi_t) d\xi_t
 *    
 *    The pdes can be approximated with three different numerical methods.
 *    These are the Euler method, Milstein method and the Runge-Kutta
 *    method.
 *
 *
 *
 */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <math.h>

#include "stockpath.h"

#include "mex.h"
#include "matrix.h"
 
static double euler_stock(const double stock_t, const double vol_t, const
		double mu, const double Dt, const double phi_stock);
static double euler_vol(const double vol_t, const double xi_t, const double p,
		const double Dt, const double phi_vol);

static double milstein_stock(const double stock_t, const double vol_t, const
		double mu, const double Dt, const double phi_stock);
static double milstein_vol(const double vol_t, const double xi_t, const double
		p, const double Dt, const double phi_vol);

static double rk_stock(const double stock_t, const double vol_t, const double
		mu, const double Dt, const double phi_stock);
static double rk_vol(const double vol_t, const double xi_t, const double p,
		const double Dt, const double phi_vol);

/* The worker routine. */
void 
stockPath(gsl_rng *rng_stock, gsl_rng *rng_vol, mwSize samples, double Dt,
      double sigma_0, double S_0, double xi_0, double mu, double p, double
      alpha, mwSize N, double r, char num_method, double
      *stock_paths, double *vol_paths, double *xi_paths) 
{
   mwSize  i;
   mwSize  j;
   mwSize  index;
   double  stock_t = S_0;
   double  vol_t = sigma_0;
   double  vol_t1;
   double  xi_t = xi_0;
   double  phi_stock;
   double  phi_vol;
   double  k_1, k_2, k_3, k_4;
   double  (*num_method_stock)(const double, const double,
         const double, const double, const double) = NULL;
   double              (*num_method_vol)(const double,  const double,
         const double, const double, const double) = NULL;

   /* Determine which functions to use. */
   if (num_method == EULER) {
      num_method_stock = &euler_stock;
      num_method_vol = &euler_vol;
   } else if (num_method == MILSTEIN) {
      num_method_stock = &milstein_stock;
      num_method_vol = &milstein_vol;
   } else if (num_method == RK) {
      num_method_stock = &rk_stock;
      num_method_vol = &rk_vol;
   } else {
      mexErrMsgTxt("Only 'euler', 'milstein' and 'rk' are implemented");
   }

   for (i = 0; i < samples; i++) {
      /* Set the initial value of the paths. */
      stock_paths[index] = S_0;
      vol_paths[index] = sigma_0;
      xi_paths[index] = xi_0;

      /* Reset the current values of the stock, volatility and xi to the start
       * of the paths. Local variables are used to increase speed. */
      stock_t = S_0;
      vol_t = sigma_0;
      xi_t = xi_0;
      /* Generate the remaining paths. */
      for (j = 1; j < N; j++) {
	      phi_stock = gsl_ran_gaussian(rng_stock, 1) * sqrt(Dt);
	      phi_vol = gsl_ran_gaussian(rng_vol, 1) * sqrt(Dt);

         /* Check if the interest rate should be included. */
         if (!(r < 0.0))
            mu = r;

         /* Compute the integral of the stock using the numerical method. */
         stock_t = num_method_stock(stock_t, vol_t, mu, Dt, phi_stock);
         /* Temporarily store the current value of the volatility to use
          * in the approximation of xi. */
         vol_t1 = vol_t;
         /* Compute the integral of the volatility using the numerical 
          * method. */
         vol_t = num_method_vol(vol_t, xi_t, p, Dt, phi_vol);

         /* Compute the approximation of the xi using Runge-Kutta fourth order
          * method. */
         k_1 = 1 / alpha * (vol_t1 - xi_t);
         k_2 = 1 / alpha * (vol_t1 + 0.5 * Dt * k_1 - xi_t);
         k_3 = 1 / alpha * (vol_t1 + 0.5 * Dt * k_2 - xi_t);
         k_4 = 1 / alpha * (vol_t1 + Dt * k_3 - xi_t);

         xi_t += Dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
         
         if (stock_t < 0)
            mexErrMsgTxt("Stock price cannot be negative!");
         if (vol_t < 0)
            mexErrMsgTxt("Vol cannot be negative!");
         if (xi_t < 0)
            mexErrMsgTxt("Xi cannot be negative!");

         /* Compute the mean using a recurring relation. */
         stock_paths[index] =  stock_t;
         vol_paths[index] =  vol_t;
         xi_paths[index] =  xi_t;
      }
   }

	return;
}

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

