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
stockPath(double rng_state, long int samples, double Dt, double sigma_0, double
      S_0, double xi_0, double mu, double p, double alpha, long int N,
      char num_method, double *stock_mn, double *vol_mn, double *xi_mn) 
{
	const gsl_rng_type *rng_type;
	gsl_rng				 *rng_stock, *rng_vol;
   long int            i;
   long int            j;
   double              stock_t = S_0;
   double              vol_t = sigma_0;
   double              stock_t_AV = S_0;
   double              vol_t_AV = sigma_0;
   double              vol_t1;
   double              xi_t = xi_0;
   double              xi_t_AV = xi_0;
   double              phi_stock;
   double              phi_vol;
   double              k_1, k_2, k_3, k_4;
   double              (*num_method_stock)(const double, const double,
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

   /* The Mersenne Twister generator should be good enough for our purpose. */
	rng_type = gsl_rng_mt19937;

   /* Create two random number generators. */
	rng_stock = gsl_rng_alloc(rng_type);
	rng_vol = gsl_rng_alloc(rng_type);
 
   /* Initialize the state of the random number generator. */
   gsl_rng_set(rng_stock, rng_state);
   gsl_rng_set(rng_vol, rng_state + 1);
  
   /* Set the initial values of the path. */
   stock_mn[0] = stock_t;
   vol_mn[0] = vol_t;
   xi_mn[0] = xi_t;

   for (i = 0; i < samples; i++) {
      /* Reset the current values of the stock, volatility and xi to the start
       * of the path. */
      stock_t = stock_t_AV = S_0;
      vol_t = vol_t_AV = sigma_0;
      xi_t = xi_t_AV = xi_0;
      /* Generate the remaining path. */
      for (j = 1; j < N; j++) {
	      phi_stock = gsl_ran_gaussian(rng_stock, 1) * sqrt(Dt);
	      phi_vol = gsl_ran_gaussian(rng_vol, 1) * sqrt(Dt);

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

         /* Compute the integral of the stock using the numerical method. */
         stock_t_AV = num_method_stock(stock_t_AV, vol_t_AV, mu, Dt, -phi_stock);
         /* Temporarily store the current value of the volatility to use
          * in the approximation of xi. */
         vol_t1 = vol_t_AV;
         /* Compute the integral of the volatility using the numerical 
          * method. */
         vol_t_AV = num_method_vol(vol_t_AV, xi_t_AV, p, Dt, -phi_vol);

         /* Compute the approximation of the xi using Runge-Kutta fourth order
          * method. */
         k_1 = 1 / alpha * (vol_t1 - xi_t_AV);
         k_2 = 1 / alpha * (vol_t1 + 0.5 * Dt * k_1 - xi_t_AV);
         k_3 = 1 / alpha * (vol_t1 + 0.5 * Dt * k_2 - xi_t_AV);
         k_4 = 1 / alpha * (vol_t1 + Dt * k_3 - xi_t_AV);

         xi_t_AV += Dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
         
         /* Compute the mean using a recurring relation. */
         stock_mn[j] +=  (0.5 * (stock_t + stock_t_AV) - stock_mn[j]) / (i + 1);
         vol_mn[j] +=  (0.5 * (vol_t + vol_t_AV) - vol_mn[j]) / (i + 1);
         xi_mn[j] +=  (0.5 * (xi_t + xi_t_AV) - xi_mn[j]) / (i + 1);
      }
   }
	gsl_rng_free(rng_stock);
	gsl_rng_free(rng_vol);

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

