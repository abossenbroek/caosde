#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <string.h>

#include "mex.h"
#include "matrix.h"
   
#define  RANDSTATE_IN prhs[0]
#define  SAMPLES_IN   prhs[1]
#define  DT_IN	       prhs[2]
#define  SIGMA0_IN    prhs[3]
#define  S_0_IN	       prhs[4]
#define  XI0_IN       prhs[5]
#define  MU_IN	       prhs[6]
#define  P_IN	       prhs[7]
#define  ALPHA_IN     prhs[8]
#define  T_IN	       prhs[9]
#define  NUMMETHOD_IN prhs[10]

#define  STOCKAVG     plhs[0]
#define  VOLAVG       plhs[1]
#define  XIAVG        plhs[2]

#define  EULER        0
#define  MILSTEIN     1
#define  RK           2
#define  NO_METHOD    4

double euler_stock(const double stock_t, const double vol_t, const double mu,
      const double Dt, const double phi_stock);
double euler_vol(const double vol_t, const double xi_t, const double p,
      const double Dt, const double phi_vol);

double milstein_stock(const double stock_t, const double vol_t, const double mu,
      const double Dt, const double phi_stock);
double milstein_vol(const double vol_t, const double xi_t, const double p,
      const double Dt, const double phi_vol);

double rk_stock(const double stock_t, const double vol_t, const double mu,
      const double Dt, const double phi_stock);
double rk_vol(const double vol_t, const double xi_t, const double p,
      const double Dt, const double phi_vol);


/* The worker routine. */
static void 
stockPath(double rng_state, long int samples, double Dt, double sigma_0, double
      S_0, double xi_0, double mu, double p, double alpha, long int N,
      char num_method, double *stock_mn, double *vol_mn, double *xi_mn) 
{
	const gsl_rng_type *rng_type;
	gsl_rng				 *rng_stock, *rng_vol;
   long int            i;
   long int            j;
   double              stock_t = S_0, stock_t_AV = S_0;
   double              vol_t = sigma_0, vol_t_AV = sigma_0;
   double              vol_t1;
   double              xi_t = xi_0, xi_t_AV = xi_0;
   double              phi_stock;
   double              phi_vol;
   double              k_1, k_2, k_3, k_4;
   double              (*num_method_stock)(const double, const double,
         const double, const double, const double);
   double              (*num_method_vol)(const double,  const double,
         const double, const double, const double);

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
         
         /* Antithetic variance is used to reduce the number of computations
          * required to converge. */

         /* Compute the integral of the stock using the numerical method. */
         stock_t_AV = num_method_stock(stock_t_AV, vol_t_AV, mu, Dt,
               -phi_stock);
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

void 
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
   long int N;
   double   *stock_mn;
   double   *vol_mn;
   double   *xi_mn;
   char     num_method = NO_METHOD;

	if (nrhs != 11) {
		mexErrMsgTxt("Eleven input arguments required.");
	} else if (nlhs != 3) {
		mexErrMsgTxt("Six output arguments required."); 
   } 

   /* Determine if all the parameters are the proper type. */
   if (!mxIsDouble(RANDSTATE_IN) || !mxIsDouble(SAMPLES_IN) ||
         !mxIsDouble(DT_IN) || !mxIsDouble(SIGMA0_IN) ||
         !mxIsDouble(S_0_IN) || !mxIsDouble(XI0_IN) ||
         !mxIsDouble(MU_IN) || !mxIsDouble(P_IN) ||
         !mxIsDouble(ALPHA_IN) || !mxIsDouble(T_IN))
      mexErrMsgTxt("All parameters have to be numbers."); 

   if (!mxIsClass(NUMMETHOD_IN, "char"))
      mexErrMsgTxt("num_method must be a char");

   /* How many steps should be simulated? */
   N = lround(*mxGetPr(T_IN) / *mxGetPr(DT_IN));
   /* Assign memory to the output parameters. */
   STOCKAVG = mxCreateDoubleMatrix(1, N, mxREAL);
   VOLAVG = mxCreateDoubleMatrix(1, N, mxREAL);
   XIAVG = mxCreateDoubleMatrix(1, N, mxREAL);

   /* Use local variables to access the memory allocated above. */
   stock_mn = mxGetPr(STOCKAVG);
   vol_mn = mxGetPr(VOLAVG);
   xi_mn = mxGetPr(XIAVG);

   /* Determine which numerical method should be used. */
   if (strncmp(mxArrayToString(NUMMETHOD_IN), "euler", strlen("euler")) == 0) 
      num_method = EULER;
   else if (strncmp(mxArrayToString(NUMMETHOD_IN), "milstein",
            strlen("milstein")) == 0) 
      num_method = MILSTEIN;
   else if (strncmp(mxArrayToString(NUMMETHOD_IN), "rk", strlen("rk")) == 0) 
      num_method = RK;
   else 
      mexErrMsgTxt("Only 'euler', 'milstein' and 'rk' are implemented");


   /* Call the worker routine. */
   stockPath(*mxGetPr(RANDSTATE_IN), lround(*mxGetPr(SAMPLES_IN)),
         *mxGetPr(DT_IN), *mxGetPr(SIGMA0_IN),
         *mxGetPr(S_0_IN), *mxGetPr(XI0_IN), *mxGetPr(MU_IN),
         *mxGetPr(P_IN), *mxGetPr(ALPHA_IN), N, num_method, stock_mn, vol_mn,
         xi_mn);
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
