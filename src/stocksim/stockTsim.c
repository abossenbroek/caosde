#include <math.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

#include "mex.h"
#include "matrix.h"

#include "numerical.h"
   
#define  RANDSTATE_IN prhs[0]
#define  SAMPLES_IN   prhs[1]
#define  DT_IN	       prhs[2]
#define  SIGMA0_IN    prhs[3]
#define  S_0_IN	    prhs[4]
#define  XI0_IN       prhs[5]
#define  MU_IN	       prhs[6]
#define  P_IN	       prhs[7]
#define  ALPHA_IN     prhs[8]
#define  T_IN	       prhs[9]
#define  NUMMETHOD_IN prhs[10]

#define  STOCK_T      plhs[0]

enum { 
   EULER,
   MILSTEIN,
   RK,
   NO_METHOD
};

/* The worker routine. */
void 
stockPath(gsl_rng *rng_stock, gsl_rng *rng_vol, mwSize samples, double Dt,
      double sigma_0, double S_0, double xi_0, double mu, double p, double
      alpha, mwSize N, char num_method, double *stock_T)
{
   mwSize  i;
   mwSize  j;
   double  phi_stock;
   double  phi_vol;
   double  k_1, k_2, k_3, k_4;
   double  stock_t;
   double  vol_t;
   double  vol_t1;
   double  xi_t;
   double  (*num_method_stock)(const double, const double,
         const double, const double, const double) = NULL;
   double  (*num_method_vol)(const double,  const double,
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
      stock_t = S_0;
      vol_t = sigma_0;
      xi_t = xi_0;

      /* Generate the rest of the path. */
      for (j = 1; j < N; j++) {
	      phi_stock = gsl_ran_gaussian(rng_stock, 1) * sqrt(Dt);
	      phi_vol = gsl_ran_gaussian(rng_vol, 1) * sqrt(Dt);

         /* Compute the integral of the stock using the numerical method. */
         stock_t = num_method_stock(stock_t, vol_t, mu, Dt, phi_stock);
         /* Compute the integral of the volatility using the numerical method.
          * Store the current value of the volatility to use for the
          * approximation of xi. */
         vol_t1 = vol_t;
         vol_t = num_method_vol(vol_t, xi_t, p, Dt, phi_vol);

         /* Compute the approximation of the xi using Runge-Kutta fourth order
          * method. */
         k_1 = 1 / alpha * (vol_t1 - xi_t);
         k_2 = 1 / alpha * (vol_t1 + 0.5 * Dt * k_1 - xi_t);
         k_3 = 1 / alpha * (vol_t1 + 0.5 * Dt * k_2 - xi_t);
         k_4 = 1 / alpha * (vol_t1 + Dt * k_3 - xi_t);

         xi_t += Dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
      }
      stock_T[i] = stock_t;
   }

	return;
}


void 
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
   mwSize   N;
   mwSize   samples;
   double   *stock_T;
   char     num_method = NO_METHOD;
   gsl_rng				 *rng_stock, *rng_vol;
	const gsl_rng_type *rng_type;

	if (nrhs != 11) {
		mexErrMsgTxt("Eleven input arguments required.");
	} else if (nlhs != 1) {
		mexErrMsgTxt("Three output arguments required."); 
   } 

   /* Determine if all the parameters are the proper type. */
   if (!mxIsDouble(RANDSTATE_IN) || !mxIsDouble(SAMPLES_IN) ||
         !mxIsDouble(DT_IN) || !mxIsDouble(SIGMA0_IN) ||
         !mxIsDouble(S_0_IN) || !mxIsDouble(XI0_IN) ||
         !mxIsDouble(MU_IN) || !mxIsDouble(P_IN) ||
         !mxIsDouble(ALPHA_IN) || !mxIsDouble(T_IN))
      mexErrMsgTxt("All parameters have to be numbers."); 

   if (!(*mxGetPr(ALPHA_IN) > 0))
      mexErrMsgTxt("alpha must be larger than zero");

   if (!mxIsClass(NUMMETHOD_IN, "char"))
      mexErrMsgTxt("num_method must be a char");

   /* How many steps should be simulated? */
   N = lround(*mxGetPr(T_IN) / *mxGetPr(DT_IN));
   samples = (mwSize)lround(*mxGetPr(SAMPLES_IN));
   /* Assign memory to the output parameters. */
   STOCK_T = mxCreateDoubleMatrix(samples, 1,  mxREAL);

   /* Use local variables to access the memory allocated above. */
   stock_T = mxGetPr(STOCK_T);

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

   /* The Mersenne Twister generator should be good enough for our purpose. */
	rng_type = gsl_rng_mt19937;

   /* Create two random number generators. */
	rng_stock = gsl_rng_alloc(rng_type);
	rng_vol = gsl_rng_alloc(rng_type);
 
   /* Initialize the state of the random number generator. */
   gsl_rng_set(rng_stock, *mxGetPr(RANDSTATE_IN));
   gsl_rng_set(rng_vol, *mxGetPr(RANDSTATE_IN) + 1);
  

   /* Call the worker routine. */
   stockPath(rng_stock, rng_vol, samples, *mxGetPr(DT_IN),
         *mxGetPr(SIGMA0_IN), *mxGetPr(S_0_IN), *mxGetPr(XI0_IN),
         *mxGetPr(MU_IN), *mxGetPr(P_IN), *mxGetPr(ALPHA_IN), N,
         num_method, stock_T);
	gsl_rng_free(rng_stock);
	gsl_rng_free(rng_vol);
}

/* vim: set et : tw=80 : spell spelllang=en: */
