#include <math.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

#include "stockpath.h"

#include "mex.h"
#include "matrix.h"
   
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

#define  STOCKPATHS     plhs[0]
#define  VOLPATHS       plhs[1]
#define  XIPATHS        plhs[2]

void 
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
   mwSize   N;
   mwSize   samples;
   double   *stock_paths;
   double   *vol_paths;
   double   *xi_paths;
   char     num_method = NO_METHOD;
   gsl_rng				 *rng_stock, *rng_vol;
	const gsl_rng_type *rng_type;

	if (nrhs != 11) {
		mexErrMsgTxt("Eleven input arguments required.");
	} else if (nlhs != 3) {
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
   STOCKPATHS = mxCreateDoubleMatrix(samples, N, mxREAL);
   VOLPATHS = mxCreateDoubleMatrix(samples, N, mxREAL);
   XIPATHS = mxCreateDoubleMatrix(samples, N, mxREAL);

   /* Use local variables to access the memory allocated above. */
   stock_paths = mxGetPr(STOCKPATHS);
   vol_paths = mxGetPr(VOLPATHS);
   xi_paths = mxGetPr(XIPATHS);

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
         *mxGetPr(MU_IN), *mxGetPr(P_IN), *mxGetPr(ALPHA_IN), N, -1,
         num_method, stock_paths, vol_paths, xi_paths);
	gsl_rng_free(rng_stock);
	gsl_rng_free(rng_vol);
}

/* vim: set et : tw=80 : spell spelllang=en: */
