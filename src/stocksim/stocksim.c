#include <math.h>
#include <string.h>

#include "stockpath.h"

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

/* vim: set et : tw=80 : spell spelllang=en: */
