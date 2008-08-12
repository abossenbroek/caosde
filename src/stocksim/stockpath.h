#ifndef _STOCKPATH_H_
#define _STOCKPATH_H_

#define  EULER        0
#define  MILSTEIN     1
#define  RK           2
#define  NO_METHOD    4

void stockPath(double rng_state, long int samples, double Dt, double sigma_0,
		double S_0, double xi_0, double mu, double p, double alpha, long int N,
		char num_method, double *stock_mn, double *vol_mn, double *xi_mn);

#endif
