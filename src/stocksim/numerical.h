#ifndef _NUMERICAL_H_
#define _NUMERICAL_H_

double euler_stock(const double stock_t, const double vol_t, const
		double mu, const double Dt, const double phi_stock);
double euler_vol(const double vol_t, const double xi_t, const double p,
		const double Dt, const double phi_vol);

double milstein_stock(const double stock_t, const double vol_t, const
		double mu, const double Dt, const double phi_stock);
double milstein_vol(const double vol_t, const double xi_t, const double
		p, const double Dt, const double phi_vol);

double rk_stock(const double stock_t, const double vol_t, const double
		mu, const double Dt, const double phi_stock);
double rk_vol(const double vol_t, const double xi_t, const double p,
		const double Dt, const double phi_vol);

#endif /* _NUMERICAL_H_ */
