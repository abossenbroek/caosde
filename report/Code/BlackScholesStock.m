function integralMonteCarlo = BlackScholesStock(interestRate, deltat, ...
	T0, T, sigma, phi, priceZero)

% Initialize the result.
integralMonteCarlo      = zeros(1, (T - T0) / deltat);
% Set T_0 = S_0
integralMonteCarlo(1)   = priceZero;

for i = 2 : (T / deltat + 1)
	% Compute the stock price realization using an integral of
	% dS_T = (r - \sigma^2 / 2)T + \sigma dW_T.
	%
	% The time index has to be adjusted since the first element of the resulting
	% vector is equal to S_0.
	integralMonteCarlo(i) = priceZero * exp((interestRate - (sigma^2) / 2) ...
		* i * deltat  + sigma * phi(i - 1) * sqrt(i * deltat));
end    

