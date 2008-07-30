function premium = BlackScholes(priceZero, volatility, T, strikePrice, ...
	interestRate, type)
% Compute the premium of an option.
% Parameters:
%	   priceZero:     the initial stock price.
%	   volatility:    the volatility of the stock.
%	   T:             the expiry date
%	   interestRate:	
%	   strikePrice:
%	   interestRate:
%	   type:          the type of the option 'put' or 'call'
% Output:
%     premium:       

d1 = (log(priceZero / strikePrice) + (interestRate + volatility^2 / 2) * T) ...
	/ (volatility * sqrt(T));

d2 = d1 - (volatility * sqrt(T));


if strcmp(type, 'call')
	premium = priceZero * normcdf(d1) - ...
		strikePrice * (exp(-interestRate * T)) * normcdf(d2, 0, 1);
elseif strcmp(type, 'put')
	premium = strikePrice * exp(-interestRate * T) * normcdf(-d2) - ... 
		priceZero * normcdf(-d1);
end

% vim: expandtab:ft=matlab
