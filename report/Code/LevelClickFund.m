function valueAtT = LevelClickFund(T0, T, dt, p, alpha, method, S0, ...
	sigma, r, initialStocks, level)

N = (T - T0) / dt;

randn('state', 1);

walk1 = randn(N);
walk2 = randn(N);

fund = Portfolio;
fund.numberOfOptions = initialStocks;
fund.numberOfStocks = initialStocks;
fund.lastStockPrice = S0;
fund.level = level;

% At this point we use (1 + level) * S_0 as the strike price of the option.
if strcmp(method, 'BS')
	% Compute the initial value of the fund.
	fund.value = initialStocks * S0 + BlackScholes(S0, sigma, 5, ...
		(1 + level) * S0, r, 'put');
	% Generate a stock path using the standard Black-Scholes model.
	stock = BlackScholesStock(r, dt, T0, T, sigma, walk1, S0);

	for i = 1 : N
		if (fund.lastStockPrice * (1 + fund.level)) < stock(i)
			% Compute the profit we can make by selling the current options. 
			currentOptionsValue = fund.numberOfOptions * ...
				BlackScholes(fund.lastStockPrice, sigma, 5 - fund.lastChange, ...
				(1 + fund.level) * fund.lastStockPrice, r, 'put');

			% Compute the price of options which will have to be bought.
			newOptionPrice = BlackScholes(stock(i), sigma, 5 - i * dt, ...
				(1 + fund.level) * stock(i), r, 'put');
		end
	end

elseif strcmp(method, 'ST')
else
	ME = MException('LevelClickFund:wrongInput', 'Invalid method');
	throw(ME);
end	






