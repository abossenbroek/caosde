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

if strcmp(method, 'BS')
	% Compute the initial value of the fund.
	fund.value = initialStocks * S0 + BlackScholes(S0, sigma, 5, ...
		S0, r, 'put');
	% Generate a stock path using the standard Black-Scholes model.
	stock = BlackScholesStock(r, dt, T0, T, sigma, walk1, S0);

	for i = 1 : N
		if (fund.lastStockPrice * (1 + fund.level)) < stock(i)
			% Compute the profit we can make by selling the current options. 
			currentOptionsValue = fund.numberOfOptions * ...
				BlackScholes(fund.lastStockPrice, sigma, 5 - fund.lastChange, ...
				fund.lastStockPrice, r, 'put');

			% Compute the price of options which will have to be bought.
			newOptionPrice = BlackScholes(stock(i), sigma, 5 - i * dt, ...
				stock(i), r, 'put');

			% Using this price we can compute the number of options which have to
			% be bought.
			quantity = (fund.numberOfStocks * stock(i) + currentOptionsValue) ...
				/ (newOptionPrice + stock(i))

			if quantity < 0
				ME = MException('LevelClickFund:FewShares', ...
					'Too few shares');
				throw(ME);
			end

			if abs(newOptionPrice * quantity - stock(i) * ...
				(fund.numberOfStocks - quantity) - currentOptionsValue) > 0.00001
				ME = MException('LevelClickFund:invalidPor', ...
					'Portfolio is incorrect');
				throw(ME);
			end

			% Adjust the number of stocks.
			fund.numberOfStocks = quantity;
			% Set the new number of options.
			fund.numberOfOptions = quantity;
			fund.lastChange = i * dt;
			fund.lastStockPrice = stock(i);
		end

		valueAtT =fund.numberOfOptions * ...
				BlackScholes(fund.lastStockPrice, sigma, 5 - fund.lastChange, ...
				fund.lastStockPrice, r, 'put') + stock(i) * fund.numberOfStocks;
 
	end

elseif strcmp(method, 'ST')
else
	ME = MException('LevelClickFund:wrongInput', 'Invalid method');
	throw(ME);
end	






