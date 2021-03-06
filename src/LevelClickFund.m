function valueAtT = LevelClickFund(T0, T, dt, ...
	sigmaPath, r, initialStocks, level, stockPath)

N = (T - T0) / dt;

if N ~= length(stockPath)
	ME = MException('LevelClickFund:FewShares', ...
		'Too points in the stock path');
	throw(ME);
end

if N ~= length(sigmaPath)
	ME = MException('LevelClickFund:FewShares', ...
		'Too points in the vol path');
	throw(ME);
end


S0 = stockPath(1);

randn('state', 1);

fund = Portfolio;
fund.numberOfOptions = initialStocks;
fund.numberOfStocks = initialStocks;
fund.lastStockPrice = S0;
fund.level = level;

% Compute the initial value of the fund.
fund.value = initialStocks * S0 + BlackScholes(S0, sigmaPath(1), 5, ...
	S0, r, 'put');

for i = 1 : N
		if (fund.lastStockPrice * (1 + fund.level)) < stockPath(i)
			% Compute the profit we can make by selling the current options. 
			currentOptionsValue = fund.numberOfOptions * ...
			BlackScholes(fund.lastStockPrice, sigmaPath(i), 5 - fund.lastChange, ...
			fund.lastStockPrice, r, 'put');

			% Compute the price of options which will have to be bought.
			newOptionPrice = BlackScholes(stockPath(i), sigmaPath(i), 5 - i * dt, ...
			stockPath(i), r, 'put');

			% Using this price we can compute the number of options which have to
			% be bought.
			quantity = (fund.numberOfStocks * stockPath(i) + currentOptionsValue) ...
				/ (newOptionPrice + stockPath(i));

			if quantity < 0
				ME = MException('LevelClickFund:FewShares', ...
					'Too few shares');
				throw(ME);
			end

			if abs(newOptionPrice * quantity - stockPath(i) * ...
					(fund.numberOfStocks - quantity) - currentOptionsValue) > 0.00001
				ME = MException('LevelClickFund:invalidPor', ...
					'Portfolio is incorrect');
				throw(ME);
			end

			% Adjust the number of stockPaths.
			fund.numberOfStocks = quantity;
			% Set the new number of options.
			fund.numberOfOptions = quantity;
			fund.lastChange = i * dt;
			fund.lastStockPrice = stockPath(i);
		end

	end
	
if stockPath(i) > fund.lastStockPrice
	valueAtT = fund.numberOfOptions * ...
		BlackScholes(fund.lastStockPrice, sigmaPath(i), 5 - fund.lastChange, ...
		fund.lastStockPrice, r, 'put') + stockPath(i) * fund.numberOfStocks;
else
	naked = fund.numberOfOptions - fund.numberOfStocks;

	if naked < 0
		ME = MException('LevelClickFund:FewShares', ...
			'Too few shares');
		throw(ME);
	end
	% Execute the options to prevent loss.
	valueAtT = fund.numberOfOptions * fund.lastStockPrice;
	% Sell the remaining stocks at lower price.
	valueAtT = valueAtT + naked * stockPath(i);
end
