function stockInt = RKStock(stock, vol, xi, mu, Dt, phiStock)

stockHat = stock + mu * stock * Dt + vol * stock * sqrt(Dt);
stockInt = stock + mu * stock * Dt + vol * stock * phiStock ...
	+ 1 / (2 * sqrt(Dt)) ...
	* (phiStock * phiStock - Dt) * (vol * stockHat - vol * stock);

