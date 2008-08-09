function stockInt = MilsteinStock(stock, vol, mu, Dt, phiStock)

stockInt = stock + mu * stock * Dt + vol * stock * phiStock ...
	+ 1 / 2 * vol * vol * stock * ((phiStock * phiStock) - Dt);

