function stockInt = EulerStock(stock, vol, xi, mu, Dt, phiStock)

stockInt = stock + mu * stock * Dt + vol * stock * phiStock;

