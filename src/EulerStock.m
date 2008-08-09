function stockInt = EulerStock(stock, vol, mu, Dt, phiStock)

stockInt = stock + mu * stock * Dt + vol * stock * phiStock;

