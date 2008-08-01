function volInt = EulerVol(stock, vol, xi, p, Dt, phiVol);

volInt = vol - (vol - xi) * Dt + p * vol * phiVol;

