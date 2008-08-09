function volInt = MilsteinVol(vol, xi, p, Dt, phiVol);

volInt = vol - (vol - xi) * Dt + p * vol * phiVol ...
	+ 1 / 2 * p * p * vol * (phiVol * phiVol - Dt);

