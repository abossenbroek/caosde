function volInt = RKVol(vol, xi, p, Dt, phiVol);

volhat = vol - (vol - xi) * Dt + p * vol * sqrt(Dt);
volInt = vol - (vol - xi) * Dt + p * vol * phiVol ...
	+ 1 / (2 * sqrt(Dt)) * (phiVol * phiVol - Dt) ...
	* (p * volhat - p * vol);

