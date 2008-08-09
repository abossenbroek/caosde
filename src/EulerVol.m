function volInt = EulerVol(vol, xi, p, Dt, phiVol);

volInt = vol - (vol - xi) * Dt + p * vol * phiVol;

