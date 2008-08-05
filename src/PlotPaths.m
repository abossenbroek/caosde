function PlotPaths...
   (randstate, samples, Dt, sigma0, S0, xi0, mu, p, alpha, T, ...
      numMethod)

[stockAvg, volAvg, xiAvg, stockPaths, volPaths, xiPaths] = NumWrapper...
   (randstate, samples, Dt, sigma0, S0, xi0, mu, p, alpha, T, ...
      numMethod);

subplot(2, 2, [1 2]), plot(0 + Dt : Dt : 1, stockAvg)
title 'Expectation of the Stock path'
xlabel 'Time'
ylabel 'Stock Price'
whitebg('white')
subplot(2, 2, 3), plot(0 + Dt : Dt : 1, volAvg)
title 'Expectation of the Volatility path'
xlabel 'Time'
ylabel 'Volatility'
subplot(2, 2, 4), plot(0 + Dt : Dt : 1, xiAvg)
title 'Xi Average Path'
xlabel 'Time'
ylabel 'Volatility'
