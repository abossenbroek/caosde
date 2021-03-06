function PlotPaths...
   (randstate, samples, Dt, sigma0, S0, xi0, mu, p, alpha, T, ...
      numMethod, printFig)

% Generate the paths.
%[stockAvg, volAvg, xiAvg, stockPaths, volPaths, xiPaths] = NumWrapper...
%   (randstate, samples, Dt, sigma0, S0, xi0, mu, p, alpha, T, ...
%      numMethod);
[stockAvg, volAvg, xiAvg] = NumWrapperRedMEM(randstate, samples, Dt, ...
	sigma0, S0, xi0, mu, p, alpha, T, numMethod);

% Create the figure.
CreateStockPlot(0 + Dt : Dt : 1, stockAvg, volAvg, xiAvg, samples, p, alpha, printFig)

