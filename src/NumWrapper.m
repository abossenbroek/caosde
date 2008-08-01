function [stockAvg, volAvg, stockPaths, volPaths, xiPaths] = NumWrapper...
   (randstate, samples, Dt, sigma0, S0, xi0, mu, p, alpha, T, ...
      numMethod)

% Set the state of the random number generator.
randn('state', randstate);

% Determine the number of steps which will have to be simulated.
N = T / Dt;

if strcmp(numMethod, 'Euler')
   numMethodStock = @(stock, vol, xi, mu, Dt, phiStock)...
      EulerStock(stock, vol, xi, mu, Dt, phiStock);
   numMethodVol   = @(stock, vol, xi, p, Dt, phiVol)...
      EulerVol(stock, vol, xi, p, Dt, phiVol);
end

% Create a matrix which will hold all the paths.
stockPaths = zeros(samples, N);
volPaths = zeros(samples, N);
xiPaths = zeros(samples, N);

% Set the initial values.
stockPaths(:, 1) = S0;
volPaths(:, 1) = sigma0;
xiPaths(:, 1) = xi0;

for i = 1 : samples
  % Generate two sets of random numbers.
  phiStock = randn(N - 1);
  phiVol = randn(N - 1);
  % Create the Brownian motion
  phiStock = phiStock .* sqrt(Dt);
  phiVol = phiVol .* sqrt(Dt);

  % First compute the volatility path since this cannot be compuated in a
  % vector operation.
  for j = 2 : N
     volt = volPaths(i, j - 1);
     xit = xiPaths(i, j - 1);
     stockt = stockPaths(i, j - 1);
     % Compute the approximation of the volatility using the Euler method.
     volPaths(i, j) = numMethodVol(stockPaths(i, j - 1), ...
         volPaths(i, j - 1), xiPaths(i, j - 1), p, Dt, phiVol(j - 1));
     % Compute the approximation of the xi using Runge-Kutta fourth order
     % method.
     k1 = 1 / alpha * (volt - xit);
     k2 = 1 / alpha * (volt + Dt / 2 * k1 - xit);
     k3 = 1 / alpha * (volt + Dt / 2 * k2 - xit);
     k4 = 1 / alpha * (volt + Dt * k3 - xit);
     xiPaths(i, j) = xit + Dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
     % Use the Backward Euler method to compute xi.
     %xiPaths(i, j) = xit + 1 / alpha * (volt - xit) * Dt;

     % Compute the approximation of the stock using the Euler method.
     stockPaths(i, j) = numMethodStock(stockPaths(i, j - 1), ...
         volPaths(i, j - 1), xiPaths(i, j - 1), mu, Dt, phiStock(j - 1));
  end

end

stockAvg = sum(stockPaths);
stockAvg = stockAvg / samples;
volAvg = sum(volPaths);
volAvg = volAvg / samples;

% vim: tabstop=2:expandtab
