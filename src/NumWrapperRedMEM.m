function [stockAvg, volAvg, xiAvg] = NumWrapperRedMEM(randstate, samples, ...
   Dt, sigma0, S0, xi0, mu, p, alpha, T, numMethod)
% This is a less memory consuming method than NumWrapper.


% Set the state of the random number generator.
randn('state', randstate);

% Determine the number of steps which will have to be simulated.
N = T / Dt;

% Determine which methods to use for the specified numerical method.
if strcmp(numMethod, 'Euler')
   numMethodStock = @(stock, vol, mu, Dt, phiStock)...
      EulerStock(stock, vol, mu, Dt, phiStock);
   numMethodVol   = @(vol, xi, p, Dt, phiVol)...
      EulerVol(vol, xi, p, Dt, phiVol);
elseif strcmp(numMethod, 'Milstein')
   numMethodStock = @(stock, vol, mu, Dt, phiStock)...
      MilsteinStock(stock, vol, xi, mu, Dt, phiStock);
   numMethodVol   = @(vol, xi, p, Dt, phiVol)...
      MilsteinVol(vol, xi, p, Dt, phiVol);
elseif strcmp(numMethod, 'RK')
   numMethodStock = @(stock, vol, mu, Dt, phiStock)...
      RKStock(stock, vol, mu, Dt, phiStock);
   numMethodVol   = @(vol, xi, p, Dt, phiVol)...
      RKVol(vol, xi, p, Dt, phiVol);
end

% Create a matrix which will hold all the paths.
stockAvg = zeros(N);
volAvg = zeros(N);
xiAvg = zeros(N);

stockAvg(1) = S0;
volAvg(1) = sigma0;
xiAvg(1) = xi0;

for i = 1 : samples
   % Set the initial values.
   stockt = S0;
   stocktAV = S0;
   volt = sigma0;
   voltAV = sigma0;
   xit = xi0;
   xitAV = xi0;

   % First compute the volatility path since this cannot be compuated in a
   % vector operation.
   for j = 2 : N
      % Generate two sets of random numbers.
      phiStock = randn(1);
      phiVol = randn(1);
      % Create the Brownian motion
      phiStock = phiStock * sqrt(Dt);
      phiVol = phiVol * sqrt(Dt);

      % Compute the approximation of the stock using a numerical method.
      stockt = numMethodStock(stockt, volt, mu, Dt, phiStock);
    
      % Compute the approximation of the volatility using a numerical method.
      % Note that the current value of the volatility has to be stored
      % since it is required for the approximation of xi_{t + 1}
      volt1 = volt;
      volt = numMethodVol(volt, xit, p, Dt);


      % Compute the approximation of the xi using Runge-Kutta fourth order
      % method.
      k1 = 1 / alpha * (volt1 - xit);
      k2 = 1 / alpha * (volt1 + Dt / 2 * k1 - xit);
      k3 = 1 / alpha * (volt1 + Dt / 2 * k2 - xit);
      k4 = 1 / alpha * (volt1 + Dt * k3 - xit);

      xit = xit + Dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

      % Antithetic variance is used to reduce the number of computations
      % required to converge.
      stocktAV = numMethodStock(stocktAV, voltAV, mu, Dt, -phiStock);
    
      % Compute the approximation of the volatility using a numerical method.
      % Note that the current value of the volatility has to be stored
      % since it is required for the approximation of xi_{t + 1}
      volt1AV = voltAV;
      voltAV = numMethodVol(voltAV, xitAV, p, Dt);

      % Compute the approximation of the xi using Runge-Kutta fourth order
      % method.
      k1 = 1 / alpha * (volt1AV - xitAV);
      k2 = 1 / alpha * (volt1AV + Dt / 2 * k1 - xitAV);
      k3 = 1 / alpha * (volt1AV + Dt / 2 * k2 - xitAV);
      k4 = 1 / alpha * (volt1AV + Dt * k3 - xitAV);

      xitAV = xitAV + Dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

      % Add these values to the average
      stockAvg(j) = stockAvg(j) + ((stockt + stocktAV) / 2 - stockAvg(j)) / i;
      volAvg(j) = volAvg(j) + ((volt + voltAV) / 2 - volAvg(j)) / i;
      xiAvg(j) = xiAvg(j) + ((xit + xitAV) / 2 - xiAvg(j)) / i;
  end
end

% vim: tabstop=2:expandtab:ft=matlab
