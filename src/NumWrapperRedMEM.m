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
      MilsteinStock(stock, vol, mu, Dt, phiStock);
   numMethodVol   = @(vol, xi, p, Dt, phiVol)...
      MilsteinVol(vol, xi, p, Dt, phiVol);
elseif strcmp(numMethod, 'RK')
   numMethodStock = @(stock, vol, mu, Dt, phiStock)...
      RKStock(stock, vol, mu, Dt, phiStock);
   numMethodVol   = @(vol, xi, p, Dt, phiVol)...
      RKVol(vol, xi, p, Dt, phiVol);
end

% Create a matrix which will hold all the paths.
stockAvg = zeros(N, 1);
volAvg = zeros(N, 1);
xiAvg = zeros(N, 1);

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
      % Generate two random numbers used for the Brownian motion.
      phiStock = randn() * sqrt(Dt);
      phiVol = randn() * sqrt(Dt);

      % Compute the approximation of the stock using a numerical method.
      stockt = numMethodStock(stockt, volt, mu, Dt, phiStock);
    
      % Compute the approximation of the volatility using a numerical method.
      % Note that the current value of the volatility has to be stored
      % since it is required for the approximation of xi_{t + 1}
      volt1 = volt;
      volt = numMethodVol(volt, xit, p, Dt, phiVol);


      % Compute the approximation of the xi using Runge-Kutta fourth order
      % method.
      k1 = 1 / alpha * (volt1 - xit);
      k2 = 1 / alpha * (volt1 + Dt / 2 * k1 - xit);
      k3 = 1 / alpha * (volt1 + Dt / 2 * k2 - xit);
      k4 = 1 / alpha * (volt1 + Dt * k3 - xit);

      xit = xit + Dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

      % Add these values to the average
      stockAvg(j) = stockAvg(j) + (stockt  - stockAvg(j)) / i;
      volAvg(j) = volAvg(j) + (volt - volAvg(j)) / i;
      xiAvg(j) = xiAvg(j) + (xit  - xiAvg(j)) / i;
  end
end

% vim: tabstop=2:expandtab:ft=matlab
