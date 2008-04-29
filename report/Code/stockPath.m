function [t, epsilon] = stockPath(T0, T, dt, N1, N2, p, alpha)
%   path  Simulates the asset price path
%
%   [S] = stockPath(T0, T, dt, N1, N2, p, alpha) computes stock trajectory
%
%
n = (T-T0)/dt;% n timesteps

t = zeros(n,1);
t(1) = T0;
for i = 2 : n
    t(i) = t(i-1)+dt;
end

S = zeros(n,1);
sigma = zeros(n,1);
epsilon = zeros(n,1);

S(1) = 50.0; sigma(1) = 0.2; epsilon(1) = 0.2; mju = 0.1;

for i = 2:n
    
    S(i) = S(i-1) + mju * S(i-1)* dt + sigma(i-1)*S(i-1)*sqrt(dt)*N1(i);

    sigma(i) = sigma(i-1) - (sigma(i-1)-epsilon(i-1)) * dt + p * sigma(i-1) * sqrt(dt) * N2(i);

    epsilon(i) = epsilon(i-1) + (1.0/alpha) * (sigma(i-1)-epsilon(i-1)) * dt;

end