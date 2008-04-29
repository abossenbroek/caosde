T = 5;%years
T0 = 0;
dt = 2^-9;
n = (T-T0)/dt;% n timesteps

t = zeros(n,1);
t(1) = T0;
for i = 2 : n
    t(i) = t(i-1)+dt;
end

N1 = randn(n,1);
N2 = randn(n,1);

%The effect of alpha or p in the stock price
p = 3.0;

alpha = 0.01;
[t,xi] = stockPath(T0, T, dt, N1, N2, p, alpha);
plot(t,xi,'b-');
hold on;
alpha = 0.4;
[t,xi] = stockPath(T0, T, dt, N1, N2, p, alpha);
plot(t,xi,'r-');
alpha = 5.0;
[t,xi] = stockPath(T0, T, dt, N1, N2, p, alpha);
plot(t,xi,'g-');
alpha = 100.0;
[t,xi] = stockPath(T0, T, dt, N1, N2, p, alpha);
plot(t,xi,'black-');
hold off;