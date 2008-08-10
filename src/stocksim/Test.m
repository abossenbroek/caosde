function Test()

for p = 0 : 0.1 : 2
	for alpha = 0.01 : 0.1 : 1
		PlotPaths(2, 10000, 0.0005, 0.2,  50, 0.2, 0.2, p, alpha, 1, 'euler', true); 
		close all;
	end
	for alpha = 1 : 2 : 30 
		PlotPaths(2, 10000, 0.0005, 0.2,  50, 0.2, 0.2, p, alpha, 1, 'euler', true); 
		close all;
	end
end

clear all;
[stockMnEuler01, volMnEuler01, xiMnEuler01] = ...
	stocksim(1, 10000, 0.01, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'euler');

[stockMnEuler001, volMnEuler001, xiMnEuler001] = ...
	stocksim(1, 10000, 0.001, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'euler');

[stockMnEuler0001, volMnEuler0001, xiMnEuler0001] = ...
	stocksim(1, 10000, 0.0001, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'euler');

[stockMnEuler00001, volMnEuler00001, xiMnEuler00001] = ...
	stocksim(1, 10000, 0.00001, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'euler');

save -V6 euler.mat
clear all;

[stockMnMilstein01, volMnMilstein01, xiMnMilstein01] = ...
	stocksim(1, 10000, 0.01, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'milstein');

[stockMnMilstein001, volMnMilstein001, xiMnMilstein001] = ...
	stocksim(1, 10000, 0.001, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'milstein');

[stockMnMilstein0001, volMnMilstein0001, xiMnMilstein0001] = ...
	stocksim(1, 10000, 0.0001, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'milstein');

[stockMnMilstein00001, volMnMilstein00001, xiMnMilstein00001] = ...
	stocksim(1, 10000, 0.00001, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'milstein');

save -V6 milstein.mat;
clear all;

[stockMnRK01, volMnRK01, xiMnRK01] = ...
	stocksim(1, 10000, 0.01, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'rk');

[stockMnRK001, volMnRK001, xiMnRK001] = ...
	stocksim(1, 10000, 0.001, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'rk');

[stockMnRK0001, volMnRK0001, xiMnRK0001] = ...
	stocksim(1, 10000, 0.0001, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'rk');

[stockMnRK00001, volMnRK00001, xiMnRK00001] = ...
	stocksim(1, 10000, 0.00001, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'rk');

save -V6 rk.mat;

