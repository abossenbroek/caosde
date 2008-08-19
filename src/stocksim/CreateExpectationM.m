function CreateExpectationE(p, alpha, seed)

disp('Compute exact solution with Milstein dt = 0.0001');
tic;
approxE1 = stockTsim(5489, 5489, 10000, 1, 0.0001, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc

save(['mil_ran5489.mat'], 'approxE1', '-V6');
clear all;

disp('Compute exact solution with Milstein dt = 0.0002');
tic;
approxE2 = stockTsim(5489, 5489, 10000, 2, 0.0002, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc

save(['mil_ran5489.mat'], 'approxE2', '-V6', '-append');
clear all;

disp('Compute exact solution with Milstein dt = 0.0005');
tic;
approxE5 = stockTsim(5489, 5489, 10000, 5, 0.0005, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc

save(['mil_ran5489.mat'], 'approxE5', '-V6', '-append');
clear all;


disp('Compute exact solution with Milstein dt = 0.0008');
tic;
approxE8 = stockTsim(5489, 5489, 10000, 8, 0.0008, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc

save(['mil_ran5489.mat'], 'approxE8', '-V6', '-append');
clear all;

disp('Compute exact solution with Milstein dt = 0.001');
tic;
approxE10 = stockTsim(5489, 5489, 10000, 10, 0.001, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc

save(['mil_ran5489.mat'], 'approxE10', '-V6', '-append');
clear all;

disp('Compute exact solution with Milstein dt = 0.002');
tic;
approxE20 = stockTsim(5489, 5489, 10000, 20, 0.002, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc

save(['mil_ran5489.mat'], 'approxE20', '-V6', '-append');
clear all;

disp('Compute exact solution with Milstein dt = 0.005');
tic;
approxE50 = stockTsim(5489, 5489, 10000, 50, 0.005, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc

save(['mil_ran5489.mat'], 'approxE50', '-V6', '-append');
clear all;


disp('Compute exact solution with Milstein dt = 0.008');
tic;
approxE80 = stockTsim(5489, 5489, 10000, 80, 0.008, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc
save(['mil_ran5489.mat'], 'approxE80', '-V6', '-append');
clear all;

disp('Compute exact solution with Milstein dt = 0.01');
tic;
approxE100 = stockTsim(5489, 5489, 10000, 100, 0.01, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc
save(['mil_ran5489.mat'], 'approxE100', '-V6', '-append');
clear all;


disp('Compute exact solution with Milstein dt = 0.02');
tic;
approxE200 = stockTsim(5489, 5489, 10000, 200, 0.02, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc
save(['mil_ran5489.mat'], 'approxE200', '-V6', '-append');
clear all;

disp('Compute exact solution with Milstein dt = 0.04');
tic;
approxE400 = stockTsim(5489, 5489, 10000, 400, 0.04, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc
save(['mil_ran5489.mat'], 'approxE400', '-V6', '-append');
clear all;

disp('Compute exact solution with Milstein dt = 0.05');
tic;
approxE500 = stockTsim(5489, 5489, 10000, 500, 0.05, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc
save(['mil_ran5489.mat'], 'approxE500', '-V6', '-append');
clear all;

disp('Compute exact solution with Milstein dt = 0.1');
tic;
approxE1000 = stockTsim(5489, 5489, 10000, 1000, 0.1, 0.2, 50, 0.2, 0.2, 0.2, 0.2, 1, 'milstein');
toc

save(['mil_ran5489.mat'], 'approxE1000', '-V6', '-append');
clear all;

