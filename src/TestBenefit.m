function TestBenefit()

disp('Generating stockpath 2');
cd stocksim
tic
[stockpath2, volpath2, xipath] = stocksim(10, 1500, 0.01, 0.2, 50, 0.2, 0.05, 0, 0.2, 5, 'euler');
toc;
cd ..
clear xipath;
disp('Compute benefit 2.0');
tic;
[n m] = size(stockpath2);
benefit2_0 = 100 * stockpath2(n, :)';
toc;
%save('benefits2_3.mat', 'benefit2_0', '-V6', '-append');
save('benefits2_3.mat', 'benefit2_0', '-V6');
clear benefit2_0;

disp('Compute benefit 2.1');
tic;
benefit2_01 = ComputeProfit(stockpath2, volpath2, 0.01, 0.05, 5, 100, 0.01);
toc;
save('benefits2_3.mat', 'benefit2_01', '-V6', '-append');
clear benefit2_01;

disp('Compute benefit 2.2');
tic;
benefit2_05 = ComputeProfit(stockpath2, volpath2, 0.01, 0.05, 5, 100, 0.05);
toc;
save('benefits2_3.mat', 'benefit2_05', '-V6', '-append');
clear benefit2_05;

disp('Compute benefit 2.3');
tic;
benefit2_1 = ComputeProfit(stockpath2, volpath2, 0.01, 0.05, 5, 100, 0.1);
toc;
save('benefits2_3.mat', 'benefit2_1', '-V6', '-append');
clear benefit2_1;

disp('Compute benefit 2.4');
tic;
benefit2_2 = ComputeProfit(stockpath2, volpath2, 0.01, 0.05, 5, 100, 0.2);
toc;
save('benefits2_3.mat', 'benefit2_2', '-V6', '-append');
clear benefit2_2;

disp('Compute benefit 2.5');
tic;
benefit2_5 = ComputeProfit(stockpath2, volpath2, 0.01, 0.05, 5, 100, 0.5);
toc;

save('benefits2_3.mat', 'benefit2_5', '-V6', '-append');

disp('Compute benefit 2.6');
tic;
benefit2_10 = ComputeProfit(stockpath2, volpath2, 0.01, 0.05, 5, 100, 1);
toc;

save('benefits2_3.mat', 'benefit2_10', '-V6', '-append');
clear benefit2_10;

disp('Compute benefit 2.7');
tic;
benefit2_15 = ComputeProfit(stockpath2, volpath2, 0.01, 0.05, 5, 100, 1.5);
toc;

save('benefits2_3.mat', 'benefit2_15', '-V6', '-append');
clear benefit2_15;

