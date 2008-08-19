function TestBenefit()

cd stocksim

disp('Generating stockpath 1');
tic
[stockpathc, volpathc, xipath] = stocksim(4537, 10000, 0.001, 0.2, 50, 0.2, 0.05, 0, 0.2, 5, 'euler');
toc;
cd ..
clear xipath;
disp('Compute benefit 1');
tic;
benefit_1 = ComputeProfit(stockpathc, volpathc, 0.001, 0.05, 5, 1000, 0.1);
toc;
save('benefits.mat', 'benefit_1', '-V6');
clear benefit_1;

disp('Compute benefit 2');
tic;
benefit_2 = ComputeProfit(stockpathc, volpathc, 0.001, 0.05, 5, 1000, 0.2);
toc;
save('benefits.mat', 'benefit_2', '-V6', '-append');
clear benefit_2;

disp('Compute benefit 3');
tic;
benefit_5 = ComputeProfit(stockpathc, volpathc, 0.001, 0.05, 5, 1000, 0.5);
toc;
save('benefits.mat', 'benefit_5', '-V6', '-append');
clear benefit_5;

disp('Compute benefit 4');
tic;
benefit_10 = ComputeProfit(stockpathc, volpathc, 0.001, 0.05, 5, 1000, 1);
toc;
save('benefits.mat', 'benefit_10', '-V6', '-append');
clear benefit_10;
disp('Compute benefit 5');
tic;
benefit_15 = ComputeProfit(stockpathc, volpathc, 0.001, 0.05, 5, 1000, 1.5);
toc;
save('benefits.mat', 'benefit_15', '-V6', '-append');
clear benefit_15;
clear stockpathc;
clear volpathc;

disp('Generating stockpath 2');
cd stocksim
tic
[stockpath2, volpath2, xipath] = stocksim(4537, 10000, 0.001, 0.2, 50, 0.2, 0.05, 0.2, 0.2, 5, 'euler');
toc;
cd ..
clear xipath;

disp('Compute benefit 2.1');
tic;
benefit2_1 = ComputeProfit(stockpath2, volpath2, 0.001, 0.05, 5, 1000, 0.1);
toc;
save('benefits.mat', 'benefit2_1', '-V6', '-append');
clear benefit2_1;

disp('Compute benefit 2.2');
tic;
benefit2_2 = ComputeProfit(stockpath2, volpath2, 0.001, 0.05, 5, 1000, 0.2);
toc;
save('benefits.mat', 'benefit2_2', '-V6', '-append');
clear benefit2_2;

disp('Compute benefit 2.3');
tic;
benefit2_5 = ComputeProfit(stockpath2, volpath2, 0.001, 0.05, 5, 1000, 0.5);
toc;
save('benefits.mat', 'benefit2_5', '-V6', '-append');
clear benefit2_5;

disp('Compute benefit 2.4');
tic;
benefit2_10 = ComputeProfit(stockpath2, volpath2, 0.001, 0.05, 5, 1000, 1);
toc;
save('benefits.mat', 'benefit2_10', '-V6', '-append');
clear benefit2_10;

disp('Compute benefit 2.5');
tic;
benefit2_15 = ComputeProfit(stockpath2, volpath2, 0.001, 0.05, 5, 1000, 1.5);
toc;

save('benefits.mat', 'benefit2_15', '-V6', '-append');
clear benefit2_15;

