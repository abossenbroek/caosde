function Verify()

[stockMnEuler, volMnEuler, xiMnEuler] = ...
	stocksim(1, 1, 0.001, 0.2,  50, 0.2, 0.1, 0, 0.1, 1, 'euler');

[stockMnMilstein, volMnMilstein, xiMnMilstein] = ...
	stocksim(1, 1, 0.001, 0.2,  50, 0.2, 0.1, 0, 0.1, 1, 'milstein');

[stockMnRK, volMnRK, xiMnRK] = ...
	stocksim(1, 1, 0.001, 0.2,  50, 0.2, 0.1, 0, 0.1, 1, 'rk');

disp('For 1, 1, 0.001, 0.2, 50, 0.2, 0.1, 0, 0.1, 1');
disp(['max difference euler <-> milstein stock ' num2str(max(stockMnEuler - stockMnMilstein)) ]);
disp(['max difference euler <-> milstein vol ' num2str(max(volMnEuler - volMnMilstein)) ]);
disp(['max difference euler <-> milstein xi ' num2str(max(xiMnEuler - xiMnMilstein)) ]);

disp(['max difference euler <-> rk stock ' num2str(max(stockMnEuler - stockMnRK)) ]);
disp(['max difference euler <-> rk vol ' num2str(max(volMnEuler - volMnRK)) ]);
disp(['max difference euler <-> rk xi ' num2str(max(xiMnEuler - xiMnRK)) ]);

disp(['max difference milstein <-> rk stock ' num2str(max(stockMnMilstein - stockMnRK)) ]);
disp(['max difference milstein <-> rk vol ' num2str(max(volMnMilstein - volMnRK)) ]);
disp(['max difference milstein <-> rk xi ' num2str(max(xiMnMilstein - xiMnRK)) ]);
disp(' ');
disp(' ');


[stockMnEuler, volMnEuler, xiMnEuler] = ...
	stocksim(1, 1, 0.001, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'euler');

[stockMnMilstein, volMnMilstein, xiMnMilstein] = ...
	stocksim(1, 1, 0.001, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'milstein');

[stockMnRK, volMnRK, xiMnRK] = ...
	stocksim(1, 1, 0.001, 0.2,  50, 0.2, 0.1, 0.1, 0.1, 1, 'rk');

disp('For 1, 1, 0.001, 0.2, 50, 0.2, 0.1, 0.1, 0.1, 1');
disp(['max difference euler <-> milstein stock ' num2str(max(stockMnEuler - stockMnMilstein)) ]);
disp(['max difference euler <-> milstein vol ' num2str(max(volMnEuler - volMnMilstein)) ]);
disp(['max difference euler <-> milstein xi ' num2str(max(xiMnEuler - xiMnMilstein)) ]);

disp(['max difference euler <-> rk stock ' num2str(max(stockMnEuler - stockMnRK)) ]);
disp(['max difference euler <-> rk vol ' num2str(max(volMnEuler - volMnRK)) ]);
disp(['max difference euler <-> rk xi ' num2str(max(xiMnEuler - xiMnRK)) ]);

disp(['max difference milstein <-> rk stock ' num2str(max(stockMnMilstein - stockMnRK)) ]);
disp(['max difference milstein <-> rk vol ' num2str(max(volMnMilstein - volMnRK)) ]);
disp(['max difference milstein <-> rk xi ' num2str(max(xiMnMilstein - xiMnRK)) ]);




