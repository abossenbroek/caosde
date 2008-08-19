function FindParams(approxError)

%approxError = [0.0002 1/10000 * sum(abs(approxE1 - approxE2))]
%approxError = [approxError ; [0.0005 1/10000 * sum(abs(approxE1 - approxE5))]];
%approxError = [approxError ; [0.0008 1/10000 * sum(abs(approxE1 - approxE8))]];
%approxError = [approxError ; [0.001 1/10000 * sum(abs(approxE1 - approxE10))]];
%approxError = [approxError ; [0.002 1/10000 * sum(abs(approxE1 - approxE20))]];
%approxError = [approxError ; [0.005 1/10000 * sum(abs(approxE1 - approxE50))]];
%approxError = [approxError ; [0.008 1/10000 * sum(abs(approxE1 - approxE80))]];
%approxError = [approxError ; [0.01 1/10000 * sum(abs(approxE1 - approxE100))]];
%approxError = [approxError ; [0.02 1/10000 * sum(abs(approxE1 - approxE200))]];
%approxError = [approxError ; [0.05 1/10000 * sum(abs(approxE1 - approxE500))]];
%approxError = [approxError ; [0.1 1/10000 * sum(abs(approxE1 - approxE1000))]];
%FindParams(approxError)

function residual = strongOrder(x)
	residual = 0;

	for i = 1 : length(approxError)
		residual = residual + abs(approxError(i, 2) - x(1) * approxError(i, 1)^x(2));
	end
end

A = [-0.1 0 ; 0 -0.1];
b = [0 ; 0];
x0 = [0.01 ; 0.42];

x = fmincon(@strongOrder, x0, A, b)

end
