function CreateStockPlot(time, stock, vol, xi, samples, p, alpha, printFig)
%CREATEFIGURE(time,stock,vol,xi)
%  time:  vector of x data
%  stock:  vector of y data
%  vol:  vector of y data
%  xi:  vector of y data

%  Auto-generated by MATLAB on 05-Aug-2008 20:32:58

% Create figure
figure1 = figure('XVisual',...
    '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'Color', 'white');

% Create axes
axes1 = axes('Parent',figure1,'Position',[0.13 0.5838 0.775 0.3412]);
box('on');
hold('all');

% Create plot
plot(time,stock,'Parent',axes1);

% Create title
title('Expectation of the Stock path');

% Create xlabel
xlabel('Time');

% Create ylabel
ylabel('Stock Price');

% Create subplot
subplot1 = subplot(2,2,3,'Parent',figure1);
box('on');
hold('all');

% Create plot
plot(time,vol,'Parent',subplot1);

% Create title
title('Expectation of the Volatility path');

% Create xlabel
xlabel('Time');

% Create ylabel
ylabel('Volatility');

% Create subplot
subplot2 = subplot(2,2,4,'Parent',figure1);
box('on');
hold('all');

% Create plot
plot(time,xi,'Parent',subplot2);

% Create title
title('Xi Average Path');

% Create xlabel
xlabel('Time');

% Create ylabel
ylabel('Xi');

% Create textbox
annotation(figure1,'textbox', [0.3810 0.01584 0.5 0.03464],...
    'Interpreter','none',...
    'String',{['Samples = ' num2str(samples) '    Alpha = '  num2str(alpha)  '    p = '  num2str(p)]},...
    'FitBoxToText','off',...
    'LineStyle','none');

if (printFig)
	print(figure1, '-depsc2', ['s' num2str(samples) '_a' num2str(alpha) '_p' num2str(p) '.eps']);
end
