clear

load('IBDEvent_P25D28.mat');
% IBDEvent = IBDEvent(:,2);
timeGap = [];
for ii = 1:size(IBDEvent, 1)
    time = IBDEvent{ii,2}(:,2);
    time = cell2mat(time);
    gap = time(2:end) - time(1:end-1);
    gap(gap > 300) = [];
    timeGap = [timeGap; gap];
end
%%
figure
tEdges = 0:8:300;
t = histogram(timeGap, tEdges, DisplayStyle = "stairs");

tEdges = (tEdges(1:end-1) + tEdges(2:end))./2;
tEdges = tEdges(3:end);
tCounts = t.Values(3:end);

% bar(tEdges, tCounts, 1, EdgeColor = "flat", FaceColor = "none");

hold on

beta0 = [50000; 10; 100];
ff = @(beta0, x)(beta0(1) / beta0(2) * exp(-x / beta0(2)) + beta0(3));
x = 0:0.1:300;
% fitOp = statset('RobustWgtFun', 'bisquare');
[beta, r, J, CovB, MSE, errInfo] = nlinfit(tEdges, tCounts, ff, beta0);%, fitOp);
[y, ~] = nlpredci(ff, tEdges, beta, r, "Jacobian", J, 'PredOpt','curve');
R = corrcoef(y, tCounts);
R2 = R(1,2).^2;
[y, delta] = nlpredci(ff, x, beta, r, "Jacobian", J, 'PredOpt','curve');
ci = nlparci(beta, r, "Jacobian", J);

plot(x, y);
xconf = [x, x(end:-1:1)];
yconf = [y + delta, y(end:-1:1) - delta(end:-1:1)];
p = fill(xconf, yconf, 'red');
p.FaceColor = [1, 0.8, 0.8];      
p.EdgeColor = 'none';
p.FaceAlpha = 0.5;

set(gca, 'fontname', 'Times New Roman', 'xgrid', 'off', 'ygrid', 'off', 'Box','off', FontSize = 10);
xlabel(gca, "\fontname{宋体}时间差 / \mus");
ylabel(gca, "\fontname{宋体}计数 / 8 \mus");
xlim([0, 300]);
legend(gca, '实验值', '拟合曲线', '拟合曲线\fontname{times new roman}95%\fontname{宋体}置信区间', fontname = '宋体', fontsize = 10);
legend(gca, 'boxoff');

% bar(tEdges,abs(tCounts - beta(3)), 1);
% plot(x, beta(1).*exp(-x ./ beta(2)));

hold off
