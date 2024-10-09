clear

load('IBDEvent_P25D28_corrected.mat');
load('IBDEvent_P25D28_1st6.mat');
load('IBDEvent_P25D28_withmuonic.mat');
load('IBDEvent_P25D28_1st6_rest700.mat');
% IBDEvent = IBDEvent(:,2);
timeGap = [];
edep = [];
for ii = 1:size(IBDEvent, 1)
    % time = IBDEvent{ii,2}(:,2);
    % time = cell2mat(time);
    time = IBDEvent{ii, 2}.time;
    gap = time(2:end) - time(1:end-1);
    gap(gap > 300) = [];
    timeGap = [timeGap; gap];
    e = IBDEvent{ii,2}.edep;
    edep = cat(3, edep, e);
end
edep = reshape(permute(edep, [2, 1, 3]), [4, size(edep, 3) * 4]);
edep = reshape(edep, [16, numel(edep) ./ 16])';
[edepSorted, eIndex] = sort(edep, 2, 'descend');
[r1st, c1st] = ind2sub([4, 4], eIndex);
edepPrompt = edepSorted(1:2:end,:);
edepDelayed = edepSorted(2:2:end,:);
%%
% load("MuonNeutronTimeGap_corrected.mat");
% load("MuonNeutronTimeGap_1st6.mat");

%%
sliceNum = 10;
nTot = size(timeGap, 1);
inc = ceil(nTot ./ sliceNum);
incT = 8;
tEdges = 0:incT:308;
jj = 1;
for ii = 0:inc:nTot
    if ii + inc < nTot
        data = timeGap(ii + 1:ii + inc);
    else
        data = timeGap(ii + 1:end);
    end
    [tCounts(jj,:), ~] = histcounts(data, tEdges);
    jj = jj + 1;
end
%%
fcolor = '#6279c1';
falpha = 0.3;
myfigure(Name = 'IBDTime');
tCounts = tCounts ./ incT;
varTCounts = var(tCounts);
sigma = sqrt(varTCounts .* sliceNum);
tEdges = (tEdges(1:end - 1) + tEdges(2:end))./2;
tCounts = sum(tCounts);

errorbar(tEdges(2:end - 1), tCounts(2:end - 1), sigma(2:end - 1), ...
    'o', Color = fcolor, MarkerFaceColor = 'auto');

tStart = 4;
tEdges = tEdges(tStart:end - 1);
tCounts = tCounts(tStart:end - 1);

wei = sigma(tStart:end - 1) ./ tCounts;
% wei = ones(size(tCounts));

% t = histogram(timeGap, tEdges, DisplayStyle = "stairs");
% tCounts = t.Values(3:end);

hold on

beta0 = [4e3; 20; 1e2];
ff = @(beta0, x)(beta0(1) / beta0(2) * exp(-x / beta0(2)) + beta0(3));
x = 0:0.1:300;
% fitOp = statset('RobustWgtFun', 'bisquare');
[beta, r, J, CovB, MSE, errInfo] = nlinfit(tEdges, tCounts, ff, beta0, Weights = wei);%, fitOp);
nlm = fitnlm(tEdges, tCounts, ff, beta0, Weights = wei);
[ypred, ypredci] = predict(nlm, x', 'Simultaneous', false);
[y, ~] = nlpredci(ff, tEdges, beta, r, "Jacobian", J, 'PredOpt','curve');
R = corrcoef(y, tCounts);
R2 = R(1,2).^2;
[y, delta] = nlpredci(ff, x, beta, r, "Jacobian", J, 'PredOpt','curve');
ci = nlparci(beta, r, "Jacobian", J, alpha =  1 - 0.6826);

plot(x, y, Color = fcolor);
xconf = [x, x(end:-1:1)];
yconf = [y + delta, y(end:-1:1) - delta(end:-1:1)];
p = fill(xconf, yconf, 'red');
p.FaceColor = fcolor;
p.EdgeColor = 'none';
p.FaceAlpha = falpha;
% p = fill(xconf, [ypredci(:,1)', fliplr(ypredci(:,2)')], 'red');
% p.FaceColor = 'none';
% p.EdgeColor = 'red';
% p.FaceAlpha = falpha;

set(gca, 'fontname', 'Times New Roman', 'xgrid', 'off', 'ygrid', 'off', 'Box','off', FontSize = 12);
xlabel(gca, "\fontname{times new roman}{\it\DeltaT} (\mus)");
ylabel(gca, "\fontname{times new roman}Events / " + num2str(incT) + " \mus");
xlim([0, 300]);
legend(gca, 'Measurement', 'Fit curve', "\fontname{times new roman}95% C.I.", ...
    fontname = 'times new roman', fontsize = 20);
legend(gca, 'boxoff');
myfigstyle(gca);

% bar(tEdges,abs(tCounts - beta(3)), 1);
% plot(x, beta(1).*exp(-x ./ beta(2)));

hold off
%%
IBDNum = cell2mat(IBDEvent(:,3)) ./ 2;
% IBDNum = char(IBDEvent(:,1));
IBDNum(169:186,:) = [];
IBDNum(end - 18:end,:) = [];
IBDNum = reshape(IBDNum, 24, []);
IBDNumDay = sum(IBDNum);
varIBDNum = var(IBDNum);
stdIBDNum = std(IBDNum);
stdIBDNumDay = sqrt(varIBDNum .* 24);
meanIBDNumDay = mean(IBDNumDay);
meanIBDNumDay = ones(size(IBDNumDay)) .* meanIBDNumDay;
stdMeanIBDNumDay = std(IBDNumDay);
myfigure(Name = "IBDNum");
hold on;

errorbar(IBDNumDay, stdIBDNumDay, 'o', Color = '#6279c1', MarkerFaceColor = 'auto');
plot(1:10, meanIBDNumDay, Color = fcolor, LineWidth = 1.5);
p = fill([1:10, 10:-1:1], [meanIBDNumDay + stdMeanIBDNumDay, meanIBDNumDay - stdMeanIBDNumDay], 'red');
p.FaceColor = fcolor;
p.EdgeColor = 'none';
p.FaceAlpha = falpha;
set(gca, 'xlim', [0, 11], 'XTick', 1:10);
hold off;

set(gca, 'fontname', 'Times New Roman', 'xgrid', 'off', 'ygrid', 'off', 'Box','off', FontSize = 12);
xlabel(gca, "\fontname{times new roman}Day");
ylabel(gca, "\fontname{times new roman}Events (day^{-1})");
legend(gca,'Events per day',  'Mean', '$1\sigma$ C.I.', 'Interpreter','latex', fontname = 'times new roman', fontsize = 20);
legend(gca, 'boxoff');
myfigstyle(gca);

%%

tCountsT = predict(nlm, tEdges') .* incT;
tCounts = tCounts' .* incT;
chi2 = sum((tCounts - tCountsT) .^ 2 ./ tCountsT);
chi2_ = sum(tCounts .^ 2 ./ tCountsT) - sum(tCounts);
dof = length(tCounts) - 1;


%%
% exportgraphics(gca, 'IBDTime.pdf', 'ContentType', 'vector');
