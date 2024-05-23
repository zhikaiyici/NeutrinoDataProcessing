clear

IBDMuonic = load('IBDEvent_P25D28_muonic.mat');
IBD = load('IBDEvent_P25D28_1st6_rest700.mat');

IBDMuonic = IBDMuonic.IBDEvent;
IBD = IBD.IBDEvent;

muonic = cell2mat(IBDMuonic(:,3)) - cell2mat(IBD(:,3));
muonic = muonic > 0;

IBDMuonic = IBDMuonic(muonic,:);
IBD = IBD(muonic,:);

IBDMuonic = cell2mat(IBDMuonic(:,2));
IBD = cell2mat(IBD(:,2));

timeGap = [];
edep = [];
for ii = 1:size(IBD, 1)
    IDMuonic = IBDMuonic(ii).ID;
    ID = IBD(ii).ID;
    [ID, ind] = setdiff(IDMuonic, ID);

    time = IBDMuonic(ii).time(ind);
    gap = time(2:end) - time(1:end-1);
    gap(gap > 300) = [];
    timeGap = [timeGap; gap];

    e = IBDMuonic(ii).edep(:,:,ind);
    edep = cat(3, edep, e);
end

sliceNum = 10;
nTot = size(timeGap, 1);
inc = ceil(nTot ./ sliceNum);
incT = 16;
tEdges = 0:incT:320;
jj = 1;
for ii = 0:inc:nTot - 1
    if ii + inc < nTot
        data = timeGap(ii + 1:ii + inc);
    else
        data = timeGap(ii + 1:end);
    end
    [tCounts(jj,:), ~] = histcounts(data, tEdges);
    jj = jj + 1;
end

fcolor = '#6279c1';
falpha = 0.3;
myfigure(Name = 'IBDTime');
varTCounts = var(tCounts);
sigma = sqrt(varTCounts .* sliceNum);
tEdges = (tEdges(1:end - 1) + tEdges(2:end))./2;
tCounts = sum(tCounts);

errorbar(tEdges(1:end - 1), tCounts(1:end - 1), sigma(1:end - 1), ...
    'o', Color = fcolor, MarkerFaceColor = 'auto');

tEdges = tEdges(1:end - 1);
tCounts = tCounts(1:end - 1);

% t = histogram(timeGap, tEdges, DisplayStyle = "stairs");
% tCounts = t.Values(3:end);

hold on

beta0 = [30; 50; 10];
ff = @(beta0, x)(beta0(1) / beta0(2) * exp(-x / beta0(2)) + beta0(3));
x = 0:0.1:300;
% fitOp = statset('RobustWgtFun', 'bisquare');
[beta, r, J, CovB, MSE, errInfo] = nlinfit(tEdges, tCounts, ff, beta0);%, fitOp);
nlm = fitnlm(tEdges, tCounts, ff, beta0);
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


