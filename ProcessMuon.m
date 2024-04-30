clear
dirName = 'muon';
dirName = 'muon_meanofLR';
dirName = 'muon_corrected';
cd(dirName);
name = ls;
cd ..
name = name(3:end,:);
fileNum = size(name, 1);
muNum = zeros(fileNum, 1);
cHit = zeros(fileNum, 16);
eHit = zeros(fileNum, 17);
hit = [];
hitarray = zeros(4);
for ii = 1:fileNum
    fileName = [dirName, '\', name(ii,:)];
    while contains(fileName, ' ')
        fileName = strrep(fileName, ' ', '');
    end
    temp = load(fileName);
    event = temp.event;
    edep = event.edep;
    edep(edep <= 200) = 0;
    edep(:,:,sum(sum(edep)) < 10000) = [];
    % edep = edep(1:3,1:3,:);
    fhit = permute(sum(sum(edep > 0)), [3, 1, 2]);
    [cHit(ii,:), eHit(ii,:)] = histcounts(fhit, 0.5:17);
    muNum(ii, 1) = length(fhit);
    hit = [hit; fhit];
    hitarray = hitarray + sum(edep > 0, 3);
end
%%
stdHit = sqrt(std(cHit(:,2:end))' .* 277);
stdMuNum = sqrt(std(muNum) .* 277);

muRateHour = mean(muNum);
stdMuRH = std(muNum);

muonRate = length(hit) ./ fileNum ./ 3600;
stdMuRate = stdMuRH ./3600;

eff = 0.7823;
stdEff = 0.0010;

muonRateReal = muonRate ./ eff;
stdMuRateReal = sqrt((stdMuRate ./ eff) .^ 2 + (stdEff ./ muonRate ./ eff .^ 2) .^ 2);

%%
fcolor = '#6279c1';
falpha = 0.5;

hfig = myfigure;
axhit = axes(hfig);
h = histogram(axhit, hit, 'Normalization', 'probability',...
    FaceColor = fcolor, FaceAlpha = falpha, EdgeColor = fcolor);
xlabel(axhit, 'Number of hit');
ylabel(axhit, 'Probability');
set(axhit, "XLim", [0, 16])

%%
hitcounts = h.BinCounts;
weight = [0, hitcounts ./ length(hit)]';
% weight = ones(1, 16)';

wnlm = load('hitfit-2-4.mat').wnlm;
wnlm1 = load('hitfit-5-7.mat').wnlm1;
wnlm2 = load('hitfit-8-11.mat').wnlm2;

[ypred, ypredci] = predict(wnlm, (2:4)', 'Prediction', 'observation', 'W', weight(2:4));
[ypred1, ypredci1] = predict(wnlm1, (5:7)', 'Prediction', 'observation', 'W', weight(5:7));
[ypred2, ypredci2] = predict(wnlm2, (8:16)', 'Prediction', 'observation', 'W', weight(8:16));

ypred = [ypred; ypred1; ypred2];
ypredci = [ypredci; ypredci1; ypredci2];
sigma = (ypredci(:,2) - ypredci(:,1)) ./ 2 ./ 1.96;

totalTrackLength = sum(ypred .* hitcounts');
s2TTL = sum((sigma .* hitcounts') .^ 2 + (ypred .* stdHit) .^ 2);
% s2TTL = sum(sigma .^2 .* hitcounts');
sigmaTTL = sqrt(s2TTL); 

meanTrackLength = totalTrackLength ./ length(hit);
stdMTL = sqrt((sigmaTTL ./ length(hit)) .^ 2 + ...
    (totalTrackLength ./ (-length(hit)).^2 .* stdMuNum) .^ 2);
%%
trackLength1 = 7.5480;
stdTL1 = 4.8551;

eff1 = 0.2177;
stdEff1 = 0.0011;

muon1Rate = muonRateReal .* eff1;
stdMu1Rate = sqrt((muonRateReal .* stdEff1).^2 + (eff1 .* stdMuRateReal).^ 2);

numHit1 = muon1Rate .* fileNum .* 3600;
stdNumHit1 = stdMu1Rate .* fileNum .* 3600;

totalTrackLength1 = trackLength1 .* numHit1;
s2TTL1 = (stdNumHit1 .* trackLength1) .^ 2 + (stdTL1 .* numHit1) .^ 2;
stdTTL1 = sqrt(s2TTL1);

%%
totalTrackLengthReal = totalTrackLength + totalTrackLength1;
stdTTLR = sqrt(s2TTL + s2TTL1);

numMuonReal = muonRateReal .* fileNum .* 3600;
meanTrackLengthReal = totalTrackLengthReal ./ numMuonReal;
stdMTLR = sqrt((stdTTLR ./ numMuonReal) .^ 2 + (stdMuRateReal .*...
    totalTrackLengthReal ./ (3600 .* fileNum .* muonRateReal .^ 2)) .^ 2);

%%
ebfig = myfigure;
axeb = axes(ebfig);
errorbar(axeb, 2:16, ypred', sigma, 'o', Color = fcolor, MarkerFaceColor = 'auto');
set(axeb, 'XLim', [0, 16]);
xlabel(axeb, 'Number of hit');
ylabel(axeb, 'Muon trajectory length / mm');

myfigstyle([axhit, axeb]);
