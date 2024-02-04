clear
dirName = 'muon';
dirName = 'muon_meanofLR';
dirName = 'muon_corrected';
cd(dirName);
name = ls;
cd ..
name = name(3:end,:);
fileNum = size(name, 1);

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
    hit = [hit; fhit];
    hitarray = hitarray + sum(edep > 0, 3);
end

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

[ypred,ypredci] = predict(wnlm, (2:4)', 'Prediction', 'observation', 'W', weight(2:4));
[ypred1,ypredci1] = predict(wnlm1, (5:7)', 'Prediction', 'observation', 'W', weight(5:7));
[ypred2,ypredci2] = predict(wnlm2, (8:16)', 'Prediction', 'observation', 'W', weight(8:16));

ypred = [ypred; ypred1; ypred2];
ypredci = [ypredci; ypredci1; ypredci2];
sigma = (ypredci(:,2) - ypredci(:,1)) ./ 2 ./ 1.96;

totaltracklength = sum(ypred .* hitcounts');
sigmattl = sqrt(sum((sigma .* hitcounts') .^ 2)); 

%%
ebfig = myfigure;
axeb = axes(ebfig);
errorbar(axeb, 2:16, ypred', sigma, 'o', Color = fcolor, MarkerFaceColor = 'auto');
set(axeb, 'XLim', [0, 16]);
xlabel(axeb, 'Number of hit');
ylabel(axeb, 'Muon trajectory length / mm');

myfigstyle([axhit, axeb]);
