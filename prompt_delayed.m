clear

load('IBDEvent_P25D28_corrected.mat');
load('IBDEvent_P25D28_1st6.mat');
load('IBDEvent_P25D28_1st6_rest700.mat');
arraySize = 4;
channelWidth = 0.05;
edep = [];
for ii = 1:size(IBDEvent, 1)
    % edep = IBDEvent{ii,2}(:,3);
    e = IBDEvent{ii,2}.edep;
    edep = cat(3, edep, e);
end
edepPrompt = edep(:,:,1:2:end) ./ 1000;
edepDelayed = edep(:,:,2:2:end) ./ 1000;

% edepPrompt = ReshapeDataMatrix(arraySize, edepPrompt);
% edepDelayed = ReshapeDataMatrix(arraySize, edepDelayed);

promptTotal = sum(sum(edepPrompt));
delayedTotal = sum(sum(edepDelayed));

%%
figure;
tiledlayout(2, 2)
nexttile([2, 1]);
histogram2(promptTotal, delayedTotal, 'BinWidth', 0.2, 'DisplayStyle', 'bar3',...
    'ShowEmptyBins', 'off', 'LineStyle', 'none', 'FaceColor', 'flat',FaceAlpha=0.8);
colorbar;
colormap(cool);
set(gca,"XGrid","off","YGrid","off","ZGrid","off","Box","off");
nexttile;
PlotSpectrum(promptTotal, channelWidth, 'e');
set(gca,"XGrid","off","YGrid","off","ZGrid","off","Box","off");
nexttile;
PlotSpectrum(delayedTotal, channelWidth, 'e');
set(gca,"XGrid","off","YGrid","off","ZGrid","off","Box","off");

pFig = figure;
axesPrompt = axes(pFig);
PlotSpectrum(promptTotal, channelWidth, 'e', axesPrompt);
dFig = figure;
axesDelayed = axes(dFig);
PlotSpectrum(delayedTotal, channelWidth, 'e', axesDelayed);
hold([axesPrompt, axesDelayed], "on");
edep = reshape(permute(edep, [2, 1, 3]), [4, size(edep, 3) * 4]);
edep = reshape(edep, [16, numel(edep) ./ 16])';
[edepSorted, eIndex] = sort(edep, 2, 'descend');
[r1st, c1st] = ind2sub(4, eIndex);
edepPrompt = edepSorted(1:2:end,:) ./ 1000;
edepDelayed = edepSorted(2:2:end,:) ./ 1000;
for ii = 1:arraySize * arraySize
    PlotSpectrum(edepPrompt(:,ii), channelWidth, 'e', axesPrompt);
    PlotSpectrum(edepDelayed(:,ii), channelWidth, 'e', axesDelayed);
end
hold([axesPrompt, axesDelayed], "off");
