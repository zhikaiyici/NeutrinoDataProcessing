clear

load('IBDEvent_P25D28.mat');
arraySize = 4;
edepPrompt = [];
edepDelayed = [];
for ii = 1:size(IBDEvent, 1)
    edep = IBDEvent{ii,2}(:,3);
    prompt = edep(1:2:end);
    prompt = cell2mat(prompt);
    edepPrompt = [edepPrompt; prompt];

    delayed = edep(2:2:end);
    delayed = cell2mat(delayed);
    edepDelayed = [edepDelayed; delayed];
end

edepPrompt = ReshapeDataMatrix(arraySize, edepPrompt);
edepDelayed = ReshapeDataMatrix(arraySize, edepDelayed);

promptTotal = sum(sum(edepPrompt));
delayedTotal = sum(sum(edepDelayed));
%%
figure;
tiledlayout(2, 2)
nexttile([2, 1]);
histogram2(promptTotal, delayedTotal, 'BinWidth', 200, 'DisplayStyle', 'bar3',...
    'ShowEmptyBins', 'off', 'LineStyle', 'none', 'FaceColor', 'flat',FaceAlpha=0.8);
colorbar;
colormap(cool);
set(gca,"XGrid","off","YGrid","off","ZGrid","off","Box","off");
nexttile;
PlotSpectrum(promptTotal, 100, 'e');
set(gca,"XGrid","off","YGrid","off","ZGrid","off","Box","off");
nexttile;
PlotSpectrum(delayedTotal, 100, 'e');
set(gca,"XGrid","off","YGrid","off","ZGrid","off","Box","off");
