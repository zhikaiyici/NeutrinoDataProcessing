clear;
muonFamily = load("MuonFamily_corrected.mat").muonFamily;
muonFamily = load("MuonFamily_1st6.mat").muonFamily;
arraySize = 4;
neutronEdep = [];
for ii = 1:size(muonFamily, 1)
    muonF = muonFamily{ii};
    for jj = 1:size(muonF, 1)
        neutronEdep = cat(3, neutronEdep, muonF(jj, 2).edep);
    end
end

tempSpec = neutronEdep ./ 1000;
tempSpec1 = reshape(permute(tempSpec, [2, 1, 3]), [arraySize, size(tempSpec, 3) * arraySize]);
tempSpec1 = reshape(tempSpec1, [arraySize * arraySize, numel(tempSpec1) ./ (arraySize * arraySize)])';
sortedSpec = sort(tempSpec1, 2, 'descend');

edepTotal = sum(sortedSpec, 2);

axesE = axes(figure);
PlotSpectrum(edepTotal, 0.01, 'e', axesE);
hold(axesE, "on");
for ii = 1:arraySize * arraySize
    temp = sortedSpec(:, ii);
    temp(temp == 0) = [];
    [c, e] = PlotSpectrum(temp, 0.01, 'e', axesE);
end
legend(axesE, '{\itE}_{total}', '{\itE}_{1st}', '{\itE}_{2nd}', '{\itE}_{3rd}', '{\itE}_{4th}', ...
    'FontName', 'Times New Roman','Box','off');
hold(axesE,'off');
axesE.XLim = [0,10];
set(axesE, 'yscale', 'log');
