clear
% dirName = 'sim_1600';
dirName = 'data';
cd(dirName);
% dataFolder = pwd;
name = ls;
cd ..
name = name(3:end,:);
fileNum = size(name, 1) / 2;
listName = name(1:fileNum,:);
timeName = name(fileNum + 1:end,:);
IBDEvent = cell(fileNum, 5);
muonEvent = {};
muonFamily = {};
timeGap = [];
fCalibD = @(ADC) ADC .* 1.734 - 194.3; % 刻度探测器 E(keV) = 1.734 * ADC - 194.3
% fCalibS = @(ADC) ADC .* 1.702 - 281.6; % 刻度模拟信号发生器 E(keV) = 1.702 * ADC - 281.6 % 1900
fCalibS = @(ADC) ADC .* 2.015 - 323.3; % 刻度模拟信号发生器 E(keV) = 2.015 * ADC - 323.3 % 1600

for ii = 1:fileNum
% parfor ii = 1:fileNum
    listFileName = [dirName, '\', listName(ii,:)];
    timeFileName = [dirName, '\', timeName(ii,:)];
    name = listFileName(10:end);
    name = strrep(name, '.txt', '');
    while contains(name, ' ')
        name = strrep(name, ' ', '');
    end
    name = name(end-11:end);
    event = Preprocess(listFileName, timeFileName, fCalibD);
    ParSave(['event_corrected\event_', name], event, '-v7.3');
    % [muonE, muonF, timeG] = MuonSelection(event);
    % ParSave(['muon\muon', name], muonE, '-v7.3');
    % timeGap = [timeGap; timeG];
    % muonEvent = [muonEvent; muonE];
    % muonFamily{ii, 1} = muonF;
    % temp = IBDSelection(event);
    % IBDEvent(ii,:) = {name, temp, size(temp, 1)};
    % IBDEvent(ii,:) = {name, temp, size(temp, 1), size(event,1), 100 .* size(temp, 1) ./2 ./ 9000};% 最后一列为甄别效率（模拟）
end

% save('IBDEvent_P25D28.mat', 'IBDEvent');
% save('MuonEvent_10MeV-v7.mat', 'muonEvent', '-v7');
% save('MuonGenerated.mat', 'muonGenerated');
% save('MuonFamily.mat', 'muonFamily');
%%
% t = cell2mat(muonGenerated(:,2));
% tGap = t(2:end)-t(1:end-1);

% muonEvent1 = muonEvent(1:length(muonEvent) / 10,:);

