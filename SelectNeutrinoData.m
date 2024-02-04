clear
dirName = 'event_corrected';
offset = size(dirName, 2);
cd(dirName);
name = ls;
cd ..
name = name(3:end,:);
fileNum = size(name, 1);

IBDEvent = cell(fileNum, 5);
timeGap = [];
muonEvent = [];
muonFamily = cell(fileNum, 1);
for ii = 1:fileNum
% parfor ii = 1:fileNum
    fileName = [dirName, '\', name(ii,:)];
    while contains(fileName, ' ')
        fileName = strrep(fileName, ' ', '');
    end
    temp = load(fileName);
    event = temp.event;
    [muonE, muonF, timeG] = MuonSelection(event);
    timeGap = [timeGap; timeG];
    % muonEvent = [muonEvent; muonE];
    muonFamily{ii, 1} = muonF;
    % ParSave(['muon_corrected\muon_', fileName(end-15:end-4)], muonE, '-v7.3');
    temp = IBDSelection(event);
    IBDNum = size(temp.ID, 1);
    eventNum = size(event.ID, 1);
    IBDEvent(ii,:) ={fileName(end-15:end-4), temp, IBDNum, eventNum, IBDNum ./ eventNum};
end
% save("IBDEvent_P25D28-s.mat", "IBDEvent", '-v7');

