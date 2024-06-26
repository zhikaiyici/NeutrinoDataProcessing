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
% for ii = 1:fileNum
parfor ii = 1:fileNum
    fileName = [dirName, '\', name(ii,:)];
    while contains(fileName, ' ')
        fileName = strrep(fileName, ' ', '');
    end
    temp = load(fileName);
    event = temp.event;
    % [muonE, muonF, timeG] = MuonSelection(event);
    % timeGap = [timeGap; timeG];
    % % muonEvent = [muonEvent; muonE];
    % muonFamily{ii, 1} = muonF;
    % ParSave(['muon_corrected\muon_', fileName(end-15:end-4)], muonE, '-v7.3');
    temp = IBDSelection(event);
    IBDNum = size(temp.ID, 1);
    eventNum = size(event.ID, 1);
    IBDEvent(ii,:) ={fileName(end-15:end-4), temp, IBDNum, eventNum, IBDNum ./ eventNum};
    % 'soft triggered' means selecting events only in terms of energy
    % temp.eventNum = eventNum;
    % temp.num = IBDNum;
    % ParSave(['event_soft_triggered\event_st_', fileName(end-15:end-4)], temp, '-v7.3');
end
save("IBDEvent_P25D28_1st6_rest700.mat", "IBDEvent", '-v7');
% save("MuonFamily_1st6.mat", "muonFamily", '-v7');
% save("MuonNeutronTimeGap_1st6.mat", "timeGap", '-v7');

