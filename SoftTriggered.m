clear

dirName = 'event_soft_triggered';
cd(dirName);
name = ls;
cd ..
name = name(3:end,:);

fileNum = size(name, 1);
stNum = zeros(fileNum, 1);
for ii = 1:fileNum
    fileName = [dirName, '\', name(ii,:)];
    while contains(fileName, ' ')
        fileName = strrep(fileName, ' ', '');
    end
    temp = load(fileName);
    event = temp.event;
    stNum(ii, 1) = event.num;
end

IBDEvent = load('IBDEvent_P25D28_1st6_rest700.mat');
IBDEvent = IBDEvent.IBDEvent;
num = [stNum, cell2mat(IBDEvent(:,4))];
%%
rate = num ./ 3600;
meanRate = mean(rate);
stdRate = std(rate);
