function event = Preprocess(listFileName, timeFileName, fCalib)
% 去除编号重复的无效事件，取Time和List数据的交集
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V2更新内容
% 减少cell2mat(num2cell, mat2cell)使用，提升运行速度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% edepData = importdata('LIST0.txt');
% edepData = importdata(listFileName);
% edepData = edepData.data;
edepData = importdata(listFileName);

% timeData = importdata('TIME0.txt');
% timeData = importdata(timeFileName);
% timeData = timeData.data;
timeData = importdata(timeFileName);

time = timeData(2:2:end);
timeID = timeData(1:2:end);
clear timeData;
[~, uniqueTimeIndex] = unique(timeID,'stable');
duplicateTimeIndex = setdiff((1:size(timeID, 1))',uniqueTimeIndex);
clear uniqueTimeIndex;

time = [timeID, time];
clear timeID;
time(duplicateTimeIndex,:) = [];

edep = edepData;
edep(1:5:end,:) = [];
edep = sqrt(edep(:,1:2:end) .* edep(:,2:2:end));
edep = fCalib(edep);
edep(edep < 200) = 0;
edepID = edepData(1:5:end, 1);
clear edepData;

edep = ReshapeDataMatrix(4,edep);

[~, uniqueEdepIndex] = unique(edepID,'stable');
duplicateEdepIndex = setdiff((1:size(edepID, 1))',uniqueEdepIndex);
clear uniqueEdepIndex;

edepID(duplicateEdepIndex,:) = [];
edep(:,:,duplicateEdepIndex) = [];

eventID = intersect(time(:,1), edepID, 'stable');
[~, dropedTimeIndex] = setdiff(time(:,1), eventID);
[~, dropedEdepIndex] = setdiff(edepID, eventID);
clear eventID;
clear edepID;

time(dropedTimeIndex, :) = [];
edep(:,:,dropedEdepIndex) = [];

edep = permute(edep, [2,1,3]);
edep = reshape(edep, 4, 4 .* size(edep, 3))';
edep = mat2cell(edep, 4 .* ones(size(edep,1) ./ 4, 1), size(edep,2));
event = [num2cell(time), edep];

% name = listFileName(10:end);
% name = strrep(name, '.txt', '');
% while contains(name, ' ')
%     name = strrep(name, ' ', '');
% end
% 
% clearvars -except event name;
% save(['event', name, '.mat'], 'event');

end
