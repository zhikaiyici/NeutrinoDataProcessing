function IBDEvent = IBDSelection(event)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Prompt              |           Delayed
%          Eth > 200 keV           |        Eth > 200 keV
%         2 ≤ Ntrigger ≤ 16        |    2 ≤ Ntrigger ≤ 16
%     2.5 MeV < Etotal ≤ 7 MeV     | 2.8 MeV < Etotal ≤ 8 MeV
%      1 MeV < E1st ≤ 6 MeV        |  0.5 MeV < E1st ≤ 6 MeV
%         E2nd < 520 keV           |       E2nd ≤ 3 MeV
% Etotal – (E1st + E2nd) < 700 keV |       E3rd ≤ 2 MeV
%                 --               |       E4th ≤ 1 MeV
% ---------------------------------------------------------------
%                          8 μs < ΔT < 300 μs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 甄别顺序：
% 1、200keV预甄别
% 2、触发数量甄别
% 4、u子及后续事件（300us内）剔除
% 5、瞬时信号和延迟信号能量甄别条件交集甄别
% 6、时间甄别
% 7、瞬时信号和延迟信号能量甄别条件分别甄别
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V3更新内容
% event由cell变为struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V2更新内容
% 1、更新了缪子后续事件判别逻辑。由初版的时间序列错位比较改为时间间隔与甄别
%    时间比较；更改循环逻辑，由初版固定10次错位比较改为不固定次数比较，对应
%    时间内无事件则循环停止，提升了运行速度。
% 2、调整了甄别顺序。先进行瞬时信号与延迟信号能量甄别交集条件的甄别，再进行
%    时间甄别，最后是瞬时信号与延迟信号特有能量甄别条件的甄别。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% edep = cell2mat(event(:,3));
% edep = ReshapeDataMatrix(4, edep);
edep = event.edep;
edep(edep < 200) = 0;
dropedEvent = (permute(sum(sum(edep > 0)), [3, 1, 2]) < 2);
% event(dropedEvent,:) = []; % 去除触发数量小于2的事件
event.ID(dropedEvent,:) = [];
event.time(dropedEvent,:) = [];
event.edep(:,:,dropedEvent) = [];


minT = 8;
maxT = 300;

% edep = cell2mat(event(:,3));
% edep = ReshapeDataMatrix(4, edep);
edep = event.edep;
% muonEvent = (permute(sum(sum(edep > 0)), [3, 1, 2]) > 2);
totalEdep = permute(sum(sum(edep)), [3, 1, 2]);
muon = (totalEdep > 10000);
% muonEvent = event(muon,:);
% time = cell2mat(event(:,2));
time = event.time;
dropedEvent = false(size(time));
tempMuon = false(size(muon));
timeGap = 1;
ii = 1;
while sum(timeGap < maxT & timeGap > 0) > 0
% for ii = 1:10
    tempMuon = [false(ii, 1); muon(1:end - ii)] | tempMuon;
    timeGap = [500 .* ones(ii,1); time(1 + ii:end) - time(1:end - ii)];
    dropedEvent = dropedEvent | (timeGap < maxT);
    dropedEvent = tempMuon & dropedEvent;
    ii = ii + 1;
end
% event(muon,:) = [];
% event(dropedEvent,:) = []; %去除u子事件之后300us内的事件
event.ID(dropedEvent,:) = [];
event.time(dropedEvent,:) = [];
event.edep(:,:,dropedEvent) = [];

% edep = cell2mat(event(:,3));
% edep = ReshapeDataMatrix(4, edep);
edep = event.edep;
totalEdep = permute(sum(sum(edep)), [3, 1, 2]);
dropedEvent = (totalEdep < 2500 | totalEdep > 8000);
% event(dropedEvent,:) = []; % 去除总能量小于2.5MeV、大于8MeV的事件
event.ID(dropedEvent,:) = [];
event.time(dropedEvent,:) = [];
event.edep(:,:,dropedEvent) = [];

% edep = cell2mat(event(:,3));
edep = event.edep;
edep = reshape(permute(edep, [2, 1, 3]), [4, size(edep, 3) * 4]);
edep = reshape(edep, [16, numel(edep) ./ 16])';
edepSorted = sort(edep, 2, 'descend');
dropedEvent = (edepSorted(:,1) < 500 | edepSorted(:,1) > 6000 |...
    edepSorted(:,2) > 3000 | edepSorted(:,3) > 2000 | edepSorted(:,4) > 1000);
% 去除0.5MeV > E1st > 6MeV、E2nd > 3MeV、E3rd > 2 MeV、E4th > 1MeV的事件
% event(dropedEvent,:) = []; 
event.ID(dropedEvent,:) = [];
event.time(dropedEvent,:) = [];
event.edep(:,:,dropedEvent) = [];

% time = cell2mat(event(:,2));
time = event.time;
timeGap = time(2:end) - time(1:end - 1);
dropedEvent = (timeGap > maxT | timeGap < minT);
index1 = [dropedEvent; dropedEvent(end)];
indes2 = [dropedEvent(1); dropedEvent];
dropedEvent = indes2 & index1;
% index = find(dropedEvent == 0);
% index = unique([index; index + 1]);
% dropedEvent = [dropedEvent; true(1)];
% dropedEvent(index) = 0;
% % event(dropedEvent,:) = []; % 去除时间差大于300us或小于8us的除相邻事件
event.ID(dropedEvent,:) = [];
event.time(dropedEvent,:) = [];
event.edep(:,:,dropedEvent) = [];

% time = cell2mat(event(:,2));
time = event.time;
timeGap = time(2:end) - time(1:end - 1);
prompt = (timeGap < 300);
% promptEvent = event(prompt,:);
promptEvent.ID = event.ID(prompt,:);
promptEvent.time = event.time(prompt,:);
promptEvent.edep = event.edep(:,:,prompt);
% edep = cell2mat(promptEvent(:,3));
edep = promptEvent.edep;
edep = reshape(permute(edep, [2, 1, 3]), [4, size(edep, 3) * 4]);
edep = reshape(edep, [16, numel(edep) ./ 16])';
edepSorted = sort(edep, 2, 'descend');
dropedEvent = (sum(edepSorted, 2) > 7000 | edepSorted(:,1) > 6000 ...
    | edepSorted(:,1) < 1000 | edepSorted(:,2) > 520);
% 去除Etotoal > 7MeV、1MeV > E1st > 6MeV、E2nd > 520keV的瞬时事件
% promptEvent(dropedEvent,:) = [];
promptEvent.ID(dropedEvent,:) = [];
promptEvent.time(dropedEvent,:) = [];
promptEvent.edep(:,:,dropedEvent) = [];

edep = promptEvent.edep;
edep = reshape(permute(edep, [2, 1, 3]), [4, size(edep, 3) * 4]);
edep = reshape(edep, [16, numel(edep) ./ 16])';
edepSorted = sort(edep, 2, 'descend');
dropedEvent = ((sum(edepSorted, 2) - edepSorted(:,1)- edepSorted(:,2)) > 700);
% 去除Etotoal - E1st - E2nd > 700keV的瞬时事件
% promptEvent(dropedEvent,:) = [];
promptEvent.ID(dropedEvent,:) = [];
promptEvent.time(dropedEvent,:) = [];
promptEvent.edep(:,:,dropedEvent) = [];

% % edep = cell2mat(promptEvent(:,3));
% edep = promptEvent.edep;
% edep = reshape(permute(edep, [2, 1, 3]), [4, size(edep, 3) * 4]);
% edep = reshape(edep, [16, numel(edep) ./ 16])';
% [e1st, iE1st] = max(edep, [], 2);
% edep(edep == e1st) = 0;
% [~, iE2nd] = max(edep, [], 2);
% [r1st, c1st] = ind2sub(4, iE1st);
% [r2nd, c2nd] = ind2sub(4, iE2nd);
% dropedEvent = (abs(r1st - r2nd) > 2 | abs(c1st - c2nd) > 2); 
% % promptEvent(dropedEvent,:) = []; % 去除E1st和E2nd不相邻的瞬时事件
% promptEvent.ID(dropedEvent,:) = [];
% promptEvent.time(dropedEvent,:) = [];
% promptEvent.edep(:,:,dropedEvent) = [];

% [~, reservedIndex] = intersect(cell2mat(event(:,1)), cell2mat(promptEvent(:,1)), 'stable');
[~, reservedIndex] = intersect(event.ID, promptEvent.ID, 'stable');
reservedEventIndex = unique([reservedIndex; reservedIndex + 1]);
% 保留Etotoal ≤ 7MeV、E1st ≤ 6MeV、E2nd ≤ 520keV的瞬时事件对应的事件对
% event = event(reservedEventIndex,:); 
event.ID = event.ID(reservedEventIndex,:);
event.time = event.time(reservedEventIndex,:);
event.edep = event.edep(:,:,reservedEventIndex);

% time = cell2mat(event(:,2));
time = event.time;
timeGap = time(2:end) - time(1:end - 1);
prompt = (timeGap < 300);
delayedEvent = event;
% delayedEvent(prompt,:) = [];
delayedEvent.ID(prompt,:) = [];
delayedEvent.time(prompt,:) = [];
delayedEvent.edep(:,:,prompt) = [];

% edep = cell2mat(delayedEvent(:,3));
edep = delayedEvent.edep;
edep = reshape(permute(edep, [2,1,3]), [4, size(edep, 3) * 4]);
edep = reshape(edep, [16, numel(edep) ./ 16])';
edepSorted = sort(edep, 2, 'descend');
dropedEvent = (sum(edepSorted, 2) < 2800);
% delayedEvent(dropedEvent,:) = []; % 去除Etotal < 2.8MeV的延迟事件
delayedEvent.ID(dropedEvent,:) = [];
delayedEvent.time(dropedEvent,:) = [];
delayedEvent.edep(:,:,dropedEvent) = [];

% [~, reservedIndex] = intersect(cell2mat(event(:,1)), cell2mat(delayedEvent(:,1)), 'stable');
[~, reservedIndex] = intersect(event.ID, delayedEvent.ID, 'stable');

reservedEventIndex = unique([reservedIndex - 1; reservedIndex]);
% IBDEvent = event(reservedEventIndex,:);
IBDEvent.ID = event.ID(reservedEventIndex,:);
IBDEvent.time = event.time(reservedEventIndex,:);
IBDEvent.edep = event.edep(:,:,reservedEventIndex);

end
