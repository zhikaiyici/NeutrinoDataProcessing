function [muonEvent, muonFamily, timeGap]= MuonSelection(event)
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
% V2更新内容
% 1、添加了μ子事件后续事件的甄别（依据延迟信号，即中子俘获信号的甄别条件）。
% 2、添加了μ子与μ生事件时间间隔计算并输出。
% 3、更新了MuonFamily事件保存结构。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% edep = cell2mat(event(:,3));
% edep = ReshapeDataMatrix(4, edep);
edep = event.edep;
edep(edep < 200) = 0;
muonG = (permute(sum(sum(edep > 0)), [3, 1, 2]) < 2);
% event(muonG,:) = []; % 去除触发数量小于2的事件
event.ID(muonG,:) = [];
event.time(muonG,:) = [];
event.edep(:,:,muonG) = [];

minT = 8;
maxT = 300;

% edep = cell2mat(event(:,3));
% edep = ReshapeDataMatrix(4, edep);
edep = event.edep;
% muonEvent = (permute(sum(sum(edep > 0)), [3, 1, 2]) > 2);
totalEdep = permute(sum(sum(edep)), [3, 1, 2]);
muon = (totalEdep > 10000); % 总能量10MeV以上认为是缪子事件
% muonEvent = event(muon,:);
muonEvent.ID = event.ID(muon,:);
muonEvent.time = event.time(muon,:);
muonEvent.edep = event.edep(:,:,muon);
% time = cell2mat(event(:,2));
time = event.time;
timeMuon = time(muon,:);
timeGap = [];
muonFamily = [];
ii = 1;
while true
    mii = [false(ii, 1); muon(1:end - ii)]; % μ子后第ii个事件（逻辑数组）
    % muonEventii = event(mii,:); % μ子后第ii个事件
    muonEventii.ID = event.ID(mii,:);
    muonEventii.time = event.time(mii,:);
    muonEventii.edep = event.edep(:,:,mii);
    timeii = time(mii,:);
    tGap = timeii - timeMuon(1:sum(mii),:); % μ子后第ii个事件与μ子事件的时间差
    mG = (tGap < maxT) & (tGap > minT); % μ子后第ii个事件与μ子事件的时间差在ΔT之内的事件（逻辑数组）
    if sum(mG) <= 0
        break;
    end
    % muonE = muonEvent(mG,:); % 产生后续事件的μ子
    muonE.ID = muonEvent.ID(mG,:);
    muonE.time = muonEvent.time(mG,:);
    muonE.edep = muonEvent.edep(:,:,mG);
    % muonGeneratedii = muonEventii(mG,:); % μ子后第ii个μ生事件
    muonGeneratedii.ID = muonEventii.ID(mG,:);
    muonGeneratedii.time = muonEventii.time(mG,:);
    muonGeneratedii.edep = muonEventii.edep(:,:,mG);

    % edep = cell2mat(muonF(:,6));
    edep = muonGeneratedii.edep;
    edep = reshape(permute(edep, [2,1,3]),[4, size(edep, 3) * 4])';
    edep = reshape(edep',[16, numel(edep) ./ 16])';
    edepSorted = sort(edep, 2, 'descend');
    totalEdep = sum(edepSorted, 2);
    dropedEvent = (totalEdep > 8000 | totalEdep < 2800 | edepSorted(:,1) > 6000 ...
        | edepSorted(:,2) > 3000 | edepSorted(:,3) > 2000 | edepSorted(:,4) > 1000);
    % 去除2.8MeV > Etotoal > 8MeV、
    % E1st > 6MeV、E2nd > 3MeV、E3rd > 2 MeV、E4th > 1MeV的μ子后第ii个事件
    % muonF(dropedEvent,:) = [];
    muonE.ID(dropedEvent,:) = [];
    muonE.time(dropedEvent,:) = [];
    muonE.edep(:,:,dropedEvent) = [];
    muonGeneratedii.ID(dropedEvent,:) = [];
    muonGeneratedii.time(dropedEvent,:) = [];
    muonGeneratedii.edep(:,:,dropedEvent) = [];
    muonF = [muonE, muonGeneratedii];
    timeG = tGap(mG,:);
    timeG(dropedEvent,:) = [];
    % if ~isempty(muonF)
    if ~isempty(muonF(1).ID)
        muonFamily = [muonFamily; muonF];
        timeGap = [timeGap; timeG]; % μ子与后续μ生事件的时间差
    end
    ii = ii + 1;
end

end
