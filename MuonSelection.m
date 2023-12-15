function [muonEvent, muonFamily, timeGap]= MuonSelectionV2(event)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          |          瞬时信号         |           延迟信号
% 甄别项目 |          甄别条件         |           甄别条件
% 预甄别   |      Eth > 200 keV        |        Eth > 200 keV
% Ntrigger |     2 ≤ Ntrigger ≤ 16     |    2 ≤ Ntrigger ≤ 16
% Etotal   | 2.5 MeV < Etotal ≤ 7 MeV  | 2.8 MeV < Etotal ≤ 8 MeV
% E1st     |  1 MeV < E1st ≤ 6.5 MeV   |  0.5 MeV < E1st ≤ 7 MeV
% E2nd     |       E2nd < 520 keV      |       E2nd ≤ 3 MeV
% E3rd     |             --            |       E3rd ≤ 2 MeV
% E4th     |             --            |       E4th ≤ 1 MeV
% -------------------------------------------------------------------------
% ΔT       |                    8 μs < ΔT < 300 μs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V2更新内容
% 1、添加了μ子事件后续事件的甄别（依据延迟信号，即中子俘获信号的甄别条件）。
% 2、添加了μ子与μ生事件时间间隔计算并输出。
% 3、更新了MuonFamily事件保存结构。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edep = cell2mat(event(:,3));
edep = ReshapeDataMatrix(4, edep);
edep(edep < 200) = 0;
muonG = (permute(sum(sum(edep > 0)), [3, 1, 2]) < 2);
event(muonG,:) = []; % 去除触发数量小于2的事件

minT = 8;
maxT = 300;

edep = cell2mat(event(:,3));
edep = ReshapeDataMatrix(4, edep);
% muonEvent = (permute(sum(sum(edep > 0)), [3, 1, 2]) > 2);
totalEdep = permute(sum(sum(edep)), [3, 1, 2]);
muon = (totalEdep > 10000); % 总能量10MeV以上认为是缪子事件
muonEvent = event(muon,:);
time = cell2mat(event(:,2));
timeMuon = time(muon,:);
timeGap = [];
muonFamily = {};
ii = 1;
while true
    mii = [false(ii, 1); muon(1:end - ii)]; % μ子后第ii个事件（逻辑数组）
    muonEventii = event(mii,:); % μ子后第ii个事件
    timeii = time(mii,:);
    tGap = timeii - timeMuon(1:sum(mii),:); % μ子后第ii个事件与μ子事件的时间差
    mG = (tGap < maxT) & (tGap > minT); % μ子后第ii个事件与μ子事件的时间差在ΔT之内的事件（逻辑数组）
    if sum(mG) <= 0
        break;
    end
    muonE = muonEvent(mG,:); % 产生后续事件的μ子
    muonGeneratedii = muonEventii(mG,:); % μ子后第ii个μ生事件
    muonF = [muonE, muonGeneratedii];

    edep = cell2mat(muonF(:,6));
    edep = reshape(edep',[16, numel(edep) ./ 16])';
    edepSorted = sort(edep, 2, 'descend');
    totalEdep = sum(edepSorted, 2);
    dropedEvent = (totalEdep > 8000 | totalEdep < 2800 | edepSorted(:,1) > 7000 ...
        | edepSorted(:,2) > 3000 | edepSorted(:,3) > 2000 | edepSorted(:,4) > 1000);
    % 去除2.8MeV > Etotoal > 8MeV、
    % E1st > 7MeV、E2nd > 3MeV、E3rd > 2 MeV、E4th > 1MeV的μ子后第ii个事件
    muonF(dropedEvent,:) = [];
    timeG = tGap(mG,:);
    timeG(dropedEvent,:) = [];
    if ~isempty(muonF)
        muonFamily{1, ii} = muonF;
        timeGap = [timeGap; timeG]; % μ子与后续μ生事件的时间差
    end
    ii = ii + 1;
end

end
