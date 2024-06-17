clear

dirName = 'event_soft_triggered';
cd(dirName);
name = ls;
cd ..
name = name(3:end,:);
fileNum = size(name, 1);
num = zeros(fileNum, 1);
for ii = 1:fileNum
    fileName = [dirName, '\', name(ii,:)];
    while contains(fileName, ' ')
        fileName = strrep(fileName, ' ', '');
    end
    temp = load(fileName);
    event = temp.event;
    num(ii, 1) = event.num;
end
%%
stRate = num ./ 3600;
meanSTRate = mean(stRate);
stdSTRate = std(stRate);
