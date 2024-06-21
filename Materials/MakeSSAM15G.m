clear();
clc;
data = xlsread('./DataFiles/SSAM15GBlueshGreenf05.xlsx');

data = [data(:,1), data(:,2)];

save('./AM15GBlueshGreenf05.mat')