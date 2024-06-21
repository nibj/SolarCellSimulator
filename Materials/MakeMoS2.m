
clear();
eps = xlsread('./Materials/DataFiles/epsMoSe2.xlsx');
VEg=[1.13];
VChi=[4.2];
icm3Nc = [7e17];
icm3Nv = [4.5e18];
cm2iVismun = [20];
cm2iVismup = [5];
epsdc = [4];
eps = [eps(:,1), (eps(:,2) + 1i*eps(:,3)).^2];
alpha = [1e-10];
istausrhn = [1e-11];
istausrhp = [1e-11];
save('./Materials/Semiconductors/MoSe2.mat')