%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear();
VEg=[0.92];%updated in OptoElec
VChi=[0];
icm3Nc = [7.8e17];
icm3Nv = [4.5e18];
cm2iVismun = [40];
cm2iVismup = [12.6];
epsdc = [14];
eps = '@(nmlambda, VEg) CIGs(nmlambda, VEg)';
alpha = [1e-10];
istausrhn = [1e-9];%updated in OptoElec
istausrhp = [1e-6];%updated in OptoElec

save('./Materials/Semiconductors/CIGS.mat')