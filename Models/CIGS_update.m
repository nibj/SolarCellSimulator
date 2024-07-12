function sim = CIGS_update(sim)
% Update parameters since badgap Eg0 has changes
% These are experimental/theoretical laws
xxeg0=sim.material.VEg0(6);
xeg0=abs((xxeg0-0.947))/0.679;

xchi0=4.5-0.6*xeg0;
sim.material.VChi0(6) = xchi0;
%Egrange=[Eg0,Eg0];

Nf1=(xeg0<=0.3).*1e11.*(5-15.33.*xeg0)+...
    (xeg0>0.3).*(1e11.*(153.83-1007.57.*xeg0+1654.74.*xeg0.^2));

sim.material.icm3Nf0(6)=-Nf1;
taun1=1/(Nf1*1e7*5e-13);
sim.material.istausrhn0(6)=taun1;

taup1=1/(Nf1*1e7*1e-15);
sim.material.istausrhp0(6)=taup1;

% also update overall thickness
ls=sim.material.nmSec(6);
sim.optical.setup.nslices(6) = ceil(ls);
sim.electrical.setup.nx_sec(3)=ceil(ls);

% After changes to Ls, must recompute some parameters
% that depend on it (because of non-dimensionalization)
sim.phys.nmLs = sim.material.nmSec *...
    sim.material.elecsec.'; % in nm
sim.phys.cmLs = cm_from_nm(sim.phys.nmLs); % in nm
sim.phys.icm3isGs = ...
    sim.phys.cm2iVismus * sim.phys.VVs * ...
    sim.phys.icm3Ns/sim.phys.cmLs^2; %need 
sim.phys.cm3isalphas = (sim.phys.cm2iVismus * ...
     sim.phys.VVs)/(sim.phys.cmLs^2 * sim.phys.icm3Ns);
sim.phys.isTaus = ...
     sim.phys.cmLs^2/(sim.phys.cm2iVismus * ...
     sim.phys.VVs);
sim.phys.cm6isCs = (sim.phys.cm2iVismus * ...
     sim.phys.VVs)/(sim.phys.cmLs^2 * sim.phys.icm3Ns^2);
sim.phys.mAicm2Js = 1000 * sim.phys.Cq *... %need
     sim.phys.VVs *...
     sim.phys.cm2iVismus * ...
     sim.phys.icm3Ns/sim.phys.cmLs; % in mA cm^-2  %need
sim.phys.icm2isGRs = sim.phys.VVs *...
     sim.phys.cm2iVismus *...
     sim.phys.icm3Ns/sim.phys.cmLs^2; % in cm^-3s^-1
sim.phys.lambda2s = sim.phys.CiVicmeps0 * ...
     sim.phys.VVth/(sim.phys.cmLs^2 * sim.phys.Cq *...
     sim.phys.icm3Ns);
end