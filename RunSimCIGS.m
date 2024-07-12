%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear();
clc;
close all;
addpath(genpath('./DDNewton'))
addpath(genpath('./Models'))
addpath(genpath('./Conversions'))
addpath(genpath('./DEA'))
addpath(genpath('./RCWA'))
addpath(genpath('./Materials'))
addpath(genpath('./FD'))
% Create the simulations
sim = DesignJunctionCIGSRunSim();
sim = DesignSimCIGSRunSim(sim);
CSun=1;

%ls=input("Please enter thickness of CIGS layer between 100 nm and 2200 nm: Ls= ");
ls=200
if ls >= 100 && ls <= 2200
    sim.material.nmSec(6)=ls;
else
    disp('Please Re-Run the code with 100<=Ls<=2200 (nm)')
    return
end

%Eg0=input("Please enter the bandgap Eg [0.947 1.626] (eV): ");
Eg0=1.1
if Eg0 >=0.947  && Eg0 <= 1.626
    sim.material.VEg0(6) = Eg0;
else
    disp('Please Re-Run the code with bandgap 0.947<=Eg<=1.626 (eV): ')
    return
end

% Since bandgap is a user input, we need to compute bandgap dependent parameters
% in the model
xxeg0=sim.material.VEg0(6);
xeg0=abs((xxeg0-0.947))/0.679;
%**** disp(xeg0);
xchi0=4.5-0.6*xeg0;
%**** disp(xchi0);
sim.material.VChi0(6) = xchi0;
Egrange=[Eg0,Eg0];

Nf1=(xeg0<=0.3).*1e11.*(5-15.33.*xeg0)+...
    (xeg0>0.3).*(1e11.*(153.83-1007.57.*xeg0+1654.74.*xeg0.^2));

sim.material.icm3Nf0(6)=-Nf1;
taun1=1/(Nf1*1e7*5e-13);
sim.material.istausrhn0(6)=taun1;

taup1=1/(Nf1*1e7*1e-15);
sim.material.istausrhp0(6)=taup1;

% also update overall thickness
sim.optical.setup.nslices(6) = ceil(ls);
% disp(sim.optical.setup.nslices(6));
% sim.optical.setup.Nx = ceil(params.nmLx/3);
sim.electrical.setup.nx_sec(3)=ceil(ls);


disp('Starting Optical Simulation')
if sim.optical.setup.optical_toggle
    tic
    % Copy input to optical model
    sim = BuildOptInput(sim);
    % Initialise optical model next line moved to BuildOptInput
    % sim.optical = RCWAinitilisation(sim.optical);
    % Run optical model
    sim.optical = RunOptical(sim.optical);
    % Transfer results to electrical model
    sim.input.icm3isG = sim.optical.results.icm3isG;
    optical_toc=toc;
end
disp(['Time for optical solve= ',num2str(optical_toc)])
disp('Optical done');

% After changes parameters, redefine constants to ensure
% correct values

sim.phys.nmLs = sim.material.nmSec *...
    sim.material.elecsec.'; % in nm
sim.phys.cmLs = cm_from_nm(sim.phys.nmLs); % in nm
sim.phys.icm3Ns = ...
    max([max(abs(sim.material.icm3ND)),1e14]);
% in cm^-3
sim.phys.VVs = sim.phys.VVth; % in V
sim.phys.cm2iVismus = 1.5; % in cm^2V^-1s^-1
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
% Other interal scaling factors
sim.phys.cmLs = cm_from_nm(sim.phys.nmLs);
% Ls in cm
sim.phys.cm2isDs = sim.phys.cm2iVismus *...
    sim.phys.VVs; %in cm^2s^-1
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
% Dimensionless (needs multiplying by eps)
sim.phys.Qscale = sim.phys.Oeta0 * ...
    sim.phys.CiVimeps0/(sim.phys.Jshbar * ...
    sim.phys.Enormconst^2);

sim.electrical.setup.start_electrical = tic;
sim = BuildElecInput(sim);% Function call

if sim.electrical.setup.FD == 0 %use HDG

    if sim.electrical.setup.electrical_toggle
        sim.electrical = RunElectrical(sim.electrical);
    end
    sim.electrical.results.total_electrical = ...
        toc(sim.electrical.setup.start_electrical);
    disp('Results computed by HDG')
else %use FD

    tic; sim=computeFD(sim); total_time_FD= toc; %runs the FD method

    %total time for electrical computation
    sim.electrical.results.total_electrical =total_time_FD;

    disp('Results computed by FD')
    disp(['Time for FD solve=',num2str(total_time_FD)])
end
eta = sim.electrical.results.Wim2Pmax/10;
Jsc = sim.electrical.results.mAicm2Jsc;
JscOpt = sim.electrical.results.JscOpt;
Voc = sim.electrical.results.VVOC;
FF = sim.electrical.results.FF;
Vmax = sim.electrical.results.VVmax;
disp('Input parameters:')
disp(['L_s=',num2str(ls),' , Eg0=',num2str(Eg0)])
disp('Results:')
disp(['eta = ',num2str(eta),'; Jsc = ',num2str(Jsc),'; Voc = ',...
    num2str(Voc),'; FF = ', num2str(FF),...
    '; Vmax = ',num2str(Vmax)]);

disp('Electrical done');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%