%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear();
clc;
close all;
addpath(genpath('./DDNewton'))
addpath(genpath('./Models'))
addpath(genpath('./Conversions'))
%addpath(genpath('./DEA'))
addpath(genpath('./RCWA'))
addpath(genpath('./Materials'))
addpath(genpath('./FD'))
% Create the simulations
sim = DesignJunctionCIGSRunSim();
sim = DesignSimCIGSRunSim(sim);

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

% Since parameters length ls and bandgap Eg0 have changed we now need to update the
% entries in the simulation data including non-dimensionalization
% ls and Eg0 are already stored in sim
sim=CIGS_update(sim);


%% Start optical simulation
disp('Starting Optical Simulation')
if sim.optical.setup.optical_toggle
    tic
    % Copy input to optical model
    sim = BuildOptInput(sim);
    % Run optical model
    sim.optical = RunOptical(sim.optical);
    % Transfer results to electrical model
    sim.input.icm3isG = sim.optical.results.icm3isG;
    optical_toc=toc;
end
disp(['Time for optical solve= ',num2str(optical_toc)])
disp('Optical done');

%% start electrical simulation
sim.electrical.setup.start_electrical = tic;
sim = BuildElecInput(sim);

if sim.electrical.setup.FD == 0 %use HDG
    if sim.electrical.setup.electrical_toggle
        sim.electrical = RunElectrical(sim.electrical);
    end
    sim.electrical.results.total_electrical = ...
        toc(sim.electrical.setup.start_electrical);
    disp('Results computed by HDG')
else %use FD
    tic; sim=computeFD(sim); total_time_FD= toc; %runs the FD method
    sim.electrical.results.total_electrical =total_time_FD;
    disp('Results computed by FD')
    disp(['Time for FD solve=',num2str(total_time_FD)])
end
disp('Electrical done');

%% Unpack and report results
eta = sim.electrical.results.Wim2Pmax/10; % converts from m to cm
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%