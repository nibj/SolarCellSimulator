%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sim] = DesignJunctionCIGSRunSim()
%JUNCTION Create a junction object.


sim.input = struct(...
    'Linear',0,...
    'nmLx', 500, ...     % Grating period
    'nmLg',100,...
    'NSun',1,...
    'VVext', 0.0, ...       % External voltage
    'icm3isG', -inf ...
    ...% Generation Rate (value, function,
    ...%file, or -inf for results from RCWA)
    );

% Design the junction
sim.material = struct(...
    'nsec', 9, ... % Number of independent
    ...%section in study
    'elecsec', [0,0,0,1,1,1,0,0,0], ...
    ...% Flags for which sections to include in the
    ...%electrical simulation
    'periodicsec', [0,0,0,0,0,0,0,0,0], ...
    ...% grating section
    'nmSec', [100,110,100,80,70,730,40,10,500], ...
    ...% Lengths of the sections
    'icm3ND', [0,0,0,1e17,5e17,-2e16,0,0,0], ...
    ...% Doping of the sections
    'icm3Nf0', [0,0,0,1e16,5e17,0,0,0,0], ...
    'icm3Nf1',[0,0,0,0,0,0,0,0,0], ...
    ...% Doping of the sections
    ... % Material properties will be
    ...%overwritten if 'material'
    ...%not set to 0. If material has a range
    ...%of possible bandgaps,
    ...%will tailor it to specified Eg.
    'material','Air&MgF2&AZO&ZnO&CdS&CIGS&Al2O3&Al2O3&Mo',...
    ...%% Loads a material from the database, 0 for
    ...%user input
    'VEg0', [0,0,0,3.3,2.4,1.62,0,0,0], ...
    ...% Baseline bandgap
    'VEg1', [0,0,0,0,0,0,0,0,0], ...   % Bandgap perturbation
    ... % Electrical Properties
    'icm3Nc0', [0,0,0,3e18,1.3e18,6.8e17,0,0,0], ...
    ...%% Conduction band DOS
    'icm3Nv0', [0,0,0,1.7e19, 9.1e19,1.5e19,0,0,0], ...
    ...%% Valence band DOS
    'icm3Nc1', 0, ... % Conduction band DOS
    'icm3Nv1', 0, ... % Valence band DOS
    ... % Electrical Properties
    'cm2iVismun0',[0,0,0,100,72,100,0,0,0], ...
    ...% Electron mobility
    'cm2iVismup0',  [0,0,0,31,20,13,0,0,0], ...
    ...% Hole mobility
    'cm2iVismun1', [0,0,0,0,0,0,0,0,0], ...  % Electron mobility
    'cm2iVismup1', [0,0,0,0,0,0,0,0,0], ...  % Hole mobility
    ...% Electrical Properties
    'VChi0', [0,0,0,4.4,4.2,4.2,0,0,0],...
    ...% Electron Affinity
    'VChi1', [0,0,0,0,0,0,0,0,0],...  % Electron Affinity
    ... % 'VChi1',0,...  % Electron Affinity
    'epsdc', [0,0,0,9,5.4,13.6,0,0,0], ...
    ...% DC relative permittivity of bulk silicon
    ... % Recombination Parameters
    'icm3isalpha', 1e-10, ...
    ...% It is infact R_B (cm^3 s^{-1})
    ...
    'istausrhn0',[0,0,0,2e-11,4e-13,0,0,0,0], ...
    ...% SRH recombination constant 1/s
    'istausrhp0',[0,0,0,1e-8, 2e-10,0,0,0,0], ...
    ... % SRH recombination constant 1/s
    'istausrhn1',[0,0,0,0,0,0,0,0,0], ... % SRH recombination constant 1/s
    'istausrhp1',[0,0,0,0,0,0,0,0,0], ... % SRH recombination constant 1/s
    ...
    'Cn', 0e-30, ... % Auger recombination constant
    'Cp', 0e-30, ... % Auger recombination constant
    'ET', [0,0,0,0.5,0.5,0.5,0,0,0], ...
    ...% Relative position of trap level in bandgap
    ... % Optical Properties
    'eps', [1,1,1,1,1,1,1,1,1] ... % Optical Permittivity
    );

sim.optical.periodic = struct(...
    ... % Takes a vector of inputs. One entry
    ...% per periodic region specified above
    'zr_flag', [1], ...
    ...%% Flag to toggle between 1 zeta, and 2 relief
    'zeta', 0.5, ...  % Graing Duty Cycle -
    ...%input from 0 (bottom) to 1 (top), output 0 to 1.
    'relief', '@(zeta) cos(2*pi*zeta)', ...
    ...% Function describing the surface relief,input
    ...%-1/2 to 1/2 (I think), output trimmed to [0,1]
    'material2', 'Mo', ...  % Grating material
    'eps2', -10.4356+0.8294i ...
    ...% Grating metal permittivity
    );

% Constants
sim.phys.Cq = 1.6021766e-19;
sim.phys.Jshbar = 1.0545718e-34;
sim.phys.VVth = 300 * 8.61733034e-5;
sim.phys.CiVimeps0 = 8.85418782e-12;
% converted im to icm
sim.phys.CiVicmeps0 = sim.phys.CiVimeps0 / 100;
% converted im to icm

% Units needed here
sim.phys.mkgis2iA2mu0=4*pi*10^(-7);
% permeability of free space
sim.phys.Oeta0 = 376.73031346177;
% Impedence of free space
sim.phys.misc = 299792458;
% speed of light
sim.phys.m2kgish = 6.6260696e-34;
% Planck constant
sim.phys.Enormconst = ...
    sqrt(1/(sim.phys.CiVimeps0 * sim.phys.misc));
% sqrt(2*sim.phys.Oeta0); % Scale E by this
%and H by 1/this
%such that incident power density is 1 W/m^2

% Scaling (These can be changed)
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

%%%%%%%%%%%%%%%%%%%%%%
% Do not edit below this box %
%%%%%%%%%%%%%%%%%%%%%%
nsec = sim.material.nsec;
if length(sim.material.nmSec) ~= nsec
    error(['Lengths not defined for each'...
        'section of junction']);
end

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

end

