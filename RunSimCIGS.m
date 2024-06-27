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

ls=input("Please enter thickness of CIGS layer between 100 nm and 2200 nm: Ls= ");
if ls >= 100 && ls <= 2200
    sim.material.nmSec(6)=ls;
else
    disp('Please Re-Run the code with 100<=Ls<=2200 (nm)')
    return
end

Eg0=input("Please enter the bandgap Eg [0.947 1.626] (eV): ");
if Eg0 >=0.947  && Eg0 <= 1.626
    sim.material.VEg0(6) = Eg0;
else
    disp('Please Re-Run the code with bandgap 0.947<=Eg<=1.626 (eV): ')
    return
end


xxeg0=sim.material.VEg0(6);
xeg0=abs((xxeg0-0.947))/0.679;
%**** disp(xeg0);
xchi0=4.5-0.6*xeg0;
%**** disp(xchi0);
sim.material.VChi0(6) = xchi0;
Egrange=[Eg0,Eg0];

Nf1=(xeg0<=0.3).*1e11.*(5-15.33.*xeg0)+...
    (xeg0>0.3).*(1e11.*(153.83-1007.57.*xeg0+1654.74.*xeg0.^2));
%**** disp(Nf1);

sim.material.icm3Nf0(6)=-Nf1;
taun1=1/(Nf1*1e7*5e-13);
%*** disp(taun1);
sim.material.istausrhn0(6)=taun1;

taup1=1/(Nf1*1e7*1e-15);
%**** disp(taup1);
sim.material.istausrhp0(6)=taup1;

sim.optical.setup.nslices(6) = ceil(ls);
% disp(sim.optical.setup.nslices(6));
% sim.optical.setup.Nx = ceil(params.nmLx/3);
sim.electrical.setup.nx_sec(3)=ceil(ls);
disp('Starting Optical Simulation')
if sim.optical.setup.optical_toggle
    tic
    % Copy input to optical model
    sim = BuildOptInput(sim);
    % Initialise optical model
    sim.optical = RCWAinitilisation(sim.optical);
    % Run optical model
    sim.optical = RunOptical(sim.optical);
    % Transfer results to electrical model
    sim.input.icm3isG = sim.optical.results.icm3isG;
    optical_toc=toc;
end
disp('Optical done');

% After optoelec changes parameters, redefine constants to ensure
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


if sim.electrical.setup.FD == 1 %use the FD method instead of HDG
    %make sure to use nanometer scaling
    sim.electrical.setup.nx_sec = sim.electrical.material.nmSec; %h=1
    %make sure that the initial damping parameter = 1
    sim.electrical.setup.dampmin = 1;
    simben = initializeDevice(sim.electrical);
    simben = initializeSim(simben);
    simben = BCond(simben);
    simben = SetInitialConditions(simben);


    G = Dataout(simben.sbp.G,simben,0)';
    %disp(G);
    n1 = Dataout(simben.sbp.n1,simben,0)';
    p1 = Dataout(simben.sbp.p1,simben,0)';
    taun = Dataout(simben.sbp.tausrhn,simben,0)';
    taup = Dataout(simben.sbp.tausrhp,simben,0)';
    mun = Dataout(simben.sbp.mun,simben,0)';
    mup = Dataout(simben.sbp.mup,simben,0)';
    nbar = Dataout(simben.sbp.ni,simben,0)';
    Nf = Dataout(simben.sbp.ND,simben,0)';
    lambda = sqrt(Dataout(simben.sbp.lambda2,simben,0))';
    Eg = Dataout(simben.sbp.Eg,simben,0)';
    Chi = Dataout(simben.sbp.Chi,simben,0)';
    Nc = Dataout(simben.sbp.Nc,simben,0)';
    Nv = Dataout(simben.sbp.Nv,simben,0)';
    Vth = simben.phys.VVth;
    N0 = 1;
    alpha = Dataout(simben.sbp.alpha,simben,0)'; %radiative recombination
    Cn = Dataout(simben.sbp.Cn,simben,0)'; %Auger recombination
    Cp = Dataout(simben.sbp.Cp,simben,0)'; %Auger recombination
    phin = Chi;
    phip = Chi + Eg;
    %Diriclet boundary conditions for n
    n0 = simben.bc.n0;
    nN = simben.bc.n1;

    %Diriclet boundary conditions for p
    p0 = simben.bc.p0;
    pN = simben.bc.p1;

    %Diriclet boundary conditions for phi
    phi0 = simben.bc.phi0;
    phiN = simben.bc.phi1;

    %Dirichlet bondary conditions for phi_n
    phin0 = polyval(simben.sbp.Chi(:,1),-1);
    phinN = polyval(simben.sbp.Chi(:,end),1);

    %Dirichlet bondary conditions for phi_p
    phip0 = polyval(simben.sbp.Chi(:,1),-1) + polyval(simben.sbp.Eg(:,1),-1);
    phipN = polyval(simben.sbp.Chi(:,end),1) + polyval(simben.sbp.Eg(:,end),1);

    Jscale = simben.phys.mAicm2Js;
    npscale = simben.phys.icm3Ns;
    Gscale = simben.phys.icm3isGs;
    CSun=sim.input.NSun;
    %sectionLength = cumsum(simben.material.nmSec(1:end-1));
    %save('DataCIGS.mat','G','phin','phip','mun','mup','n1','p1','taun',...
    %'taup','Nf','lambda','phi0','phiN','n0','nN','p0','pN',...
    %'phin0','phinN','phip0','phipN','nbar','Jscale','npscale','Gscale','Vth');
else
end

if sim.electrical.setup.FD == 0 %use HDG

    if sim.electrical.setup.electrical_toggle
        sim.electrical = RunElectrical(sim.electrical);
    end
    sim.electrical.results.total_electrical = ...
        toc(sim.electrical.setup.start_electrical);

    eta = sim.electrical.results.Wim2Pmax/10;
    Jsc = sim.electrical.results.mAicm2Jsc;
    JscOpt = sim.electrical.results.JscOpt;
    Voc = sim.electrical.results.VVOC;
    FF = sim.electrical.results.FF;
    Vmax = sim.electrical.results.VVmax;



else %use FD

    tic; computeFD; total_time= toc; %runs the FD method


    sim.electrical.results.Wim2Pmax = eta*10; %max power density
    sim.electrical.results.mAicm2Jsc = Jvec(1); %short-circuit current

    if length(sim.input.icm3isG) == 1 %optical short-circuit current
        sim.electrical.results.JscOpt = 1000 * cm_from_nm(simben.setup.nmLz)...
            * sim.phys.Cq * sim.input.icm3isG;
    else
        GFn = str2func(sim.input.icm3isG);
        z = 0:0.01:1;
        sim.electrical.results.JscOpt = 1000 * cm_from_nm(simben.setup.nmLz)...
            * sim.phys.Cq * sum(GFn(z))/length(z);
    end

    %open-circuit voltage
    sim.electrical.results.VVOC = Vext(end-1);

    %fill-factor
    sim.electrical.results.FF = CSun*sim.electrical.results.Wim2Pmax/...
        (10*sim.electrical.results.VVOC*sim.electrical.results.mAicm2Jsc);

    %maximum power point voltage
    sim.electrical.results.VVmax = mPpoint(1);

    %total time for electrical computation
    sim.electrical.results.total_electrical =total_time;

    eta = sim.electrical.results.Wim2Pmax/10;
    Jsc = sim.electrical.results.mAicm2Jsc;
    JscOpt = sim.electrical.results.JscOpt;
    Voc = sim.electrical.results.VVOC;
    FF = sim.electrical.results.FF;
    Vmax = sim.electrical.results.VVmax;

    disp(['Time for optical solve= ',num2str(optical_toc),', Time for FD solve=',num2str(total_time)])
    disp(['L_s=',num2str(ls),' , Eg0=',num2str(Eg0)])
    disp(['eta = ',num2str(eta),'; Jsc = ',num2str(Jsc),'; Voc = ',...
        num2str(Voc),'; FF = ', num2str(FF),...
        '; Vmax = ',num2str(Vmax)]);
end
disp('Electrical done');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%