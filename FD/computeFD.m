function sim=computeFD(sim)
% initialize data for FD and run
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

    Vext = 0/Vth; %initial Vext = 0

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
    %'phin0','phinN','phip0','phipN','nbar','Jscale','npscale','Gscale','Vth')
    Run; JVcurve; 



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

% MakePlot;