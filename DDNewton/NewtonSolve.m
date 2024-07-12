function [newsim, abschange, relchange] = NewtonSolve(sim)

% Set up the sim.problem
newsim = sim;
nx = sim.setup.nx;
pdeg1 = sim.setup.pdeg1;
np = sim.setup.np;

% Create the iteration matrixes

%newsim.global.U = (tiprod(newsim.global.Psi, 3, newsim.coeffs.p).')*newsim.coeffs.n;

newsim.global.D_taun = zeros(np, np);
for i = 1:nx
    index1 = (i-1)*pdeg1+1;
    index2 = i*pdeg1;
    newsim.global.D_taun(index1:index2,index1:index2) = sim.global.taunr(i) * sim.local.DR - sim.global.taunl(i)* sim.local.DL;
end

newsim.global.D_taup = zeros(np, np);
for i = 1:nx
    index1 = (i-1)*pdeg1+1;
    index2 = i*pdeg1;
    newsim.global.D_taup(index1:index2,index1:index2) = sim.global.taupr(i) * sim.local.DR - sim.global.taupl(i)* sim.local.DL;
end

newsim.global.D_tauphi = zeros(np, np);
for i = 1:nx
    index1 = (i-1)*pdeg1+1;
    index2 = i*pdeg1;
    newsim.global.D_tauphi(index1:index2,index1:index2) = sim.global.tauphir(i) * sim.local.DR - sim.global.tauphil(i)* sim.local.DL;
end

% Calculate the jump at a node gamma  -tau_r @ + tau_l
newsim.global.S_taun_n = diag(sim.global.DeltaChi * sim.global.taunl(2:end) - sim.global.taunr(1:end-1));
newsim.global.S_taup_p = diag(sim.global.DeltaChiEg * sim.global.taupl(2:end) - sim.global.taupr(1:end-1));
newsim.global.S_taun = diag(sim.global.taunl(2:end) - sim.global.taunr(1:end-1));
newsim.global.S_taup = diag(sim.global.taupl(2:end) - sim.global.taupr(1:end-1));

newsim.global.S_tauphi_n = diag(sim.global.DeltaChi * sim.global.tauphil(2:end) - sim.global.tauphir(1:end-1));
newsim.global.S_tauphi_p = diag(sim.global.DeltaChiEg * sim.global.tauphil(2:end) - sim.global.tauphir(1:end-1));
newsim.global.S_tauphi = diag(sim.global.tauphil(2:end) - sim.global.tauphir(1:end-1));

% Jump at a node gamma
newsim.global.S_taunpsi_n = zeros(nx-1, np);
newsim.global.S_taunpsi_p = zeros(nx-1, np);
newsim.global.S_taunpsi = zeros(nx-1, np);
for i = 1:nx-1
    index1 = (i-1) * pdeg1 + 1;
    index2 = i * pdeg1;
    
    index3 = i * pdeg1 + 1;
    index4 = (i+1) * pdeg1;
    
    newsim.global.S_taunpsi_n(i,index1:index2) = - newsim.global.taunr(i) * newsim.local.SR;
    newsim.global.S_taunpsi_n(i,index3:index4) = sim.global.DeltaChi(i,i) * newsim.global.taunl(i+1) * newsim.local.SL;
    
    newsim.global.S_taunpsi_p(i,index1:index2) = - newsim.global.taunr(i) * newsim.local.SR;
    newsim.global.S_taunpsi_p(i,index3:index4) = sim.global.DeltaChiEg(i,i) * newsim.global.taunl(i+1) * newsim.local.SL;
    
     newsim.global.S_taunpsi(i,index1:index2) = - newsim.global.taunr(i) * newsim.local.SR;
    newsim.global.S_taunpsi(i,index3:index4) = newsim.global.taunl(i+1) * newsim.local.SL;
end



newsim.global.S_tauppsi_n = zeros(nx-1, np);
newsim.global.S_tauppsi_p = zeros(nx-1, np);
newsim.global.S_tauppsi = zeros(nx-1, np);
for i = 1:nx-1
    index1 = (i-1) * pdeg1 + 1;
    index2 = i * pdeg1;
    
    index3 = i * pdeg1 + 1;
    index4 = (i+1) * pdeg1;
    
    newsim.global.S_tauppsi_n(i,index1:index2) = - newsim.global.taupr(i) * newsim.local.SR;
    newsim.global.S_tauppsi_n(i,index3:index4) =  sim.global.DeltaChi(i,i) * newsim.global.taupl(i+1) * newsim.local.SL;
       
    newsim.global.S_tauppsi_p(i,index1:index2) = - newsim.global.taupr(i) * newsim.local.SR;
    newsim.global.S_tauppsi_p(i,index3:index4) =  sim.global.DeltaChiEg(i,i) * newsim.global.taupl(i+1) * newsim.local.SL;
       
    newsim.global.S_tauppsi(i,index1:index2) = - newsim.global.taupr(i) * newsim.local.SR;
    newsim.global.S_tauppsi(i,index3:index4) =  newsim.global.taupl(i+1) * newsim.local.SL;
end

newsim.global.S_tauphipsi_n = zeros(nx-1, np);
newsim.global.S_tauphipsi_p = zeros(nx-1, np);
newsim.global.S_tauphipsi = zeros(nx-1, np);
for i = 1:nx-1
    index1 = (i-1) * pdeg1 + 1;
    index2 = i * pdeg1;
    
    index3 = i * pdeg1 + 1;
    index4 = (i+1) * pdeg1;
    
    newsim.global.S_tauphipsi_n(i,index1:index2) = - newsim.global.tauphir(i) * newsim.local.SR;
    newsim.global.S_tauphipsi_n(i,index3:index4) = sim.global.DeltaChi(i,i) * newsim.global.tauphil(i+1) * newsim.local.SL;
        
    newsim.global.S_tauphipsi_p(i,index1:index2) = - newsim.global.tauphir(i) * newsim.local.SR;
    newsim.global.S_tauphipsi_p(i,index3:index4) = sim.global.DeltaChiEg(i,i) * newsim.global.tauphil(i+1) * newsim.local.SL;
        
    newsim.global.S_tauphipsi(i,index1:index2) = - newsim.global.tauphir(i) * newsim.local.SR;
    newsim.global.S_tauphipsi(i,index3:index4) = newsim.global.tauphil(i+1) * newsim.local.SL;
end


% Create simulation matrixes
newsim = Func(newsim);

[J, newsim] = Jacobian( newsim, 0);
newsim.global.J = J;

newsim.setup.counter = newsim.setup.counter + 1;

% if sim.setup.test_J
%     newsim = Test(newsim);
% end

sim=newsim;

x= [newsim.coeffs.n; newsim.coeffs.Jn; newsim.coeffs.n_hat; newsim.coeffs.p; newsim.coeffs.Jp; newsim.coeffs.p_hat; newsim.coeffs.phi; newsim.coeffs.E; newsim.coeffs.phi_hat];

minnp = 0.1*min([sim.bc.n0,sim.bc.n1,sim.bc.p0,sim.bc.p1]);

dx = -newsim.global.J\newsim.global.F;

xnew = x + dx;
newsim = globalcoeff2coeffs(xnew, newsim);

abschange = sum(abs(dx));
relchange = sum(abs(dx./(abs(x)+1)));
end

