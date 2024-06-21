function [ sim ] = SetInitialConditions( sim)
% Set the initial conditions for the simulation

damping = sim.setup.damping;

% Initialise rho fixed
ND = damping * sim.sbp.ND;
%sim.coeffs.ND = sbp2coeffs(sim.sbp.ND, sim);

% Calculate initial values for n and p
sim.sbp.n_te = 0*ND;
sim.sbp.p_te = 0*ND;

for i = 1:sim.setup.nx
    ni2 = polyval(sim.sbp.ni2(:,i),0);
    ND2 = polyval(damping^2 * sim.sbp.ND2(:,i),0);
    NDh = polyval(ND(:,i),0);
    if NDh>0
        sim.sbp.n_te(end,i) = (NDh +  sqrt(ND2 + 4*ni2))/2;
        sim.sbp.p_te(end,i) = ni2/2/sim.sbp.n_te(end,i);
    else
        sim.sbp.p_te(end,i) = -(NDh -  sqrt(ND2 + 4*ni2))/2;
        sim.sbp.n_te(end,i) = ni2/2/sim.sbp.p_te(end,i);
    end
end
figure(100)
semilogy(sim.sbp.p_te(end,:))
hold on
semilogy(sim.sbp.n_te(end,:))
title('Initial n and p')
drawnow

sim = BCond(sim);

sim.setup.phi_bi = sim.bc.phi1 - sim.bc.phi0;

% Calculate the initial field as linear between phi0 and phi1
sim.setup.dphi = sim.setup.phi_bi / sim.setup.Lz;
phi_te = @(x) sim.bc.phi0 + sim.setup.dphi * x * sim.setup.Lz;
sim.sbp.phi_te = func2sbp(phi_te, sim);

% Set coefficient lists and hat variables
sim.coeffs.n = sbp2coeffs(sim.sbp.n_te, sim);
sim.coeffs.n_hat = node_av_global(sim.sbp.n_te);

sim.coeffs.p = sbp2coeffs(sim.sbp.p_te, sim);
sim.coeffs.p_hat = node_av_global(sim.sbp.p_te);

sim.coeffs.phi = sbp2coeffs(sim.sbp.phi_te, sim);
sim.coeffs.phi_hat = node_av_global(sim.sbp.phi_te);

% Find intial efield and phi
phi0 = func2sbp(@(x) sim.setup.dphi + 0*x, sim);
E0 = zeros(sim.setup.pdeg1, sim.setup.nx);
for i = 1:sim.setup.nx
    product = polytimes(phi0(:,i),sim.sbp.lambda2(:,i));
    E0(:,i) =  - product(end-sim.setup.pdeg:end);
end

sim.coeffs.E = sbp2coeffs(E0, sim);

% Find initial currents
Jn0 = zeros(sim.setup.pdeg1, sim.setup.nx);
for i = 1:sim.setup.nx
    product = polytimes(sim.sbp.mun(:,i), polytimes(phi0(:,i),sim.sbp.n_te(:,i)));
    Jn0(:,i) =  product(end-sim.setup.pdeg:end);
end
sim.coeffs.Jn = sbp2coeffs(Jn0, sim);

Jp0 = zeros(sim.setup.pdeg1, sim.setup.nx);
for i = 1:sim.setup.nx
    product = polytimes(sim.sbp.mup(:,i), polytimes(phi0(:,i),sim.sbp.p_te(:,i)));
    Jp0(:,i) =  product(end-sim.setup.pdeg:end);
end
sim.coeffs.Jp = sbp2coeffs(Jp0, sim);

end
