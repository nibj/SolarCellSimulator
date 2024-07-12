function sim = BuildOptInput(sim)
%
% intialize optical parameters for RCWA

sim.optical.input = sim.input;
sim.optical.phys = sim.phys;
sim.optical.material = sim.material;

% Initialise optical model
sim.optical = RCWAinitilisation(sim.optical);
end

