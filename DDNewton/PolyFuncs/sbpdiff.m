function [ out ] = sbpdiff(p, d, sim)
%sbpdiff(p, d, sim) Differentiates an sbp variable, returning a new sbp
% function

out = zeros(sim.setup.pdeg1, sim.setup.nx);

for i = 1:sim.setup.nx
    temp = pdiff(p(:,i),d)*2/sim.setup.mesh(i);  % the bug is gone
    out(:,i) = polymatch(temp, sim.setup.pdeg1);
end



end

