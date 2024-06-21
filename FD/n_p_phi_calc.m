function [psin,psip,phi,n,p,R,dRpsin,dRpsip,dRdphi] = n_p_phi_calc(x,nt,...
    pt,ni,taun,taup,n1,p1,N,phin,phip,damp,alpha,Cn,Cp)

%n_p_phi_calc - calculates n,p,phi and the recombination rate

psin = x(1:N);
psip = x(N+1:2*N);
phi = x(2*N+1:3*N);

n = exp(psin+phin+phi+nt*damp);
p = exp(psip-phip-phi+pt*damp);

%define R(n,p) and dR/dpsin
Rt = n.*p - ni.^2;
Rb = taun.*(p+p1)+taup.*(n+n1);
R = Rt./Rb + alpha.*(n.*p-ni.^2)+...
    Cn.*n.*(n.*p-ni.^2)+Cp.*p.*(n.*p-ni.^2); %SRH + radiative + Auger

dRdpsint = n.*p.*Rb-Rt.*taup.*n;
dRdpsinb = Rb.^2;
dRpsin = dRdpsint./dRdpsinb + alpha.*n.*p + ...
    Cn.*(n.*(n.*p-ni.^2)+n.*n.*p)+Cp.*n.*p.*p; %SRH + radiative + Auger

dRdpsipt = n.*p.*Rb-Rt.*taun.*p;
dRdpsipb = Rb.^2;
dRpsip = dRdpsipt./dRdpsipb + alpha.*n.*p + ...
    Cn.*(n.*n.*p)+Cp.*(p.*(n.*p-ni.^2)+p.*p.*n); %SRH + radiative + Auger 

dRdphi = (Rt.*(taun.*p-taup.*n))./(Rb.^2) + ...
    Cn.*n.*(n.*p-ni.^2)-Cp.*p.*(n.*p-ni.^2); %SRH + radiative + Auger 
% Note: dRad/dphi=0


end

