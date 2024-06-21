%script that takes the soution with Vext = 0 and computes JV curve
x = xnew;
Jvec = zeros(1,1);
Vext = zeros(1,1);
J = Jn+Jp;
Jvec(1) = Jscale*(mean(J));
Vext(1) = 0;

if mean(J) > 0 %use positive external voltage
    backwards = 1;
else
    backwards = -1; %use negative external voltage
end
    
for i = 1:200
    Vext(i+1) = 0 + backwards*0.01*i; %V (0-2V or -2-0V)
    phiNnew = phiN + Vext(i+1)/Vth; %non-dimensionalized
    psinN = log(nN/1) - phinN-phiNnew-nt; %update psinN
    psipN = log(pN/1) + phipN+phiNnew-pt; %update psipN
    j = 1;
while Fl2 >= 10^-5 || j<=15
    [~,~,~,Jnn,Jpn,~,~,~,Fl2,xn] = Newton_solve(xnew,G,mun,mup,...
            phin,phip,nbar,taun,taup,lambda,n1,p1,Nf,psin0,psinN,...
            psip0,psipN,phi0,phiNnew,nt,pt,N,h,damp,Linear,alpha,Cn,Cp);
    xnew = xn;
    l2error(j) = Fl2;
    if j >=15
        break;
    else
    end
    j = j+1;
end
J = Jnn+Jpn;
Jvec(i+1) = Jscale*(mean(J));

if Jscale*(mean(Jnn+Jpn)) <= 0 && backwards == 1 %positive J
    break;
elseif Jscale*(mean(Jnn+Jpn)) >= 0 && backwards == -1 %negative J
    break;
else
end
    
    J = Jnn+Jpn;
    verbose = 1;
    if verbose == 1
    clc;
    disp('')
    disp('Using FD method')
    disp('************************************')
    disp(['Damping parameter: ',num2str(damp)]);
    disp(['Relative error: ',num2str(l2error(end))])
    disp(['Jn+Jp: ',num2str(Jscale*(mean(J)))])
    disp('************************************')
    disp('')
    else
    end
end
P = Jvec.*Vext;
[Pmax, indm] = max(P,[],'all','linear'); %maximum power density mW/cm^2
Pin = 1000; %W/m^2
eta = ((Pmax*10/Pin)*100)/CSun; %efficiency 
Vmax = 1000*Vext(indm);
mPpoint = [Vmax,Pmax]; %maximum power point
disp(['efficiency: ',num2str(eta),' %']); %display the efficiency 