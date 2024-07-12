mplot = 1;
N=size(G,1);
h=ones(N,1)/N; %an array for h (for nanometer scaling)
x = ones(3*N,1); % x = [psin,psip,phi]
%x = rand(3*N,1);
nt = 1;
pt = 1;

Linear = sim.input.Linear;

damp = 0.001;
alp0D = (log(1/p0)+1-phip0*damp-phi0*damp-nt*damp)/log(damp);
alpND = (log(1/pN)+1-phipN*damp-phiN*damp-nt*damp)/log(damp);
aln0D = (log(1/n0)+1+phin0*damp+phi0*damp-pt*damp)/log(damp);
alnND = (log(1/nN)+1+phinN*damp+phiN*damp-pt*damp)/log(damp);

taun_save=taun;
taup_save=taup;

for i=1:1000
    close all;
    damp = i*0.001;
    taun=taun_save/damp^2;
    taup=taup_save/damp^2;
    nND = nN*(damp^alnND);
    p0D = p0*(damp^alp0D);
    pND = pN*(damp^alpND);
    n0D = n0*(damp^aln0D);
    phin0D = damp*phin0;
    phinND = damp*phinN;
    phip0D = damp*phip0;
    phipND = damp*phipN;
    
    phinD = damp*phin;
    phipD = damp*phip;
    
    phiND = damp*(phiN+Vext);
    phi0D = damp*phi0;
    

    NfD = damp^2*Nf;
    
    psin0 = log(n0D/1) - phin0D-phi0D-damp*nt;
    psinN = log(nND/1) - phinND-phiND-damp*nt;
    
    psip0 = log(p0D/1) + phip0D+phi0D-damp*pt;
    psipN = log(pND/1) + phipND+phiND-damp*pt;
    
    Fl2 = 1;
    i = 1;
    l2error = zeros(1,1);
    
    while Fl2 >= 10^-10
        
        [psin,psip,phi,Jn,Jp,n,p,R,Fl2,xnew] = Newton_solve(x,G,mun,mup,...
            phinD,phipD,nbar,taun,taup,lambda,n1,p1,NfD,psin0,psinN,...
            psip0,psipN,phi0D,phiND,nt,pt,N,h,damp,Linear,alpha,Cn,Cp);
        x = xnew;
        l2error(i) = Fl2;
        i = i + 1;
        
        if i >=15
            break;  
        else
        end
        %nt = geomean(n);
        %pt = geomean(p);

    end
    J = Jn+Jp;
    verbose = 1;
    if verbose == 1
    %clc;
    disp('')
    disp('Using FD method')
    disp('************************************')
    disp(['Damping parameter: ',num2str(damp)]);
    disp(['nt,pt: ',num2str(nt),' ',num2str(pt)]);
    disp(['Relative error: ',num2str(l2error(end))])
    disp(['Jn+Jp: ',num2str(Jscale*(mean(J)))])
    disp('************************************')
    disp('')
    else
    end
end



