function [psin,psip,phi,Jn,Jp,n,p,R,Fl2,x] = Newton_solve(x,G,mun,mup,...
    phin,phip,ni,taun,taup,lambda,n1,p1,Nf,psin0,psinN,psip0,psipN,phi0,...
    phiN,nt,pt,N,h,damp,Linear,alpha,Cn,Cp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         define n,p and R(n,p)                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[psin,psip,phi,n,p,R,dRpsin,dRpsip,dRdphi] = n_p_phi_calc(x,nt,...
    pt,ni,taun,taup,n1,p1,N,phin,phip,damp,alpha,Cn,Cp);

Roff = 0;
if Linear || Roff
    dRpsin = 0*dRpsin;
    dRpsip = 0*dRpsip;
    dRdphi = 0*dRdphi;
    R = R*0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         define DF11                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = -G+R; df = dRpsin; u = mun.*n; du = mun.*n;

Hv1 = h(1:N-1)./u(1:N-1); 
Hv2 = h(2:N)./u(2:N);
H = 2./(Hv1+Hv2);

w = psin;
w0 = psin0;
wN = psinN;
dHr = h(2:N).*(H.^2)./(2*u(2:N));
dHl = h(1:N-1).*(H.^2)./(2*u(1:N-1));

tag = 'n';
[F1,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
[DF11,~] = define_jacobiandiag(df,u,w,h,N,w0,wN,H,dHl,dHr,du,1,...
    tag,phip,phin,phi,pt,mup,damp,nt,mun,f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         define DF22                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = -G+R; df = dRpsip; u = mup.*p; du = mup.*p;

Hv1 = h(1:N-1)./u(1:N-1);
Hv2 = h(2:N)./u(2:N);
H = 2./(Hv1+Hv2);

w = psip;
w0 = psip0;
wN = psipN;
dHr = h(2:N).*(H.^2)./(2*u(2:N));
dHl = h(1:N-1).*(H.^2)./(2*u(1:N-1));

tag = 'p';
[F2,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
[DF22,~] = define_jacobiandiag(df,u,w,h,N,w0,wN,H,dHl,dHr,du,1,tag,...
    phip,phin,phi,pt,mup,damp,nt,mun,f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         define DF33                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = n-p-Nf; df = n+p; u = lambda.^2; du = zeros(N,1);

if Linear
    f = 0*f;
    df=0*df;
end

Hv1 = h(1:N-1)./u(1:N-1);
Hv2 = h(2:N)./u(2:N);
H = 2./(Hv1+Hv2);

w = phi;
w0 = phi0;
wN = phiN;

dHr = zeros(N-1,1);
dHl = zeros(N-1,1);

tag = 'phi';
[F3,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
[DF33,~] = define_jacobiandiag(df,u,w,h,N,w0,wN,H,dHl,dHr,du,1,tag,...
    phip,phin,phi,pt,mup,damp,nt,mun,f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         define DF13                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = -G+R; df = dRdphi; u = mun.*n; du = mun.*n;

Hv1 = h(1:N-1)./u(1:N-1);
Hv2 = h(2:N)./u(2:N);
H = 2./(Hv1+Hv2);

w = psin;
w0 = psin0;
wN = psinN;
dHr = h(2:N).*(H.^2)./(2*u(2:N));
dHl = h(1:N-1).*(H.^2)./(2*u(1:N-1));

tag = 'n';
[DF13,~] = define_jacobiandiag(df,u,w,h,N,w0,wN,H,dHl,dHr,du,0,tag,...
    phip,phin,phi,pt,mup,damp,nt,mun,f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         define DF23                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = -G+R; df = dRdphi; u = mup.*p; du = -mup.*p;
Hv1 = h(1:N-1)./u(1:N-1);
Hv2 = h(2:N)./u(2:N);
H = 2./(Hv1+Hv2);

w = psip;
w0 = psip0;
wN = psipN;
dHr = -h(2:N).*(H.^2)./(2*u(2:N));
dHl = -h(1:N-1).*(H.^2)./(2*u(1:N-1));
tag = 'p';
[DF23,~] = define_jacobiandiag(df,u,w,h,N,w0,wN,H,dHl,dHr,du,0,tag,...
    phip,phin,phi,pt,mup,damp,nt,mun,f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         define DF12,DF21,DF31,DF32                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = 0;
DF12 = -(3/4)*spdiags(dRpsip,d,N,N);
DF21 = -(3/4)*spdiags(dRpsin,d,N,N); 
DF31 = spdiags(-n,d,N,N);
DF32 = spdiags(p,d,N,N);

if Linear
    DF31 = 0*DF31;
    DF32 = 0*DF32;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         define DF and solve                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DF = [DF11 DF12 DF13;...
      DF21 DF22 DF23;...
      DF31 DF32 DF33];
F = [F1;F2;F3];

testJacobian = 0;

if testJacobian == 1 %test the Jacobian
    [test11,test22,test33,test12,test13,test23,test31,test32,test21]=...
        approx_jacobian(x,G,mun,mup,...
    phin,phip,ni,taun,taup,lambda,n1,p1,Nf,psin0,psinN,psip0,psipN,phi0,...
    phiN,nt,pt,N,h,damp,F1,F2,F3,DF11,DF12,DF13,DF21,DF22,DF23,DF31,...
    DF32,DF33);   
else
end
%d=0;
%B = [Fd1;Fd2;Fd3];
%Binv = 1./B;
%P = spdiags(B,d,3*N,3*N);
%Pinv = spdiags(Binv,d,3*N,3*N);
precond = 0;
if precond ==1
    s = zeros(3*N,1);
    y = zeros(3*N,1);
    p = symrcm(DF);
    mat(p,p) = DF(p,p)*Pinv(p,p);
    y = mat\F;
    s = -P\y;
    x = x + s;
else
    s = zeros(3*N,1);
    p= symrcm(DF);
    s(p) = -DF(p,p)\F(p);
    x = x + s;
end
Fl2 = norm(F,2)/norm(x,2); %relative l2 norm of residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         post-processing                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[psin,psip,phi,n,p,R,~,~,~] = n_p_phi_calc(x,nt,...
    pt,ni,taun,taup,n1,p1,N,phin,phip,damp,alpha,Cn,Cp);

if Linear || Roff
    R = 0*R;
end
Jn = zeros(N+1,1);
Jp = zeros(N+1,1);

u = mun.*n;
Hv1 = h(1:N-1)./u(1:N-1);
Hv2 = h(2:N)./u(2:N);
H = 2./(Hv1+Hv2);

Jn(2:N) = H.*(psin(2:end)-psin(1:end-1));
%Jn(1) = 2*u(1)*(psin(1)-psin0)/h(1);
%Jn(N+1) = 2*u(N)*(psinN-psin(N))/h(N);
Jn(1) = -2*exp(nt)*mun(1)*exp(phin(1)+phi(1))*(exp(psin0)-...
    exp(psin(1)))/h(1)-h(1)*(R(1)-G(1))/4;
Jn(N+1) = 2*exp(nt)*mun(N)*exp(phin(N)+phi(N))*(exp(psinN)...
    -exp(psin(N)))/h(N)-h(N)*(R(N)-G(N))/4;
u = mup.*p;
Hv1 = h(1:N-1)./u(1:N-1);
Hv2 = h(2:N)./u(2:N);
H = 2./(Hv1+Hv2);

Jp(2:N) = -H.*(psip(2:end)-psip(1:end-1));
%Jp(1) = -2*u(1)*(psip(1)-psip0)/h(1);
Jp(1) = 2*exp(pt)*mup(1)*exp(-(phip(1)+phi(1)))*(exp(psip0)-...
    exp(psip(1)))/h(1)-h(1)*(G(1)-R(1))/4;
%Jp(N+1) = -2*u(N)*(psipN-psip(N))/h(N);
Jp(N+1) = -2*exp(pt)*mup(N)*exp(-(phip(N)+phi(N)))*(exp(psipN)...
    -exp(psip(N)))/h(N)-h(N)*(G(N)-R(N))/4;

%modification at a bandgap jump
jumpIndex = []; %empty array
i = 1;
for j=1:N-1 %find 
    phipLeft = phip(j); phipRight = phip(j+1);
    if abs(phipRight-phipLeft) > 10
        jumpIndex(i) = j+1;  % node index where jump occurs (1,2,...,N+1)
        i = i + 1;           % eqn index will be j-1 and j
    else 
    end    
end

f = G-R; %for p
if isempty(jumpIndex)~=1 %if there is a jump
        for j = 1:size(jumpIndex,2)
            iL = jumpIndex(j)-1; %left index
            iR = jumpIndex(j);   %right index
            ind = jumpIndex(j);
            Jp(ind) = (2*exp(pt*damp)*(exp(psip(iL))-exp(psip(iR))) + ...
                h(iL)^2*exp(phip(iL)+phi(iL))/(4*mup(iL))*f(iL) -...
                h(iR)^2*exp(phip(iR)+phi(iR))/(4*mup(iR))*f(iR))/...
                (h(iL)*exp(phip(iL)+phi(iL))/(mup(iL)) + ...
                h(iR)*exp(phip(iR)+phi(iR))/(mup(iR)));
        end
end

f = -G+R; %for n
if isempty(jumpIndex)~=1 %if there is a jump
        for j = 1:size(jumpIndex,2)
            iL = jumpIndex(j)-1; %left index
            iR = jumpIndex(j);   %right index
            ind = jumpIndex(j);
            Jn(ind) = (2*exp(nt*damp)*(exp(psin(iL))-exp(psin(iR))) + ...
                h(iL)^2*exp(-(phin(iL)+phi(iL)))/(4*mun(iL))*f(iL) -...
                h(iR)^2*exp(-(phin(iR)+phi(iR)))/(4*mun(iR))*f(iR))/...
                (h(iL)*exp(-(phin(iL)+phi(iL)))/(mun(iL)) + ...
                h(iR)*exp(-(phin(iR)+phi(iR)))/(mun(iR)));
        end
end

end  

