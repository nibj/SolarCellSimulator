function [test11,test22,test33,test12,test13,test23,test31,test32,test21] = approx_jacobian(x,G,mun,mup,...
    phin,phip,ni,taun,taup,lambda,n1,p1,Nf,psin0,psinN,psip0,psipN,phi0,...
    phiN,nt,pt,N,h,damp,F1,F2,F3,DF11,DF12,DF13,DF21,DF22,DF23,DF31,...
    DF32,DF33)
% DF11 approximation
vec = zeros(3*N,1);
Japprox1 = zeros(N,N);
for j = 1:N
    vec = zeros(3*N,1);
    vec(j) = 1;
    xh = x + vec*sqrt(eps);
    [psin,psip,phi,n,p,R,dRpsin,dRpsip,dRdphi] = n_p_phi_calc(xh,nt,...
    pt,ni,taun,taup,n1,p1,N,phin,phip,damp);
    f = -G+R; df = dRpsin; u = mun.*n; du = mun.*n;
    
    Hv1 = h(1:N-1)./u(1:N-1); %changed
    Hv2 = h(2:N)./u(2:N);
    H = 2./(Hv1+Hv2);
    
    w = psin;
    w0 = psin0;
    wN = psinN;
    dHr = h(2:N).*(H.^2)./(2*u(2:N));
    dHl = h(1:N-1).*(H.^2)./(2*u(1:N-1));
    
    tag = 'n';
    [Fh,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
    
    val = (Fh-F1)/sqrt(eps);
    Japprox1(:,j) = val;
end
test11 = DF11-Japprox1;

% DF22 approximation
vec = zeros(3*N,1);
Japprox2 = zeros(N,N);
for j = 1:N
    vec = zeros(3*N,1);
    vec(N+j) = 1;
    xh = x + vec*sqrt(eps);
    [psin,psip,phi,n,p,R,dRpsin,dRpsip,dRdphi] = n_p_phi_calc(xh,nt,...
    pt,ni,taun,taup,n1,p1,N,phin,phip,damp);
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
    [Fh,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
    
    val = (Fh-F2)/sqrt(eps);
    Japprox2(:,j) = val;
end
test22 = DF22-Japprox2;

% DF33 approximation
vec = zeros(3*N,1);
Japprox3 = zeros(N,N);
for j = 1:N
    vec = zeros(3*N,1);
    vec(2*N+j) = 1;
    xh = x + vec*sqrt(eps);
    [psin,psip,phi,n,p,R,dRpsin,dRpsip,dRdphi] = n_p_phi_calc(xh,nt,...
    pt,ni,taun,taup,n1,p1,N,phin,phip,damp);
    f = n-p-Nf; df = n+p; u = lambda.^2; du = zeros(N,1);

    Hv1 = h(1:N-1)./u(1:N-1);
    Hv2 = h(2:N)./u(2:N);
    H = 2./(Hv1+Hv2);
    
    w = phi;
    w0 = phi0;
    wN = phiN;
    
    dHr = zeros(N-1,1);
    dHl = zeros(N-1,1);
    
    tag = 'phi';
    [Fh,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
    
    val = (Fh-F3)/sqrt(eps);
    Japprox3(:,j) = val;
end
test33 = DF33-Japprox3;


% DF12 approximation
vec = zeros(3*N,1);
Japprox1 = zeros(N,N);
for j = 1:N
    vec = zeros(3*N,1);
    vec(N+j) = 1;
    xh = x + vec*sqrt(eps);
    [psin,psip,phi,n,p,R,dRpsin,dRpsip,dRdphi] = n_p_phi_calc(xh,nt,...
    pt,ni,taun,taup,n1,p1,N,phin,phip,damp);
    f = -G+R; df = dRpsin; u = mun.*n; du = mun.*n;
    
    Hv1 = h(1:N-1)./u(1:N-1); %changed
    Hv2 = h(2:N)./u(2:N);
    H = 2./(Hv1+Hv2);
    
    w = psin;
    w0 = psin0;
    wN = psinN;
    dHr = h(2:N).*(H.^2)./(2*u(2:N));
    dHl = h(1:N-1).*(H.^2)./(2*u(1:N-1));
    
    tag = 'n';
    [Fh,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
    
    val = (Fh-F1)/sqrt(eps);
    Japprox1(:,j) = val;
end
test12 = DF12-Japprox1;

% DF13 approximation
vec = zeros(3*N,1);
Japprox1 = zeros(N,N);
for j = 1:N
    vec = zeros(3*N,1);
    vec(2*N+j) = 1;
    xh = x + vec*sqrt(eps);
    [psin,psip,phi,n,p,R,dRpsin,dRpsip,dRdphi] = n_p_phi_calc(xh,nt,...
        pt,ni,taun,taup,n1,p1,N,phin,phip,damp);
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
    [Fh,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
    
    val = (Fh-F1)/sqrt(eps);
    Japprox1(:,j) = val;
end
test13 = DF13-Japprox1;


% DF21 approximation
vec = zeros(3*N,1);
Japprox2 = zeros(N,N);
for j = 1:N
    vec = zeros(3*N,1);
    vec(j) = 1;
    xh = x + vec*sqrt(eps);
    [psin,psip,phi,n,p,R,dRpsin,dRpsip,dRdphi] = n_p_phi_calc(xh,nt,...
    pt,ni,taun,taup,n1,p1,N,phin,phip,damp);
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
    [Fh,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
    
    val = (Fh-F2)/sqrt(eps);
    Japprox2(:,j) = val;
end
test21 = DF21-Japprox2;

% DF23 approximation
vec = zeros(3*N,1);
Japprox2 = zeros(N,N);
for j = 1:N
    vec = zeros(3*N,1);
    vec(2*N+j) = 1;
    xh = x + vec*sqrt(eps);
    [psin,psip,phi,n,p,R,dRpsin,dRpsip,dRdphi] = n_p_phi_calc(xh,nt,...
        pt,ni,taun,taup,n1,p1,N,phin,phip,damp);
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
    [Fh,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
    
    val = (Fh-F2)/sqrt(eps);
    Japprox2(:,j) = val;
end
test23 = DF23-Japprox2;

% DF31 approximation
vec = zeros(3*N,1);
Japprox3 = zeros(N,N);
for j = 1:N
    vec = zeros(3*N,1);
    vec(j) = 1;
    xh = x + vec*sqrt(eps);
    [psin,psip,phi,n,p,R,dRpsin,dRpsip,dRdphi] = n_p_phi_calc(xh,nt,...
    pt,ni,taun,taup,n1,p1,N,phin,phip,damp);
    f = n-p-Nf; df = n+p; u = lambda.^2; du = zeros(N,1);

    Hv1 = h(1:N-1)./u(1:N-1);
    Hv2 = h(2:N)./u(2:N);
    H = 2./(Hv1+Hv2);
    
    w = phi;
    w0 = phi0;
    wN = phiN;
    
    dHr = zeros(N-1,1);
    dHl = zeros(N-1,1);
    
    tag = 'phi';
    [Fh,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
    
    val = (Fh-F3)/sqrt(eps);
    Japprox3(:,j) = val;
end
test31 = DF31-Japprox3;


% DF32 approximation
vec = zeros(3*N,1);
Japprox3 = zeros(N,N);
for j = 1:N
    vec = zeros(3*N,1);
    vec(N+j) = 1;
    xh = x + vec*sqrt(eps);
    [psin,psip,phi,n,p,R,dRpsin,dRpsip,dRdphi] = n_p_phi_calc(xh,nt,...
    pt,ni,taun,taup,n1,p1,N,phin,phip,damp);
    f = n-p-Nf; df = n+p; u = lambda.^2; du = zeros(N,1);

    Hv1 = h(1:N-1)./u(1:N-1);
    Hv2 = h(2:N)./u(2:N);
    H = 2./(Hv1+Hv2);
    
    w = phi;
    w0 = phi0;
    wN = phiN;
    
    dHr = zeros(N-1,1);
    dHl = zeros(N-1,1);
    
    tag = 'phi';
    [Fh,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
    
    val = (Fh-F3)/sqrt(eps);
    Japprox3(:,j) = val;
end
test32 = DF32-Japprox3;

% DF21 approximation
vec = zeros(3*N,1);
Japprox2 = zeros(N,N);
for j = 1:N
    vec = zeros(3*N,1);
    vec(j) = 1;
    xh = x + vec*sqrt(eps);
    [psin,psip,phi,n,p,R,dRpsin,dRpsip,dRdphi] = n_p_phi_calc(xh,nt,...
    pt,ni,taun,taup,n1,p1,N,phin,phip,damp);
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
    [Fh,~] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun);
    
    val = (Fh-F2)/sqrt(eps);
    Japprox2(:,j) = val;
end
test21 = DF21-Japprox2;

end