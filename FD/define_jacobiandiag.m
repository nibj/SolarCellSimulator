function [DF,Fd] = define_jacobiandiag(df,u,w,h,N,w0,wN,H,dHl,dHr,du,pf,tag,...
    phip,phin,phi,pt,mup,damp,nt,mun,f)
%define_func - defines the diagonal blocks of the jacobian DF 
% 10/9/21
%
%   Inputs: df - (Nx1)  vector for df wrt w_{j+1/2}
%         : u - (Nx1)   vector for mu_n n, mu_p p or lambda^2
%         : w - (Nx1)   vector for the unknown variable psi_n, psi_p or phi
%         : h - (Nx1)   vector for the mesh-size 
%         : N - postive integer for number of intervals
%         : w0 - boundary condition for w(-1)
%         : wN - boundary conditions for w(1)
%         : dHl - (N-1x1) - derivative of H_{j}^{u} wrt w_{j-1/2}
%         : dHr - (N-1x1) - derivative of H_{j}^{u} wrt w_{j+1/2}
%         : du  - (Nx1) vector for the derivative of u wrt w_{j+1/2}
%         : pf  - flag for the blocks of DF associated to phi (0 or 1)
%  Outputs: DF - (NxN)   matrix that is a tridiagonal block of Jacobian
%         : Fd - (Nx1) vector containing diagonal entries 

hv  = ones(N,1).*h.^-1;

Fd  = zeros(N,1);
Fdm = zeros(N-1,1);
Fdp = zeros(N-1,1);

% j = 0 (left boundary term)
Fd(1) = hv(1)*(dHl(1)*w(2)-(dHl(1)+2*du(1)*hv(1))*w(1)+2*du(1)*w0*hv(1)-pf*H(1)-...
    pf*2*u(1)*hv(1))- df(1);

if tag == 'p' %modified boundary term
Fd(1) = hv(1)*(dHl(1)*w(2)-dHl(1)*w(1)-pf*H(1)-... 
    pf*(2)*u(1)*hv(1) - (1-pf)*(2*u(1)*hv(1)*(exp(w0-w(1))-1)))- 3*df(1)/4;   
end

if tag == 'n' %modified boundary term
Fd(1) = hv(1)*(dHl(1)*w(2)-dHl(1)*w(1)-pf*H(1)-... 
    pf*(2)*u(1)*hv(1) - (1-pf)*(2*u(1)*hv(1)*(exp(w0-w(1))-1)))- 3*df(1)/4;   
end

Fd(N) = hv(N)*(dHr(N-1)*w(N-1)-(dHr(N-1)+2*du(N)*hv(N))*w(N)+2*du(N)*wN*hv(N)-...
    pf*H(N-1)-pf*2*u(N)*hv(N)) - df(N);

if tag == 'n' %modified boundary term
Fd(N) = hv(N)*(dHr(N-1)*w(N-1)-dHr(N-1)*w(N)-pf*H(N-1)...
    -pf*2*u(N)*hv(N) - (1-pf)*(2*u(N)*hv(N)*(exp(wN-w(N))-1))) - 3*df(N)/4;
end

if tag == 'p' %modified boundary term
Fd(N) = hv(N)*(dHr(N-1)*w(N-1)-dHr(N-1)*w(N)-pf*H(N-1)...
    -pf*2*u(N)*hv(N) - (1-pf)*(2*u(N)*hv(N)*(exp(wN-w(N))-1))) - 3*df(N)/4;
end
% j = N (right boundary term)


% j = 1,2,..,N-1 (diagonal terms)
Fd(2:N-1) = hv(2:N-1).*(dHl(2:N-1).*w(3:N)-(dHl(2:N-1)+...
    dHr(1:N-2)).*w(2:N-1)+dHr(1:N-2).*w(1:N-2)-pf*H(2:N-1)-pf*H(1:N-2))...
    - df(2:N-1);

% j = 0 (left boundary term)
Fdp(1) = hv(1)*(dHr(1)*w(2)-dHr(1)*w(1)+pf*H(1));
% j = 2,3,..,N-1 
Fdp(2:N-1)=hv(2:N-1).*(dHr(2:N-1).*w(3:N)+pf*H(2:N-1)-...
    (dHr(2:N-1)).*w(2:N-1));
Fdp = [0;Fdp];


% j = N-1 (right boundary term)
Fdm(N-1) = hv(N)*(dHl(N-1)*w(N-1)-dHl(N-1)*w(N)+pf*H(N-1));
% j = 1,2,..,N-2 
Fdm(1:N-2)=hv(2:N-1).*(pf*H(1:N-2)-(dHl(1:N-2)).*w(2:N-1)+...
    dHl(1:N-2).*w(1:N-2));

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

if tag == 'p'
    %replace equations on either side of jump for p
    if isempty(jumpIndex)~=1 %if there is a jump
        for j = 1:size(jumpIndex,2)
            iL = jumpIndex(j)-1; %left index
            iR = jumpIndex(j);   %right index
            bottom = h(iL)*exp(phip(iL)+phi(iL))/(mup(iL)) + ...
                    h(iR)*exp(phip(iR)+phi(iR))/(mup(iR));
            top = 2*exp(pt*damp)*exp(w(iL)-w(iR))+...
                f(iL)*h(iL)^2*exp(phip(iL)+phi(iL))/(4*mup(iL))-...
                f(iR)*h(iR)^2*exp(phip(iR)+phi(iR))/(4*mup(iR));
            if pf == 1  %wrt psip
                dJu = (2*exp(pt*damp+w(iR))+...
                    h(iR)^2*df(iR)*exp(phip(iR)+phi(iR))/(4*mup(iR)))/...
                    bottom;
                dJl = -(2*exp(pt*damp+w(iL))+...
                    h(iL)^2*df(iL)*exp(phip(iL)+phi(iL))/(4*mup(iL)))/...
                    bottom;
            elseif pf == 0 %wrt phi
                dJu = (bottom*h(iR)^2*exp(phip(iR)+phi(iR))/(4*mup(iR))*...
                    (df(iR)+f(iR))+top*exp(phip(iR)+phi(iR))*h(iR)/mup(iR))/...
                    bottom^2;
                dJl = (-bottom*h(iL)^2*exp(phip(iL)+phi(iL))/(4*mup(iL))*...
                    (df(iL)+f(iL))+top*exp(phip(iL)+phi(iL))*h(iL)/mup(iL))/...
                    bottom^2;
            end
            %F(iL) = hv(iL)*(J-H(iL-1)*w(iL)+H(iL-1)*w(iL-1))-f(iL);
            %F(iR) = hv(iR)*(H(iR)*(w(iR+1)-w(iR))-J)-f(iR);
            %diagonal terms
            Fd(iL) = hv(iL).*(dJl+dHr(iL-1).*(w(iL-1)-w(iL))-pf*H(iL-1))...
                - df(iL);
            Fd(iR) = hv(iR).*(dHl(iR).*(w(iR+1)-w(iR))-pf*H(iR)-...
                dJu)- df(iR);
            Fdp(iR) = hv(iR)*dJu; 
            Fdm(iL) = -hv(iL)*dJl;
        end
    end
else
end

if tag == 'n'
    %replace equations on either side of jump for n
    if isempty(jumpIndex)~=1 %if there is a jump
        for j = 1:size(jumpIndex,2)
            iL = jumpIndex(j)-1; %left index
            iR = jumpIndex(j);   %right index
            bottom = h(iL)*exp(-(phin(iL)+phi(iL)))/(mun(iL)) + ...
                    h(iR)*exp(-(phin(iR)+phi(iR)))/(mun(iR));
            top = 2*exp(nt*damp)*exp(w(iL)-w(iR))+...
                f(iL)*h(iL)^2*exp(-(phin(iL)+phi(iL)))/(4*mun(iL))-...
                f(iR)*h(iR)^2*exp(-(phin(iR)+phi(iR)))/(4*mun(iR));
            if pf == 1  %wrt psip
                dJu = -(2*exp(nt*damp+w(iR))+...
                    h(iR)^2*df(iR)*exp(-(phin(iR)+phi(iR)))/(4*mun(iR)))/...
                    bottom;
                dJl = (2*exp(nt*damp+w(iL))+...
                    h(iL)^2*df(iL)*exp(-(phin(iL)+phi(iL)))/(4*mun(iL)))/...
                    bottom;
            elseif pf == 0 %wrt phi
                dJu = -(bottom*h(iR)^2*exp(-(phin(iR)+phi(iR)))/(4*mun(iR))*...
                    (df(iR)-f(iR))-top*exp(-(phin(iR)+phi(iR)))*h(iR)/mun(iR))/...
                    bottom^2;
                dJl = -(-bottom*h(iL)^2*exp(-(phin(iL)+phi(iL)))/(4*mun(iL))*...
                    (df(iL)-f(iL))+top*exp(-(phin(iL)+phi(iL)))*h(iL)/mun(iL))/...
                    bottom^2;
            end
            dJu = -(2*exp(nt*damp+w(iR))+...
                    h(iR)^2*df(iR)*exp(-(phin(iR)+phi(iR)))/(4*mun(iR)))/...
                    bottom;
                dJl = (2*exp(nt*damp+w(iL))+...
                    h(iL)^2*df(iL)*exp(-(phin(iL)+phi(iL)))/(4*mun(iL)))/...
                    bottom;
            %F(iL) = hv(iL)*(J-H(iL-1)*w(iL)+H(iL-1)*w(iL-1))-f(iL);
            %F(iR) = hv(iR)*(H(iR)*(w(iR+1)-w(iR))-J)-f(iR);
            %diagonal terms
            Fd(iL) = hv(iL).*(dJl+dHr(iL-1).*(w(iL-1)-w(iL))-pf*H(iL-1))...
                - df(iL);
            Fd(iR) = hv(iR).*(dHl(iR).*(w(iR+1)-w(iR))-pf*H(iR)-...
                dJu)- df(iR);
            Fdp(iR) = hv(iR)*dJu; 
            Fdm(iL) = -hv(iL)*dJl;
        end
    end
else
end

Fdm = [Fdm;0];
B = [Fdm Fd Fdp];
d = [-1 0 1];
DF = spdiags(B,d,N,N);

end

