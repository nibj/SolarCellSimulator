function [F,Fl2] = define_func(f,u,w,h,N,w0,wN,H,tag,phip,phin,phi,pt,mup,damp,nt,mun)
%define_func - defines the FD system to solve
% 1/6/21
%
%   Inputs: f - (Nx1)   vector for the right-hand side values at x_{j+1/2}
%         : u - (Nx1)   vector for mu_n n, mu_p p or lambda^2
%         : w - (Nx1)   vector for the unknown variable psi_n, psi_p or phi
%         : h - (Nx1)   vector for the mesh-size 
%         : N - postive integer for number of intervals
%         : w0 - boundary condition for w(-1)
%         : wN - boundary conditions for w(1)
%         : H - (N-1x1) vector for the haronic average of u at the interior
%         : tag - tag used to specify modified method
%  Outputs: F - (Nx1)   vector for the N equations

hv = ones(N,1).*h.^-1;
F = zeros(N,1);

% j = 0 (left boundary term)
F(1) = hv(1)*(H(1)*w(2)-(H(1)+2*u(1)*hv(1))*w(1)+2*u(1)*w0*hv(1)) - f(1);
% j = N-1 (right boundary term)
F(N) = hv(N)*(H(N-1)*w(N-1)-(H(N-1)+2*u(N)*hv(N))*w(N)+2*u(N)*wN*hv(N)) - f(N);

if tag == 'p' %modified boundary term
    F(1) = hv(1)*(H(1)*w(2)-H(1)*w(1)+2*exp(pt*damp)*mup(1)*hv(1)*exp(-(phi(1)+...
        phip(1)))*(exp(w0)-exp(w(1)))) - 3*f(1)/4; 
end

if tag == 'n' %modified boundary term
    F(1) = hv(1)*(H(1)*w(2)-H(1)*w(1)+2*exp(nt*damp)*mun(1)*hv(1)*exp(phi(1)+...
        phin(1))*(exp(w0)-exp(w(1)))) - 3*f(1)/4; 
end

if tag == 'n' %modified boundary term
    F(N) = hv(N)*(H(N-1)*w(N-1)-H(N-1)*w(N)+2*exp(nt*damp)*mun(N)*hv(N)*exp(phi(N)+...
        phin(N))*(exp(wN)-exp(w(N)))) - 3*f(N)/4;
end

if tag == 'p' %modified boundary term
    F(N) = hv(N)*(H(N-1)*w(N-1)-H(N-1)*w(N)+2*exp(pt*damp)*mup(N)*hv(N)*exp(-(phi(N)+...
        phip(N)))*(exp(wN)-exp(w(N)))) - 3*f(N)/4;
end

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

% j = 1,2,..,N-2 
F(2:N-1) = hv(2:N-1).*(H(2:N-1).*w(3:N)-(H(2:N-1)+H(1:N-2)).*w(2:N-1)...
    +H(1:N-2).*w(1:N-2)) - f(2:N-1);

if tag == 'p'
    %replace equations on either side of jump for p
    if isempty(jumpIndex)~=1 %if there is a jump
        for j = 1:size(jumpIndex,2)
            iL = jumpIndex(j)-1; %left index
            iR = jumpIndex(j);   %right index
            J = -(2*exp(pt*damp)*(exp(w(iL))-exp(w(iR))) + ...
                h(iL)^2*exp(phip(iL)+phi(iL))/(4*mup(iL))*f(iL) -...
                h(iR)^2*exp(phip(iR)+phi(iR))/(4*mup(iR))*f(iR))/...
                (h(iL)*exp(phip(iL)+phi(iL))/(mup(iL)) + ...
                h(iR)*exp(phip(iR)+phi(iR))/(mup(iR)));
            F(iL) = hv(iL)*(J-H(iL-1)*w(iL)+H(iL-1)*w(iL-1))-f(iL);
            F(iR) = hv(iR)*(H(iR)*(w(iR+1)-w(iR))-J)-f(iR);
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
            J = (2*exp(nt*damp)*(exp(w(iL))-exp(w(iR))) + ...
                h(iL)^2*exp(-(phin(iL)+phi(iL)))/(4*mun(iL))*f(iL) -...
                h(iR)^2*exp(-(phin(iR)+phi(iR)))/(4*mun(iR))*f(iR))/...
                (h(iL)*exp(-(phin(iL)+phi(iL)))/(mun(iL)) + ...
                h(iR)*exp(-(phin(iR)+phi(iR)))/(mun(iR)));
            F(iL) = hv(iL)*(J-H(iL-1)*w(iL)+H(iL-1)*w(iL-1))-f(iL);
            F(iR) = hv(iR)*(H(iR)*(w(iR+1)-w(iR))-J)-f(iR);
        end
    end
else
end

Fl2 = norm(F,2); %norm of the residual F in the l2 norm

%F = sparse(F); remove
end

