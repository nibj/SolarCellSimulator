%%%plot
close all
h=ones(N,1)/N;
mplot = 1;
xmid = 2*(cumsum(h)-h/2)-1;
xnode = [-1;2*(cumsum(h))-1];
if mplot ==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         plot n and p                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nout=npscale*n;
pout=npscale*p;
csvwrite("nout.csv",nout);
csvwrite("pout.csv",pout);
    figure(1)
    semilogy(xmid,npscale*n,'LineWidth',2); hold on;
    semilogy(xmid,npscale*p,'LineWidth',2)
    title('$n,p$','fontsize',24,'interpreter','latex')
    legend('$n$','$p$','fontsize',24,'interpreter','latex',...
        'AutoUpdate','off')
    xlabel('$z/L_{z}$','fontsize',24,'interpreter','latex');
    ylabel('Charge-carrier Density (cm$^{-3}$)','fontsize',24,...
        'interpreter','latex');
    hold on 
    %scatter(-1,npscale*n0,'filled');
    %scatter(-1,npscale*p0,'filled');
    %scatter(1,npscale*nN,'filled');
    %scatter(1,npscale*pN,'filled');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         plot Jn, Jp and Jn+Jp                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(2)
    plot(xnode,Jscale*(Jn+Jp),'LineWidth',2)
    hold on
    plot(xnode,Jscale*Jn,'LineWidth',2); hold on; 
    plot(xnode,Jscale*Jp,'LineWidth',2)
    ylim([Jscale*min([Jn;Jp]);1.5*Jscale*max([Jn;Jp])])
    title('$J_{n}+J_{p}$','fontsize',24,'interpreter','latex')
    xlabel('$z/L_{z}$','fontsize',24,'interpreter','latex');
    ylabel('J (mA cm$^{-2}$)','fontsize',24,'interpreter','latex')
    legend('$J_{n}+J_{p}$','$J_n$','$J_p$','fontsize',24,'interpreter',...
        'latex','AutoUpdate','off')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         plot phi                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(4)
    plot(xmid,Vth*phi,'LineWidth',2);
    title('$\phi$','fontsize',24,'interpreter','latex')
    xlabel('$z/L_{z}$','fontsize',24,'interpreter','latex');
    hold on
    %scatter(-1,phi0,'filled');
    %scatter(1,phiN,'filled');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         plot psin and psip                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(5)
    plot(xmid,psin,'LineWidth',2); hold on; plot(xmid,psip,'LineWidth',2);
    title('$\psi_{n},\psi_{p}$','fontsize',24,'interpreter','latex')
    legend('$\psi_{n}$','$\psi_{p}$','fontsize',24,'interpreter',...
        'latex','AutoUpdate','off')
    xlabel('$z/L_{z}$','fontsize',24,'Interpreter','latex');
    hold on
    %scatter(-1,psin0,'filled');
    %scatter(-1,psip0,'filled');
    %scatter(1,psinN,'filled');
    %scatter(1,psipN,'filled');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         plot R and G                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gavs=Gscale*G;
Ravs=Gscale*R;
csvwrite("Xav.csv",xmid);
csvwrite("Gav.csv",Gavs);
csvwrite("Rav.csv",Ravs);
    figure(6)
    
    semilogy(xmid,Gscale*R,'--','LineWidth',2); hold on; 
    semilogy(xmid,Gscale*G,'LineWidth',2);
    legend('$R$','$G$','fontsize',24,'interpreter',...
        'latex','AutoUpdate','off')
    title('$G$,$R$','fontsize',24,'interpreter','latex')
    ylabel('Rate (cm$^{-3}$s$^{-1}$)','fontsize',24,'interpreter','latex');
    xlabel('$z/L_{z}$','fontsize',24,'interpreter','latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         plot J and JV curve                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pav=Jvec.*Vext;
csvwrite("Vext.csv",(1000*Vext)');
csvwrite("JV.csv",Jvec');
csvwrite("PV.csv",Pav');

    figure(7)
    plot(1000*Vext,Jvec.*Vext,'LineWidth',2);
    hold on
    plot(1000*Vext,Jvec,'LineWidth',2)
    hold on
    %plot maximum power point
    scatter(mPpoint(1),mPpoint(2),'*','LineWidth',3); 
    yyaxis left
    ylabel('Jn+Jp (mA cm$^-2$)','fontsize',24,'Interpreter','latex')
    yyaxis right
    ylim([min(Jvec)-2,max(Jvec)+2])
    ylabel('P (mW cm$^-2$)','fontsize',24,'interpreter','latex',...
        'color','k')
    xlabel('Vext (mV)','interpreter','latex','fontsize',24)
    set(gca,'ycolor','k')
else
      
end
