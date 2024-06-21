%this loads Tom's solution and produces plots comparing to the FD 
%solution 
close all;
testname = "CIGS";
N = size(G,1);
os = 1; %'windows=2'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         plot phi                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(1)
    plot(xmid,phi,'LineWidth',2);
    title('$\phi$','fontsize',24,'interpreter','latex')
    xlabel('$z/L_{z}$','interpreter','latex');
    hold on
    %scatter(-1,phi0,'filled');
    %scatter(1,phiN,'filled');
    
    if os == 2
        Tomphi = csvread(strcat("Tom's solutions\",testname,"\phidata.csv"));
    elseif os == 1
        Tomphi = csvread(strcat("Tom's solutions/",testname,"/phidata.csv"));
    end
    Tomphi1 = Tomphi(:,1);
    Tomphi2 = Tomphi(:,2);
    Tlength = size(Tomphi,1);
    Lx = Tomphi1(end);
    plot(2*(Tomphi1-Lx/2)/Lx,Tomphi2,'--','LineWidth',2) %plot on the interval (-1,1)

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         plot n and p                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(2)
    semilogy(xmid,npscale*n,'LineWidth',2); hold on;
    semilogy(xmid,npscale*p,'LineWidth',2)
    title('$n,p$','fontsize',24,'interpreter','latex')
    legend('$n$','$p$','fontsize',24,'interpreter','latex',...
        'AutoUpdate','off')
    xlabel('$z/L_{z}$','interpreter','latex');
    ylabel('Charge-carrier Density (cm$^{-3}$)','interpreter','latex');
    hold on 
    %scatter(-1,npscale*n0,'filled');
    %scatter(-1,npscale*p0,'filled');
    %scatter(1,npscale*nN,'filled');
    %scatter(1,npscale*pN,'filled');
    
    if os == 2
        Tomn = csvread(strcat("Tom's solutions\",testname,"\ndata.csv"));
    elseif os == 1
        Tomn = csvread(strcat("Tom's solutions/",testname,"/ndata.csv"));
    end
    Tomn1 = Tomn(:,1);
    Tomn2 = npscale*Tomn(:,2);
    Tlength = size(Tomn,1);
    Lx = Tomn1(end);
    semilogy(2*(Tomn1-Lx/2)/Lx,Tomn2,'--','LineWidth',2) %plot on the interval (-1,1)
    hold on
    
    if os == 2
        Tomp = csvread(strcat("Tom's solutions\",testname,"\pdata.csv"));
    elseif os == 1
            Tomp = csvread(strcat("Tom's solutions/",testname,"/pdata.csv"));
    end
    
    Tomp1 = Tomp(:,1);
    Tomp2 = npscale*Tomp(:,2);
    Tlength = size(Tomp,1);
    Lx = Tomp1(end);
    semilogy(2*(Tomp1-Lx/2)/Lx,Tomp2,'--','LineWidth',2) %plot on the interval (-1,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         plot R and G                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(3)
    semilogy(xmid,Gscale*R,'LineWidth',2); hold on; 
    semilogy(xmid,Gscale*G,'LineWidth',2);
    legend('$R$','$G$','fontsize',24,'interpreter',...
        'latex','AutoUpdate','off')
    title('$G$,$R$','fontsize',24,'interpreter','latex')
    ylabel('Rate (cm$^{-3}$s$^{-1}$)','interpreter','latex');
    xlabel('$z/L_{z}$','interpreter','latex');


    %load R
    if os == 2
        TomR = csvread(strcat("Tom's solutions\",testname,"\RSRHdata.csv"));
    elseif os == 1
        TomR = csvread(strcat("Tom's solutions/",testname,"/RSRHdata.csv"));
    end
    TomR1 = TomR(:,1);
    TomR2 = Gscale*TomR(:,2);
    Tlength = size(TomR,1);
    Lx = TomR1(end);
    figure(3)
    semilogy(2*(TomR1-Lx/2)/Lx,TomR2,'--','LineWidth',2) %plot on the interval (-1,1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         plot Jn, Jp and Jn+Jp                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(4)
    plot(xnode,Jscale*(Jn+Jp),'LineWidth',2)
    hold on
    plot(xnode,Jscale*Jn,'LineWidth',2); hold on; 
    plot(xnode,Jscale*Jp,'LineWidth',2)
    delta = mean(Jscale*abs(Jn+Jp));
    ylim([mean(Jscale*(Jn+Jp))-delta;mean(Jscale*(Jn+Jp))+delta])
    title('$J_{n}+J_{p}$','fontsize',24,'interpreter','latex')
    xlabel('$z/L_{z}$','interpreter','latex');
    ylabel('J (mA cm$^{-2}$)','interpreter','latex')
    legend('$J_{n}+J_{p}$','$J_n$','$J_p$','fontsize',24,'interpreter',...
        'latex','AutoUpdate','off')
    %load Jn and Jp
    if os == 2
        TomJn = csvread(strcat("Tom's solutions\",testname,"\Jndata.csv"));
    elseif os == 1
        TomJn = csvread(strcat("Tom's solutions/",testname,"/Jndata.csv"));
    end
    TomJn1 = TomJn(:,1);
    TomJn2 = Jscale*TomJn(:,2);
    
    if os == 2
        TomJp = csvread(strcat("Tom's solutions\",testname,"\Jpdata.csv"));
    elseif os == 1
        TomJp = csvread(strcat("Tom's solutions/",testname,"/Jpdata.csv"));
    end
    
    TomJp1 = TomJp(:,1);
    TomJp2 = Jscale*TomJp(:,2);
    figure(4)
    plot(2*(TomJn1-Lx/2)/Lx,TomJn2,'--','LineWidth',2)
    hold on
    plot(2*(TomJp1-Lx/2)/Lx,TomJp2,'--','LineWidth',2)
    hold on
    plot(2*(TomJp1-Lx/2)/Lx,TomJn2+TomJp2,'--','LineWidth',2);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         plot psin and psip                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(5)
    plot(xmid,psin,'LineWidth',2); hold on; plot(xmid,psip,'LineWidth',2);
    title('$\psi_{n},\psi_{p}$','fontsize',24,'interpreter','latex')
    legend('$\psi_{n}$','$\psi_{p}$','fontsize',24,'interpreter',...
        'latex','AutoUpdate','off')
    xlabel('$z/L_{z}$','Interpreter','latex');
    hold on
    %scatter(-1,psin0,'filled');
    %scatter(-1,psip0,'filled');
    %scatter(1,psinN,'filled');
    %scatter(1,psipN,'filled');
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         plot J and JV curve                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(7)
    plot(1000*Vext,Jvec.*Vext,'LineWidth',2);
    hold on
    plot(1000*Vext,Jvec,'LineWidth',2)
    ylim([0,max(Jvec)+2])
    hold on
    scatter(mPpoint(1),mPpoint(2),'*','LineWidth',2); %plot maximum power point
    yyaxis left
    ylabel('Jn+Jp (mA cm$^-2$)','Interpreter','latex')
    
    if os == 2
        TomV = csvread(strcat("Tom's solutions\",testname,"\Vdata.csv"));
        TomJ = csvread(strcat("Tom's solutions\",testname,"\Jdata.csv"));
    elseif os == 1
        TomV = csvread(strcat("Tom's solutions/",testname,"/Vdata.csv"));
        TomJ = csvread(strcat("Tom's solutions/",testname,"/Jdata.csv"));
    end
    
    hold on
    plot(1000*TomV,TomJ,'--','LineWidth',2)
    hold on
    plot(1000*TomV,TomV.*TomJ,'--','LineWidth',2)
    yyaxis right
    ylim([0,max(Jvec)+2])
    ylabel('P (mW cm$^-2$)','interpreter','latex')
    xlabel('Vext (mV)')
    disp(['efficiency: ',num2str(eta)]);

    