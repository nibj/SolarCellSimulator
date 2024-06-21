function [dataout,h] = Dataout(datain,sim,outflag)
%Dataout- Takes the data from Tom's code, and outputs data so it can be 
%read by new electical solver
%
% Inputs:   datain - data formatted in Tom's code (ex: sim.sbp.G)
%           sim 
%           outflag - 0 for midpoint values, 1 for endpoint values
% Outputs:  dataout - data formatted for new electical solver

scale = 1; %nm scaling 1,1/2,1/4 etc
if outflag == 0
    xsamp = linspace(-1+scale,1-scale,1/scale);
    dataout = zeros(1,(1/scale)*sum(sim.setup.nx_sec));
    for i = 1:size(datain,2)
        sigvals = polyval(datain(:,i), xsamp); %choosing 0 should compute in the middle
        ind1 = 1+(i-1)/scale;
        ind2 = i/scale;
        dataout(ind1:ind2) = sigvals;
    end
elseif outflag == 1
    dataout = zeros(1,(1/scale)*sum(sim.setup.nx_sec)+1);
    xsamp = linspace(-1,1,(1/scale)+1);
    xsamp = xsamp(1:end-1); %prevent double counting
    for i = 1:size(datain,2)
        sigvals = polyval(datain(:,i), xsamp);%choosing -1 should compute on left
        ind1 = 1+(i-1)/scale;
        ind2 = i/scale;
        dataout(ind1:ind2) = sigvals;
    end
        dataout(end) = polyval(datain(:,end), 1); %right endpoint for right most endpoint
elseif outflag == 2 %variable h
    datainEg = sim.sbp.Eg;
    xsamp = linspace(-1+scale,1-scale,1/scale);
    dataout = zeros(1,(1/scale)*sum(sim.setup.nx_sec));
    for i = 1:size(datainEg,2)
        sigvals = polyval(datainEg(:,i), xsamp); %choosing 0 should compute in the middle
        ind1 = 1+(i-1)/scale;
        ind2 = i/scale;
        dataout(ind1:ind2) = sigvals;
    end
    %N = size(dataout,2)-4*(1/scale); %works
    N = size(dataout,2);
    %h = ones(1:N); %sum has to be N

    %create a mesh

    alpha = 2.5;
    M = floor(alpha/log10(3/(2+10^-alpha)));
    %solve N = alpha/log10(3/(2+10^-alpha)) for alpha
    f = @(x) (3/(2+10^-x))-10^(x/M);
    fprime = @(x) (3*log(10)*10^-x)/((2+10^-x))^2-(log(10)/M)*10^(x/M);
    x = alpha;
    for i = 1:1000
        xnew = x -f(x)/fprime(x);
        x = xnew;
    end
    alpha = x;
    delta = 3/(2+10^-alpha);

    %check
    j = M:-1:1;
    vec = (scale)*(1/delta).^j;
    val = sum(vec); %should be equal to 2h
    fix = M/2;
    
    jmp_ind = zeros(1,1);
    j = 1;
    %check for a jump in the bandgap Eg
    for i = 1:size(dataout,2)-1
        val1 = dataout(i);
        val2 = dataout(i+1);
        if abs(val1-val2)>=1
            jmp_ind(j)=i;
            j = j+1;
        end
        
    end
    %construct h
    h = [];
    if jmp_ind(1)==0 %no jumps
        jndnew = [0 N];
    else
        jndnew = [0 jmp_ind N];
    end
    for i = 1:size(jndnew,2)-1
        onevec = scale*ones(jndnew(i+1)-2-(jndnew(i)+2),1)';
        h = [ h vec onevec flip(vec)];
    end
   
    sampvals = cumsum(h)-h/2; %sample in the middle
    for i = 1:size(h,2)
        x = sampvals(i);
        k = ceil(x);
        n = k-1;
        xsamp = 2*x - (1+2*n);
        if k> size(datain,2)
            k = size(datain,2);
        else
        end
        sigvals = polyval(datain(:,k),xsamp);
        dataout(i) = sigvals;
    end
    
    h = (1/scale)*h'/N;
    
    
%     h = scale*[vec ones(1,N-4) flip(vec)];
%     h = scale*[vec ones(1,93) flip(vec) vec ones(1,101) flip(vec)];
%     %h = ones(1,N);
%     sampvals = cumsum(h)-h/2; %sample in the middle
%     for i = 1:size(h,2)
%         x = sampvals(i);
%         k = ceil(x);
%         n = k-1;
%         xsamp = 2*x - (1+2*n);
%         if k> size(datain,2)
%             k = size(datain,2);
%         else
%         end
%         sigvals = polyval(datain(:,k),xsamp);
%         dataout(i) = sigvals;
%     end
%      h = [vec ones(1,N-4) flip(vec)]/N;
%      h = [vec ones(1,93) flip(vec) vec ones(1,101) flip(vec)]/N;
%      h = h';
end

end

