
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  PID-based search algorithm (PSA) source codes version 1.0
%
%  Developed in:	MATLAB 9.13 (R2022b)
%
%  Programmer:		Yuansheng Gao
%
%  Original paper:	Yuansheng Gao,
%                   PID-based search algorithm: A novel metaheuristic
%                   algorithm based on PID algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TargetX,TargetF,ConvergenceCurve] = PSA (N,T,lb,ub,nvars,fun)

    %% Initialization
    ConvergenceCurve = zeros(1,T);
    t = 1;
    Kp = 1;
    Ki = 0.5;
    Kd = 1.2;
    
    LogT = log(T);
    % Create the population
    [x, lb, ub]=Initialization(N,nvars,ub,lb);
    lbExtended = repmat (lb,[N,1]);
    ubExtended = repmat (ub,[N,1]);
    
    % Create the fitness vector
    f = zeros(N,1);
    for i = 1:N
        f(i,:) = fun (x(i,:));
    end
    
    %% main loop
    while t<=T
        % Step 1: Guide selection mechanism
        if t == 1
            [TargetF, index]=min(f);
            TargetX = x(index,:);
            newTargetX = TargetX;
            newTargetF = TargetF;
    
        else
            for i = 1:N
                f(i,:) = fun (x(i,:));
            end
            if min(f)<TargetF
                [newTargetF, index]=min(f);
                newTargetX = x(index,:);
            end
        end
        ConvergenceCurve(t) = newTargetF;
    
        % Step 2: Search operations
        if t == 1
            Ek = TargetX - x;
            Ek_1 = Ek;
            Ek_2 = Ek;
        else
            Ek_2 = Ek_1;
            Ek_1 = Ek + newTargetX - TargetX;
            Ek = newTargetX - x;
            TargetF = newTargetX;
            TargetX = newTargetF;
        end
    
        a = (log(T-t+2)/LogT)^2;
        out0 = (cos(1-t/T) + a*rand(N,nvars).*PSA_levyFlight(N,nvars)).*Ek;
        pid = rand(N,1)*Kp.*(Ek-Ek_1) + rand(N,1)*Ki.*Ek + rand(N,1)*Kd.*(Ek - 2*Ek_1 + Ek_2);
    
        % Step 3: Update mechanism
        r = rand(N,1)*cos(t/T);
        x = x + r.*pid + (1-r).*out0;
    
        lbViolated = x < lbExtended;
        ubViolated = x > ubExtended;
    
        x (lbViolated) = lbExtended (lbViolated);
        x (ubViolated) = ubExtended (ubViolated);
    
        % Next generation until termination criterion
        t = t + 1;
    end
end


function o=PSA_levyFlight(n,d)
beta=1.5;
sigma=(gamma(1+beta).*sin(pi*beta/2)./(gamma((1+beta)/2).*beta.*2.^((beta-1)/2))).^(1/beta);
u=randn(n,d)*sigma;
v=randn(n,d);
step=u./abs(v).^(1/beta);
o=step;
end


% This function is used to initialize the population
function [x, new_lb, new_ub]=Initialization(PopSize,nvars,ub,lb)

    num= size(ub,2); % Number of boundaries
    new_lb = lb;
    new_ub = ub;
    
    % If the boundaries of all variables are equal and user enter a signle
    % number for both ub and lb
    if num==1
        x=rand(PopSize,nvars).*(ub-lb)+lb;
        new_lb = lb*ones(1,nvars);
        new_ub = ub*ones(1,nvars);
    end
    
    % If each variable has a different lb and ub
    if num>1
        for i=1:nvars
            ubi=ub(i);
            lbi=lb(i);
            x(:,i)=rand(PopSize,1).*(ubi-lbi)+lbi;
        end
    end
end