



%%% Designed and Developed by Mohammad Dehghani %%%


function[Best_score,Best_pos,DTBO_curve]=DTBO(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness)

lowerbound=ones(1,dimension).*(lowerbound);                              % Lower limit for variables
upperbound=ones(1,dimension).*(upperbound);                              % Upper limit for variables

%%
for i=1:dimension
    X(:,i) = lowerbound(i)+rand(SearchAgents,1).*(upperbound(i) - lowerbound(i));                          % Initial population
end

for i =1:SearchAgents
    L=X(i,:);
    fit(i)=fitness(L);
end
%%

for t=1:Max_iterations
    %% update the best member
    [best , blocation]=min(fit);
    if t==1
        Xbest=X(blocation,:);                                           % Optimal location
        fbest=best;                                           % The optimization objective function
    elseif best<fbest
        fbest=best;
        Xbest=X(blocation,:);
    end
    
    %%
    XF=[X fit'];
    XFsort=sortrows(XF,dimension+1);
    X=XFsort(:,1:dimension);
    fit=(  XFsort(:,1+dimension) )';
    N_DI=1+round(0.2*SearchAgents*(1-t/Max_iterations));
    
    DI=XFsort(1:N_DI,1:dimension);
    F_DI=(  XFsort(1:N_DI,1+dimension)  )';
    
    %% update DTBO population
    
    for i=1:SearchAgents
        
        %% Phase 1: training by the driving instructor (exploration)
        k_i=randperm(N_DI,1);
        DI_ki=DI(k_i);
        F_DI_ki=F_DI(k_i);
        I=round(1+rand(1,1));
        if F_DI_ki< fit (i)
            X_P1=X(i,:)+rand(1,1) .* (DI_ki-I.*X(i,:)); % Eq. (5)
        else
            X_P1=X(i,:)+rand(1,1) .* (1.*X(i,:)-I.*DI_ki); % Eq. (5)
        end
        X_P1= max(X_P1,lowerbound);X_P1 = min(X_P1,upperbound);
        
        % Update X_i based on Eq(6)
        F_P1 = fitness(X_P1);
        if F_P1 <= fit (i)
            X(i,:) = X_P1;
            fit (i)=F_P1;
        end
        
        %% END Phase 1: training by the driving instructor (exploration)
        
        %% Phase 2: learner driver patterning from instructor skills (exploration)
        PIndex=0.01+0.9*(1-t/Max_iterations);
        X_P2=(PIndex).* X(i,:)+(1-PIndex) .* (Xbest); % Eq. (7)
        X_P2= max(X_P2,lowerbound);X_P2 = min(X_P2,upperbound);
        
        % Update X_i based on Eq(8)
        F_P2 = fitness(X_P2);
        if F_P2 <= fit (i)
            X(i,:) = X_P2;
            fit (i)=F_P2;
        end
        %% END Phase 2: learner driver patterning from instructor skills (exploration)
        
        %% Phase 3: personal practice (exploitation)
        R=0.05;
        X_P3= X(i,:)+ (1-2*rand(1,dimension))*R*(1-t/Max_iterations).*X(i,:);% Eq.(9)
        X_P3= max(X_P3,lowerbound);X_P3 = min(X_P3,upperbound);
        
        % Update X_i based on Eq(10)
        F_P3 = fitness(X_P3);
        if F_P3 <= fit (i)
            X(i,:) = X_P3;
            fit (i)=F_P3;
        end
        %% END Phase 3: personal practice (exploitation)
        
    end% END for i=1:SearchAgents
    
    best_so_far(t)=fbest;
    
end% END for t=1:Max_iterations
Best_score=fbest;
Best_pos=Xbest;
DTBO_curve=best_so_far;
end

