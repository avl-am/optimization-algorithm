



%%% Designed and Developed by Eva Trojovská and Mohammad Dehghani %%%


function[Best_score,Best_pos,CBOA_curve]=CBOA(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness)

lowerbound=ones(1,dimension).*(lowerbound);                              % Lower limit for variables
upperbound=ones(1,dimension).*(upperbound);                              % Upper limit for variables

%% INITIALIZATION
for i=1:dimension
    X(:,i) = lowerbound(i)+rand(SearchAgents,1).*(upperbound(i) - lowerbound(i));                          % Initial population
end

for i =1:SearchAgents
    L=X(i,:);
    fit(i)=fitness(L);
end
%%

for t=1:Max_iterations
    %% update the best condidate solution
    [Fbest , Blocation]=min(fit);
    if t==1
        Xbest=X(Blocation,:);                                           % Optimal location
        fbest=Fbest;                                           % The optimization objective function
    elseif Fbest<fbest
        fbest=Fbest;
        Xbest=X(Blocation,:);
    end
    
    %% UPDATE sorted X matrix based on objective function value using Eqs. (4) and (5).
    XF=[X fit'];
    XFsort=sortrows(XF,dimension+1);
    X=XFsort(:,1:dimension);
    fit=(XFsort(:,dimension+1))';
    N_C=max (1,round(0.2*SearchAgents*(1-t/Max_iterations)));
    %%
    
    %% PHASE 1: Chef instructors update process.
    for i=1:N_C
        %% start S1
        I=(1+rand); % I is a random number from set {1,2}.
        BC=X(1,:);
        X_C_S1=X(i,:)+rand(1,dimension).*(BC-I.*X(i,:)); % Eq(6)
        X_C_S1 = max(X_C_S1,lowerbound);X_C_S1 = min(X_C_S1,upperbound);
        
        % update position based on Eq (7)
        L=X_C_S1;
        fit_C_S1=fitness(L);
        if fit_C_S1<fit(i)
            X(i,:) = X_C_S1;
            fit(i) = fit_C_S1;
        end
        %% end S1
        
        %% start S2
        lo_local=lowerbound./t; %Eq(8)
        hi_local=upperbound./t; %Eq(9)
        X_C_S2=X(i,:)+lo_local+rand(1,1).*(hi_local-lo_local);% Eq(10)
        X_C_S2 = max(X_C_S2,lowerbound);X_C_S2 = min(X_C_S2,upperbound);
        
        % update position based on Eq (11)
        L=X_C_S2;
        fit_C_S2=fitness(L);
        if fit_C_S2<fit(i)
            X(i,:) = X_C_S2;
            fit(i) = fit_C_S2;
        end
        %% start S2
    end
    %% end Chef instructors update process.
    %%
    
    %% start Phase 2: Cooking students update process
    for i=1+N_C:SearchAgents
        %% S1
        % Random selection of Chef Instructor
        k=randperm(N_C,1);
        CI=X(k,:);
        F_CHEF=fit(k);
        I=(1+rand);
        
        X_S_S1=X(i,:)+rand(1,1).*(CI-I.*X(i,:));% Eq(12)
        X_S_S1 = max(X_S_S1,lowerbound);X_S_S1 = min(X_S_S1,upperbound);
        
        % update position based on Eq (13)
        L=X_S_S1;
        fit_S_S1=fitness(L);
        if fit_S_S1<fit(i)
            X(i,:) = X_S_S1;
            fit(i) = fit_S_S1;
        end
        %% S2
        % based on Eq (14)
        k=randperm(dimension,1);
        X_S_S2=X(i,:);
        X_S_S2(1,k)=CI(1,k);
        X_S_S2 = max(X_S_S2,lowerbound);X_S_S2 = min(X_S_S2,upperbound);
        
        % update position based on Eq (15)
        L=X_S_S2;
        fit_S_S2=fitness(L);
        if fit_S_S2<fit(i)
            X(i,:) = X_S_S2;
            fit(i) = fit_S_S2;
        end
        %% S3
        % based on Eq(16)
        k=randperm(dimension,1);
        X_S_S3=X(i,:);
        X_S_S3(1,k)=X_S_S3(1,k)+lo_local(k)+rand*(hi_local(k)-lo_local(k));
        X_S_S3 = max(X_S_S3,lowerbound);X_S_S3 = min(X_S_S3,upperbound);
        
        % update position based on Eq (17)
        L=X_S_S3;
        fit_S_S3=fitness(L);
        if fit_S_S3<fit(i)
            X(i,:) = X_S_S3;
            fit(i) = fit_S_S3;
        end
        
        %% start Phase 2: Cooking students update process
        %%
    end
    %%
    
    best_so_far(t)=fbest;
    average(t) = mean (fit);
    
end
Best_score=fbest;
Best_pos=Xbest;
CBOA_curve=best_so_far;
end

