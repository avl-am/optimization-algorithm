



%%% Designed and Developed by Pavel Trojovský and Mohammad Dehghani %%%


function[Best_score,Best_pos,OOA_curve]=OOA(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness)
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

for t=1:Max_iterations  % algorithm iteration
    
    %%  update: BEST proposed solution
    [Fbest , blocation]=min(fit);
    
    if t==1
        xbest=X(blocation,:);                                           % Optimal location
        fbest=Fbest;                                           % The optimization objective function
    elseif Fbest<fbest
        fbest=Fbest;
        xbest=X(blocation,:);
    end
    %%
    %%
    for i=1:SearchAgents
        %% Phase 1: : POSITION IDENTIFICATION AND HUNTING THE FISH (EXPLORATION)
        fish_position=find(fit<fit(i));% Eq(4)
        if size(fish_position,2)==0
            selected_fish=xbest;
        else
            if rand <0.5
                selected_fish=xbest;
            else
                k=randperm(size(fish_position,2),1);
                selected_fish=X(fish_position(k));
            end
        end
        %
        I=round(1+rand);
        X_new_P1=X(i,:)+rand(1,1).*(selected_fish-I.*X(i,:));%Eq(5)
        X_new_P1 = max(X_new_P1,lowerbound);X_new_P1 = min(X_new_P1,upperbound);
        
        % update position based on Eq (6)
        L=X_new_P1;
        fit_new_P1=fitness(L);
        if fit_new_P1<fit(i)
            X(i,:) = X_new_P1;
            fit(i) = fit_new_P1;
        end
        %% END Phase 1
        
        %%
        %% PHASE 2: CARRYING THE FISH TO THE SUITABLE POSITION (EXPLOITATION)
        X_new_P1=X(i,:)+(lowerbound+rand*(upperbound-lowerbound))/t;%Eq(7)
        X_new_P1 = max(X_new_P1,lowerbound);X_new_P1 = min(X_new_P1,upperbound);
        % update position based on Eq (8)
        L=X_new_P1;
        fit_new_P1=fitness(L);
        if fit_new_P1<fit(i)
            X(i,:) = X_new_P1;
            fit(i) = fit_new_P1;
        end
        %% END Phase 2
        %%
    end
    %%
    
    best_so_far(t)=fbest;
    average(t) = mean (fit);
    
end
Best_score=fbest;
Best_pos=xbest;
OOA_curve=best_so_far;
end

