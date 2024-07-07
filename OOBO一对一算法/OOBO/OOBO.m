



%%% Designed and Developed by Mohammad Dehghani %%%


function[Best_score,Best_pos,OOBO_curve]=OOBO(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness)

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
    %% update the global best (fbest)
    [best , location]=min(fit); 
    if t==1
        Xbest=X(location,:);                                           % Optimal location
        fbest=best;                                           % The optimization objective function
    elseif best<fbest
        fbest=best;
        Xbest=X(location,:);
    end
    
    %% update agents position
    K=1:SearchAgents; %Eq(9)
    for i=1:SearchAgents
        k=randperm(SearchAgents-i+1,1);
        k=K(k);
        
        if k==i
            k=k+1;
            if k>SearchAgents
                k=k-2;
            end
        end
        K(find(K==k))=[];
        X_k=X(k,:);
        I=round(1+rand); %Eq(12)
        
        if fit(i)>fit(k)
            X_new=X(i,:)+ rand(1,dimension).*(X_k-I.* X(i,:)); %Eq(11)
        else
            X_new=X(i,:)+ rand(1,dimension).*(X(i,:)-X_k); %Eq(11)
        end
        X_new= max(X_new,lowerbound);X_new = min(X_new,upperbound);
        
        % based on Eq(13)
        ff = fitness(X_new);
        if ff <= fit (i)
            X(i,:) = X_new;
            fit (i)=ff;
        end
    end
    
    best_so_far(t)=fbest;
    average(t) = mean (fit);

end
    Best_score=fbest;
    Best_pos=Xbest;
    OOBO_curve=best_so_far;
end

