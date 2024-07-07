% 灰狼优化算法
function [Alpha_score,Alpha_pos,Convergence_curve]=GWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
    % initialize alpha, beta, and delta_pos
    Alpha_pos=zeros(1,dim);
    Alpha_score=inf; %change this to -inf for maximization problems
    
    Beta_pos=zeros(1,dim);
    Beta_score=inf; %change this to -inf for maximization problems
    
    Delta_pos=zeros(1,dim);
    Delta_score=inf; %change this to -inf for maximization problems
    
    %Initialize the positions of search agents
    Positions=initialization(SearchAgents_no,dim,ub,lb);
    
    Convergence_curve=zeros(1,Max_iter);
    
    l=0;% Loop counter
    
    % Main loop
while l<Max_iter
    for i=1:size(Positions,1)  
       % 返回超出搜索空间边界的搜索值
        Flag4ub=Positions(i,:) > ub;
        Flag4lb=Positions(i,:) < lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
        % 计算每个搜索的目标函数值
        fitness=fobj(Positions(i,:));
        
        % 更新位置、Update alpha
        if fitness<Alpha_score 
            Alpha_score=fitness; 
            Alpha_pos=Positions(i,:);
        end
        
        if fitness > Alpha_score && fitness<Beta_score 
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score 
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
        end
    end
    
%     a∈[0,2]
    a=2-l*((2)/Max_iter); % a decreases linearly fron 2 to 0
    
    % 更新位置
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
              % r1 is a random number in [0,1]         
            r1=rand(); 
             % r2 is a random number in [0,1]
            r2=rand();
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
    end
    l=l+1;    
    Convergence_curve(l)=Alpha_score;
end




%%函数声明，输入参数包括搜索代理数量（SearchAgents_no）、问题维度（dim）、上界（ub）和下界（lb）。输出为搜索代理的位置矩阵
function Positions=initialization(SearchAgents_no,dim,ub,lb)
% 获取上界矩阵的列数，即边界数量，判断变量边界的情况
    Boundary_no= size(ub,2); % numnber of boundaries
    
    % If the boundaries of all variables are equal and user enter a signle
    % number for both ub and lb
    if Boundary_no==1
        Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
    end
    
    % If each variable has a different lb and ub
    if Boundary_no>1
        for i=1:dim
            ub_i=ub(i);
            lb_i=lb(i);
            Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
        end
    end


