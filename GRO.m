%  Gold Rush Optimizer (GRO)  
%
%  Source code version 1.0                                                                      
%                                                                                                     
%  Developed in MATLAB R2013b(v 8.2)                                                                   
%                                                                                                     
%  Author and programmer: Kamran Zolfi                                                          
%                                                                                                     
%  e-Mail: zolfi@iauln.ac.ir   
%          kamran.zolfi@gmail.com    
%  main article:
%           Gold Rush optimizer. A new population-based metaheuristic algorithm
%           DOI: 10.37190/ord230108
%           journal: Operations Research and Decisions
%           Vol. 33, No. 1 (2023)                                                                                       
%%
% You can define the cost function in a seperate file and load its handle to fobj 
% The initial parameters of the problem are:
% fobj = @CostFunction
% dim = dimension size
% Max_iter = maximum number of iterations
% N = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single numbers

% To run GRO: [best_score,best_pos,Convergence_curve]=GRO(dim,N,Max_iter,lb,ub,fobj )
%______________________________________________________________________________________________

function [best_score,best_pos,Convergence_curve]=GRO(N,Max_iter,lb,ub,dim,fobj )

lb = lb .* ones(1, dim);
ub = ub .* ones(1, dim);

%% GRO parameter initialization
sigma_initial = 2; 
sigma_final = 1 / Max_iter ;

% Initialize best position X* (global best)
best_pos=zeros(1, dim);
best_score=inf; %change this to -inf for maximization problems

%Initialize the gold prospectors?population Xi, i = 1, 2, . . . , N
Positions=initialization(N, dim, lb, ub);
Fit = inf(1,N);

%Initialize the gold prospectors?new positions Xnewi = Xi , i = 1, 2, . . . , N
X_NEW = Positions;
Fit_NEW = Fit;

Convergence_curve=zeros(1, Max_iter);
Convergence_curve(1) = min(Fit);
iter = 1;% Loop counter

%% Main loop
while iter <= Max_iter

    for i= 1:N   
        %Calculate fitness of current search agent at new position XNewi
        Fit_NEW(i) =  fobj(X_NEW(i,:));
        
        %Update position of current search agent Xi according to Equation (13)
        if Fit_NEW(i) < Fit(i)
            Fit(i) = Fit_NEW(i);
            Positions(i,:) = X_NEW(i,:);
        end
      
        %Update best search agent X*
        if Fit(i) < best_score
            % new gold mine is found
            best_score = Fit(i); 
            best_pos   = Positions(i,:);
        end

    end
    
   %Update l1, l2 by Equation (7)   
   l2 =  ((Max_iter - iter)/(Max_iter-1) )^2 * (sigma_initial - sigma_final) + sigma_final;
   l1 =  ((Max_iter - iter)/(Max_iter-1) )^1 * (sigma_initial - sigma_final) + sigma_final;
  
     
    %calculate the next position of current search agent XNewi with one of
    %... the migration, mining or collaboration methods
    for i = 1:size(Positions,1)

        coworkers = randperm(N-1,2);
        diggers = 1:N;
        diggers(i) = [];
        coworkers = diggers(coworkers);

        digger1 = coworkers(1);    %random prospector g1
        digger2 = coworkers(2);    %random prospector g2
        
        m = rand;
        %collaboration
        if m <  1/3
            for d  = 1:dim
                r1 = rand;                                         % r1 is a random number in [0,1]
                D3 = Positions(digger2,d) -  Positions(digger1,d); % Equation (11)
                X_NEW(i,d) = Positions(i,d) +  r1 * D3;            % Equation (12) 
            end
        %mining method
        elseif m < 2/3
            for d = 1:dim
                r1 = rand;                                          % r1 is a random number in [0,1]
                A2 = 2*l2*r1 - l2 ;                                 % Equation (10)               
                D2 = Positions(i,d) - Positions(digger1,d) ;        % Equation (8)
                X_NEW(i,d) = Positions(digger1,d) + A2*D2;          % Equation (9)                              
            end
        %migartion method    
        else
            for d = 1:dim
                r1 = rand; % r1 is a random number in [0,1]
                r2 = rand; % r2 is a random number in [0,1]
                C1 = 2 * r2;                                        % Equation (6)
                A1 = 1 + l1 * (r1 - 1/2);                           % Equation (5)
                D1 = C1 * best_pos(d) - Positions(i,d) ;            % Equation (3)
                X_NEW(i,d) = Positions(i,d) + A1 * D1;              % Equation (4)                      
            end
        end
            
        %Domain control
        X_NEW(i,:) = boundConstraint(X_NEW(i,:),Positions(i,:), lb , ub);
       
    end
    
    Convergence_curve(iter) = best_score;  
    iter = iter+1;
   
end


end


function Positions=initialization(SearchAgents_no,dim,ub,lb)

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
end

function newPos = boundConstraint(newPos, oldPos , lb, ub)

[NP, ~] = size(newPos);  % the Population size and the problem's dimension
%% check the lower bound
xl = repmat(lb, NP, 1);
pos = newPos < xl;
newPos(pos) =    oldPos(pos)  ;

%% check the upper bound
xu = repmat(ub, NP, 1);
pos = newPos > xu;
newPos(pos) =    oldPos(pos)  ;

end

