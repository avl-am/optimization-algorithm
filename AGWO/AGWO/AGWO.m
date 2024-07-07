%___________________________________________________________________%
%  Grey Wold Optimizer (GWO) source codes version 1.0               %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, A. Lewis             %
%               Grey Wolf Optimizer, Advances in Engineering        %
%               Software , in press,                                %
%               DOI: 10.1016/j.advengsoft.2013.12.007               %
%                                                                   %
%___________________________________________________________________%

% Grey Wolf Optimizer
function [Alpha_score,Alpha_pos,Convergence_curve,Trajectories,fitness_history, position_history]=AGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems
fitness_history=zeros(1,Max_iter);
Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems
Trajectories=zeros(1,Max_iter);
Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems
%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
Convergence_curve=zeros(1,Max_iter);
position_history=zeros(SearchAgents_no,Max_iter,dim);
Trajectories=zeros(SearchAgents_no,Max_iter);
l=1;

% Loop counter

%____________________________________________________________________________________________
%for i=1:size(Positions,1)
%            position_history(i,1,:)=Positions(i,:);
%            Trajectories(:,1)=Positions(:,1);
%            fitness_history(i,1)=fitness(1,i);
%end
%_____________________________________________________________________________________________

% Main loop
while l<Max_iter+1

    for i=1:size(Positions,1)  

       % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));

        % Update Alpha, Beta, and Delta
        if fitness<Alpha_score 
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness<Beta_score 
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        
         if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score 
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
         end
 
              position_history(i,l,:)=Positions(i,:);
              Trajectories(:,l)=Positions(:,1);
              fitness_history(i,l)=fobj(Positions(i,:));
         
    end
    
    % QF=l^((2*rand()-1)/(1-Max_iter)^2);
     %G2=2*rand()-1; % Eq. (16)
     %G1=2*(1-(l/Max_iter));  % Eq. (17)
     a=cos((pi/2)*((l/Max_iter).^4));
   % a=2-(cos(rand())*l/Max_iter);% a decreases linearly fron 2 to 0
    % Update the Position of search agents including omegas

    for i=1:size(Positions,1)
        

        for j=1:size(Positions,2)   
        %________________________________________________________________________________
            r1 = rand();             % r1ÊÇ[0,1]ÖÐµÄËæ»úÊý
            r2 = rand();             % r2ÊÇ[0,1]ÖÐµÄËæ»úÊý
            A1 = 2*a*r1-a;        % ¹«Ê½£¨4£©
            C1 = 2*r2;            % ¹«Ê½£¨5£©
            D_alpha = abs(C1*Alpha_pos(j)-Positions(i, j));  % ¹«Ê½£¨6£©-µÚÒ»²¿·Ö
 
            r1 = rand();
            r2 = rand();
            A2 = 2*a*r1-a;         % ¹«Ê½£¨4£©
            C2 = 2*r2;             % ¹«Ê½£¨5£©
            D_beta = abs(C2*Beta_pos(j)-Positions(i, j));   % ¹«Ê½£¨6£©-µÚ¶þ²¿·Ö
            
            r1=rand();
            r2=rand(); 
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
           % Positions_old = Positions(i, :);
   
        %________________________________________________________________________________         
            if rand() < 0.8
                Positions(i, j) = 0.5*(Alpha_pos(j)-A1*D_alpha)+0.3*(Beta_pos(j)-A2*D_beta)+0.2*(Delta_pos(j)-A3*D_delta);
            else
               %Positions(i, j)=QF*Alpha_pos(j)-(G2*Positions(i,j)*rand)-G1*Levy+rand*G2; % Eq. (14)
              % rnew = rand();
               %p = rand();
                if rand <0.5
                  Positions(i,j)=Alpha_pos(j)*(1-l/Max_iter)+(mean(Positions(i,j))-Alpha_pos(j))*rand(); % Eq. (3) and Eq. (4)
         
                else
                %-------------------------------------------------------------------------------------
                    beta = 2*rand();
                    sigma_u = ((gamma(1+beta)*sin(pi*beta/2))/(gamma((1+beta)/2)*beta*2^(0.5*(beta-1))))^(1/beta);
                    u = normrnd(0, sigma_u);
                    v = normrnd(0, 1);
                    alpha_levi = 0.01*u/abs(v)^(-beta)*(Positions(i, j)-Alpha_pos(j));
                    Positions(i, j) =Alpha_pos(j)*alpha_levi;
               
                end  
            end
        end
    end
        disp(['In iteration #', num2str(l), ' , target''s objective = ', num2str(Alpha_score)])
        Convergence_curve(l)=Alpha_score;
        l=l+1;
end
end


