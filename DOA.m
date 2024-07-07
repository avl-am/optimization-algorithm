%_________________________________________________________________________%
%  Dingo Optimization Algorithm (DOA) source code                         %
%                                                                         %
%  Developed in MATLAB 9.4.0.813654 (R2018a)                              %
%                                                                         %
%  Author: Dr. Hernan Peraza-Vazquez                                      %
%          MTA. Gustavo Echavarria-Castillo                               %
%                                                                         %
%  e-mail:  hperaza@ipn.mx        gechavarriac1700@alumno.ipn.mx          %
%                                                                         %
%  Programmer:  Dr. Hernan Peraza-Vazquez                                 %
%  Main paper:                                                            %
%  A Bio-Inspired Method for Engineering Design Optimization Inspired by  %
%  Dingoes Hunting Strategies.                                            %
%  Mathematical Problems in Engineering. (2021). Hindawi.                 %                                                      %
%  DOI:   doi.org/10.1155/2021/9107547                                    %
%_________________________________________________________________________%
function [vMin,theBestVct,Convergence_curve]=DOA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
P= 0.5;  % Hunting or Scavenger?  rate.  See section 3.0.4, P and Q parameters analysis
Q= 0.7;  % Group attack or persecution? 
beta1= -2 + 4* rand();  % -2 < beta < 2     Used in Eq. 2, 
beta2= -1 + 2* rand();  % -1 < beta2 < 1    Used in Eq. 2,3, and 4
naIni= 2; % minimum number of dingoes that will attack
naEnd= SearchAgents_no /naIni; % maximum number of dingoes that will attack
na= round(naIni + (naEnd-naIni) * rand()); % number of dingoes that will attack, used in Attack.m Section 2.2.1: Group attack
Positions=initialization(SearchAgents_no,dim,ub,lb);
 for i=1:size(Positions,1)
      Fitness(i)=fobj(Positions(i,:)); % get fitness     
 end
[vMin minIdx]= min(Fitness);  % the min fitness value vMin and the position minIdx
theBestVct= Positions(minIdx,:);  % the best vector
[vMax maxIdx]= max(Fitness); % the max fitness value vMax and the position maxIdx
Convergence_curve=zeros(1,Max_iter);
Convergence_curve(1)= vMin;
survival= survival_rate(Fitness,vMin,vMax);  % Section 2.2.4 Dingoes'survival rates
t=0;% Loop counter
% Main loop
for t=1:Max_iter       
   for r=1:SearchAgents_no
      sumatory=0;
    if rand() < P  % If Hunting?
           sumatory= Attack( SearchAgents_no, na, Positions, r );     % Section 2.2.1, Strategy 1: Part of Eq.2   
           if rand() < Q  % If group attack?                
                 v(r,:)=  beta1 * sumatory-theBestVct; % Strategy 1: Eq.2
           else  %  Persecution
               r1= round(1+ (SearchAgents_no-1)* rand()); % 
               v(r,:)= theBestVct + beta1*(exp(beta2))*((Positions(r1,:)-Positions(r,:))); % Section 2.2.2, Strategy 2:  Eq.3
           end  
    else % Scavenger
        r1= round(1+ (SearchAgents_no-1)* rand());
        v(r,:)=   (exp(beta2)* Positions(r1,:)-((-1)^getBinary)*Positions(r,:))/2; % Section 2.2.3, Strategy 3: Eq.4
    end
    if survival(r) <= 0.3  % Section 2.2.4, Algorithm 3 - Survival procedure
         band=1; 
         while band 
           r1= round(1+ (SearchAgents_no-1)* rand());
           r2= round(1+ (SearchAgents_no-1)* rand());
           if r1 ~= r2 
               band=0;
           end
         end
              v(r,:)=   theBestVct + (Positions(r1,:)-((-1)^getBinary)*Positions(r2,:))/2;  % Section 2.2.4, Strategy 4: Eq.6
    end 
     % Return back the search agents that go beyond the boundaries of the search space .  
        Flag4ub=v(r,:)>ub;
        Flag4lb=v(r,:)<lb;
        v(r,:)=(v(r,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    % Evaluate new solutions
     Fnew= fobj(v(r,:));
     % Update if the solution improves
     if Fnew <= Fitness(r)
        Positions(r,:)= v(r,:);
        Fitness(r)= Fnew;
     end
     if Fnew <= vMin
         theBestVct= v(r,:);
         vMin= Fnew;
     end 
   end
   Convergence_curve(t+1)= vMin; 
   [vMax maxIdx]= max(Fitness);
    survival= survival_rate( Fitness, vMin, vMax); % Section 2.2.4 Dingoes'survival rates
 end
%_____________________________________________________End DOA Algorithm]%
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

function [ o ] =  survival_rate(  fit, min, max )
    for i=1:size(fit,2)
         o(i)= (max-fit(i))/(max-min);
    end
end

function [ val] = getBinary( )
if rand() < 0.5
     val= 0;
else
     val=1;
end
end

function [ sumatory ] = Attack( SearchAgents_no, na, Positions, r )
sumatory=0;
vAttack= vectorAttack( SearchAgents_no, na );
     for j=1:size(vAttack,2)
           sumatory= sumatory + Positions(vAttack(j),:)- Positions(r,:);                       
     end 
     sumatory=sumatory/na;
end

function [ vAttack ] = vectorAttack( SearchAgents_no,na )
c=1; 
vAttack=[];
 while(c<=na)
    idx =round( 1+ (SearchAgents_no-1) * rand());
    if ~findrep(idx, vAttack)
        vAttack(c) = idx;
        c=c+1;
    end
 end
end

function [ band ] = findrep( val, vector )
% return 1= repeated  0= not repeated
band= 0;
for i=1:size(vector, 2)
    if val== vector(i)
        band=1;
        break;
    end
end
end