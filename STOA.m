%################################################################################################%
% Designed and Developed by Dr. Gaurav Dhiman (http://dhimangaurav.com/)                         %
% Doi: https://doi.org/10.1016/j.engappai.2019.03.021                                            %
% Link: https://www.sciencedirect.com/science/article/abs/pii/S0952197619300715                  %
% Title: STOA: A bio-inspired based optimization algorithm for industrial engineering problems   %
% Authors: Gaurav Dhiman and Amandeep Kaur                                                       %
%################################################################################################%


function[Score,Position,Convergence]=STOA(Search_Agents,Max_iterations,Lb,Ub,dimension,objective)

Position=zeros(1,dimension);
Score=inf; 

Positions=init(Search_Agents,dimension,Ub,Lb);

Convergence=zeros(1,Max_iterations);

l=0;

while l<Max_iterations
    for i=1:size(Positions,1)  
        
        Flag4Ub=Positions(i,:)>Ub;
        Flag4Lb=Positions(i,:)<Lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4Ub+Flag4Lb)))+Ub.*Flag4Ub+Lb.*Flag4Lb;               
        
        fitness=objective(Positions(i,:));
        
        if fitness<Score 
            Score=fitness; 
            Position=Positions(i,:);
        end
        

    end
    
    
    Sa=2-l*((2)/Max_iterations); 
    
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
            r1=0.5*rand(); 
            r2=0.5*rand();             
            X1=2*Sa*r1-Sa;  
            b=1;             
            ll=(Sa-1)*rand()+1;         
            D_alphs=Sa*Positions(i,j)+X1*((Position(j)-Positions(i,j)));                   
            res=D_alphs*exp(b.*ll).*cos(ll.*2*pi)+Position(j);
            Positions(i,j)=res;            
        end
    end
    l=l+1;    
    Convergence(l)=Score;
end
end

function Pos=init(SearchAgents,dimension,upperbound,lowerbound)

Boundary= size(upperbound,2); 
if Boundary==1
    Pos=rand(SearchAgents,dimension).*(upperbound-lowerbound)+lowerbound;
end

if Boundary>1
    for i=1:dimension
        ub_i=upperbound(i);
        lb_i=lowerbound(i);
        Pos(:,i)=rand(SearchAgents,1).*(ub_i-lb_i)+lb_i;
    end
end
end


