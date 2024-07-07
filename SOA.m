
%%% Designed and Developed by Dr. Gaurav Dhiman (http://dhimangaurav.com/) %%%


function[Score,Position,Convergence]=SOA(Search_Agents,Max_iterations,Lower_bound,Upper_bound,dimension,objective)

Position=zeros(1,dimension);
Score=inf; 

Positions=init(Search_Agents,dimension,Upper_bound,Lower_bound);

Convergence=zeros(1,Max_iterations);

l=0;

while l<Max_iterations
    for i=1:size(Positions,1)  
        
        Flag4Upper_bound=Positions(i,:)>Upper_bound;
        Flag4Lower_bound=Positions(i,:)<Lower_bound;
        Positions(i,:)=(Positions(i,:).*(~(Flag4Upper_bound+Flag4Lower_bound)))+Upper_bound.*Flag4Upper_bound+Lower_bound.*Flag4Lower_bound;               
        
        fitness=objective(Positions(i,:));
        
        if fitness<Score 
            Score=fitness; 
            Position=Positions(i,:);
        end
        

    end
    
    
    Fc=2-l*((2)/Max_iterations); 
    
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
                       
            r1=rand(); 
            r2=rand(); 
            
            A1=2*Fc*r1-Fc; 
            C1=2*r2; 
            b=1;             
            ll=(Fc-1)*rand()+1;  
       
            D_alphs=Fc*Positions(i,j)+A1*((Position(j)-Positions(i,j)));                   
            X1=D_alphs*exp(b.*ll).*cos(ll.*2*pi)+Position(j);
            Positions(i,j)=X1;
            
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

