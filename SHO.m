
function [Best_hyena_score,Best_hyena_pos,Convergence_curve]=SHO(N,Max_iterations,lowerbound,upperbound,dimension,fitness)

hyena_pos=init(N,dimension,upperbound,lowerbound);

Convergence_curve=zeros(1,Max_iterations);

Iteration=1;


while Iteration<Max_iterations
    
    
    for i=1:size(hyena_pos,1)
        

        H_ub=hyena_pos(i,:)>upperbound;
        H_lb=hyena_pos(i,:)<lowerbound;
        hyena_pos(i,:)=(hyena_pos(i,:).*(~(H_ub+H_lb)))+upperbound.*H_ub+lowerbound.*H_lb;  
        hyena_fitness(1,i)=fitness(hyena_pos(i,:));  
        
    end
       
    if Iteration==1
        [fitness_sorted FS]=sort(hyena_fitness);
        sorted_population=hyena_pos(FS,:);
        best_hyenas=sorted_population;
        best_hyena_fitness=fitness_sorted;
        
    else
        double_population=[pre_population;best_hyenas];
        double_fitness=[pre_fitness best_hyena_fitness];
        [double_fitness_sorted FS]=sort(double_fitness);
        double_sorted_population=double_population(FS,:);
        fitness_sorted=double_fitness_sorted(1:N);
        sorted_population=double_sorted_population(1:N,:);
        best_hyenas=sorted_population;
        best_hyena_fitness=fitness_sorted;
    end
    
    NOH=noh(best_hyena_fitness);
    
   
    Best_hyena_score=fitness_sorted(1);
    Best_hyena_pos=sorted_population(1,:);
    pre_population=hyena_pos;
    pre_fitness=hyena_fitness;
    
    a=5-Iteration*((5)/Max_iterations);
    HYE=0;
    CV=0;
    for i=1:size(hyena_pos,1)
        
        for j=1:size(hyena_pos,2)
             for k=1:NOH
            HYE=0;
             r1=rand();
             r2=rand(); 
            Var1=2*a*r1-a; 
            Var2=2*r2; 
            distance_to_hyena=abs(Var2*sorted_population(k)-hyena_pos(i,j));
            HYE=sorted_population(k)-Var1*distance_to_hyena;
            CV=CV+HYE;        
            distance_to_hyena=0;
             end
        
        hyena_pos(i,j)=(CV/(NOH+1));
        CV=0;
        end
        
    end
    
    Convergence_curve(Iteration)=Best_hyena_score;
    
     Iteration=Iteration+1; 
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

function X=noh(best_hyena_fitness)
min = 0.5;
max = 1;
count=0;
M=(max-min).*rand(1,1) + min;
M=M+best_hyena_fitness(1);

for i=2:numel(best_hyena_fitness)
    if M>=best_hyena_fitness(i)
        count=count+1;
    end
    
end
 
X=count;
clear count
clear M
clear px
count=0;
end





