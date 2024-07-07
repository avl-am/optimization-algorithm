function [Destination_fitness,Destination_position,Convergence_curve]=SCHO(N,Max_iteration,lb,ub,dim,fobj)

Destination_position=zeros(1,dim);
Destination_fitness=inf;
Destination_position_second=zeros(1,dim);
Convergence_curve=zeros(1,Max_iteration);
Position_sort = zeros(N,dim);
%Initialize SCHO parameters
u=0.388;
m=0.45;
n=0.5;
p=10;
q=9;
Alpha=4.6;
Beta=1.55;
BS=floor(Max_iteration/Beta);
ct=3.6;
T=floor(Max_iteration/ct);
BSi=0;
BSi_temp=0;
ub_2=ub;
lb_2=lb;
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);
Objective_values = zeros(1,size(X,1));
% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    Objective_values(1,i)=fobj(X(i,:));
    if Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
end
Convergence_curve(1)=Destination_fitness;
t=2; 
%Main loop
while t<=Max_iteration    
    for i=1:size(X,1) % in i-th solution
        for j=1:size(X,2) % in j-th dimension
            %update A by using Eq. (17)
            cosh2=(exp(t/Max_iteration)+exp(-t/Max_iteration))/2;
            sinh2=(exp(t/Max_iteration)-exp(-t/Max_iteration))/2;
            r1=rand();
            A=(p-q*(t/Max_iteration)^(cosh2/(sinh2)))*r1; 
            % enter the bounded search strategy
            if t==BSi
ub_2=Destination_position(j)+(1-t/Max_iteration)*abs(Destination_position(j)-Destination_position_second(j));
lb_2=Destination_position(j)-(1-t/Max_iteration)*abs(Destination_position(j)-Destination_position_second(j));
                if ub_2>ub
                    ub_2=ub;
                end
                if lb_2<lb
                    lb_2=lb;
                end
                X=initialization(N,dim,ub_2,lb_2);                
                BSi_temp=BSi;
                BSi=0;
            end
        % the first phase of exploration and exploitation    
        if t<=T%3.6-3.62  
            r2=rand();
            r3=rand();
            a1=3*(-1.3*t/Max_iteration+m);
            r4=rand();
            r5=rand();
            if A>1
                sinh=(exp(r3)-exp(-r3))/2;
                cosh=(exp(r3)+exp(-r3))/2;
                W1=r2*a1*(cosh+u*sinh-1);
                if r5<=0.5
                    X(i,j)=Destination_position(j)+r4*W1*X(i,j);
                else
                    X(i,j)=Destination_position(j)-r4*W1*X(i,j);  
                end                
            else
                sinh=(exp(r3)-exp(-r3))/2;
                cosh=(exp(r3)+exp(-r3))/2;
                W3=r2*a1*(cosh+u*sinh);
                if r5<=0.5
                    X(i,j)=Destination_position(j)+r4*W3*X(i,j);
                else
                    X(i,j)=Destination_position(j)-r4*W3*X(i,j);  
                end
            end
        else
            % the second phase of exploration and exploitation
            r2=rand();
            r3=rand();
            a2=2*(-t/Max_iteration+n);
            W2=r2*a2;
            r4=rand();
            r5=rand();
            if A<1
                sinh=(exp(r3)-exp(-r3))/2;
                cosh=(exp(r3)+exp(-r3))/2;
                X(i,j)= X(i,j)+(r5*sinh/cosh*abs(W2*Destination_position(j)-X(i,j)));
            else
                if r4<=0.5
                    X(i,j)=X(i,j)+(abs(0.003*W2*Destination_position(j)-X(i,j)));
                else
                    X(i,j)=X(i,j)+(-abs(0.003*W2*Destination_position(j)-X(i,j)));  
                end
            end 
        end
        end
        BSi=BSi_temp;
    end
    for i=1:size(X,1)         
        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=X(i,:)>ub_2;
        Flag4lb=X(i,:)<lb_2;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+(ub_2+lb_2)/2.*Flag4ub+lb_2.*Flag4lb;        
        % Calculate the objective values
        Objective_values(1,i)=fobj(X(i,:));
%         % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    end
    %find the second solution
    if t==BS
        BSi=BS+1;
        BS=BS+floor((Max_iteration-BS)/Alpha);
        temp = zeros(1,dim);
        temp2 = zeros(N,dim);
 %sorting
        for i=1:(size(X,1)-1)
          for j=1:(size(X,1)-1-i)
            if Objective_values(1,j) > Objective_values(1,j+1)
              temp(1,j) = Objective_values(1,j);
              Objective_values(1,j) = Objective_values(1,j+1);
              Objective_values(1,j+1) = temp(1,j);
              temp2(j,:) = Position_sort(j,:);
              Position_sort(j,:) = Position_sort(j+1,:);
              Position_sort(j+1,:) = temp2(j,:);   
                end
            end
        end
        Destination_position_second=Position_sort(2,:);%the second solution
    end
    Convergence_curve(t)=Destination_fitness;
    t=t+1;
end


function X=initialization(nP,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    X=rand(nP,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(nP,1).*(ub_i-lb_i)+lb_i;
    end
end