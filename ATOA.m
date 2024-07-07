
function [Best_FF,Best_P,Conv_curve]=ATOAcossin(N,M_Iter,LB,UB,Dim,F_obj)
%display('ATOA Working');
%Two variables to keep the positions and the fitness value of the best-obtained solution

Best_P=zeros(1,Dim);
Best_FF=inf;
Conv_curve=zeros(1,M_Iter);

%Initialize the positions of solution
X=initialization(N,Dim,UB,LB);
Xnew=X;
Ffun=zeros(1,size(X,1));% (fitness values)
Ffun_new=zeros(1,size(Xnew,1));% (fitness values)

MOP_Max=1;
MOP_Min=0.2;
C_Iter=1;
Alpha=5;
Mu=0.499;


for i=1:size(X,1)
    Ffun(1,i)=F_obj(X(i,:));  %Calculate the fitness values of solutions
    if Ffun(1,i)<Best_FF
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
    end
end
    
    

while C_Iter<M_Iter+1  %Main loop
    MOP=1-((C_Iter)^(1/Alpha)/(M_Iter)^(1/Alpha));   % Probability Ratio 
    MOA=MOP_Min+C_Iter*((MOP_Max-MOP_Min)/M_Iter); %Accelerated function
   
    %Update the Position of solutions
    for i=1:size(X,1)   % if each of the UB and LB has a just value 
        for j=1:size(X,2)
           r1=rand();
            if (size(LB,2)==1)
                if r1<MOA
                    r2=rand();
                    if r2>0.5
                        Xnew(i,j)=Best_P(1,j)/(MOP+eps)*((cos(UB-LB))*Mu+cos(LB));
                    else
                        Xnew(i,j)=Best_P(1,j)*MOP*((cos(UB-LB))*Mu+cos(LB));
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,j)=Best_P(1,j)-MOP*((sin(UB-LB))*Mu+sin(LB));
                    else
                        Xnew(i,j)=Best_P(1,j)+MOP*((sin(UB-LB))*Mu+sin(LB));
                    end
                end               
            end
            
           
            if (size(LB,2)~=1)   % if each of the UB and LB has more than one value 
                r1=rand();
                if r1<MOA
                    r2=rand();
                    if r2>0.5
                        Xnew(i,j)=Best_P(1,j)/(MOP+eps)*(cos((UB(j)-LB(j)))*Mu+cos(LB(j)));
                    else
                        Xnew(i,j)=Best_P(1,j)*MOP*(cos((UB(j)-LB(j)))*Mu+cos(LB(j)));
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,j)=Best_P(1,j)-MOP*(sin((UB(j)-LB(j)))*Mu+sin(LB(j)));
                    else
                        Xnew(i,j)=Best_P(1,j)+MOP*(sin((UB(j)-LB(j)))*Mu+sin(LB(j)));
                    end
                end               
            end
            
        end
        
        Flag_UB=Xnew(i,:)>UB; % check if they exceed (up) the boundaries
        Flag_LB=Xnew(i,:)<LB; % check if they exceed (down) the boundaries
        Xnew(i,:)=(Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB;
 
        Ffun_new(1,i)=F_obj(Xnew(i,:));  % calculate Fitness function 
        if Ffun_new(1,i)<Ffun(1,i)
            X(i,:)=Xnew(i,:);
            Ffun(1,i)=Ffun_new(1,i);
        end
        if Ffun(1,i)<Best_FF
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
        end
       
    end
    

    %Update the convergence curve
    Conv_curve(C_Iter)=Best_FF;
    
    %Print the best solution details after every 50 iterations
    if mod(C_Iter,50)==0
        %display(['At iteration ', num2str(C_Iter), ' the best solution fitness is ', num2str(Best_FF)]);
    end
     
    C_Iter=C_Iter+1;  % incremental iteration
   
end
end

function X=initialization(N,Dim,UB,LB)

B_no= size(UB,2); % numnber of boundaries

if B_no==1
    X=rand(N,Dim).*(UB-LB)+LB;
end

% If each variable has a different lb and ub
if B_no>1
    for i=1:Dim
        Ub_i=UB(i);
        Lb_i=LB(i);
        X(:,i)=rand(N,1).*(Ub_i-Lb_i)+Lb_i;
    end
end

end

function F = tracklsq(pid)
         % Track the output of optsim to a signal of 1
        
         % Variables a1 and a2 are shared with RUNTRACKLSQ
         Kp = pid(1);
         Ki = pid(2);
         %Kd = pid(3);
         L = pid(3);
        % Kd = pid(2);
         
         %sprintf('The value of interation Kp= %3.0f, Ki= %3.0f', pid(1),pid(2)); 
         % Compute function value
         simopt = simset('solver','ode45','SrcWorkspace','Current','DstWorkspace','Current');  % Initialize sim options
         yout=sim('arunpressure',[0 500],simopt);
         F=ITAE;
         e=yout-1;
         sys_overshoot=max(yout)-1;
         %          
%          [tout,xout,yout] = sim('optiModel',[0 3],simopt);
%          e=yout-1 ;  % compute the error 
%          sys_overshoot=max(yout)-1; % compute the overshoot
%          
%       alpha=10;beta=10;
%       F=e2*beta+sys_overshoot*alpha;
         
    end



