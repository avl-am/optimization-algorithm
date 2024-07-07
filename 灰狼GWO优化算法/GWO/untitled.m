% 主程序 GWO
clear
close all
clc
 
%%完整代码见微信公众号：电力系统与算法之美
 
%输入关键字：灰狼算法
 
SearchAgents_no = 30 ; %　种群规模
dim = 10 ; % 粒子维度
Max_iter = 1000 ; %　迭代次数
ub = 5 ;
lb = -5 ;
 
%% 初始化三匹头狼的位置
Alpha_pos=zeros(1,dim);
Alpha_score=inf; 
 
Beta_pos=zeros(1,dim);
Beta_score=inf; 
 
Delta_pos=zeros(1,dim);
Delta_score=inf; 
 
 
 
Convergence_curve = zeros(Max_iter,1);
 
%% 开始循环
for l=1:Max_iter
    for i=1:size(Positions,1)  
        
       %% 返回超出搜索空间边界的搜索代理
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
        %% 计算每个搜索代理的目标函数
        fitness=sum(Positions(i,:).^2);
        
        %% 更新 Alpha, Beta, and Delta
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
    end
    
    
    a=2-l*((2)/Max_iter); % a decreases linearly fron 2 to 0
    
    %% 更新搜索代理的位置，包括omegas
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
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
  
    Convergence_curve(l)=Alpha_score;
    disp(['Iteration = ' num2str(l)  ', Evaluations = ' num2str(Alpha_score)]);
 
end
%========可视化==============
figure('unit','normalize','Position',[0.3,0.35,0.4,0.35],'color',[1 1 1],'toolbar','none')
%% 目标空间
subplot(1,2,1);
x = -5:0.1:5;y=x;
L=length(x);
f=zeros(L,L);
for i=1:L
    for j=1:L
       f(i,j) = x(i)^2+y(j)^2;
    end
end
surfc(x,y,f,'LineStyle','none');
xlabel('x_1');
ylabel('x_2');
zlabel('F')
title('Objective space')
%% 狼群算法 
subplot(1,2,2);
semilogy(Convergence_curve,'Color','r','linewidth',1.5)
title('Convergence_curve')
xlabel('Iteration');
ylabel('Best score obtained so far');
 
axis tight
grid on
box on
legend('GWO')
display(['The best solution obtained by GWO is : ', num2str(Alpha_pos)]);
display(['The best optimal value of the objective funciton found by GWO is : ', num2str(Alpha_score)]);
 