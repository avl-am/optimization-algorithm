
% You can simply define your cost in a seperate file and load its handle to fobj 
% 需要初始化的参数: 
%__________________________________________
% fobj = @目标函数
% dim = 变量的数量
% Max_iteration = 最大迭代次数
% SearchAgents_no = 灰狼数量
% lb=[lb1,lb2,...,lbn] 变量n的下界
% ub=[ub1,ub2,...,ubn] 变量n的上界

% 如果所有变量的下界都相等，可以把lb和ub定义为两个单独的数字

% To run GWO: [Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%__________________________________________

clear all;
clc;

SearchAgents_no=30; % Number of search agents

Function_name='F1'; % 测试函数名称 f1-f24

Max_iteration=500; % 最大迭代次数

% 加载所选基准函数的详细信息
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

[Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

PSO_cg_curve=PSO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); % run PSO to compare to results

figure('Position',[500 500 660 290])
%Draw search space
subplot(1,2,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])

%Draw objective space
subplot(1,2,2);
semilogy(GWO_cg_curve,'Color','r')
hold on
semilogy(PSO_cg_curve,'Color','b')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('GWO','PSO')

display(['The best solution obtained by GWO is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by GWO is : ', num2str(Best_score)]);

        



