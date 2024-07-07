clear all 
clc
Function_name='F15';

% 初始化
SearchAgents_no=30; % 种群数量
 % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
Max_iteration=1000; % 最大迭代次数

%加载函数细节
% lb为下界；ub为上届；dim为狼的寻值范围
[lb,ub,dim,fobj] = Get_F (Function_name);
[Best_score,Best_pos,GWO_cg_curve]= GWO (SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
[aBest_score,aBest_pos,aGWO_cg_curve]= AGWO (SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%绘图，固定搜索区间
% figure('Position',[500 500 660 290])

%Draw search space
subplot(1,2,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])

%绘制搜索区间
subplot(1,2,2);
%绘制半对数坐标图
semilogy(GWO_cg_curve,'Color','r')
hold on
semilogy(aGWO_cg_curve,'Color','y')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

%自动调整坐标显示范围
axis tight
% 网格开
grid on
%绘制边框
box on
%图例
legend('GWO','agwo')

display(['The best solution obtained by GWO is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by GWO is : ', num2str(Best_score)]);
