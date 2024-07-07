


%% 
% OOBO.
% OOBO: A New One to One Based Optimizer for Solving Various Optimization Problems

% " Optimizer"
%%
clc
clear
close all
%%
Fun_name='F2'; 
SearchAgents=30;  
Max_iterations=1000;
%%

[lowerbound,upperbound,dimension,fitness]=fun_info(Fun_name);
[Best_score,Best_pos,OOBO_curve]=OOBO(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);

%%
plots=semilogx(OOBO_curve,'Color','g');
set(plots,'linewidth',2)
hold on
title('Objective space')
xlabel('Iterations');
ylabel('Best score');

axis tight
grid on
box on
legend('OOBO')


display(['The best solution obtained by OOBO is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by OOBO is : ', num2str(Best_score)]);

        