clear all 
clc

Solution_no=20;  
F_name='F3';    
M_Iter=1000;    
tic;
RT=2;
% N=20;
[LB,UB,Dim,F_obj]=Get_F(F_name); 
[Best_FF,Best_P,conv]=GWCA(Solution_no,M_Iter,LB,UB,Dim,F_obj);  

figure('Position',[454   445   694   297]);
subplot(1,3,1);
func_plot(F_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([F_name,'( x_1 , x_2 )'])

subplot(1,2,2);
semilogy(conv,'Color','r','LineWidth',1)
title('Convergence curve')
xlabel('Iteration#');
ylabel('Best fitness function');
axis tight

display(['The best-obtained solution by GWCA is : ', num2str(Best_P)]);
display(['The best optimal values of the objective funciton found by GWCA is : ', num2str(Best_FF)]);    