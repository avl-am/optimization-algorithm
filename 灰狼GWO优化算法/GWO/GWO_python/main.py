# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 15:19:27 2023

@author: ypzhao
"""
from function import Get_F
from GWO import GWO
from func_plot import func_plot 
import numpy as np
import matplotlib.pyplot as plt

#main函数
SearchAgents_no = 30  # #寻值的狼的数量
Function_name = 'F4'  # Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
Max_iteration = 50  # 最大迭代次数

# lb为下界；ub为上届；dim为狼的寻值范围
lb, ub, dim, fobj = Get_F(Function_name)

# 加载函数细节
# x = GWO(fobj, lb, ub, dim, SearchAgents_no, Max_iteration)
[Best_score,Best_pos,GWO_cg_curve]= GWO (SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

# Draw search space
fig = plt.figure(figsize=(10, 5))
ax1 = fig.add_subplot(121, projection='3d')
func_plot(Function_name)
ax1.set_title("Parameter space")
ax1.set_xlabel('x_1')
ax1.set_ylabel('x_2')
ax1.set_zlabel(Function_name + '( x_1 , x_2 )')


Best_score, Best_pos, GWO_cg_curve = GWO(SearchAgents_no, Max_iteration, lb, ub, dim, fobj)
# Draw objective space
ax2 = fig.add_subplot(122)
ax2.semilogy(GWO_cg_curve, color='r')
ax2.set_title("Objective space")
ax2.set_xlabel('Iteration')
ax2.set_ylabel('Best score obtained so far')

plt.tight_layout()
plt.grid(True)
plt.legend(['GWO'])

print('The best solution obtained by GWO is:', Best_pos)
print('The best optimal value of the objective funciton found by GWO is:', Best_score)
