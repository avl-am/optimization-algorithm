# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 15:16:30 2023

@author: ypzhao
"""

import numpy as np
import matplotlib.pyplot as plt
from function import Get_F

def func_plot(func_name):
    lb, ub, dim, fobj = Get_F(func_name)
    switch_func = {
        'F1': (-100, 100),
        'F2': (-10, 10),
        'F3': (-100, 100),
        'F4': (-100, 100),
        'F5': (-5, 5),
        'F6': (-100, 100),
        'F7': (-1, 1),
        'F8': (-500, 500),
        'F9': (-5, 5),
        'F10': (-500, 500),
        'F11': (-500, 500),
        'F12': (-10, 10),
        'F13': (-5, 5),
        'F14': (-100, 100),
        'F15': (-5, 5),
        'F16': (-5, 5),
        'F17': (-5, 5),
        'F18': (-5, 5),
        'F19': (-5, 5),
        'F20': (-5, 5),
        'F21': (-5, 5),
        'F22': (-5, 5),
        'F23': (-5, 5)
    }

    x = np.arange(switch_func[func_name][0], switch_func[func_name][1], 0.1)
    y = np.copy(x)

    L = len(x)
    f = np.zeros((L, L))

    # for i in range(L):
    #     for j in range(L):
    #         if func_name != 'F15' and func_name != 'F19' and func_name != 'F20' and func_name != 'F21' and func_name != 'F22' and func_name != 'F23':
    #             f[i, j] = fobj([x[i], y[j]])
    #         if func_name == 'F15':
    #             f[i, j] = fobj([x[i], y[j], 0, 0])
    #         if func_name == 'F19':
    #             f[i, j] = fobj([x[i], y[j], 0])
    #         if func_name == 'F20':
    #             f[i, j] = fobj([x[i], y[j], 0, 0, 0, 0])
    #         if func_name == 'F21' or func_name == 'F22' or func_name == 'F23':
    #             f[i, j] = fobj([x[i], y[j], 0, 0])

    X, Y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, f, cmap='viridis', edgecolor='none')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
