# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 14:59:32 2023

@author: ypzhao
"""

import numpy as np
def initialization(SearchAgents_no, dim, ub, lb):
    Boundary_no = np.size(ub)  # number of boundaries

    # If the boundaries of all variables are equal and user enter a single
    # number for both ub and lb
    if Boundary_no == 1:
        Positions = np.random.rand(SearchAgents_no) * (ub - lb) + lb

    # If each variable has a different lb and ub
    if Boundary_no > 1:
        Positions = np.zeros((SearchAgents_no, dim))
        for i in range(dim):
            if isinstance(ub, int):
                ub = [ub] * dim  # 将整数转换为长度为 dim 的列表
                ub_i = ub[i]
                lb_i = lb[i]
            Positions[:, i] = np.random.rand(SearchAgents_no) * (ub_i - lb_i) + lb_i

    return Positions
