# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 18:37:38 2019

@author: Jens
"""

def timesteps(u, a, dx):
#   Calculate dt based on CFL limit
#   dt = C*dx/(u+a)    
#
#
    C = 0.75
    dt = C*dx/(u+a)
    dt_min = min(dt)
    
    return dt_min
    