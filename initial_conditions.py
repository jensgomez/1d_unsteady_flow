# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 14:07:26 2019

@author: Jens
"""

from scipy import array

def initial_conditions( ):
#    
#   Initial conditions
#   v     = [m/s]    
#   rho   = [kg/m**3]
#   p     = [N/m**2]
#   x     = m

    x0   = array([0.10160, 0.12700, 0.15240])
    u0   = array([188.280, 235.300, 282.530]) 
    p0   = array([111.770E5, 110.780E5, 109.540E5])
    rho0 = array([14.0640, 13.9590, 13.8270])

    return [x0, u0, p0, rho0]
    