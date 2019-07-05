# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 00:15:54 2019

@author: Jens
"""

from scipy import sqrt


def speed_of_sound( p, dens ):
# --------------------------------------------
#   Inputs
#   p    -- pressure 
#           units: N/m**2
#    
#   dens -- density
#           units: kg/m**3    
#   
#   Outputs
#   a    -- speed of sound
#           units: m/s    
# --------------------------------------------        
    
#   Equation 19.114 in Zucrow, Hoffman Gas Dynamics Volume 2    
    a = sqrt( 1.2*p / dens)
    
    return a