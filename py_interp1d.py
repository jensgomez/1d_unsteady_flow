# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 00:15:54 2019

@author: Jens
"""

from scipy import interp

def interp1d( xi, xn, yn ):
# --------------------------------------------
#   Inputs
#   x1 -- value to interpolate
#   xn -- 1D array representing x-values
#   yn -- 1D array representing y-values
    
#   Outputs
#   yi -- interpolated value    
# --------------------------------------------    
    yi = interp( xi, xn, yn )
    
    return yi
