# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 20:17:16 2018

@author: Jens
"""

from scipy import sqrt


class unsteady:
    def __init__(self):       
        self.x = []
        self.p = []
        self.u = []
        self.a = 0.0
        self.ti = 0.0
        self.pi = 0.0
        
        
    def calculate ( self ):
        node = 0
        
        
    def inputs ( self, x, u, p ):
        self.x = x
        self.u = u
        self.p = p
        
    def initialcond ( self, ti, pi ):
        self.ti = ti
        self.pi = pi

#   Linear 1-D interpolation    
    def interpolate( x1, x2, x3, y1, y2 ):
        
        mu = ( y1 - y2 )/( x1 - x2 )
        bu = y2 - ( mu * x2 )
        
        y3 = mu*x3 + bu
        
        return y3


#   Calculate speed of sound        
    def calcspeedofsound ( self, p ):
#       Gas constant 
#       J/(kg-K)       
        R = 320.0        
        
        gamma_gas = 1.20 
        
#       Equation 13.30(a)        
#       m/s
        a = sqrt((gamma_gas*R*self.ti)) * (p/self.pi)**((gamma_gas-1.0)/(2.0*gamma_gas))
        
        echo = 0
        if ( echo == 1 ):
            
            print(" R     = {:8.3f}".format(R))
            print(" γ     = {:8.3f}".format(gamma_gas))
            print(" P     = {:8.3f}".format(p))
            print(" Tinit = {:8.3f}".format(self.ti))            
            print(" Pinit = {:8.3f}".format(self.pi))
            print("")
            print(" a     = {:8.3f}".format(a))
            print("")
            print("")
            
            
        return a


#   Calculate density
    def calcdensity ( self, p ):
#       Gas constant 
#       J/(kg-K)       
        R = 320.0        
        
        gamma_gas = 1.20 
        
        rho = self.pi*10**5/(R*self.ti)*(p/self.pi)**(1/gamma_gas)
        
        echo = 0
        if ( echo == 1 ):
            
            print(" R     = {:8.3f}".format(R))
            print(" γ     = {:8.3f}".format(gamma_gas))
            print(" P     = {:8.3f}".format(p))
            print(" Tinit = {:8.3f}".format(self.ti))                        
            print(" Pinit = {:8.3f}".format(self.pi))
            print("")
            print(" dens  = {:8.3f}".format(rho))
            print("")
            print("")            
        
        
        return rho 
    
    def predictor ( self ):
        
        
        


# Inputs

# m
x = [0.10, 0.125, 0.150]

# m/s
u = [188.30, 235.30, 282.50]

# N/m**2 (10E+05)
p = [111.77, 110.78, 109.54]


# Initial conditions
# K
ti = 3330.0

# N/m**2 (10E+05)
#
pi = 690 

# milli-second
dt = 0.020



test = unsteady()
test.inputs(x, u, p)
test.initialcond(ti, pi)
test.calcspeedofsound( p[0] )
test.calcdensity( p[0] )




