# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 14:16:15 2019

@author: Jens
"""
from import_header import initial_conditions as _IC
from import_header import py_speed_of_sound as _speedofsound
from import_header import py_CFL
from import_header import py_interp1d 

from scipy import array


class interior:
    def __init__(self):
        self.x   = array([])
        self.p   = array([])
        self.u   = array([])
        self.rho = array([])
        
    def initialcond(self, data):       
        self.x0    = array(data[0])
        self.u0    = array(data[1])
        self.p0    = array(data[2])
        self.rho0  = array(data[3])
        
        self.n = len(self.x0)
#       Calculate speed of sound        
        self.a0 = _speedofsound.speed_of_sound(self.p0, self.rho0)
        
#       Calculate minimum time step based on CFL
        self.dt = py_CFL.timesteps(self.u0, self.a0, self.x0)

#       Set initial condition values to domain
        self.u   = self.u0
        self.p   = self.p0
        self.x   = self.x0
        self.rho = self.rho0
        
#               
        self.dt  = 0.016604E-3
        
        echo = 0
        if ( echo == 1 ):            
            print("X0     = ", self.x0)
            print("U0     = ", self.u0)
            print("P0     = ", self.p0)
            print("RHO0   = ", self.rho0)
            print("A0     = ", self.a0)
            print("DT     = ", self.dt)

        
    def iterate_func(self, x0, i):
#       Iterate to find x
        tol = 1E-6
        diff = 1
        xinitial = x0
        xold     = x0
        x1       = xinitial
        k_counter = 1
        
        if ( x0 < self.x[i-1] or x0 > self.x[i]):
            print("ERROR. X0 NOT IN DOMAIN OF (I-1) AND (I).\n")            
       
        echo = 1        
        if ( echo == 1 ):
            print("Inputs: ")
            print("Grid point   X(i-1) = {0:10.8E}".format(self.x[i-1]))
            print("Intial value X      = {0:10.8E}".format(xinitial))            
            print("Grid point   X(i+1) = {0:10.8E}\n\n".format(self.x[i]))            
            print("Calculations: ")
                        
        
        while ( diff > tol ): 
                        
            u   = py_interp1d.interp1d(x1, (self.x[i-1], self.x[i]), (  self.u[i-1], self.u[i]))
            p   = py_interp1d.interp1d(x1, (self.x[i-1], self.x[i]), (  self.p[i-1], self.p[i]))
            rho = py_interp1d.interp1d(x1, (self.x[i-1], self.x[i]), (self.rho[i-1], self.rho[i]))            
            a   = _speedofsound.speed_of_sound(p, rho)
                      
            lamb = 1.0/(u + a)
            
            x1 = self.x[i] - (self.dt/lamb)
            
            diff = abs(x1 - xold)
            xold = x1
            
            k_counter += 1           

            
            if ( k_counter > 20 ):
                print("ERROR. CONVERGENCE CRITERIA FAILED.\n")
                break
            
            if ( echo == 1 ):
                print("Counter = {0}".format(k_counter))
                print("U          = {0:10.8E}".format(u))
                print("P          = {0:10.8E}".format(p))
                print("RHO        = {0:10.8E}".format(rho))
                print("A          = {0:10.8E}".format(a))
                print("X          = {0:10.8E}".format(x1))
                print("XOLD       = {0:10.8E}".format(xold))
                print("DIFF       = {0:10.8E}\n\n".format(diff))
                
        if ( echo == 1 ):
            print("FINAL X    = {0:10.8E}\n\n".format(x1))
                
            
        return (x1)

    def solve(self):
        """
            Solver for interior points
            Points are layed out spatially 
            
  
            Unknown point:                4

            Known points:        5,       6,       7
            Interpolated points:   1,  3,      2    
        """         
        
#       Locate point 1
        u_p   = self.u[0]
        p_p   = self.p[0]
        rho_p = self.rho[0]
        
        a_p   = _speedofsound.speed_of_sound(p_p, rho_p)
        
        lambda_p = 1.0/(u_p + a_p)
        
        
#       Locate point based on x4
        x4 = self.x[1]
        x1 = x4 - self.dt/lambda_p
        
                
#       Re-calculate a_p and iterate until x1 is below tolerance        
        x1   = self.iterate_func(x1, 1)
        
#       Interpolate based on x1
        u1   = py_interp1d.interp1d(x1, (self.x[0], self.x[1]), (  self.u[0], self.u[1]))
        p1   = py_interp1d.interp1d(x1, (self.x[0], self.x[1]), (  self.p[0], self.p[1]))
        rho1 = py_interp1d.interp1d(x1, (self.x[0], self.x[1]), (self.rho[0], self.rho[1]))        
        

            
            
        
            
            
            
            
            
            
            
        
        
        
                    

                

        
        
            
    
            
        
i1 = interior()
i1.initialcond(_IC.initial_conditions())
i1.solve()


