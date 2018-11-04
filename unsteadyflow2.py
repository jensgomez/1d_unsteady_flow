# -*- coding: utf-8 -*-
"""
Created on Sat Nov  3 13:30:50 2018

@author: 212714696
"""

from scipy import sqrt
from scipy.linalg import solve
from scipy import array
from scipy import zeros

# Class to calculate unsteady flow of interior points
# for 1-D unsteady flow

class interior:
    
    def __init__ ( self ):

#       Initialize Velocity, pressure and position arrays
        self.u = []
        self.p = []
        self.x = []
       
#       Time step in milliseconds
        self.dt = 0.020
        
#       Max time to solve        
        self.tmax = 0.10       
        self.tcount = int ( self.tmax / self.dt )
        

#   Function to gather inputs, excluding initial conditions        
    def inputs ( self, u, p, x ):
        
        self.u = u
        self.p = p
        self.x = x
        
        self.len_u = len(u)
        self.len_p = len(p)
        self.len_x = len(x)
        
        self.u_arr = zeros(( self.tcount, self.len_u ))
        self.p_arr = zeros(( self.tcount, self.len_p ))
        self.x_arr = zeros(( self.tcount, self.len_x ))
        
        print(self.u_arr)
        
        
        
    def initialcond ( self, p, t ):
        self.pi = p
        self.ti = t
        
        
#   1D interpolation based on for x3 using: [x1, x2], [y1, y2]
#   y3 = f(x3,[x1, x2], [y1, y2])        
    def interp1d ( self, x1, x2, x3, y1, y2 ):
        
        mu = ( y1 - y2 )/( x1 - x2 )
        bu = y2 - ( mu * x2 )        
        y3 = mu*x3 + bu
        
        return y3        


#   Falculate speed of sound        
#   a = a(p)        
    def calc_speedofsound ( self, p ):
        
#       Gas constant 
#       J/(kg-K)       
        R = 320.0        
        
        gamma_gas = 1.20 
        
#       Equation 13.30(a)        
#       m/s

        a = sqrt( (gamma_gas*R*self.ti) ) * ( p/self.pi )**(( gamma_gas - 1.0 )/( 2.0 * gamma_gas ))
        
        
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
#   rho = rho(p)        
    def calc_density ( self, p ):
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
    
#   Calculate Q = rho*a 
#   Equation 13.42 in Gas Dynamics Volume 1; Zucrow, Hoffman          
    def calc_q ( self, p ):
        
        rho = self.calc_density( p )
        a = self.calc_speedofsound( p )       
        q = rho * a
        
        return q

#   Calculate T 
#   For C+:
#   T = p + q*u        
#   
#   For C-:
#   T = p - q*u
#        
#   Equation 13.40 in Gas Dynamics Volume 1; Zucrow, Hoffman              
    def calc_t ( self, p, u, opt ):

        q = self.calc_q ( p )        
        
#       For C+
        if ( opt == 1 ): 
            t = p*(10**5) + q*u
                        
#       For C-            
        elif ( opt == 2 ):
            t = p*(10**5) - q*u

        
        echo = 0
        if ( echo == 1 ):
            print(" *** INPUTS *** ")
            print("P = ", p)
            print("U = ", u)
            print("OPT = ", opt)
            print(" ")
            print("Q = ", q)
            print("T = ", t)
            print(" ")

        return t   



#   Main function to step through iterations 
    def calc_flow ( self ):
       
        if ( len(self.x) != len(self.p) or len(self.x) != len(self.u) ):
            print("ERROR")
            exit 
        
        t = 0.0 + self.dt 
        
        while ( t < self.tmax ):
            t = t + self.dt 

#           Locate Point 1
#           For the initial guess, use u[0], p[0]            
            self.u_pl = self.u[0]
            self.p_pl = self.p[0]
            
#           Use iterator function to calculate new values for C+
            [x1, u_pl, p_pl] = self.cp_iter( self.u[0], self.p[0], self.x[0], 0 )
            
#           Set u(+) and p(+) from results of C+
            self.u_pl = u_pl
            self.p_pl = p_pl
            
#           Calculate Q+ and T+
            self.q_pl = self.calc_q( self.p_pl )
            self.t_pl = self.calc_t( self.p_pl, self.u_pl, 1 )
            
#           Locate Point 2
#           For the initial guess, use u[0], p[0]            
            self.u_mn = self.u[2]
            self.p_mn = self.p[2]
            
#           Use iterator function to calculate new values for C-
            [x2, u_mn, p_mn] = self.cm_iter( self.u[2], self.p[2], self.x[0], 1 ) 
            
#           Set u(-) and p(-) from results of C-
            self.u_mn = u_mn
            self.p_mn = p_mn
            
#           Calculate Q- and T-
            self.q_mn = self.calc_q( self.p_mn )
            self.t_mn = self.calc_t( self.p_mn, self.u_mn, 2 )
            
                        
        a = [ [1.0, self.q_pl], [ 1.0, -self.q_mn  ]]
        b = [ [self.t_pl], [self.t_mn] ]
        [ p_fin, u_fin ] = solve( a, b )
        
        self.p_fin = float(p_fin)
        self.u_fin = float(u_fin)
        
        
            
#   Function to iterate the C+ characterisitic line
    def cp_iter ( self, u, p, x, i ):
        
#       Perform iterations until under tolerance        
        diff = 1.0 
        tol  = 1E-03
        
#       Iteration counter
        kcount = 0        
        x1 = 0.0
        
        while ( diff > tol ):
            
            kcount += 1
            
            old_x1 = x1
        
#           Calculate speed of sound            
            a = self.calc_speedofsound( p )
            
            lmda = 1.0/( a + u )
            lmda *= 1000
            x1 = self.x[i+1] - (self.dt/lmda)
                        
#           Determine what two points x1 lies between
            n = len(self.x)
            for i in range(0, n):
                if ( x1 > self.x[i] and x1 < self.x[i+1] ):
                    j = i
                    break
            
#           Determine new u and p based on x1
#           interp1d inputs: x1, x2, x3, y1, y,2
            u_new = self.interp1d( self.x[j], self.x[j+1], x1, self.u[j], self.u[j+1] )
            p_new = self.interp1d( self.x[j], self.x[j+1], x1, self.p[j], self.p[j+1] )
                        
#           Determine if x1 converged 
            diff = abs( old_x1 - x1 )
            
            echo = 0
            if ( echo == 1 ):
                
                if ( kcount == 1 ):
                    print(" *** INPUTS *** ")
                    print(" X              = {:6.5f} M     ".format(x))
                    print(" U              = {:8.4f} M/SEC ".format(u))
                    print(" P              = {:8.4f} N/M**2 *(10**5) ".format(p))
                    print(" ")
                    print(" ")
                    
                print( " ITERATION      = {:2} ".format(kcount))
                print( " SPEED OF SOUND = {:6.3f} M/SEC ".format(a))
                print( " λ(+)           = {:6.4f} MSEC/M ".format(lmda))
                print( " XNEW           = {:6.5f} M     ".format(x1))
                print( " UNEW           = {:6.3f} M/SEC ".format(u_new))
                print( " PNEW           = {:6.3f} N/M**2 *(10**5) ".format(p_new))
                print( " ")
                print( " ")    

        return [ x1, u_new, p_new ]


#   Function to iterate the C- characterisitic line
    def cm_iter ( self, u, p, x, i ):
        
#       Perform iterations until under tolerance        
        diff = 1.0 
        tol  = 1E-03
        
#       Iteration counter
        kcount = 0        
        x1 = 0.0
        
        while ( diff > tol ):
            
            kcount += 1
            
            old_x1 = x1
        
#           Calculate speed of sound            
            a = self.calc_speedofsound( p )
            
            lmda = 1.0/( u - a )
            lmda *= 1000
            x1 = self.x[i] - (self.dt/lmda)
            
            
#           Determine what two points x1 lies between
            n = len(self.x)
            for i in range(0, n):
                if ( x1 > self.x[i] and x1 < self.x[i+1] ):
                    j = i
                    break
            
#           Determine new u and p based on x1
#           interp1d inputs: x1, x2, x3, y1, y,2
            u_new = self.interp1d( self.x[j], self.x[j+1], x1, self.u[j], self.u[j+1] )
            p_new = self.interp1d( self.x[j], self.x[j+1], x1, self.p[j], self.p[j+1] )
                        
#           Determine if x1 converged 
            diff = abs( old_x1 - x1 )
            
            echo = 0
            if ( echo == 1 ):
                
                if ( kcount == 1 ):
                    print(" *** INPUTS *** ")
                    print(" X              = {:6.5f} M     ".format(x))
                    print(" U              = {:8.4f} M/SEC ".format(u))
                    print(" P              = {:8.4f} N/M**2 *(10**5) ".format(p))
                    print(" ")
                    print(" ")
                    
                print( " ITERATION      = {:2} ".format(kcount))
                print( " SPEED OF SOUND = {:6.3f} M/SEC ".format(a))
                print( " λ(-)           = {:6.4f} MSEC/M ".format(lmda))
                print( " XNEW           = {:6.5f} M     ".format(x1))
                print( " UNEW           = {:6.3f} M/SEC ".format(u_new))
                print( " PNEW           = {:6.3f} N/M**2 *(10**5) ".format(p_new))
                print( " ")
                print( " ")    

        return [ x1, u_new, p_new ]


    def output( self ):
        
        f = open( "outputUnsteady.out", "w")
        
        f.write(" *** INITIAL VALUES *** \n")       
        for i in range(0, len(self.x)):
            f.write(" X({:2} ) = {:5.4f} M\n".format(i, self.x[i]))
        f.write("\n")            
        for i in range(0, len(self.p)):    
            f.write(" P({:2} ) = {:8.4f} N/M**2 * (10**5)\n".format(i, self.p[i]))   
        f.write("\n")            
        for i in range(0, len(self.u)):
            f.write(" U({:2} ) = {:8.4f} M/SEC\n".format(i, self.u[i]))

        f.write("\n\n")            
        f.write(" *** INITIAL CONDITIONS *** \n")
        f.write(" P(INITIAL) = {:8.4f} N/M**2 * (10**5)\n".format(self.pi))
        f.write(" T(INITIAL) = {:8.4f} K\n\n".format(self.ti))
        
        f.write(" \n\n\n ")
        f.write(" *** TRANSIENT SOLUTION *** \n")
        f.write("  T    = {:5.4f} SEC\n ".format(self.dt))
        f.write(" P    = {:8.4f} N/M**2 * (10**5)\n ".format(self.p_fin/1E+05))
        f.write(" U    = {:8.4f} M/SEC ".format(self.u_fin))        
        
        
        
        
        

      
def main( ):
    
    x = [ 0.1000, 0.1250, 0.1500 ]
    u = [ 188.30, 235.30, 282.50 ]
    p = [ 111.77, 110.78, 109.54 ]
    
    pi = 690
    ti = 3330
    
    
    test = interior( )
    test.inputs( u, p, x )
    test.initialcond( pi, ti )
    test.calc_flow( )
    
    
    
    test.output( )
    

if __name__ == '__main__':
    main( ) 


    
    
    
    