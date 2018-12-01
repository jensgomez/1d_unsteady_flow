# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 20:35:05 2018

@author: Jens
"""
from scipy import sqrt
from scipy.linalg import solve
from scipy import array
from scipy import zeros

class unsteady:

    def __init__ ( self ):

#       Initialize Velocity, pressure and position arrays
        self.u   = []
        self.p   = []
        self.x   = []
        self.rho = []
       
#       Time step in milliseconds
        self.dt = 0.016604
        
#       Max time to solve and number of time steps
        self.tmax = 0.10       
        self.tcount = int ( self.tmax / self.dt )    
        
        
#   Function to gather inputs, excluding initial conditions        
    def inputs ( self, u, p, x, rho ):
        
        self.u   = u
        self.p   = p
        self.x   = x
        self.rho = rho 
        
        if ( len(self.u) == len(self.p) and len(self.p) == len(self.x) and len(self.p) == len(self.rho) ):
            n = len(self.p)
        
        else:
            print(" ERROR! ")
            exit 
        
        self.u_arr   = zeros(( self.tcount, n ))
        self.p_arr   = zeros(( self.tcount, n ))
        self.rho_arr = zeros(( self.tcount, n ))
        self.x_arr   = zeros( n )

        
#       Copy initial values into 1st row for each array        
        for i in range( 0, n ):
            self.u_arr[0][i] = self.u[i]
            self.p_arr[0][i] = self.p[i]
            self.rho_arr[0][i] = self.rho[i]
            self.x_arr[i] = self.x[i]
        
    
    def initialcond ( self, p, t ):
        self.pi = p
        self.ti = t
        
#   1D interpolation based on for x3 using: [x1, x2], [y1, y2]
#   y3 = f(x3,[x1, x2], [y1, y2])        
    def interp1d ( self, x1, x2, x3, y1, y2 ):
        
        mu = ( y1 - y2 )/( x1 - x2 )
        bu = y2 - ( mu * x2 )        
        y3 = float ( mu*x3 + bu )
        
        return y3        


#   Falculate speed of sound        
#   a = a(p)        
    def calc_speedofsound ( self, p, rho ):
        
        gamma_gas = 1.20 
        a = sqrt( gamma_gas*p /rho )
        
        echo = 0
        if ( echo == 1 ):
            
            print(" RHO   = {:8.3f}".format(rho))
            print(" γ     = {:8.3f}".format(gamma_gas))
            print(" P     = {:8.3f}".format(p))
            print("")
            print(" A     = {:8.3f}".format(a))
            print("")
            print("")
            
            
        return a  
    
#   Calculate Q = rho*a 
#   Equation 13.42 in Gas Dynamics Volume 1; Zucrow, Hoffman          
    def calc_q ( self, p, rho ):
        
        a = self.calc_speedofsound( p, rho )       
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
    def calc_t ( self, p, u, rho, opt ):

        q = self.calc_q ( p, rho )
        
#       For C+
        if ( opt == 1 ): 
            t = p + q*u
            return t   
                        
#       For C-            
        elif ( opt == -1 ):
            t = p - q*u
            return t   
            
            
#       For C0
        elif ( opt == 0 ):
            a = self.calc_speedofsound( p, rho )
            a2 = a**2
            t = p - a2*rho
            
            return [ t, a2 ]

        
        echo = 0
        if ( echo == 1 ):
            print(" *** INPUTS *** ")
            print("P = ", p)
            print("U = ", u)
            print("OPT = ", opt)
            print(" ")
            print("Q = ", q)
            print("T = ", t)
            if ( opt == 0 ):
                print("A**2 = ", a2)
            print(" ")

    def calc_flow ( self ):

#	   Start time at ( 0 + dt )
    	t = 0.0 + self.dt 
	

#	   Start time counter 
    	tcounter = 0
	
    	while ( t < self.tmax ):
            t = t + self.dt
            n = len(self.u)
            
            for i in range( 0, 1 ):
                
                if ( i == 0 ):
                    
                    self.closed_bc( i, tcounter )                    
                    
                elif ( i > 0 and i < (n-1) ):
                    print (" INTERIOR ")
                    
                elif ( i == (n-1) ):
                    print (" RIGHT BOUNDARY") 
                
            tcounter += 1
            t = self.tmax + 1
		


    def closed_bc ( self, i, j ):
        
#       Determining values        
        u0 = self.u_arr[j][i]
        p0 = self.p_arr[j][i]
        rho0 = self.rho_arr[j][i]
        x0 = self.x_arr[i]
        
        u1 = self.u_arr[j][i+1]
        p1 = self.p_arr[j][i+1]
        rho1 = self.rho_arr[j][i+1]
        x1 = self.x_arr[i+1]
        
        p = p1
        rho = rho1

        diff = 1.0
        tol  = 1e-3
        k = 0
        xold = 1
        
        while ( diff > tol ):
            
#           Calcule speed of sound
            a = self.calc_speedofsound( p, rho )
        
#           Calculate characterisitic 
            lmda = 1.0/(u1 - a)
            lmda = lmda * 1000
        
#           Find point 2
            x2 = x0 - self.dt/(lmda)
        
#           Interpolate values for point 2
            u2 = self.interp1d( x0, x1, x2, u0, u1 )
            p2 = self.interp1d( x0, x1, x2, p0, p1 )
            rho2 = self.interp1d( x0, x1, x2, rho0, rho1 ) 
            
#           Update p and rho per interpolated values for while-loop            
            p = p2
            rho = rho2
            k = k + 1

            if ( k >= 1 ):
                xold = x2

#           Check for tolerance 
            if ( k >= 1 ):
                xold = x2 
                diff = abs( x2 - xold )            
                
                
#       Calculate Q, A, T for C-
        q  = rho2 * a        
        t  = p2 - q*u2
        
#       Calculate Q, A, T for C0
        a0  = self.calc_speedofsound(p0, rho0)        
        a02 = a0**2    
        q0  = rho0 * a0
        t0  = p0 - a02*rho0
                       
        a = array( [[1.0, x0],[1.0, -a02]] )
        b = array( [[t], [t0]] )
                
        [pnew, rhonew] = solve( a, b )
        
        echo = 0
        if ( echo == 1 ):
            print ( "INPUtS: ")
            print ( "U0 = ", u0 )
            print ( "P0 = ", p0)
            print ( "X0 = ", x0 )
            print ( "RHO0 = ", rho0 )        
            print ( " " )
            print ( "U1 = ", u1 )
            print ( "P1 = ", p1)
            print ( "X1 = ", x1 )
            print ( "RHO1 = ", rho1 )
            print ( "CALCULATED ")
            print ( "")
            print ( "A = ", a)
            print ( "λ = ", lmda)
            print ( "x2 = ", x2)
            print ( "u2 = ", u2)
            print ( "p2 = ", p2)
            print ( "rho2 = ", rho2)
            
        return ([pnew, rhonew ])
            
    def interior ( self, i, j ):
        


			
"""
#   Main function to step through iterations 
    def calc_flow ( self ):
       
        if ( len(self.x) != len(self.p) or len(self.x) != len(self.u) ):
            print("ERROR")
            exit 
        
#       Start time at 0 + dt
        t = 0.0 + self.dt 

#       Start tcounter at 1 so that it copies into time+1               
        tcounter = 1
        
        while ( t < self.tmax ):
            t = t + self.dt  
            
#           
#           Limit range from 0 to (n-1) because 
#           point 1 is: i
#           point 2 is (i + 2)           
#           Don't want to hit boundary
#
#           for i in range( 0, self.n-2 )

            n = len(self.x_arr)
                           
            for i in range( 0, 1 ):
               
#               Locate Point 1
#               For the initial guess, use u[i], p[i], rho[i]                                      
#               Use iterator function to calculate new values for C+                
                [x1, u_p1, p_p1, rho_p1] = self.cp_cm_c0_iter ( i,  tcounter, 1 )
                
                
#               Locate Point 2
#               For the initial guess, use u[i+2], p[i+2], rho[i+2]                        
#               Use iterator function to calculate new values for C-
                [x2, u_p2, p_p2, rho_p2] = self.cp_cm_c0_iter ( i, tcounter, -1 )
                

#               Locate Point 3
#               Use the initial guess, use u[i+1], p[i+1], rho[i+1]  
#               Use iterator function to calculate new values for C0                
                [x3, u_p3, p_p3, rho_p3] = self.cp_cm_c0_iter ( i, tcounter,  0 )                 

                               
#               Calculate Q+/Q- and T+/T- for Points 1, 2
                q_pl = self.calc_q ( p_p1, rho_p1 )
                t_pl = self.calc_t ( p_p1, u_p1, rho_p1, 1 )
                
                q_mn = self.calc_q ( p_p2, rho_p2 )
                t_mn = self.calc_t ( p_p2, u_p2, rho_p2, -1 )
                
                
#               Determine T0 for point 3
                [t0, a0] = self.calc_t ( p_p3, u_p3, rho_p3, 0 )
                
# ---------------------------------------------------------------------             
#               Determine values at next time interval  
#               
#               p4 + Q+ * u4   = T+
#               p4 - Q- * u4   = T-
#               p4 - A0 * rho4 = T0         
#                     
#               a = [ ( 1, q+, 0), ( 1, -(q-), 0), (1, 0, -a0) ]
#               b = [ (T+), (T-), (T0) ]
# --------------------------------------------------------------------- 

                a = array( [[ 1.0,  q_pl,  0.0], 
                            [ 1.0, -q_mn,  0.0], 
                            [ 1.0,  0.0,  -a0 ]] )
    
                b = array( [[t_pl], [t_mn], [t0]] )
                
                [p_fin, u_fin, rho_fin] = solve( a, b )


#               Update arrays          
#               Update to (i+1)       
                self.p_arr[tcounter][i]   = p_fin
                self.u_arr[tcounter][i]   = u_fin
                self.rho_arr[tcounter][i] = rho_fin
                
                print ( "TEST" )
                
#           Increase tcount                
            tcounter += 1

#   Function to iterate C+, C-, C0 characterisitic
#
#   Determine is C+ or C- from sign variable
#   C+ = 1
#   C- = -1
#   C0 = 0
    def cp_cm_c0_iter ( self, node, dtstep, sign ):
        
#       Perform iterations until under tolerance        
        diff = 1.0
        tol  = 1.0E-06
        
#       Iteration counter
        kcount = 0
#        x1 = 0.0        
        
#       Determine scalar values from arrays 
        if ( sign == 1 ): 
            x0 = self.x_arr[node]
            x1 = self.x_arr[node+1]           
            
            u0 = self.u_arr[dtstep-1][node]
            u1 = self.u_arr[dtstep-1][node+1]
            
            p0 = self.p_arr[dtstep-1][node]
            p1 = self.p_arr[dtstep-1][node+1]
            
            rho0 = self.rho_arr[dtstep-1][node]
            rho1 = self.rho_arr[dtstep-1][node+1]
            
        elif ( sign == -1 ):
            
            node = ( node + 2 )
            
            x0 = self.x_arr[node]
            x1 = self.x_arr[node-1]           
            
            u0 = self.u_arr[dtstep-1][node]
            u1 = self.u_arr[dtstep-1][node-1]
            
            p0 = self.p_arr[dtstep-1][node]
            p1 = self.p_arr[dtstep-1][node-1]
            
            rho0 = self.rho_arr[dtstep-1][node]
            rho1 = self.rho_arr[dtstep-1][node-1]            

        
        elif ( sign == 0 ):
            
            node = ( node + 1 )
            
            x0 = self.x_arr[node-1]
            x1 = self.x_arr[node]
            
            u0 = self.u_arr[dtstep-1][node-1]
            u1 = self.u_arr[dtstep-1][node]
            
            p0 = self.p_arr[dtstep-1][node-1]
            p1 = self.p_arr[dtstep-1][node]
            
            rho0 = self.rho_arr[dtstep-1][node-1]
            rho1 = self.rho_arr[dtstep-1][node]                    
            
        
        xnew = 0
        while ( diff > tol ):
            
            kcount += 1                
            old_x1 = xnew
            
#           Calculate speed of sound
            a = self.calc_speedofsound( p0, rho0 )

#           Calculate lambda based on C+,  C- or C0
#           and convert from msec/m to sec/m
            if ( sign == 1 ):
                lmda = 1.0/( u0 + a )
                
            elif ( sign == -1 ):
                lmda = 1.0/( u0 - a )
                
            elif ( sign == 0 ):
                lmda = 1.0 / u0

            lmda *= 1000
            
#           Calculate new x1
            xnew = x1 - (self.dt/lmda) 
            
#           Determine new u, p, rho based on x1
            u_new   = self.interp1d ( x0, x1, xnew, u0, u1 )            
            p_new   = self.interp1d ( x0, x1, xnew, p0, p1 )           
            rho_new = self.interp1d ( x0, x1, xnew, rho0, rho1 )   
                
#               For re-calculation of speed of sound
            rho0 = rho_new
            p0   = p_new
                            
#           Determine if the new calculated x-value converged
            diff = abs( old_x1 - xnew )
            
            if ( kcount > 100 ):
                print(" CAN NOT CONVERGENCE. " )
                break
            
        echo = 0
        if ( echo == 1 ):              
            
            print(" *** INPUTS *** ")
            print(" X              = {:6.5f} M       ".format(x0))
            print(" U              = {:8.4f} M/SEC   ".format(u0))
            print(" P              = {:8.4f} N/M**2  ".format(p0))
            print(" RHO            = {:7.3f} KG/M**3 ".format(rho0))         

            if ( sign == 1 ):
                print(" CHARACTERISTIC: C+ ")
            elif ( sign == -1 ):
                print(" CHARACTERISTIC: C- ")               
            elif ( sign == 0 ):
                print(" CHARACTERISTIC: C0 ")
                

            print(" ")                      
            print( " ITERATION      = {:2} ".format(kcount))
            print( " SPEED OF SOUND = {:6.3f} M/SEC   ".format(a))
            print( " λ              = {:6.4f} MSEC/M  ".format(lmda))
            print( " XNEW           = {:6.5f} M       ".format(xnew))
            print( " UNEW           = {:6.3f} M/SEC   ".format(u_new))
            print( " PNEW           = {:6.3f} N/M**2  ".format(p_new))
            print( " RHO NEW        = {:7.3f} KG/M**3 ".format(rho_new))                
            print( " ")
            print( " ")    

        return [ xnew, u_new, p_new, rho_new ]



    def output( self ):
        
        f = open( "outputUnsteady.out", "w")
        
        f.write(" *** INITIAL VALUES *** \n")       
        for i in range(0, len(self.x_arr)):
            f.write(" X({:2} ) = {:5.4f} M\n".format(i, self.x_arr[i]))
        f.write("\n")            
        for i in range(0, len(self.x_arr)):    
            f.write(" P({:2} ) = {:8.4f} N/M**2 \n".format(i, self.p_arr[0][i]))   
        f.write("\n")            
        for i in range(0, len(self.x_arr)):
            f.write(" U({:2} ) = {:8.4f} M/SEC\n".format(i, self.u_arr[0][i]))
        f.write("\n")            
        for i in range(0, len(self.x_arr)):
            f.write(" RHO({:2} ) = {:8.4f} KG/M**3\n".format(i, self.rho_arr[0][i]))            
            
        f.write(" \n\n\n ")
        f.write(" *** TRANSIENT SOLUTION *** \n")
        
        n = len( self.x_arr )
        for i in range( 1, self.tcount ):        
            f.write("  T    = {:8.6f} SEC\n ".format(self.dt * i))
            
            for j in range ( 0, n ):
                f.write(" P( {:3} )    = {:8.4f} N/M**2 \n ".format(j, self.p_arr[i][j]))
            
            
                
            #f.write(" P{:3}    = {:8.4f} N/M**2 \n ".format(j, self.p_arr[i][j]))
            
            
            
            for j in range ( 0, len(self.x_arr) ):
                f.write(" P{:3}    = {:8.4f} N/M**2 \n ".format(j, self.p_arr[j][i]))
                f.write(" U{:3}    = {:8.4f} M/SEC ".format(j, self.u_arr[j][i]))        
                f.write(" RHO{:3}    = {:8.4f} M/SEC ".format(j, self.rho_arr[j][i])) 

        
"""        
        
        

      
def main( ):
    
    x   = array( [ 0.0000, 0.02450, 0.1016, 0.1270, 0.1524 ] )
    u   = array( [ 0.0000, 47.1600, 188.28, 235.30, 282.53 ] )
    p   = array( [ 113.49, 113.42, 111.77, 110.78, 109.54 ] )
    p   = p*(10**5)
    rho = array( [ 14.251, 14.239, 14.064, 13.959, 13.827 ] )
    
    pi = 690
    ti = 3330
    
    
    test = unsteady( )
    test.inputs( u, p, x, rho )
    test.initialcond( pi, ti )
    test.calc_flow( )
    #test.output()
    
    
    
#    test.output( )
    

if __name__ == '__main__':
    main( ) 


    
    
    
        


            
    