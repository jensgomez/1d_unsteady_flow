module mod_flow

    use mod_dp
    use mod_solvers

! -------------------------------------------------------------------------- 80
!   Units:
!   x   = m
!   u   = m/s
!   p   = Pa
!   rho = kg/m**3
! -------------------------------------------------------------------------- 80

    implicit none

    integer :: n

	real ( kind = dp ), parameter :: dt = 1.00E-5
	real ( kind = dp ), parameter :: tmax = 0.500_dp	
	integer                       :: tstep = int ( tmax / dt )
	
	
    real ( kind = dp ), dimension (:), allocatable :: xinit
    real ( kind = dp ), dimension (:), allocatable :: pinit
    real ( kind = dp ), dimension (:), allocatable :: uinit
    real ( kind = dp ), dimension (:), allocatable :: rhoinit
	
	real ( kind = dp ), dimension (:), allocatable :: t
	
	

    contains

        subroutine unsteady ( )

            implicit none
			
			integer :: i 
			
			
			
            call initialize ( )
			
			
			do i = 1, tstep 
			
				
			
			
			
			end do 
			
			
			
					
			
			
			
			
			
			
			
			
			
			call outputfile ( ) 
			
			



        end subroutine unsteady

        subroutine initialize ( )
! -------------------------------------------------------------------------- 80
!  Subroutine to initialize data
!
!  Inputs:  N/A
!
!  Outputs: N/A
!
! -------------------------------------------------------------------------- 80
            implicit none

!           Local parameters
			integer :: i 
            integer :: echo
            integer :: readinfo

!
!           Readinfo variable to read from source code or
!           from file
!           readinfo = 0 (from source code)
!           readinfo = 0 (from input file)
!
            readinfo = 0
            if ( readinfo == 0 ) then

!               Initialize arrays
                allocate ( xinit(5) )
                allocate ( pinit(5) )
                allocate ( uinit(5) )
                allocate ( rhoinit(5) )

                xinit   = (/ 0.0000_dp, 0.02450_dp, 0.1016_dp, 0.1270_dp, 0.1524_dp /)
                uinit   = (/ 0.0000_dp, 47.1600_dp, 188.28_dp, 235.30_dp, 282.53_dp /)
                pinit   = (/ 113.49_dp,  113.42_dp, 111.77_dp, 110.78_dp, 109.54_dp /)
                rhoinit = (/ 14.251_dp,  14.239_dp, 14.064_dp, 13.959_dp, 13.827_dp /)

!               Update P
                pinit   = pinit * ( 10**5 )

!               Determine number of spatial steps
                n = size(xinit)

            end if
			
			allocate ( t(tstep) )
			t( 1 : n ) = 0.0_dp 

			


            echo = 0
            if ( echo == 1 ) then

                write(*,*) ' X = ', xinit
                write(*,*) ' P = ', pinit
                write(*,*) ' U = ', uinit
                write(*,*) ' RHO = ', rhoinit


            end if



        end subroutine initialize

		
        subroutine outputfile ( )
! -------------------------------------------------------------------------- 80
!  Subroutine to create output file 
!
!  Inputs:  N/A
!
!  Outputs: N/A
!
! -------------------------------------------------------------------------- 80	

			implicit none 
			
			
			integer :: i 
			
			open ( unit = 20, file = 'data.out',         status = 'replace')
			open ( unit = 21, file = 'node1.out',        status = 'replace')
			open ( unit = 22, file = 'node2.out',        status = 'replace')
			open ( unit = 23, file = 'node3.out',        status = 'replace')
			open ( unit = 24, file = 'node4.out',        status = 'replace')
			open ( unit = 25, file = 'node5.out',        status = 'replace')
			
			


		end subroutine outputfile




end module mod_flow
