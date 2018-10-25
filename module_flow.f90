module mod_flow

    use mod_dp
    use mod_solvers

    implicit none

!   Module declarations
!   Allocatable arrays
    real ( dp ), allocatable :: u(:)
    real ( dp ), allocatable :: p(:)
    real ( dp ), allocatable :: x(:)
    real ( dp ), allocatable :: dens(:)
    real ( dp ), allocatable :: a(:)
    real ( dp ), allocatable :: q(:)

    integer :: n


!   Initial conditions
!
!   Units: k
    real ( dp ) :: ti = 3330.0_dp

!   Units: N/m**2 * (10E+05)
    real ( dp ) :: pinit = 690E+5/1.0E+05



!   Physical constants
!
!   Units: N/A
    real ( dp ), parameter :: gammagas = 1.20_dp

!   Units: J/(kg-k)
    real ( dp ), parameter :: rgas = 320.0_dp




    contains

    subroutine unsteady ( )
! -------------------------------------------------------------------------- 80
! Name: unsteady
!
! Abstract: Subroutine to handle "high level" unsteady flow information
!
! Inputs:
!   N/A
!
! Outputs:
!   N/A
!
!
! Description:
!
!	All lower level subroutines accessed through unsteady subroutine
!
! -------------------------------------------------------------------------- 80
        implicit none

        n = 3
        call allocatearray1d ( )

!       Units: m
        x = (/ 0.1000_dp, 0.1250_dp, 0.1500_dp /)

!       Units: m/s
        u = (/ 188.30_dp, 235.30_dp, 282.50_dp /)

!       Units: N*m**2 * (10E+05)
        p = (/ 111.77_dp, 110.78_dp, 109.54_dp /)

        call interiorpoint ( 1 )

        call writedata ( )

    end subroutine unsteady



    subroutine allocatearray1d ( )
! -------------------------------------------------------------------------- 80
! Name: allocatearray1d
!
! Abstract: Allocates array of n-elements to zero
!
! Inputs:
!   N/A
!
! Outputs:
!   N/A
!
!
! Description:
!
!	Uses intrinsic Fortran allocation command to allocate all allocatable 1-D
!   arrays within mod_data to zero.
!
! -------------------------------------------------------------------------- 80
        implicit none

!       Local declaration
        integer :: echo

!       Allocate and initialize arrays to zero
        allocate ( u(n) )
        allocate ( p(n) )
        allocate ( x(n) )
        allocate ( dens(n) )
        allocate ( a(n) )
        allocate ( q(n) )
        u(1:n) = 0.0_dp
        p(1:n) = 0.0_dp
        x(1:n) = 0.0_dp
        dens(1:n) = 0.0_dp
        a(1:n) = 0.0_dp
        q(1:n) = 0.0_dp

        echo = 0

        if ( echo == 1 ) then

            write(*,*) " U = ", u
            write(*,*) " "
            write(*,*) " P = ", p
            write(*,*) " "
            write(*,*) " X = ", x
            write(*,*) " "
            write(*,*) " DENS = ", dens
            write(*,*) " "

        end if

    end subroutine allocatearray1d


    subroutine interiorpoint ( i )
! -------------------------------------------------------------------------- 80
! Name: interiorpoint
!
! Abstract: Perform calculations for points that are not on the boundary
!
! Inputs:
!   i, integer, represents i-th node
!
! Outputs:
!   N/A
!
!
! Description:
!
!	This subroutine performs the required calculations to determine the
!   propagation of the wave on interior (non-boundary) nodes.
!
! -------------------------------------------------------------------------- 80

        implicit none

!       Argument declaration
        integer, intent ( in ) :: i

!       Local declaration
        real ( dp ) :: lambda_plus
        real ( dp ) :: lambda_minus
        real ( dp ) :: x1
        real ( dp ) :: u1
        real ( dp ) :: p1

!       Time step
!       Units: milli-second
        real ( dp ), parameter :: dt = 0.020_dp


!       Calculate speed of sound at node-i
        call calcspeedofsound ( p(i), a(i) )

!       Calculate lambda plus using the value of u(i) and a(i)
!       Units: seconds/meter
        lambda_plus = 1.0_dp / ( u(i) + a(i) )

!       Convert from seconds/meter to milli-seconds/meter
        lambda_plus = lambda_plus * (1E+03)

!       Find x1 using equation 13.59 from Gas Dynamics Vol. 1, Zucrow & Hoffman
        x1 = x(i+1) - dt/lambda_plus

!       Find new values of u and p at x1
        call lin1dinterp ( x(i), x(i+1), u(i), u(i+1), x1, u1)
        call lin1dinterp ( x(i), x(i+1), p(i), p(i+1), x1, p1)

!       calculate the density using p1
        call calcdensity ( p1, dens(i) )

!       Calculate the Acoustic impedance
        call calcq ( i )






    end subroutine interiorpoint


    subroutine calcspeedofsound ( plocal, alocal )
! -------------------------------------------------------------------------- 80
! Name: calcspeedofsound
!
! Abstract: Calculate speed of sound
!
! Inputs:
!   p, real, the pressure of the fluid
!            units: N/m**2
!
! Outputs:
!   a, real, the speed of sound
!            units: m/s**2
!
! Description:
!
!	This subroutine calculates the speed of sound of the fluid for a given
!   pressure using equation 13.30 from Gas Dynamics Volume 1, Zucrow & Hoffman
!
! -------------------------------------------------------------------------- 80
        implicit none

!       Argument declaration
        real ( dp ), intent ( in ) :: plocal
        real ( dp ), intent ( inout ) :: alocal

!       Local declaration
        integer :: echo

!       Formula 13.30 (a) from Gas Dynamics Volume 1, Zucrow & Hoffman
        alocal = sqrt( gammagas * rgas * ti ) * ( plocal / pinit )**(( gammagas - 1.0_dp )/( 2.0_dp * gammagas ))

        echo = 0
        if ( echo == 1 ) then

            write(*,*) " PLOCAL       = ", plocal
            write(*,*) " GAMMA GAS    = ", gammagas
            write(*,*) " R GAS        = ", rgas
            write(*,*) " INITIAL T    = ", ti
            write(*,*) " INITIAL P    = ", pinit
            write(*,*) " "
            write(*,*) " ALOCAL = ", alocal

        end if





    end subroutine calcspeedofsound


    subroutine calcdensity ( plocal, denslocal)
! -------------------------------------------------------------------------- 80
! Name: calcdensity
!
! Abstract: Calculate the fluid density
!
! Inputs:
!   p, real, the pressure of the fluid
!            units: N/m**2
!
! Outputs:
!   a, real, the fluid density
!            units: kg/m**3
!
! Description:
!
!	This subroutine calculates the fluid density using the local pressure
!   using equation 13.29 from Gas Dynamics Volume 1, Zucrow & Hoffman
!
! -------------------------------------------------------------------------- 80
        implicit none

!       Argument declaration
        real ( dp ), intent ( in ) :: plocal
        real ( dp ), intent ( inout ) :: denslocal

!       Local declaration
        integer :: echo

!       Formula 13.29 (a) from Gas Dynamics Volume 1, Zucrow & Hoffman
        denslocal = pinit*(1E+5)/( rgas * ti )*( plocal / pinit )**( 1.0_dp / gammagas )

        echo = 0
        if ( echo == 1 ) then

            write(*,*) " PLOCAL       = ", plocal
            write(*,*) " GAMMA GAS    = ", gammagas
            write(*,*) " R GAS        = ", rgas
            write(*,*) " INITIAL T    = ", ti
            write(*,*) " INITIAL P    = ", pinit
            write(*,*) " "
            write(*,*) " DENS. LOCAL = ", denslocal

        end if

    end subroutine calcdensity


    subroutine calcq ( i )
! -------------------------------------------------------------------------- 80
! Name: calcq
!
! Abstract: Calculate the acoustic impedance
!
! Inputs:
!   i, integer
!
! Outputs:
!   N/A
!
!
! Description:
!
!	This subroutine calculates the acoustic impedance
!   using equation 13.42 from Gas Dynamics Volume 1, Zucrow & Hoffman
!
! -------------------------------------------------------------------------- 80
        implicit none

!       Argument declaration
        integer, intent ( in ) :: i

!       Local declaration
        integer :: echo

!       Formula 13.42 (a) from Gas Dynamics Volume 1, Zucrow & Hoffman
        q(i) = dens(i) * a(i)

        echo = 1
        if ( echo == 1 ) then

            write(*,*) " DENSITY        = ", dens(i)
            write(*,*) " SPEED OF SOUND = ", a(i)
            write(*,*) " "
            write(*,*) " Q = ", q(i)

        end if

    end subroutine calcq





    subroutine writedata ( )
! -------------------------------------------------------------------------- 80
! Name: writedata
!
! Abstract: Write ASCII data to file
!
! Inputs:
!   N/A
!
! Outputs:
!   N/A
!
!
! Description:
!
!	Write ASCII data to local machine.
!
! -------------------------------------------------------------------------- 80
        implicit none

!       Local declaration
        integer               :: i
        character ( len = 9 ) :: format1


!       Open file streams
        open ( unit = 20, file = "output_position.out",    status = "replace" )
        open ( unit = 21, file = "output_velocity.out",    status = "replace" )
        open ( unit = 22, file = "output_pressure.out",    status = "replace" )
        open ( unit = 23, file = "output_spdsound.out",    status = "replace" )

        format1 = "(g14.6)"

        do i = 1, n

            write ( 20, format1 ) x(i)
            write ( 21, format1 ) u(i)
            write ( 22, format1 ) p(i)
            write ( 23, format1 ) a(i)

        end do


!       Close streams
        close ( 20 )
        close ( 21 )
        close ( 22 )
        close ( 23 )


    end subroutine writedata




end module mod_flow

