module mod_flow

    use mod_dp
    use mod_solvers

    implicit none

    real ( dp ), allocatable :: u(:)
    real ( dp ), allocatable :: p(:)
    real ( dp ), allocatable :: x(:)
    real ( dp ), allocatable :: dens(:)

    real ( dp ) :: x4
    real ( dp ) :: u4
    real ( dp ) :: p4

    integer     :: n

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
        call allocatearray1d ( n )

!       Units: m
        x = (/ 0.1000_dp, 0.1250_dp, 0.1500_dp /)

!       Units: m/s
        u = (/ 188.30_dp, 235.30_dp, 282.50_dp /)

!       Units: N*m**2 * (10**5)
        p = (/ 111.77_dp, 110.78_dp, 109.54_dp /)


        x4 = 0.10180_dp

        call lin1dinterp ( x(1), x(2), u(1), u(2), x4, u4 )
        call lin1dinterp ( x(1), x(2), p(1), p(2), x4, p4 )

        write(*,*) " X4 = ", x4
        write(*,*) " U4 = ", u4
        write(*,*) " P4 = ", p4





        call writedata ( )

    end subroutine unsteady



    subroutine allocatearray1d ( n )
! -------------------------------------------------------------------------- 80
! Name: allocatearray1d
!
! Abstract: Allocates array of n-elements to zero
!
! Inputs:
!   n, integer,
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

!       Argument declaration
        integer, intent( in ) :: n

!       Local declaration
        integer :: echo

!       Allocate and initialize arrays to zero
        allocate ( u(n) )
        allocate ( p(n) )
        allocate ( x(n) )
        allocate ( dens(n) )
        u(1:n) = 0.0_dp
        p(1:n) = 0.0_dp
        x(1:n) = 0.0_dp
        dens(1:n) = 0.0_dp

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
!	Write ASCII data to local machine. !
!
! -------------------------------------------------------------------------- 80
        implicit none

        integer               :: i
        character ( len = 9 ) :: format1


!       Open file streams
        open ( unit = 20, file = "output_position.out",    status = "replace" )
        open ( unit = 21, file = "output_velocity.out",    status = "replace" )
        open ( unit = 22, file = "output_pressure.out",    status = "replace" )

        format1 = "(g14.6)"

        do i = 1, n

            write ( 20, format1 ) x(i)
            write ( 21, format1 ) u(i)
            write ( 22, format1 ) p(i)

        end do


!       Close streams
        close ( 20 )
        close ( 21 )
        close ( 22 )


    end subroutine writedata



end module mod_flow

