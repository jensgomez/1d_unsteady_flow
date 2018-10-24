module mod_solvers

    use mod_dp

    implicit none

    contains

    subroutine lin1dinterp( x1, x2, y1, y2, x3, y3 )
! -------------------------------------------------------------------------- 80
! Name: lin1dinterp
!
! Abstract: Subroutine performs linear, 1D interpolation for given data
!
! Inputs:
!   x1,    real, x1 coordinate
!   x2,    real, x2 coordinate
!   y1,    real, y1 coordinate
!   y2,    real, y2 coordinate
!   x3,    real, x3 coordinate
!
! Outputs:
!   y3,    real, interpolated result
!
!
! Description:
!
!	Performs linear, 1D interpolation for given data using following formula:
!
!   y3 = mu*x3 + bu
!
!	mu = (y1-y2)/(x1-x2)
!
!	bu = y2 - (mu*x2)
!
! -------------------------------------------------------------------------- 80
        implicit none

        real ( dp ), intent ( in ) :: x1
        real ( dp ), intent ( in ) :: x2
        real ( dp ), intent ( in ) :: x3
        real ( dp ), intent ( in ) :: y1
        real ( dp ), intent ( in ) :: y2
        real ( dp ), intent ( inout ) :: y3

        real ( dp ) :: mu
        real ( dp ) :: bu

        integer :: echo

        mu = ( y1 - y2 )/( x1 - x2 )
        bu = y2 - ( mu * x2 )

        y3 = mu*x3 + bu

        echo = 0
        if ( echo == 1 ) then

            write(*,*) " INPUTS: "
            write(*,*) " X1 = ", x1
            write(*,*) " X2 = ", x2
            write(*,*) " Y1 = ", y1
            write(*,*) " Y2 = ", y2
            write(*,*) " "
            write(*,*) " X3 = ", x3
            write(*,*) " "
            write(*,*) " "
            write(*,*) " OUTPUTS: "
            write(*,*) " Y3 = ", y3
            write(*,*) " "

        end if

    end subroutine lin1dinterp







end module mod_solvers
