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
!   Performs linear, 1D interpolation for given data using following formula:
!
!   y3 = mu*x3 + bu
!
!   mu = (y1-y2)/(x1-x2)
!
!   bu = y2 - (mu*x2)
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
    
    
            subroutine sor ( a, b, x, n )
! -------------------------------------------------------------------------- 80
!       This subroutine is for the Successive over-relaxation method
!       The algorithm is taken from: Numerical Methods for Engineers
!       and Scientists by Joe D. Hoffman
!       See Section 1.7.4
! -------------------------------------------------------------------------- 80


            implicit none
            ! Subroutine input/output variable declarations
            integer ( kind = 4 ), intent(in) :: n 						! n is the matrix dimension
            real    ( kind = 8 ), intent(in) :: a(n,n)						! a is the input matrix
            real    ( kind = 8 ), intent(in) :: b(n)						! b is the constraints
            real    ( kind = 8 ), intent(inout) :: x(n)						! x is the solution array

            ! Local variable declarations
            integer ( kind = 4 ) :: i
            integer ( kind = 4 ) :: j
            integer ( kind = 4 ) :: k

            ! w is the relaxation factor. If A is symmetric and positive-definte,
            ! then the following is needed for convergence: 0 < w < 2
            ! if w = 1.0, the SOR method becomes Gauss-Seidel
            real    ( kind = 8 ) :: w


            real    ( kind = 8 ) :: sum1(n)
            real    ( kind = 8 ) :: tolerance
            real    ( kind = 8 ) :: eps
            real    ( kind = 8 ) :: old
            real    ( kind = 8 ) :: new


            ! Initializations
            do i=1,n
               x(i) = 1.0
               sum1(i) = 1.0
            end do

            w = 1.05
            tolerance = 1e-5
            eps = 1
            k = 1

            do while (tolerance < eps)
                old = x(3)
                do i=1,n
                    sum1(i) = 0.0
                    do j=1,n
                        if (i /= j) then
                            sum1(i) = sum1(i) + a(i,j)*x(j)
                        end if
                    end do
                    x(i) = x(i) + w*((b(i) - sum1(i))/a(i,i) - x(i))
        !			x(i) = (1-w)*x(i) + (w/a(i,i))*(b(i) - sum1(i))
                end do
                new = x(3)
                eps = abs(new - old)
                k = k + 1
            end do

            10 format (1x, "Number of iterations til convergence: ", i5)
            20 format (1x, f15.8)

            write(*,10) k
            print*, ''
            do i=1,n
                write(*,20) x(i)
            end do

            return
        end Subroutine sor







end module mod_solvers
