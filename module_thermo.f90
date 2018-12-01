module mod_thermo

    implicit none

    contains


    subroutine p_sat ( t, psat )

! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the saturation pressure as a function
!   of temperature (valid for 89.965 deg C < T < 373.253 deg C).
!   Correlations taken from Appendix I
!
!   t,    input,  real ( kind = 8 ), the temperature
!                                    units: deg F
!
!   psat, output, real ( kind = 8 ), the saturation pressure for the given temp.
!                                    units: psia
!

        implicit none

        real ( kind = 8 ), intent ( in )    :: t
        real ( kind = 8 ), intent ( inout ) :: psat


!
!       Variables used in this subroutine
        real ( kind = 8 ) :: tc
        real ( kind = 8 ) :: psat_mpa


!       Convert temperature from deg F to deg C
        tc = (5.0D+0/9.0D+0)*(t - 32.0D+0)

        if ( tc < 89.965D+0 .or. tc > 373.253) then

            write ( *, * ) ' THE PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '


        else if ( tc >= 89.965D+0 .and. tc < 139.781D+0 ) then


            psat_mpa = ((tc+57.0D+0)/236.2315D+0)**5.602972D+0


        else if ( tc >= 139.781D+0 .and. tc < 203.662D+0 ) then

            psat_mpa = ((tc+28.0D+0)/207.9248D+0)**4.778504D+0


        else if ( tc >= 203.662D+0 .and. tc < 299.407D+0 ) then

            psat_mpa = ((tc+5.0D+0)/185.0779D+0)**4.304376D+0


        else if ( tc >= 299.407D+0 .and. tc < 355.636D+0 ) then

            psat_mpa = ((tc+5.0D+0)/185.0779D+0)**4.304376D+0


        else if ( tc >= 299.407D+0 .and. tc < 355.636D+0 ) then

            psat_mpa = ((tc+16.0D+0)/195.1819D+0)**4.460843D+0


        else if ( tc >= 355.563D+0 .and. tc <= 373.253D+0 ) then

            psat_mpa = ((tc+50.0D+0)/227.2963D+0)**4.960785D+0


        end if


!       Convert from MPa to psia
        psat = psat_mpa * 145.038D+0


    end subroutine p_sat


    subroutine t_sat ( p, tsat )

! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the saturation temperature as a function
!   of pressure (valid for 0.070 MPa < P < 21.85 MPa).
!   Correlations taken from Appendix II
!
!   p,    input,  real ( kind = 8 ), the pressure
!                                    units: psia
!
!   tsat, output, real ( kind = 8 ), the saturation temp for the given pres.
!                                    units: deg F
!
!


        implicit none

        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( inout ) :: tsat

!       Local variables
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: tsat_si

!       Convert pressures from psia to MPa
        p_si = p * 0.00689476D+0


!       Correlations from Appendix II
        if ( p_si < 0.070D+0 .and. p_si > 21.85D+0 ) then

            write ( *, * ) ' THE PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '


        else if ( p_si >= 0.070D+0 .and. p_si < 0.359D+0 ) then

            tsat_si = 236.2315D+0 * p_si**(0.1784767D+0) - 57.0D+0


        else if ( p_si >= 0.359D+0 .and. p_si <= 1.676D+0 ) then

            tsat_si = 207.9248D+0 * p_si**(0.2092705D+0) - 28.0D+0


        else if ( p_si > 1.676D+0 .and. p_si <= 8.511D+0 ) then

            tsat_si = 185.0779D+0 * p_si**(0.2323217D+0) - 5.0D+0


        else if ( p_si > 8.511D+0 .and. p_si < 17.69D+0 ) then

            tsat_si = 195.1819D+0 * p_si**(0.2241729D+0) - 16.0D+0


        else if ( p_si >= 17.69D+0 .and. p_si <= 21.85D+0 ) then

            tsat_si = 227.2963D+0 * p_si**(0.201581D+0) - 50.0D+0

        end if


!       Convert the saturation temperature into deg F
        tsat = tsat_si*(9.0D+0/5.0D+0) + 32.0D+0

    end subroutine t_sat


    subroutine liq_dens ( p, t, rho )

! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the liquid density as a function
!   of temperature and pressure
!
!   This subroutine is valid for: 91.79 deg C < T < 373.253 deg C &
!                                  saturation < P < 21.5 MPa
!
!
!      p,  input, real ( kind = 8 ), the pressure
!                                    units: psia
!
!
!      t,  input, real ( kind = 8 ), the temperature
!                                    units: deg F
!
!    rho, output, real ( kind = 8 ), the calculated density
!                                    units: lbm/ft**3
!
!

        implicit none

        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: rho


!
!       Variables used locally in subroutine
        real ( kind = 8 ) :: psat
        real ( kind = 8 ) :: psat_si
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: t_si
        real ( kind = 8 ) :: rho_si
        real ( kind = 8 ) :: rho_t


!       Calculate the saturation pressure
        call p_sat ( t, psat )


!       Convert pressures from psia to MPa
        psat_si = psat * 0.00689476D+0
        p_si = p * 0.00689476D+0

!       Convert temperature from deg F to deg C
        t_si = (5.0D+0/9.0D+0)*(t - 32.0D+0)


!       Correlations for rho_t
        if ( psat_si < 0.075D+0 .or. psat_si > 21.500D+0 ) then

            write ( *, * ) ' THE SATURATION PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '


        else if ( psat_si >= 0.075D+0 .and. psat_si <= 1.00D+0 ) then

            rho_t = ( 1.2746977E-04 * psat_si**(0.4644339D+0) + 0.001D+0 )**-1.0D+0


        else if ( psat_si > 1.00D+0 .and. psat_si <= 3.88D+0 ) then

            rho_t = ( 1.0476071E-04 * psat_si**(0.5651090D+0) + 0.001022D+0 )**-1.0D+0


        else if ( psat_si > 3.880D+0 .and. psat_si <= 8.840D+0 ) then

            rho_t = ( 3.2836717E-05 * psat_si + 1.12174735E-03 )**-1.0D+0


        else if ( psat_si > 8.840D+0 .and. psat_si <= 14.463D+0 ) then

            rho_t = ( 3.3551046E-04 * exp(5.8403566E-02*psat_si) + 0.00085D+0 )**-1.0D+0


        else if ( psat_si > 14.463D+0 .and. psat_si < 18.052D+0 ) then

            rho_t = ( 3.1014626E-08 * psat_si**(3.284754D+0) + 0.001430D+0 )**-1.0D+0


        else if ( psat_si >= 18.052D+0 .and. psat_si < 20.204D+0 ) then

            rho_t = ( 1.5490787E-11 * psat_si**(5.7205D+0) + 0.001605D+0 )**-1.0D+0


        else if ( psat_si >= 20.204D+0 .and. psat_si <= 21.500D+0 ) then

            rho_t = ( 4.1035988E-24 * psat_si**(15.03329D+0) + 0.00189D+0 )**-1.0D+0


        end if


!       Density equation from Section 4.1
        rho_si = rho_t + ((170.0D+0/(375.0D+0 - t_si)) - 0.20D+0)*( p_si - psat_si )

!       Convert from kg/m**3 to lb/ft**3
        rho = rho_si * 0.062428D+0



    end subroutine liq_dens


    subroutine liq_enth ( p, t, h )

! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the liquid enthalpy as a function
!   of temperature and pressure
!
!   This subroutine is valid for: 91.79 deg C < T < 373.253 deg C &
!                                  saturation < P < 21.5 MPa
!
!
!      p,  input, real ( kind = 8 ), the pressure
!                                    units: psia
!
!
!      t,  input, real ( kind = 8 ), the temperature
!                                    units: deg F
!
!    rho, output, real ( kind = 8 ), the calculated enthalpy
!                                    units: BTU/lbm
!
!

        implicit none

        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: h


!       Variables used locally in this subroutine
        real ( kind = 8 ) :: psat
        real ( kind = 8 ) :: psat_si
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: t_si
        real ( kind = 8 ) :: hs


!       Calculate saturation pressure
        call p_sat ( t, psat )

!       Convert pressures from psia to MPa
        psat_si = psat * 0.00689476D+0
        p_si = p * 0.00689476D+0

!       Convert temperature from deg F to deg C
        t_si = (5.0D+0/9.0D+0)*(t - 32.0D+0)


!       Correlations for hs
        if ( psat_si < 0.075D+0 .and. psat_si > 21.70D+0 ) then

            write ( *, * ) ' THE SATURATION PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '

        else if ( psat_si > 0.075D+0 .and. psat_si < 0.942D+0 ) then

            hs = 912.1779D+0 * psat_si**(0.2061637D+0) - 150.0D+0


        else if ( psat_si >= 0.942D+0 .and. psat_si < 4.020D+0 ) then

            hs = 638.0621D+0 * psat_si**(0.2963197D+0) + 125.0D+0


        else if ( psat_si >= 4.020D+0 .and. psat_si < 9.964D+0 ) then

            hs = 373.7665D+0 * psat_si**(0.4235532D+0) + 415.0D+0


        else if ( psat_si >= 9.964D+0 .and. psat_si < 16.673D+0 ) then

            hs = 75.38673D+0 * psat_si**(0.8282384D+0) + 900.0D+0


        else if ( psat_si >= 16.673D+0 .and. psat_si < 20.396D+0 ) then

            hs = 0.1150827D+0 * psat_si**(2.711412D+0) + 1440.0D+0


        else if ( psat_si >= 20.396D+0 .and. psat_si <= 21.700D+0 ) then

            hs = 9.1417257E-14 * psat_si**(11.47287D+0) + 1752.0D+0

        end if


!       Enthalpy equation from Section 4.2
        h = hs + (1.4D+0 - ( 169.0D+0 / (369.0D0 - t_si)))*(p_si - psat_si)


!       Convert from kJ/kg to BTU/lbm
        h = h * (0.4299D+0)


    end subroutine liq_enth




    subroutine liq_cp ( psat, p, t, cp )


! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the liquid specific heat as a function
!   of temperature and pressure
!
!   This subroutine is valid for: 91.79 deg C < T < 373.253 deg C &
!                                  saturation < P < 21.5 MPa
!
!
!   psat,  input, real ( kind = 8 ), the saturation pressure
!                                    units: psia
!
!
!      p,  input, real ( kind = 8 ), the pressure
!                                    units: psia
!
!
!      t,  input, real ( kind = 8 ), the temperature
!                                    units: deg F
!
!    rho, output, real ( kind = 8 ), the calculated specific heat
!                                    units: BTU/(lbm-degF)
!
        implicit none

        real ( kind = 8 ), intent ( in ) :: psat
        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: cp


!       Variables used locally in this subroutine
        real ( kind = 8 ) :: psat_si
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: t_si
        real ( kind = 8 ) :: cp_si
        real ( kind = 8 ) :: cpf


!       Convert pressures from psia to MPa
        psat_si = psat * 0.00689476D+0
        p_si = p * 0.00689476D+0

!       Convert temperature from deg F to deg C
        t_si = (5.0D+0/9.0D+0)*(t - 32.0D+0)


!       Correlations for cpf
        if ( psat_si < 0.030D+0 .and. psat_si > 21.50D+0 ) then

            write ( *, * ) ' THE SATURATION PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '


        else if ( psat_si > 0.030D+0 .and. psat_si <= 0.671D+0 ) then

            cpf = 0.247763D+0 * psat_si**(0.5704026D+0) + 4.150D+0


        else if ( psat_si > 0.671D+0 .and. psat_si <= 2.606D+0 ) then

            cpf = 0.1795305D+0 * psat_si**(0.8957323D+0) + 4.223D+0


        else if ( psat_si > 2.606D+0 .and. psat_si <= 6.489D+0 ) then

            cpf = 0.09359843D+0 * psat_si**(1.239114D+0) + 4.340D+0


        else if ( psat_si > 6.489D+0 .and. psat_si <= 11.009D+0 ) then

            cpf = 0.01068888D+0 * psat_si**(2.11376D+0) + 4.740D+0


        else if ( psat_si > 11.009D+0 .and. psat_si <= 14.946D+0 ) then

            cpf = 1.333058E-04 * psat_si**(3.707294D+0) + 5.480D+0


        else if ( psat_si > 14.946D+0 .and. psat_si <= 18.079D+0 ) then

            cpf = 6.635658E-03 * ( psat_si - 10.0D+0 )**3.223323D+0 + 7.350D+0


        else if ( psat_si > 18.079D+0 .and. psat_si <= 21.50D+0 ) then

            cpf = 4.6844786E-06 * exp( 0.7396875D+0 * psat_si) + 10.020D+0


        end if



!       Calculate specific heat per equation in Section 4.4
        cp_si = cpf + ( p_si - psat_si )*( 0.0018D+0 - 76.D+0/( 364.0D+0 - t_si)**1.8D+0)



!       Convert from kJ/(kg-k) to BTU/(lbm-degF)
        cp = cp_si * 0.2388D+0


    end subroutine liq_cp



    subroutine vap_vol ( p, t, v, rho )
! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the vapor specific volume as a function
!   of temperature and pressure
!
!   This subroutine is valid for: saturation deg C < T < 450.0 deg C &
!                                  0.085 MPa < P < 21.5 MPa
!
!
!      p,  input, real ( kind = 8 ), the pressure
!                                    units: psia
!
!
!      t,  input, real ( kind = 8 ), the temperature
!                                    units: deg F
!
!
!      v, output, real ( kind = 8 ), the calculated specific volume
!                                    units: ft**3/kg
!
!    rho, output, real ( kind = 8 ), the calculated density
!                                    units: lb/ft**3
!
!


        implicit none

        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: v
        real ( kind = 8 ), intent ( inout ) :: rho


!       Variables used locally in subroutine
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: tsat
        real ( kind = 8 ) :: tsat_si
        real ( kind = 8 ) :: t_si
        real ( kind = 8 ) :: v_si
        real ( kind = 8 ) :: vg

!       Calculate saturation temperature
        call t_sat ( p, tsat )

!       Convert pressure from psia to MPa
        p_si = p * 0.00689476D+0

!       Convert temperatures from deg F to deg C
        tsat_si = (5.0D+0/9.0D+0)*(tsat - 32.0D+0)
        t_si = (5.0D+0/9.0D+0)*(t - 32.0D+0)


        if ( p_si < 0.085D+0 .and. p_si > 21.50D+0 ) then

            write ( *, * ) ' THE PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '


        else if ( p_si > 0.085D+0 .and. p_si < 1.112D+0 ) then

            vg = ( 5.126076D+0 * p_si**(0.9475862D+0) + 0.012D+0 )**-1.0D+0


        else if ( p_si >= 1.112D+0 .and. p_si < 3.932D+0 ) then

            vg = ( 4.630832D+0 * p_si**(1.038819D+0) + 0.520D+0 )**-1.0D+0


        else if ( p_si >= 3.932D+0 .and. p_si < 8.996D+0 ) then

            vg = ( 2.868721D+0 * p_si**(1.252148D+0) + 3.800D+0 )**-1.0D+0


        else if ( p_si >= 8.996D+0 .and. p_si < 14.628D+0 ) then

            vg = ( 0.5497653D+0 * p_si**(1.831182D+0) + 18.111 )**-1.0D+0


        else if ( p_si >= 14.628D+0 .and. p_si <= 18.210D+0 ) then

            vg = ( 8.5781572E-03 * p_si**(3.176484D+0) + 50.0D+0 )**-1.0D+0


        else if ( p_si > 18.210D+0 .and. p_si <= 20.253D+0 ) then

            vg = ( 3.5587113E-06 * p_si**(5.660939D+0) + 88.0D+0 )**-1.0D+0


        else if ( p_si > 20.253D+0 .and. p_si <= 21.50D+0 ) then

            vg = ( 3.558734E-16 * p_si**(13.03774D+0) + 138.0D+0 )**-1.0D+0

        end if


!       Formula to calculate specific volume per Section 5.1
        v_si = vg + (t_si - tsat_si) * ( 0.000466D+0/p_si - &
        ((0.12D+0/(t_si+100.0D+0) - 0.000106) * p_si**(0.1D+0))/ &
        sqrt ( 1.96E-08 * (t_si + 8.0D+0)**4.0D+0 - p_si**2.0D+0 ))


!       Convert from m**3/kg to ft**3/lbm
        v = v_si * 16.0185D+0

!       Calculate vapor density in lbm/ft**3
        rho = 1.0D+0/v



    end subroutine vap_vol


    subroutine vap_enth ( p, tsat, t, h )
! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the vapor enthalpy as a function
!   of temperature and pressure
!
!   This subroutine is valid for: saturation deg C < T < 450.0 deg C &
!                                  0.085 MPa < P < 21.5 MPa
!
!
!      p,  input, real ( kind = 8 ), the pressure
!                                    units: psia
!
!   tsat,  input, real ( kind = 8 ), the saturation temperature
!                                    units: deg F
!
!
!      t,  input, real ( kind = 8 ), the temperature
!                                    units: deg F
!
!
!      h, output, real ( kind = 8 ), the calculated enthalpy
!                                    units: BTU/lbm
!
!

        implicit none

        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: tsat
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: h


!       Variables used locally in subroutine
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: tsat_si
        real ( kind = 8 ) :: t_si
        real ( kind = 8 ) :: h_si
        real ( kind = 8 ) :: hg


!       Convert pressure from psia to MPa
        p_si = p * 0.00689476D+0

!       Convert temperatures from deg F to deg C
        tsat_si = (5.0D+0/9.0D+0)*(tsat - 32.0D+0)
        t_si = (5.0D+0/9.0D+0)*(t - 32.0D+0)



!       Correlations for hg from Section 5.2
        if ( p_si < 0.075D+0 .and. p_si > 21.55D+0 ) then

            write ( *, * ) ' THE PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '


        else if ( p_si > 0.075D+0 .and. p_si <= 0.348D+0 ) then

            hg = -4.0381938E-06 * ( 3.0D+0 - p_si )**15.72364D+0 &
            + 2750.0D+0


        else if ( p_si > 0.348D+0 .and. p_si <= 1.248D+0 ) then

            hg = -0.5767304 * exp( -1.66153D+0 * ( p_si - 3.2D+0 )) &
            + 2800.0D+0


        else if ( p_si > 1.248D+0 .and. p_si < 2.955D+0 ) then

            hg = -7.835986D+0 * ( p_si - 3.001D+0 )**2.0D+0 &
            - 2.934312D+0 * ( p_si - 3.001D+0 ) + 2803.71D+0


        else if ( p_si >= 2.955D+0 .and. p_si <= 6.522D+0 ) then

            hg = -1.347244D+0 * ( p_si - 2.999D+0 )**2.0D+0 &
            - 2.326913D+0 * ( p_si - 2.999D+0 ) + 2803.35D+0


        else if ( p_si > 6.522D+0 .and. p_si < 16.497D+0 ) then

            hg = -0.9219176D+0 * ( p_si - 9.00D+0 )**2.0D+0 &
            - 16.388835D+0 * ( p_si - 9.00D+0) + 2742.03D+0


        else if ( p_si >= 16.497D+0 .and. p_si < 20.193D+0 ) then

            hg = -3.532177D+0 * ( p_si - 8.00D+0 )**2.0D+0 &
            + 29.81305D+0 * ( p_si - 8.00D+0 ) + 2565.0D+0


        else if ( p_si >= 20.193D+0 .and. p_si <= 21.550D+0 ) then

            hg = -22.92521D+0 * ( p_si - 18.0D+0 )**2.0D+0 &
            + 44.23671D+0 * ( p_si - 18.0D+0 ) + 2415.01D+0


        end if


!       Formula to calculate enthalpy from Section 5.2
        h_si = hg + ( t_si - tsat_si ) * (( 0.28D+0 * exp( -0.008D+0 * ( t_si - 162.0D+0 ))) &
        + ( - 100.0D+0/t_si - 2.225D+0) &
        + ( 4.5D+0 * p_si / sqrt( 7.4529E-06 * t_si**3 - p_si**2 )))



!       Convert enthalpy from kJ/kg to BTU/lb
        h = h_si * 0.4299D+0


    end subroutine vap_enth

    FUNCTION mu0(tau)

        IMPLICIT NONE
        !
        !> Constant coefficients "I"
        !
        INTEGER(KIND=4), PARAMETER, DIMENSION(1:4) :: I = &
            (/ 0, 1, 2, 3 /)
        !
        !> Constant coefficients "H"
        !
        REAL(KIND=8), PARAMETER, DIMENSION(1:4) :: H = &
            (/ 1.67752, 2.20462, 0.6366564, -0.241605 /)
        !
        ! Arguments
        !
        REAL(KIND=8), INTENT(IN) :: tau
        REAL(KIND=8) :: mu0
        !
        ! Compute the viscosity in the dilute-gas limit.
        !
        mu0 = 100.0D+00/(sqrt(tau)*SUM(H*(tau**I)))

    END FUNCTION mu0
    !-------------------------------------------------------------------------
    !
    !> This function computes the contribution to viscosity due to finite
    !> density.
    !
    !> @param[in] del Dimensionless density parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The contribution to viscosity due to finite density.
    !
    !-------------------------------------------------------------------------
    FUNCTION mu1(del, tau)

        IMPLICIT NONE
        !
        !> Constant coefficients "H"
        !
        REAL(KIND=8), DIMENSION(1:6,1:7) :: H
        !
        ! Arguments
        !
        REAL(KIND=8), INTENT(IN) :: del, tau
        REAL(KIND=8) :: mu1
        !
        ! Variables
        !
        INTEGER(KIND=4) :: i, j
        REAL(KIND=8) :: sum, tau1
        !
        ! Compute the contribution to viscosity due to finite density
        !
        sum  = 0.0D+0
        H(1,:) = (/ 5.20094D-01, 2.22531D-01, -2.81378D-01, 1.61913D-01, -3.25372D-02, 0.0D+00, 0.0D+00  /)
        H(2,:) = (/ 8.50895D-02, 9.99115D-01, -9.06851D-01, 2.57399D-01, 0.0D+00, 0.0D+00, 0.0D+00       /)
        H(3,:) = (/ -1.08374D+00, 1.88797D+00, -7.72479D-01, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00          /)
        H(4,:) = (/ -2.89555D-01, 1.26613D+00, -4.89837D-01, 0.0D+00, 6.98452D-02, 0.0D+00, -4.35673D-03 /)
        H(5,:) = (/ 0.0D+00, 0.0D+00, -2.57040D-01, 0.0D+00, 0.0D+00, 8.72102D-03, 0.0D+00               /)
        H(6,:) = (/ 0.0D+00, 1.20573D-01, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, -5.93264D-04               /)

        DO i = 1, 6
            tau1 = (tau-1.0D+0)**(i-1)
            DO j = 1, 7
                IF (H(i,j) == 0.0D+0) CYCLE
                sum = sum + H(i,j)*tau1*((del-1.0D+0)**(j-1))
            END DO
        END DO

        mu1 = EXP(del*sum)

    END FUNCTION mu1
    !-------------------------------------------------------------------------
    !
    !> This function computes the viscosity as a function of density and
    !> temperature
    !
    !> @param[in] density     The steam density (Kg/m3)
    !> @param[in] temperature The steam temperature (K)
    !
    !> @return The viscosity in [Pa.sec]
    !
    !-------------------------------------------------------------------------
    subroutine viscosity( p, t, mu )

        IMPLICIT NONE
        !
        !> Star viscosity for the computation of viscosity(density,temperature)
        !> in [Pa.sec]
        !
        REAL(KIND=8), PARAMETER :: VISCOSITY_MUSTAR = 1.0D-06
        !
        ! Arguments
        !
        REAL(KIND=8), INTENT(IN) :: p, t
        real ( kind = 8 ), intent ( out ) :: mu
        !
        ! Variables
        !
        REAL(KIND=8) :: del, tau

        real ( kind = 8 ) :: tkel
        real ( kind = 8 ) :: psat
        real ( kind = 8 ) :: rho
        real ( kind = 8 ) :: rho_kgm3


!       Convert input temperature from deg F to K
        tkel = (5.0D+0/9.0D+0)*( t + 459.67D+0 )

!       Following code used to calculate the density
        call p_sat ( t, psat )
        call liq_dens ( p, t, rho )

!       Convert density from lb/ft**3 to kg/m**3
        rho_kgm3 = rho * 16.0184634D+0

        !
        ! Compute the viscosity
        !
        del = rho_kgm3 / 322.0D+00
        tau = 647.096D+00 / tkel
        mu = VISCOSITY_MUSTAR*mu0(tau)*mu1(del,tau)

!       Convert from Pa*sec to lb/(ft-sec)
        mu = mu * 0.672D+0


    end subroutine viscosity


    subroutine viscosity2phase( p, t, mu_liq, mu_vap )

        IMPLICIT NONE
        !
        !> Star viscosity for the computation of viscosity(density,temperature)
        !> in [Pa.sec]
        !
        REAL(KIND=8), PARAMETER :: VISCOSITY_MUSTAR = 1.0D-06
        !
        ! Arguments
        !
        REAL(KIND=8), INTENT(IN) :: p, t
        real ( kind = 8 ), intent ( out ) :: mu_liq, mu_vap
        !
        ! Variables
        !
        REAL(KIND=8) :: del, tau

        real ( kind = 8 ) :: tkel
        real ( kind = 8 ) :: psat
        real ( kind = 8 ) :: rho_liq
        real ( kind = 8 ) :: rho_kgm3
        real ( kind = 8 ) :: rho_vap
        real ( kind = 8 ) :: dummy_vol
        real ( kind = 8 ) :: rho_liq_si
        real ( kind = 8 ) :: rho_vap_si
        real ( kind = 8 ) :: del_liq
        real ( kind = 8 ) :: tau_liq
        real ( kind = 8 ) :: del_vap
        real ( kind = 8 ) :: tau_vap


!       Convert input temperature from deg F to K
        tkel = (5.0D+0/9.0D+0)*( t + 459.67D+0 )

!       Following code used to calculate the density
        call liq_dens ( p, t, rho_liq )
        call vap_vol ( p ,t, dummy_vol, rho_vap )


!       Convert density from lb/ft**3 to kg/m**3
        rho_liq_si = rho_liq * 16.0184634D+0
        rho_vap_si = rho_vap * 16.0184634D+0

        !
        ! Compute the viscosity
        !
        del_liq = rho_liq_si / 322.0D+00
        tau_liq = 647.096D+00 / tkel
        mu_liq = VISCOSITY_MUSTAR*mu0(tau_liq)*mu1(del_liq,tau_liq)

        del_vap = rho_vap_si / 322.0D+00
        tau_vap = 647.096D+00 / tkel
        mu_vap = VISCOSITY_MUSTAR*mu0(tau_vap)*mu1(del_vap,tau_vap)

!       Convert from Pa*sec to lb/(ft-sec)
        mu_liq = mu_liq * 0.672D+0
        mu_vap = mu_vap * 0.672D+0


    end subroutine viscosity2phase




    subroutine surtension ( t, sigma )
! -------------------------------------------------------------------------- 80
!
! Discussion:
!   This subroutine calculates liquid water surface tension as a function
!   of temperature
!
! Arguments:
!
!       t,  input,     real ( kind = 8 ), this is the temperature of the coolant
!                                         units: deg F
!
!   sigma, output,     real ( kind = 8 ), the calculated surface tension of the coolant
!                                         units: lb/sec**2
!
!
!
        implicit none

        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: sigma


!       Local variables
        real ( kind = 8 ) :: tkel
        real ( kind = 8 ) :: tc = 647.15D+0
        real ( kind = 8 ) :: B = 235.8E-03
        real ( kind = 8 ) :: b2 = -0.625D+0
        real ( kind = 8 ) :: mu = 1.256D+0

!       Convert input temperature to K from deg F
        tkel = (5.0D+0/9.0D+0)*(t + 459.67D+0)

!       Calculate surface tension
        sigma = B * ((tc-tkel)/tc)**mu * (1.0D+0 + b2*((tc - tkel)/tc))

!       Convert the surface tension from N/m to lb/sec**2
        sigma = sigma * 2.205D+0


    end subroutine surtension



    subroutine saturation_p ( t, psat )
! -------------------------------------------------------------------------- 80
!
! Discussion:
!   This subroutine calculates the saturation pressure as a function of
!   temperature per Table A-1 of:
!
!   Extensions to the approximation functions for fast calculation of
!   saturated water properties ( Garland, Wilson, Bartak, Cizek, Stasny, Zentrich )
!
! Arguments:
!
!        t,   input,    real ( kind = 8 ), the input temperature
!                                          units: deg F
!
!
!     psat,  output,    real ( kind = 8 ), the calculated saturation pressure
!                                          units: psi
!
!

        implicit none

    !   Input arguments
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: psat

    !   Local arguments
        real ( kind = 8 ) :: a
        real ( kind = 8 ) :: b
        real ( kind = 8 ) :: c
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: t_si


    !   Convert input temperature to deg C
        t_si = (5.0D+0/9.0D+0)*(t - 32.0D+0)


    !   Table A-1
        if ( t_si >= 17.511D+0 .and. t_si <= 56.275D+0 ) then

            a = 99.2D+0
            b = 270.1210D+0
            c = 7.4063650D+0

        else if ( t_si >= 56.275D+0 .and. t_si < 90.880D+0 ) then

            a = 78.2D+0
            b = 254.6831D+0
            c = 6.4058216D+0

        else if ( t_si >= 90.880D+0 .and. t_si <139.781D+0 ) then

            a = 57.0D+0
            b = 236.2315D+0
            c = 5.6029720D+0

        else if ( t_si >= 139.781D+0 .and. t_si < 203.662D+0 ) then

            a = 28.0D+0
            b = 207.9248D+0
            c = 4.7785040D+0

        else if ( t_si >= 203.662D+0 .and. t_si < 299.407D+0 ) then

            a = 5.0D+0
            b = 185.0779D+0
            c = 4.3043760D+0

        else if ( t_si >= 299.407D+0 .and. t_si < 355.636D+0 ) then

            a = 16.0D+0
            b = 195.1819D+0
            c = 4.4608430D+0

        else if ( t_si >= 355.563D+0 .and. t_si < 373.00D+0 ) then

            a = 50.0D+0
            b = 227.2963D+0
            c = 4.9607850D+0

        end if

    !   Equation from Table A-1
        p_si = ((t_si + a )/b)**c


    !   Convert to psia
        psat = p_si * 145.038D+0



    end subroutine saturation_p


    subroutine saturation_t ( p, tsat )

    ! -------------------------------------------------------------------------- 80
    !
    ! Discussion:
    !   This subroutine calculates the saturation temperature as a function of
    !   pressure per Table A-2 of:
    !
    !   Extensions to the approximation functions for fast calculation of
    !   saturated water properties ( Garland, Wilson, Bartak, Cizek, Stasny, Zentrich )
    !
    ! Arguments:
    !
    !        p,   input,    real ( kind = 8 ), the input pressure
    !                                          units: psi
    !
    !
    !     tsat,  output,    real ( kind = 8 ), the calculated saturation temperature for
    !                                          the input pressure
    !                                          units: deg F
    !
    !


        implicit none

    !   Input arguments
        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( inout ) :: tsat

    !   Local arguments
        real ( kind = 8 ) :: a
        real ( kind = 8 ) :: b
        real ( kind = 8 ) :: c
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: t_si


    !   Convert from psi to MPa
        p_si = p * 0.00689476D+0


    !   From Table A-2
        if ( p_si >= 0.002D+0 .and. p_si < 0.01672D+0 ) then

            a = 270.1210D+0
            b = 0.135019D+0
            c = -99.2D+0

        else if ( p_si >= 0.01672D+0 .and. p_si < 0.07250D+0 ) then

            a = 254.6831D+0
            b = 0.156108D+0
            c = -78.2D+0

        else if ( p_si >= 0.07250D+0 .and. p_si <= 0.359D+0 ) then

            a = 236.2315D+0
            b = 0.1784767D+0
            c = -57.0D+0

        else if ( p_si > 0.3590D+0 .and. p_si <= 1.6760D+0 ) then

            a = 207.9248D+0
            b = 0.2092705D+0
            c = -28.0D+0

        else if ( p_si > 1.676D+0 .and. p_si <= 8.51100D+0 ) then

            a = 185.0779D+0
            b = 0.2323217D+0
            c = -5.0D+0

        else if ( p_si > 8.5110D+0 .and. p_si < 17.690D+0 ) then

            a = 195.1819D+0
            b = 0.2241729D+0
            c = -16.0D+0

        else if ( p_si >= 17.690D+0 .and. p_si <= 21.50D+0 ) then

            a = 227.2963D+0
            b = 0.2241729D+0
            c = -50.0D+0

        end if


    !   Equation from Table A-2
        t_si = a*p_si**b + c


    !   Convert the saturation temperature into deg F
        tsat = t_si*(9.0D+0/5.0D+0) + 32.0D+0



    end subroutine saturation_t



    subroutine saturation_h ( p, t, hliq, hvap )

! -------------------------------------------------------------------------- 80
!
! Discussion:
!   This subroutine calculates the saturation enthalpy values
!   as a function of temperature and pressure per
!   Table A-4 and Table A-8 of Ref. 1
!
! Reference:
!   1. Extensions to the approximation functions for fast calculation of
!      saturated water properties ( Garland, Wilson, Bartak, Cizek, Stasny, Zentrich )
!
! Arguments:
!
!     p,  input,   real ( kind = 8 ), the input pressure
!                                     units: psi
!
!
!     t,  input,   real ( kind = 8 ), the input temperature
!                                     units: deg F
!
!
!  hliq, output,   real ( kind = 8 ), the specific enthalpy for the liquid phase
!                                     at saturation (see: Table A-4 of Ref. 1.)
!                                     units: BTU/lbm
!
!
!  hvap, output,   real ( kind = 8 ), the specific enthalpy for the vapor phase
!                                     at saturation (see: Table A-4 of Ref. 1)
!                                     units: BTU/lbm
!
!
        implicit none

    !   Argument variables
        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: hliq
        real ( kind = 8 ), intent ( inout ) :: hvap

    !   Local variables
        real ( kind = 8 ) :: a
        real ( kind = 8 ) :: b
        real ( kind = 8 ) :: c
        real ( kind = 8 ) :: d
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: t_si
        real ( kind = 8 ) :: hliq_si
        real ( kind = 8 ) :: hvap_si
        real ( kind = 8 ) :: psat
        real ( kind = 8 ) :: psat_si


    !   Calculate the saturated pressure
        call saturation_p ( t, psat )


    !   Convert saturation pressure and input pressure from psi to MPa
        psat_si = psat * 0.00689476D+0
        p_si = p * 0.00689476D+0


    !   Calculate the liquid phase enthalpy using Table A-4 from Ref. 1
        if ( p_si >= 0.0020D+0 .and. p_si < 0.0173D+0 ) then

            a = 1128.7770D+0
            b = 0.1351960D+0
            c = -413.72D+0

        else if ( p_si >= 0.0173 .and. p_si < 0.1028D+0 ) then

            a = 1050.7085D+0
            b = 0.1617970D+0
            c = -306.50D+0

        else if ( p_si >= 0.1028D+0 .and. p_si <= 0.9420D+0 ) then

            a = 912.1779D+0
            b = 0.2061637D+0
            c = -150.0D+0

        else if ( p_si >= 0.9420D+0 .and. p_si < 4.020D+0 ) then

            a = 638.0621D+0
            b = 0.2963192D+0
            c = 125.0D+0

        else if ( p_si >= 4.020D+0 .and. p_si < 9.964D+0 ) then

            a = 373.7665D+0
            b = 0.4235532D+0
            c = 415.0D+0

        else if ( p_si >= 9.964D+0 .and. p_si < 16.673D+0 ) then

            a = 75.38673D+0
            b = 0.8282384D+0
            c = 900.0D+0

        else if ( p_si >= 16.673D+0 .and. p_si < 20.396D+0 ) then

            a = 0.1150827D+0
            b = 2.711412D+0
            c = 1440.0D+0

        else if ( p_si >= 20.396D+0 .and. p_si <= 21.50D+0 ) then

            a = 9.1417257E-14
            b = 11.47287D+0
            c = 1752.00D+0

        end if

    !   Calculate liquid phase saturated enthalpy equation in Table A-4 of Ref. 1
        hliq_si = a*psat_si**b + c


    !   Convert enthalpy from KJ/kg to BTU/lb
        hliq = hliq_si * 0.4299D+0

    !   Calculate the vapor phase enthalpy using Table A-8 from Ref. 1
        if ( p_si >= 0.002D+0 .and. p_si < 0.1379D+0 ) then

            hvap_si = 529.44008D+0 * p_si**0.108652D+0 + 2263.5D+0

        else if ( p_si >= 0.1379D+0 .and. p_si < 0.348D+0 ) then

            hvap_si = -4.0381938E-6 * (3.0D+0 - p_si)**15.72364D+0 + 2750.0D+0

        else if ( p_si > 0.348D+0 .and. p_si <= 1.248D+0 ) then

            hvap_si = -0.5767304D+0 * exp( -1.66153D+0*(p_si - 3.2D+0)) + 2800.0D+0

        else if ( p_si > 1.2488D+0 .and. p_si < 2.955D+0 ) then

            a = -7.835986D+0
            b = -3.001D+0
            c = -2.934312D+0
            d = 2803.71D+0

            hvap_si = a*( p_si + b )**2.0D+0 + c*( p_si + b ) + d

        else if ( p_si >= 2.955D+0 .and. p_si <= 6.522D+0 ) then

            a = -1.347244D+0
            b = -2.999D+0
            c = -2.326913D+0
            d = 2803.35D+0

            hvap_si = a*( p_si + b )**2.0D+0 + c*( p_si + b ) + d


        else if ( p_si > 6.522D+0 .and. p_si < 16.497D+0 ) then

            a = -0.9219176D+0
            b = -9.0D+0
            c = -16.38835D+0
            d = 2742.03D+0

            hvap_si = a*( p_si + b )**2.0D+0 + c*( p_si + b ) + d


    !   Table A-8 seems to have a typo here. Table A-8 states the
    !   upper limit as 40.193. I believe that this should 20.193
        else if ( p_si >= 16.497 .and. p_si < 20.193D+0 ) then

            a = -3.532177D+0
            b = -8.0D+0
            c = 29.81305D+0
            d = 2565.00D+0

            hvap_si = a*( p_si + b )**2.0D+0 + c*( p_si + b ) + d


        else if ( p_si >= 20.193D+0 .and. p_si <= 21.50D+0 ) then

            a = -22.92521D+0
            b = -18.0D+0
            c = 44.23671D+0
            d = 2415.01D+0

            hvap_si = a*( p_si + b )**2.0D+0 + c*( p_si + b ) + d

        end if


    !   Convert enthalpy from KJ/kg to BTU/lb
        hvap = hvap_si * 0.4299D+0


    end subroutine saturation_h



    subroutine saturation_rho ( p, t, rholiq, rhovap )

! -------------------------------------------------------------------------- 80
!
! Discussion:
!   This subroutine calculates the saturation density values
!   as a function of temperature and pressure per
!   Table A-3 and Table A-7 of Ref. 1
!
! Reference:
!   1. Extensions to the approximation functions for fast calculation of
!      saturated water properties ( Garland, Wilson, Bartak, Cizek, Stasny, Zentrich )
!
! Arguments:
!
!     p,  input,   real ( kind = 8 ), the input pressure
!                                     units: psi
!
!
!     t,  input,   real ( kind = 8 ), the input temperature
!                                     units: deg F
!
!
!  rholiq, output,   real ( kind = 8 ), the specific density for the liquid phase
!                                       at saturation (see: Table A-4 of Ref. 1.)
!                                       units: lbm/ft**3
!
!
!  rhovap, output,   real ( kind = 8 ), the specific enthalpy for the vapor phase
!                                     at saturation (see: Table A-4 of Ref. 1)
!                                     units: lbm/ft**3
!
!

        implicit none

!       Input arguments
        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: rholiq
        real ( kind = 8 ), intent ( inout ) :: rhovap


!       Local arguments
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: psat
        real ( kind = 8 ) :: psat_si
        real ( kind = 8 ) :: rholiq_si
        real ( kind = 8 ) :: rhovap_si
        real ( kind = 8 ) :: volvap_si
        real ( kind = 8 ) :: a, b, c

!       Determine saturation pressure
        call saturation_p ( t, psat )

!       Convert pressure (input and saturation) from psi to MPa
        p_si = p * 0.00689476D+0
        psat_si = psat * 0.00689476D+0


!       Calculate the density of the liquid phase at saturation using Table A-3
        if ( p_si >= 0.002D+0 .and. p_si < 0.01468D+0 ) then

            a = 1.9118E-04
            b = 0.546472D+0
            c = 0.0009947D+0

            rholiq_si = ( a*psat_si**b + c )**-1.0D+0


        else if ( p_si >= 0.01468 .and. p_si < 0.275D+0 ) then

            a = 1.380934E-04
            b = 0.388715D+0
            c = 0.000987D+0

            rholiq_si = ( a*psat_si**b + c )**-1.0D+0


        else if ( p_si <= 0.275D+0 .and. p_si <= 1.00D+0 ) then

            a = 1.2746977E-4
            b = 0.4644339D+0
            c = 0.00010D+0

            rholiq_si = ( a*psat_si**b + c )**-1.0D+0


        else if ( p_si > 1.00D+0 .and. p_si <= 3.880D+0 ) then

            a = 1.0476071E-4
            b = 0.5651090D+0
            c = 0.0010220D+0

            rholiq_si = ( a*psat_si**b + c )**-1.0D+0


        else if ( p_si > 3.880D+0 .and. p_si <= 8.840D+0 ) then

            a = 3.2836717E-5
            b = 1.00D+0
            c = 0.00112174735D+0

            rholiq_si = ( a*psat_si**b + c )**-1.0D+0


        else if ( p_si > 8.840D+0 .and. p_si <= 14.463D+0 ) then

            rholiq_si = ( 3.3551046E-4 * exp(5.8403566E-2*psat_si) + 0.00085D+0 )**-1.0D+0


        else if ( p_si > 14.463D+0 .and. p_si < 18.052D+0 ) then

            a = 3.1014626E-8
            b = 3.2847540D+0
            c = 0.001430D+0

            rholiq_si = ( a*psat_si**b + c )**-1.0D+0


        else if ( p_si >= 18.052D+0 .and. p_si < 20.204D+0 ) then

            a = 1.5490787E-11
            b = 5.72050D+0
            c = 0.0016050D+0

            rholiq_si = ( a*psat_si**b + c )**-1.0D+0


        else if ( p_si >= 20.204D+0 .and. p_si <= 21.50D+0 ) then

            a = 4.1035988E-24
            b = 15.033290D+0
            c = 0.001890D+0

            rholiq_si = ( a*psat_si**b + c )**-1.0D+0

        end if


!       Convert from kg/m**3 to lbm/ft**3
        rholiq = rholiq_si * 0.062428D+0


!       Calculating the specific volume for the vapor phase at saturation using
!       Table A-7 of Ref. 1
        if ( p_si >= 0.002D+0 .and. p_si < 0.2139D+0 ) then

            a = 5.0981616D+0
            b = 0.936226D+0
            c = -0.00025D+0


        else if ( p_si >= 0.2139D+0 .and. p_si < 1.1120D+0 ) then

            a = 5.126076D+0
            b = 0.9475862D+0
            c = 0.00120D+0


        else if ( p_si >= 1.112D+0 .and. p_si < 3.932D+0 ) then

            a = 4.630832D+0
            b = 1.038819D+0
            c = 0.520D+0


        else if ( p_si >= 3.932D+0 .and. p_si < 8.996D+0 ) then

            a = 2.868721D+0
            b = 1.252148D+0
            c = 3.800D+0


        else if ( p_si >= 8.996D+0 .and. p_si < 14.628D+0 ) then

            a = 0.5497653D+0
            b = 1.831182D+0
            c = 18.111D+0


        else if ( p_si >= 14.628D+0 .and. p_si <= 18.210D+0 ) then

            a = 8.5791582E-3
            b = 3.176484D+0
            c = 50.00D+0


        else if ( p_si > 18.210D+0 .and. p_si <= 20.253 ) then

            a = 3.5587113E-6
            b = 5.660939D+0
            c = 88.00D+0


        else if ( p_si > 20.253D+0 .and. p_si <= 21.50D+0 ) then

            a = 3.558734E-11
            b = 13.03774D+0
            c = 138.00D+0

        end if

        volvap_si = ( a*p_si**b + c )**-1.0D+0

!       Convert from specific volume to density
        rhovap_si = 1.0D+0/volvap_si


!       Convert from kg/m**3 to lbm/ft**3
        rhovap = rhovap_si * 0.062428D+0


    end subroutine saturation_rho


end module mod_thermo




