SUBROUTINE RTX(M, X, NX, NL, R, T)
!===============================================================================
! PURPOSE:
!   To compute Rlm(x), Tlm(x) polynomials for all l = 0:L and a given m.
!   R, T = 0 are returned for m > l.
!
! INPUT:
!   M    I(1)    Order of the Fourier harmonic, m = 0...M, M < NL
!   X    D(NX)   Array of abscissas, -1.0 <= x <= 1.0
!   NX   I(1)    Length of X
!   NL   I(1)    Maximum order of the polynomials + 1, L1 > max(M, 2)
!
! OUTPUT:
!   R   D(NX, NL)   Polynomial Rlm(x) for the given m and all l = 0:NL-1
!   T   D(NX, NL)   Polynomial Tlm(x) for the given m and all l = 0:NL-1
!
! TREE:
!   -
!
! COMMENTS:
!   The matrix polynomial is defined following [1]:
!
!                      | Qlm(x)  0       0       0     |
!                      | 0       Rlm(x) -Tlm(x)  0     |
!             Plm(x) = | 0      -Tlm(x)  Rlm(x)  0     |,
!                      | 0       0       0       Qlm(x)|
!
!   where
!
!             Rlm(x) = -1/2*i^m*(Plm{+2} + Plm{-2}),
!             Tlm(x) = -1/2*i^m*(Plm{+2} - Plm{-2}),
!
!   are real functions of -1 <= x <= 1 and i = sqrt(-1).
!
!             Qlm(x) = sqrt[(l-m)!/(l+m)!]*Plm,
!             Plm = (1-x^2)^(m/2)(d/dmu)^m{Pl(x)},
!
!   where Pl(x) are the ordinary Legendre polynomials of an order l. Note
!   that (-1)^m is not used in Qlm. Generalized Legendre polynomials are [2]
!
!             Plm{n} = Alm{n}*(1 - x)^(-(n - m)/2)*(1 + x)^(-(n + m)/2)*
!                      [d/dmu]^(l-n)[(1 - x)^(l-m)*(1 + x)^(l+m)],
!
!   where n = -2, -0, +0, + 2, and
!
!   Alm{n} = [(-1)^(l-m)*i^(n-m)/2^l/(l-m)!]*[(l-m)!(l+n)!/(l+m)!/(l-n)!]^1/2.
!
!   This definition coincides with that widely used in literature [3-5].
!
! REFERENCES:
!   1. Siewert CE, JQSRT (2000), 64, pp.227-254
!   2. Gelfand IM et al., 1963: Representations of the rotation and Lorentz
!      groups and their applications. Oxford: Pergamon Press.
!   3. Hovenier JW et al., 2004: Transfer of Polarized Light in Planetary
!      Atmosphere. Basic Concepts and Practical Methods, Dordrecht: Kluwer
!      Academic Publishers.
!   4. Rozanov VV and Kokhanovsky AA, Atm.Res., 2006, 79, 241.
!   5. Y.Ota et al., JQSRT, 2010, 111, 878.
!===============================================================================
!
	IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: M, NX, NL
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NX), INTENT(IN) :: X
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NX, NL), INTENT(OUT) :: R, T
!
! LOCAL VARIABLES
    INTEGER &
        IL,  & ! Loop index over l, IL = L0:NL, L0=M1+1
        IL1, & ! IL-1
        IL2, & ! IL-2
        L0,  & ! Order (+1, prog.) of the first non-zero R&T for a given M
        M1     ! M1 = M+1
    REAL*8 &
        C0, & ! Coefficient for various purposes
        C1, & ! Coefficients
        C2, & !     in recurrence
        C3    !         relations for R & T
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NX) :: &
        X2,  & ! X*X
        SX,  & ! sqrt(1 - X2)
        SX2    ! Function of X for different purposes
!===============================================================================
!
	X2 = X*X
!
    DO IL = 1, 2
        R(:, IL) = 0.0D0
        T(:, IL) = 0.0D0
    END DO ! IL = 1, 2
!
    IF (M == 0) THEN
        T(:, 3) = 0.0D0 ! T(:, 1:2) = 0.0D0 - see above
        R(:, 3) = &     ! R(:, 1:2) = 0.0D0 - see above; sqrt(3/8) = 0.612
            0.61237243569579452454932101867647D0*(1.0D0 - X2)
        DO IL = 4, NL ! L0 = M+2
            C0 = 1.0D0/DSQRT( (IL - 3.0D0)*(IL + 1.0D0) )
            C1 = (2.0D0*IL - 3.0D0)*C0
            C3 = DSQRT( IL*(IL - 4.0D0) )*C0
            T(:, IL) = 0.0D0
            R(:, IL) = C1*X*R(:, IL-1) - C3*R(:, IL-2)
        END DO ! IL = 4, NL
    ELSEIF (M == 1) THEN
        T(:, 3) = -0.5D0*DSQRT( 1.0D0 - X2 ) ! T(:, 1:2) = 0.0D0 - see above
        R(:, 3) = T(:, 3)*X                  ! R(:, 1:2) = 0.0D0 - see above
        DO IL = 4, NL ! L0 = M+2
            IL1 = IL-1
            IL2 = IL-2
!
            C0 = IL1/DSQRT( (IL - 3.0D0)*IL2*IL*(IL + 1.0D0) )
            C1 = (2.0D0*IL - 3.0D0)*C0
            C2 = (4.0D0*IL - 6.0D0)*C0/(1.0D0*IL2*IL1)
            C3 = DSQRT( (1.0D0 - 1.0D0/(IL2*IL2))*IL*(IL - 4.0D0) )*C0
!
            R(:, IL) = C1*X*R(:, IL1) - C2*T(:, IL1) - C3*R(:, IL2)
            T(:, IL) = C1*X*T(:, IL1) - C2*R(:, IL1) - C3*T(:, IL2)
        END DO ! IL = 4, NL
    ELSEIF (M == 2) THEN
        R(:, 3)  = 0.25D0 + 0.25D0*X2 ! R(:, 1:2) = 0.0D0 - see above
        T(:, 3)  = 0.5D0*X            ! T(:, 1:2) = 0.0D0 - see above
!
        DO IL = 4, NL ! L0 = M+2
            IL1 = IL-1
            IL2 = IL-2
!
            C0 = IL1/((IL + 1.0D0)*(IL - 3.0D0))
            C1 = (2.0D0*IL -  3.0D0)*C0
            C2 = (8.0D0*IL - 12.0D0)*C0/(1.0D0*IL1*IL2)
            C3 = (IL2 - 4.0D0/IL2)*C0
!
            R(:, IL) = C1*X*R(:, IL1) - C2*T(:, IL1) - C3*R(:, IL2)
            T(:, IL) = C1*X*T(:, IL1) - C2*R(:, IL1) - C3*T(:, IL2)
        END DO ! IL = 4, NL
    ELSEIF (M == 3) THEN
        R(:, 3) = 0.0D0 ! R(:, 1:2) = 0.0D0 - see above
        T(:, 3) = 0.0D0 ! T(:, 1:2) = 0.0D0 - see above
        R(:, 4) = &    ! sqrt(3/32)=0.306...
                  0.30618621784789726227466050933824D0*DSQRT( 1.0D0 - X2 )
        T(:, 4)  = 2.0D0*R(:, 4)*X
        R(:, 4)  = R(:, 4)*(1.0D0 + X2)
!
        DO IL = 5, NL ! L0 = M1+1
            IL1 = IL-1
            IL2 = IL-2
!
            C0 = IL1/DSQRT( (IL + 2.0D0)*(IL - 4.0D0)* &
                            (IL - 3.0D0)*(IL + 1.0D0) )
            C1 = ( 2.0D0*IL -  3.0D0)*C0
            C2 = (12.0D0*IL - 18.0D0)*C0/((IL - 2.0D0)*(IL - 1.0D0))
            C3 = DSQRT( (1.0D0 - 9.0D0/(1.0D0*IL2*IL2))* &
                        (1.0D0*IL2*IL2 - 4.0D0) )*C0
!
            R(:, IL) = C1*X*R(:, IL1) - C2*T(:, IL1) - C3*R(:, IL2)
            T(:, IL) = C1*X*T(:, IL1) - C2*R(:, IL1) - C3*T(:, IL2)
        END DO ! IL = 5, NL
	ELSE ! M > 3
	    M1 = M+1
        L0 = M+2
!
        DO IL = 3, M
            R(:, IL) = 0.0D0 ! R(:, 1:2) = 0.0D0 - see above
            T(:, IL) = 0.0D0 ! T(:, 1:2) = 0.0D0 - see above
        END DO ! IL = 3, M
!
        SX = 1.0D0 - X2
        SX2 = SX**(M/2-1)
        IF (MOD( M-2, 2 ) == 1) SX2 = SX2*DSQRT( SX )
!
        C0 = DSQRT( 0.03125D0*M*(M - 1.0D0) ) ! sqrt(1/32)=0.03125
        DO IL = 3, M
            C0 = 0.5D0*C0*DSQRT( 1.0D0*M/IL + 1.0D0 )
        END DO ! IL = 3, M
!
        R(:, M1)  = C0*SX2
        T(:, M1)  = 2.0D0*X*R(:, M1)
        R(:, M1)  = (1.0D0 + X2)*R(:, M1)
!
        DO IL = L0, NL
            IL1 = IL-1
            IL2 = IL-2
!
            C0 = IL1/DSQRT( 1.0D0*(IL1+M)*(IL1-M)*(IL - 3.0D0)*(IL + 1.0D0) )
            C1 = (2.0D0*IL - 3.0D0)*C0
            C2 = 2.0D0*M*(2.0D0*IL - 3.0D0)*C0/((IL - 2.0D0)*IL1)
            C3 = DSQRT( (1.0D0 - 1.0D0*M*M/(IL2*IL2))*(IL2*IL2 - 4.0D0) )*C0
!
            R(:, IL) = C1*X*R(:, IL1) - C2*T(:, IL1) - C3*R(:, IL2)
            T(:, IL) = C1*X*T(:, IL1) - C2*R(:, IL1) - C3*T(:, IL2)
        END DO ! IL = L0, NL
    END IF ! M == 0
!
END SUBROUTINE RTX
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 12Jun15 - (IL-3)*(IL+1) and similar operators are redefined as (IL - 3.0D0)*
!           (IL + 1.0D0) to avoid integer overflow: IL*IL > max_int
!
! 11Apr15 - Created from QRTm and tested against it for 1000 values of x,
!           m = 0:10, 13, 99, 256, & 1000. All l = 0:1023. 
!===============================================================================