SUBROUTINE QX0(M, X0, NL, Q)
!===============================================================================
! PURPOSE:
!   To compute the Qlm(x0) polynomials at a single point x0 for all l = 0:L
!   (L = 1:NL in the code) & a given Fourier order m = 0, 1, 2 ... . Zero values
!   are returned for all l < m.
!
! INPUT:
!   M    I(1)   Order of the Fourier harmonic, m = 0, 1, 2 ...
!   X0   D(1)   Abscissa, -1.0 <= x0 <= 1.0
!   NL   I(1)   Number of polynomials (orders) to be computed
!
! OUTPUT:
!   Q   D(NL)   Polynomials for the given M and all l = 0:NL-1
!
! TREE:
!   -
!
! COMMENTS:
!   The following polynomials are computed
!
!   Qlm(x) = sqrt[(l-m)!/(l+m)!]*Plm(x),                                     (1)
!   Plm(x) = (1-x^2)^(m/2)(dPl(x)/dx)^m,                                     (2)
!
!   using recurrence relation. Pl(x) are the ordinary Legendre polynomials of an
!   order l. Note, unlike in [2] (-1)^m is not used in Qlm. Refer to [1-4] for
!   details as well as to the QRTm.f90, QPOL.f90, QPOL1.f90 subroutines and the
!   developer's drafts.
!
! REFERENCES:
!   1. Gelfand IM et al., 1963: Representations of the rotation and Lorentz
!      groups and their applications. Oxford: Pergamon Press.
!   2. Hovenier JW et al., 2004: Transfer of Polarized Light in Planetary
!      Atmosphere. Basic Concepts and Practical Methods, Dordrecht: Kluwer
!      Academic Publishers.
!   3. http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html
!   4. http://www.mathworks.com/help/matlab/ref/legendre.html
!===============================================================================
!
	IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: M, NL
    REAL*8, INTENT(IN) :: X0
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NL), INTENT(OUT) :: Q
!
! LOCAL VARIABLES
    INTEGER &
        IL, & ! Loop index over order of polynomial, IL = M1:NL
        L0, & ! Initial order of in case of M > 0
        M1, & ! M1 = M+1 for convenience
        M2    ! M*2 for convenience
    REAL*8 &
        C0,  & ! Coefficient for polynomials with M > 0
        IL2, & ! DBLE(IL*IL)
        MM,  & ! DBLE(M*M)
        SX,  & ! [sqrt(1 - X2)]^M
        X2     ! X0*X0
!===============================================================================
!
    IF (M == 0) THEN
        Q(1) = 1.0D0
        Q(2) = X0
        DO IL = 3, NL
            Q(IL) = &
             ((2.0D0*IL - 3.0D0)*X0*Q(IL-1) - (IL - 2.0D0)*Q(IL-2))/(IL - 1.0D0)
        END DO ! IL = 3, IL1
    ELSE ! M > 0
        M1 = M+1
        M2 = M*2
        L0 = M+2
!
        X2 = X0*X0
        IF (M == 1) THEN
            SX = DSQRT( 1.0D0 - X2 )
        ELSE IF (M == 2) THEN
            SX = 1.0D0 - X2
        ELSE IF (M == 3) THEN
            SX = (1.0D0 - X2)*DSQRT( 1.0D0 - X2 )
        ELSE ! M > 3
            SX = (1.0D0 - X2)**(M/2)
            IF (MOD(M, 2) == 1) SX = SX*DSQRT( 1.0D0 - X2 )
        END IF ! M = 1
!
        DO IL = 1, M
            Q(IL) = 0.0D0
        END DO ! IL = 1, M
!
        C0 = 1.0D0
        DO IL = 2, M2, 2
            C0 = C0*(1.0D0 - 1.0D0/IL)
        END DO ! IL = 2, M2, 2
        C0 = DSQRT(C0)
!
        Q(M1) = C0*SX
        DO IL = L0, NL
            IL2 = 1.0D0*IL*IL
            MM  = 1.0D0*M*M
            Q(IL) = ( (2.0D0*IL - 3.0D0)*X0*Q(IL-1) -                          &
                      DSQRT( IL2 - 4.0D0*IL + 4.0D0 - MM )                     &
                      *Q(IL-2) )/DSQRT(IL2 - 2.0D0*IL + 1.0D0 - MM)                      
        END DO ! IL = L0, NL
    END IF ! M = 0
!
END SUBROUTINE QX0
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 12Jun15 - (IL-3)*(IL+1) and similar operators are redefined as (IL - 3.0D0)*
!           (IL + 1.0D0) to avoid integer overflow: IL*IL > max_int
!
! 05Apr15 - This subroutine was built using QX and tested against QPOL1 for 1000
!           values of x0 = [-1:1], NL = 512, M = 0:255.
!===============================================================================