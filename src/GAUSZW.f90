SUBROUTINE GAUSZW(X1, X2, N, X, W)
!===============================================================================
! PURPOSE:
!   To compute the abscissas X and weights W for the Gauss-Legendre quadrature
!   in a range [X1, X2]
!
! INPUT:
!   X1   R(1)   Interval boundary 1
!   X2   R(1)   Interval boundary 2
!   N    I(1)   The order of the Legendre polynomial
!
! OUTPUT:
!   X   D(N)    Zeros (abscissas) of the Legendre polynomial of an order N
!   W   D(N)    Corresponding weighting coefficients (Christoffel numbers)
!
! TREE:
!   -
!
! COMMENTS:
!   Note that computation time for N ~ 1000-10000 depends on how small the
!   small number YEPS is. Example of the output data for N = 8 and [X1, X2]=
!   [-1,+1] is (note the symmetry relation):
!                        X(i)          W(i)
!             i=1:  -0.96028986    0.10122854
!             i=2:  -0.79666648    0.22238103
!             i=3:  -0.52553241    0.31370665
!             i=4:  -0.18343464    0.36268378
!             i=5:   0.18343464    0.36268378
!             i=6:   0.52553241    0.31370665
!             i=7:   0.79666648    0.22238103
!             i=8:   0.96028986    0.10122854
!
!   Gauss numerical integration is
!
!             int(f(x), -1..1) ~= sum(wj*f(xj), j = 1..n).
!
!   Double-Gauss scheme for an even N is
!
!             int(f(x), -1..1) = int(f(x), -1..0) + int(f(x), 0..1) =
!   ~= sum(w2j*f(x2j), j=1..n/2, x2j<0) + sum(w2j*f(x2j), j=1..N/2, x2j>0).
!
!   The relations between Gauss and Double-Gauss weights and zeros are
!
!             w2(+j) = 0.5wj, x2(+j) = 0.5(xj + 1) > 0,
!             w2(-j) = w2(j), x2(-j) = 0.5(xj - 1) < 0,
!
!   where xj and wj are computed using Gauss scheme (i.e. THIS SUBROUTINE) and 
!   N/2 knots. Refer to [2].
!
! REFERENCES:
!   1. Taken from SCIATRAN3.1 package (GAULEG) with minor changes.
!      The original text of GAULEG is given below
!   2. Sykes J. Approximate integration of the equation of transfer,
!      Mon.Not.Roy.Astron.Soc, 1951, 11, 377.
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    REAL*8, PARAMETER :: &
        YEPS = 3.0D-14,  &        ! A small number
        PI = 3.1415926535897938D0 ! The number PI
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: N
    REAL*8, INTENT(IN) :: X1, X2
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(N), INTENT(OUT) :: X, W
!
! LOCAL VARIABLES
    INTEGER M, I, J
    REAL*8  YXM, YXL, YZ, YP1, YP2, YP3, YPP, YZ1
!===============================================================================
!
    M = (N + 1)/2
    YXM = 0.5D0*(X2 + X1)
    YXL = 0.5D0*(X2 - X1)
    DO I = 1, M
        YZ = DCOS( PI*(I - 0.25D0)/(N + 0.5D0) )
        DO WHILE (.TRUE.)
            YP1 = 1.0D0
            YP2 = 0.0D0
            DO J = 1, N
                YP3 = YP2
                YP2 = YP1
                YP1 = ( (2.0D0*J - 1.0D0)*YZ*YP2 - (J - 1.0D0)*YP3 )/J
            END DO ! J = 1, N
            YPP = N*(YZ*YP1 - YP2)/(YZ*YZ - 1.0D0)
            YZ1 = YZ
            YZ  = YZ1 - YP1/YPP
            IF ( DABS(YZ - YZ1) < YEPS ) EXIT ! From WHILE (.TRUE.)
        END DO	! WHILE (.TRUE.)
        X(I)     = YXM - YZ*YXL
        X(N+1-I) = YXM + YXL*YZ
        W(I)     = 2.0D0*YXL / ( (1.0D0 - YZ*YZ)*YPP*YPP )
        W(N+1-I) = W(I)
    END DO ! I = 1, M
!
END SUBROUTINE GAUSZW
!===============================================================================
! 02Oct14 - Some cosmetic changes in definition of variables according to SORD.
!
! 06Jan14 - The number PI now has 16 digits after the decimal point.
!
! 06Apr12 - Compared with GAULEG(SCIATRAN), QGAUSN(DISORT), qgausn(Pstar), 
!           legzo(MVDOM), GAUSS_LEGENDRE_QUADRATURE(RT3) for different N.
!           No significant error observed.
!
! 06Apr12 - Compared with GAUSS subroutine taken from ocean.phase.f, online at
!           http://www.giss.nasa.gov/staff/mmishchenko/brf/
!           The relative error, %, for zeros/weights and two PCs are
!            N       8           64        256        1024       10000
!           IVF 2D-14/6D-12 6D-14/2D-10 7D-14/3D-09 1D-13/1D-7 9D-13/1D-4
!           SVF 1D-14/6D-12 2D-14/2D-10 4D-14/3D-09 4D-14/1D-7 9D-13/1D-4
!           Computation time for N = 10000 is 2.5s on Intel Visual Fortran
!           11.1.067, Visual Studio 2008, Core2Duo, Windows 7, 64bit, 2.1GHz,
!           4Gb RAM and 3.2s on Silverfrost Plato (SVF) 4.4.0, Intel Core2Duo,
!           3.16GHz, 3.21Gb RAM, Windows XP, 32 bit.
!
! 06Apr12 - Chandrasekhar S. Radiative Transfer, Dover, New York, 1950, p.62,
!           Table III: N = 8, |muj| and wj are
!           |muj| = 0.1834346, 0.5255324, 0.7966665, 0.9602899;
!             wj  = 0.3626838, 0.3137066, 0.2223810, 0.1012285.
!           GAUSZW gives on output (one additional digit is shown in [.]):
!           |muj| = 0.1834346[4], 0.5255324[1], 0.7966664[8], 0.9602898[6];
!             wj  = 0.3626837[8], 0.3137066[5], 0.2223810[3], 0.1012285[4].
!===============================================================================
!
! ORIGINAL TEXT of the subroutine taken from SCIATRAN 3.1 package:
!
!      Subroutine GAULEG ( X1, X2, X, W, N )
!
!C     Taken from Numerical Recipes chapter 4
!C     Returns Gauss Legendre abscissa and weights for range [X1,X2]
!
!      INTEGER                 N
!      DOUBLE PRECISION        X(N), W(N), X1, X2
!
!C     Local variables
!
!      DOUBLE PRECISION        YEPS
!      PARAMETER               ( YEPS = 3.0D-14 )
!      INTEGER                 M, I, J
!      DOUBLE PRECISION        YXM, YXL, YZ, YP1, YP2, YP3, YPP, YZ1
!
!C     *************************************
!
!      M = (N+1)/2
!      YXM = 0.5D0 * (X2+X1)
!      YXL = 0.5D0 * (X2-X1)
!      DO 10 I = 1, M
!         YZ = DCOS( 3.141592654D0 * (I-0.25D0) / (N+0.5D0) )
! 1       CONTINUE
!         YP1 = 1.0D0
!         YP2 = 0.0D0
!         DO 20 J = 1, N
!            YP3 = YP2
!            YP2 = YP1
!            YP1 = ( (2.0D0*J-1.0D0) * YZ * YP2 - (J-1.0D0) * YP3 ) / J
! 20      CONTINUE
!         YPP = N * (YZ*YP1-YP2) / (YZ*YZ-1.0D0)
!         YZ1 = YZ
!         YZ = YZ1 - YP1 / YPP
!         IF ( DABS(YZ-YZ1) .GT. YEPS ) GOTO 1
!         X(I)     = YXM - YZ * YXL
!         X(N+1-I) = YXM + YXL * YZ
!         W(I)     = 2.0D0 * YXL / ( (1.0D0-YZ*YZ) * YPP * YPP )
!         W(N+1-I) = W(I)
! 10   CONTINUE
!
!      END
!===============================================================================