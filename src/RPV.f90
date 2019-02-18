SUBROUTINE RPV(MSRF, R0, RC, K, ASC, MUI, MUR, AZI, NA, BRDF)
!===============================================================================
! PURPOSE:
!   To compute the RPV or MRPV Bidirectional Reflectance Factor (BRF).
!
! INPUT:
!   MSRF   I(1)    RPV (MODE = 1) or MRPV (MODE = 2)
!   R0     D(1)    Overall reflectance, R0 = [0, 1]
!   RC     D(1)    Hot-spot parameter, RC = [0, 1]
!   K      D(1)    Minnaert function parameter, K > 0. See examples below.
!   ASC    D(1)    Henyey-Greenstein parameter (average scattering cosine),
!                  ASC = (-1, 1). See examples below.
!   MUI    D(1)    Cosine of angle of incidence, MUI > 0
!   MUR    D(1)    Cosine of angle of reflection, MUR < 0
!   AZI    D(NA)   Relative azimuth, radians. AZI = [0, 2PI] (symmetry!)
!   NA     I(1)    Number of relative azimuths
!
! OUTPUT:
!   BRDF   D(NA)   Bidirectional Reflectance Distribution Function. BRDF >= 0
!
! TREE:
!   -
!
! COMMENTS:
!   RPV [1] and MRPV [2] were adapted from the code SHARM [3]: BCondition.cpp,
!   the RPV subroutine. The following definition is used
!
!   mu = cos(VZA) > 0, mu0 = cos(SZA) > 0                                    (1)
!
!   Note, min(mu) = min(mu0) = minimum_value [3]
!
!   BRF is
!
!   BRF(SZA, VZA, AZA) = ro*M(k)*F(a)*H(rc)                                  (2)
!
!   where modified Minnaert function is
!
!   M(k) = [mu*mu0*(mu + mu0)]^(k-1)                                         (3)
!
!   Henyey-Greenstein function is (RPV MODE = 1)
!
!   F(a) = (1 - a^2)/(1 - 2*a*cos(x) + a^2)^(2/3)                            (4)
!
!   or for MRPV (MODE = 2)
!
!   F(a) = exp(a * cos(x)),                                                  (5)
!
!   cos(x) = -mu*mu' + smu*smu'*cos(fi - fi')                                (6)
!
!   and
!
!   H(rc) = (1 + (1-rc)/(1+G))                                               (7)
!
!   Finally,
!
!   G = [tan2(SZA) + tan2(VZA) + 2*tan(SZA)*tan(VZA)*cos(AZA)]               (8)
!
!   Unlike in [4, 5], hot spot is in the direction of backscattering AZA = 2PI.
!   Different coordinate system was used in [4, 5] (note the signs in  Eqs.
!   (6), (8)).
!
!   Lambert law is fulfilled numerically by setting ASC = 0, K = 1, RC = 1 [5].
!
!   3-parameters version of RPV, RC = R0 [5], is often used (e.g. in SHARM [3]).
!
!   SHARM manual gives the following examples of the input parameters:
!
!                            Visible Band              Near-IR
!                            R0=RC   ASC    K          R0=RC   ASC    K
!   1  Spruce                0.008  -0.308  0.554      0.050  -0.201  0.581
!   2  Sparse erectophile    0.064  -0.001  1.207      0.278  -0.006  0.725
!   3  Tropical forest       0.012  -0.169  0.651      0.303  -0.034  0.729
!   4  Plowed field          0.072  -0.257  0.668      0.077  -0.252  0.678
!   5  Grasses               0.014  -0.169  0.810      0.242  -0.032  0.637
!   6  Broad leaf crops      0.012  -0.281  0.742      0.204  -0.089  0.658
!   7  Savannah              0.010  -0.287  0.463      0.219  -0.050  0.673
!   8  Leaf forest           0.022  -0.228  0.633      0.285  -0.060  0.745
!   9  Conifers              0.018  -0.282  0.364      0.235  -0.095  0.758
!   10 Hardwood forest       0.028  -0.175  0.768      0.066  -0.141  0.735
!      winter
!   11 Loam soil             0.147  -0.096  0.839      0.195  -0.097  0.850
!   12 Irrigated wheat       0.027  -0.078  0.382      0.306  -0.008  0.606
!
! REFERENCES:
!    1. Rahman H et al, JGR, 1993, V98, ND11, P20791
!    2. Martonchik JV et al, IEEE Trans Geosci Rem Sens, 1998, V36, P1266
!    3. Lyapustin A, Appl Opt, 2005, V44, N36, P7764
!    4. Soleheim I et al, Rem Sens Env, 2000, V73, P78
!    5. Lavergne T et al, Rem Sens Env, 2007, V107, P362
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    REAL*8, PARAMETER :: &
        MINC = 0.03D0   ! Min value for cosines
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NA, MSRF
    REAL*8, INTENT(IN) :: MUI, MUR, R0, RC, K, ASC
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NA), INTENT(IN) :: AZI
!
! OUTPUT VARIABLES
    REAL*8, DIMENSION(NA), INTENT(OUT) :: BRDF
!
! LOCAL VARIABLES
    REAL*8 &
        ASC2, & ! ASC*ASC
        CI,   & ! mu0 > 0
        CR,   & ! mu > 0
        M,    & ! Eq.(3)
        SI,   & ! sin(SZA)
        SR,   & ! sin(VZA)
        TI,   & ! tan(SZA)
        TR      ! tan(VZA)
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NA) :: &
        CA, & ! cos(AZI)
        CX, & ! Eq.(6)
        F,  & ! Eq.(4) or Eq.(5)
        G,  & ! Eq.(8)
        H     ! Eq.(7)
!===============================================================================
!
    CA = DCOS(AZI)
    CI =  MUI
    CR = -MUR
!
!   Following [3]
    IF (CI < MINC) CI = MINC
    IF (CR < MINC) CR = MINC
!
    SI = DSQRT(1.0D0 - CI*CI)
    SR = DSQRT(1.0D0 - CR*CR)
    TI = SI/CI
    TR = SR/CR
    CX = -CI*CR + SI*SR*CA
!
    M = (CI*CR*(CI + CR))**(K-1.0D0)
!
    IF (MSRF == 1) THEN
!       RPV
        ASC2 = ASC*ASC
        F = (1.0D0 - ASC2)/(DSQRT(1.0D0 - 2.0D0*ASC*CX + ASC2))**3
    ELSEIF (MSRF == 2) THEN
!       MRPV
        F = DEXP(ASC*CX)
    END IF ! MODE == 1
!
    G = DSQRT(TI*TI + TR*TR + 2.0D0*TI*TR*CA) ! SHARM uses ABS here
    H = 1.0D0 + (1.0D0 - RC)/(1.0D0 + G)
!
    BRDF = R0*M*F*H
!
END SUBROUTINE RPV
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 16Nov14 - Created on the basis of RPV1 and tested against it for R0, K, ASC =
!           0.219 0.463 -0.287, respectively, SZA = VZA = 0:1:75, AZI = 0:1:360.
!           RC=R0. Ok.
!===============================================================================