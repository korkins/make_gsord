SUBROUTINE RTLS(MSRF, KL, KV, KG, MUI, MUR, AZI, NA, BRDF)
!===============================================================================
! PURPOSE:
!   To compute the RTLS or MRTLS surface reflection for a pair of incident-
!   reflected beams, and a set of several azimuths
!
! INPUT:
!   MSRF   I(1)    RTLS (MSRF = 1) or MRTLS (MSRF=2)
!   KL     D(1)    Lambertian scattering component, KL = [0, 1]
!   KV     D(1)    Volume scattering component, KV = [-1, 1] (?)
!   KG     D(1)    Geometric-optical scattering component, KG = [-1, 1] (?)
!   MUI    D(1)    Cosine of angle of incidence, MUI > 0
!   MUR    D(1)    Cosine of angle of reflection, MUR < 0
!   AZI    D(NA)   Relative azimuth, radians. AZI = [0, 2PI] (symmetry!)
!   NA     I(1)    Number of relative azimuths
!
! OUTPUT:
!   BRDF   D(NA)   Bidirectional Reflectance Distribution Function. BRDF >= 0.
!
! TREE:
!   -
!
! COMMENTS:
!   An example of possible LSRT parameters (SHARM):
!   KL      KV      KG
!   0.086   0.014   0.017
!   0.175   0.108   0.041
!   0.227   0.093   0.056
!   0.330   0.053   0.066
!
!   The RTLS code [1, 2] was adapted from the code SHARM [3]: BCondition.cpp,
!   the LSRT subroutine. The following definition is used
!
!   mu = cos(VZA) > 0, mu0 = cos(SZA) > 0, fi - azimuth                      (1)
!
!   Note, min(mu) = min(mu0) = minimum_value [3]. BRDF is
!
!   BRDF(SZA, VZA, AZA) = KL + KG*FG(mu0, mu, fi) + KV*FV(mu, mu, fi)        (2)
!
!   FV = [(pi/2 - x)*cos(x) + sin(x)]/[mu0 + mu] - pi/4,                     (3)
!
!   or without '-pi/4' for modified MRTLS (MSRF = 2). 'Scattering' angle x is
!
!   cos(x) = mu*mu0 + smu*smu0*cos(180-AZA)                                  (4)
!
!   Note that sign at mu*mu0 here corresponds to the one used in SHARM's LSRT
!   but differs from [3, Eq.(L-3)].
!
!   where AZA in RT methods differes from AZA for definition of surface
!   reflection. Note, in Eq.(3) not-primed x, Eq.(4), is used. Standart RTLS
!
!   FG = O(mu0, mu, fi) - 1/mu' - 1/mu0' + 0.5(1+cos(x'))/mu'/mu0',          (5)
!
!   Modified RTLS (MRTLS - see SHARM -> BCondition.cpp -> LSRT)
!
!   FG = 1 - exp(-a(1+cos(x)/2)) + exp(-a(1/mu'+1/mu0'-O)),                  (6)
!
!   The O-function is
!
!   O = 1/pi*(t - sin(t)cos(t))(1/mu' + 1/mu0')                              (7)
!
!   and
!
!   cos(t) = h/b sqrt( G'^2 + [tan(SZA')tan(VZA')sin(AZA')]^2 )/...
!            ...(1/mu' + 1/mu0')                                             (8)
!
!   with constraint that |cos(t)| <= 1. Finally,
!
!   cos(x') = -cos(VZA')*cos(SZA') + sin(VZA')*sin(SZA')*cos(AZA')           (9)
!
!   G = [tan2(SZA') + tan2(VZA') - 2*tan(SZA')*tan(VZA')*cos(180-AZA')]     (10)
!
!   Note that '-' in Eq.(10) corresponds to the one used in SHARM, but differs
!   from [3, Eq.(L-4)].
!
!   The primed values are obtained using the following rule
!
!   tan(ANG') = (b/r)*tan(ANG)                                              (11)
!
!   The ratio of structural parameters is fixed following [2, P.981; 3, P.7770]
!
!   h/b = 2, b/r = 1                                                        (12)
!
!   Tables 1 and 2 from [1] give examples of input parameters. Note that f2
!   [1, Eq.(8)] corresponds to Eq.(3) [3, Eq.(A16)] exactly except for a factor
!   4/3pi. Thus, the values of k2 = [0.001...1.33] give a good idea of typical
!   magnitude of the KV parameter. NOTE: k2 may exceed 1 in some cases!
!
!   But f1 [1, Eq.(2)] differs from Eq.(5) [3, Eq.(A17)]. So, the values of
!   k1 = [0...0.0073] gives only the idea for ratio between KG/KV as compared
!   to k1/k2. From Tables 1 and 2 [1], k1/k2 << 1.
!
!   It is noted in [3] after Eq.(A19) that FV and FG take both positive and
!   negative values.
!
!   CAUTION must be exercised when this model is used at zenith angles larger
!   then 80 deg., when the BRDF may become negative or large [Ref.31 in 3].
!   Note that MINC = 0.03 (see parameter below) corresponds to 88 deg.
!
!   *** IMPORTANT: MINC depends on KV and KG and does not guarantee from
!   *** negative values in all possible cases. Probably it would be better to
!   *** replace negative values with the last positive one.
!
! REFERENCES:
!    1. Roujean J-L, et al, JGR, 1992, V97, D18, P20455
!    2. Lucht W, Schaaf, IEEE Trans Geosci Rem Sens, 2000, V38, N2, P977
!    3. Lyapustin A, Appl Opt, 2005, V44, N36, P7764
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    REAL*8, PARAMETER :: &
        ALF = 0.2D0,                 & ! Eq.(6)
        BR = 1.0D0,                  & ! b/r in Eq.(12)
        MINC = 0.03D0,               & ! Min value of cosine. Depends on KV & KG
        HB = 2.0D0,                  & ! h/b in Eq.(12)
        IP  = 3.1830988618379067D-1, & ! 1/PI
        PI2 = 1.5707963267948966D0,  & ! PI/2
        PI4 = 7.8539816339744831D-1    ! PI/4
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NA, MSRF
    REAL*8, INTENT(IN) :: MUI, MUR, KL, KV, KG
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NA), INTENT(IN) :: AZI
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NA), INTENT(OUT) :: BRDF
!
! LOCAL VARIABLES
    REAL*8 &
        CI,  & ! mu0 = cos(SZA) > 0
        CR,  & ! mu = cos(VZA) > 0
        IC,  & ! 1/CI + 1/CR
        SI,  & ! sin(SZA)
        SR,  & ! sin(VZA)
        TI,  & ! tan(SZA)
        TI2, & ! TI*TI
        TR,  & ! tan(VZA)
        TR2    ! TR*TR
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NA) :: &
        CA, & ! cos(AZI)
        CT, & ! Eq.(7)
        CX, & ! Cos of 'scattering angle'
        FG, & ! Geometrical-optics kernel function
        FV, & ! Volumetric scattering kernel function
        G2, & ! G*G, Eq.(8)
        O,  & ! Eq.(6)
        SX, & ! sin(X)
        T,  & ! acos(CT)
        X     ! 'Scattering' angle
!===============================================================================
!
!   AZ(RTLS) = 180 - AZ(RTE), Eq.(4)
    CA = -DCOS(AZI)
    CI =  MUI
    CR = -MUR
!
!   Following [3]
    IF (CI < MINC) CI = MINC
    IF (CR < MINC) CR = MINC
!
    SI = DSQRT(1.0D0 - CI*CI)
    TI = SI/CI
    SR = DSQRT(1.0D0 - CR*CR)
    TR = SR/CR
    CX = CI*CR + SI*SR*CA
    WHERE (CX > 1.0D0) CX = 1.0D0   ! SHARM, BCondition.cpp -> LSRT
    WHERE (CX < -1.0D0) CX = -1.0D0 ! SHARM, BCondition.cpp -> LSRT
    SX = DSQRT(1.0D0 - CX*CX)
    X = DACOS(CX)
!
!   Volumetric scattering kernel
    FV = ((PI2 - X)*CX + SX)/(CI + CR) - PI4
    IF (MSRF == 2) FV = FV + PI4 ! MRTLS
!
!   Primed values: Eq.(10)
    TI = BR*TI
    TI2 = TI*TI
    CI = 1.0D0/DSQRT(TI2 + 1.0D0) ! CI > 0
    SI = DSQRT(1.0D0 - CI*CI)
    TR = BR*TR
    TR2 = TR*TR
    CR = 1.0D0/DSQRT(TR2 + 1.0D0)
    SR = DSQRT(1.0D0 - CR*CR)
    CX = CI*CR + SI*SR*CA
    WHERE (CX > 1.0D0) CX = 1.0D0   ! SHARM, BCondition.cpp -> LSRT
    WHERE (CX < -1.0D0) CX = -1.0D0 ! SHARM, BCondition.cpp -> LSRT
!
!   Eq.(7)
    IC = 1.0D0/CI + 1.0D0/CR
    G2 = DABS(TI2 + TR2 - 2.0D0*TI*TR*CA) ! SHARM, BCondition.cpp -> LSRT
    CT = HB*DSQRT(G2 + TI2*TR2*(1.0D0 - CA*CA))/IC
    WHERE (CT > 1.0D0) CT = 1.0D0   ! SHARM, BCondition.cpp -> LSRT
    WHERE (CT < -1.0D0) CT = -1.0D0 ! SHARM, BCondition.cpp -> LSRT
!
!   Eq.(6)
    T = DACOS(CT)
    O = IP*(T - DSQRT(1.0D0 - CT*CT)*CT)*IC
!
!   Eq.(5)
    IF (MSRF == 1) THEN
!       RTLS
        FG = O - IC +0.5D0*(1.0D0 + CX)/CI/CR
    ELSE ! IF MSRF = 1
!       MRTLS
        FG = 1.0D0 - DEXP(-0.5D0*ALF*(1.0D0 + CX)) + DEXP(-ALF*(IC - O))
    END IF ! MSRF = 1
!
!   Eq.(2)
    BRDF = KL + KG*FG + KV*FV
!
END SUBROUTINE RTLS
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 27Oct14 - Renamed: BRF -> BRDF. Checked against the old version.
!
! 11Oct14 - Renamed: TYP -> MSRF. Checked against the old version.
!
! 08Sep14 - Created from RTLS1 and tested against it for TYP = 1 & 2,
!           [KL KV KG] = [0.539625, 0.330065, 0.045027], SZA=VZA=0:1:89,
!           AZA=0:1:360.
!           Maximum absolute difference ~E-13, relative difference ~ E-13%.
!           Absolute value of RTLS at the point of max error is ~40 (TYP=1).
!           RTLS is about 40% faster than RTLS1. Time(RTLS)~0.7s.
!           Both  were tested.
!           Note that RTLS was called with MUI = cos(SZA) & MUR = -cos(VZA).
!===============================================================================