SUBROUTINE NODES(NM, ISRF, UDA, AZA, NUD, NUP, NDN, NG1, NG2, NMU, NAZ, NGA, &
                 MU, MUD, ZMU, WMU, AZI, ZAZ, WCMA, WSMA)
!===============================================================================
! PURPOSE:
!   To create arrays with Gaussian and user defined zenith nodes and Gaussian
!   azimuth nodes for numerical integration over corresponding angles. This
!   subroutine also converts view azimuth from degrees to radians.
!
! INPUT:
!   NM     I(1)     Total number of the Fourier terms: m = 0:NM-1
!   ISRF   I(1)     Surface index
!   UDA    D(NUD)   User defined zeniths; TOA = [0, 90), if any, comes first
!   AZA    D(NAZ)   View azimuths, degrees
!   NUD    I(1)     Number of user defined zeniths, NUD = NUP+NDN > 0
!   NUP    I(1)     Number of user defined ascending directions
!   NDN    I(1)     Number of user defined descending directions
!   NG1    I(1)     Number of zenith Gauss nodes per hemisphere
!   NG2    I(1)     Number of zenith Gauss nodes per sphere, NG2 = 2*NG1
!   NMU    I(1)     Number of Gauss and used defined nodes together
!   NAZ    I(1)     Number of view azimuths
!   NGA    I(1)     Number of azimuth Gauss nodes
!
! OUTPUT:
!   MU    D(NMU)       Array of all zenith nodes; negative nodes come first
!   MUD   D(NUD)       cos(UDA); negative values, if any, come first
!   ZMU   D(NG2)       Zenith Gauss nodes, ZMU(i) < ZMU(i+1)
!   WMU   D(NG2)       Corresponding Gauss weights
!   AZI   D(NZA)       View azimuth, radians
!   ZAZ   D(NGA)       Azimuth Gauss nodes, in radians, to compute the m-moments
!   WCMA  D(NGA, NM)   Weighted cos(m*ZAZ) for integration over azimuth
!   WSMA  D(NGA, NM)   Similar to WCMA except for the sin-function; ISRF > 5
!
! TREE:
!   NODES
!       |
!       +-GAUSZW (Compute Gauss weights and nodes)
!
! COMMENTS:
!   Either NUP or NDN may be zero, but not both at the same time.
!   UDA = [0, 90) - TOA, UDA = (90, 180] - BOA. UDA = 90, horizon, is excluded.
!
! REFERENCES:
!   -
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    REAL*8, PARAMETER :: PI  = 3.1415926535897932D0   ! The number pi
    REAL*8, PARAMETER :: D2R = 1.74532925199432958D-2 ! pi/180
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NM, ISRF, NUD, NUP, NDN, NG1, NG2, NMU, NAZ, NGA
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NUD), INTENT(IN) :: UDA
    REAL*8, DIMENSION(NAZ), INTENT(IN) :: AZA
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NMU), INTENT(OUT) :: MU
    REAL*8, DIMENSION(NUD), INTENT(OUT) :: MUD
    REAL*8, DIMENSION(NG2), INTENT(OUT) :: ZMU, WMU
    REAL*8, DIMENSION(NAZ), INTENT(OUT) :: AZI
    REAL*8, DIMENSION(NGA), INTENT(OUT) :: ZAZ
    REAL*8, DIMENSION(NGA, NM), INTENT(OUT) :: WCMA, WSMA
!
! LOCAL VARIABLES
    INTEGER &
        IM ! Loop index over azimuth orders, IM = 1, 2, ... NM
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NG1) :: &
        WN, & ! Gauss weights corresponding to negative zenith nodes
        WP, & ! Gauss weights corresponding to positive zenith nodes
        ZN, & ! Negative (ascending) zenith Gauss nodes
        ZP    ! Positive (descending) zenith Gauss nodes
    REAL*8, DIMENSION(NGA) :: &
        WAZ   ! Azimuth weights
!===============================================================================
!
    AZI =  AZA*D2R
    MUD = -DCOS(UDA*D2R)
!
!   Azimuth nodes and weights: Gauss only
    CALL GAUSZW(0.0D0, 1.0D0, NGA, ZAZ, WAZ)
    ZAZ = PI*ZAZ
!   Precompute arrays for azimuth expansion of surface
    IF (ISRF > 1) THEN
!       cos-terms are used for scalar and vector surfaces
        WCMA(:, 1) = WAZ
        IF (NM > 1) WCMA(:, 2) = WAZ*DCOS(ZAZ)
        IF (NM > 2) THEN
            DO IM = 3, NM; WCMA(:, IM) = WAZ*DCOS((IM-1)*ZAZ); END DO
        END IF ! NM > 2
!       sin-terms are used for vector surfaces only
        IF (ISRF > 5) THEN ! Vector surfaces
            WSMA(:, 1) = 0.0D0
            IF (NM > 1) WSMA(:, 2) = WAZ*DSIN(ZAZ)
            IF (NM > 2) THEN
                DO IM = 3, NM; WSMA(:, IM) = WAZ*DSIN((IM-1)*ZAZ); END DO
            END IF ! NM > 2
        END IF ! ISRF > 5
    END IF ! ISRF > 1
!
!   Zenith nodes and weights: Gauss and user defined, if any
    CALL GAUSZW(0.0D0, 1.0D0, NG1, ZP, WP)
    ZN = -ZP
    WN =  WP
    ZMU(:NG1) = ZN
    WMU(:NG1) = WN
    ZMU(NG1+1:NG2) = ZP
    WMU(NG1+1:NG2) = WP
!
!***********************************************
!   Single Gauss quad
!    CALL GAUSZW(-1.0D0, 1.0D0, NG2, ZMU, WMU)
!    ZN = ZMU(1:NG1)
!    WN = WMU(1:NG1)
!    ZP = ZMU(NG1+1:NG2)
!    WP = WMU(NG1+1:NG2)
!***********************************************
!
!   Complete set of zenith nodes: both Gaussian and user defined
    IF (NUP > 0) MU(:NUP) = MUD(:NUP)
    MU(NUP+1:NUP+NG1) = ZN
    MU(NUP+NG1+1:NUP+NG2) = ZP
    IF (NDN > 0) MU(NUP+NG2+1:) = MUD(NUP+1:)
!
END SUBROUTINE NODES
!===============================================================================
! 01May16 - NGX is replaced with excplicit NG1+1 & NUG with NUP+NG1
!
! 30Apr16 - DCMA & DSMA are removed
!
! 22Apr16 - Minor changes in comments
!
! 09Dec15 - Error fixed:
!               IF (ISRF > 5) THEN ! Vector surfaces
!                   WSMA(:, 1) = 0.0D0
!                   WSMA(:, 2) = WAZ*DSIN(ZAZ) <<<<<<<<<<<<<<<<< Before
!
!               For ISRF > 5 & NM = 1 WSMA(:, 2) is undefined
!
!               IF (ISRF > 5) THEN ! Vector surfaces
!                   WSMA(:, 1) = 0.0D0
!                   IF (NM > 1) WSMA(:, 2) = WAZ*DSIN(ZAZ) <<<<< Now
!
! 24Sep15 - MUWG & MU0 are removed from output; SZA & NSZ are removed from input
!
! 28Apr15 - Input/output is modified to provide mu0, AZI, MUD from input data.
!           Tested for Rayleigh atmosphere over RTLS surface
!
! 24Apr15 - NODES now precomputes cos&sin terms for Gauss and view azimuths &
!           all Fourier orders, excluding m = 0 (the very 1st one)
!
! 17Apr15 - Negative zeros & weights are not sorted in ascending order. Tested
!           for Rayleigh atmosphere over black surface
!
! 14Apr15 - Azimuth nodes are redefined: [0:2pi] -> [0:pi]. Azimuthal
!           symmetry (scalar case) is used. Tested for RTLS. Not sure about
!           polarized surfaces
!
! 24Sep14 - Tested as part of SORD_IP for NT ocean model (see SORD_IP.F90).
!           Perfect agreement. Also tested for the black surface ISRF = 0,
!           Lambert ISRF = 1, and RTLS ISRF = 4. Ok
!===============================================================================