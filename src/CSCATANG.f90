SUBROUTINE CSCATANG(MU0, MU, NM, AZI, NA, MUS)
!===============================================================================
! PURPOSE:
!   To compute cosine of scattering angle for a given solar zenith, view zeniths
!   and relative azimuths
!
! INPUT:
!   MU0   D(1)    Cosine of Solar zenith angle, in RT computations mu0 > 0
!   MU    D(NM)   Cosines of zenith angles, mu > 0 for transmittance
!   NM    I(1)    Number of zeniths
!   AZI   D(NA)   Relative azimuths in radians
!   NA    I(1)    Number of azimuths
!
! OUTPUT:
!   MUS   D(NM, NA)   Cos of scattering angle
!
! TREE:
!   -
!
! COMMENTS:
!   -
!
! REFERENCES:
!   -
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NM, NA
    REAL*8, INTENT(IN) :: MU0
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NM), INTENT(IN) :: MU
    REAL*8, DIMENSION(NA), INTENT(IN) :: AZI
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NM, NA), INTENT(OUT) :: MUS
!
! LOCAL VARIABLES
    INTEGER &
        IA ! Loop index over azimuth
    REAL*8 &
        SMU0 ! sqrt(1-mu0*mu0)
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NA) :: &
        CAZ ! Cosine of azimuth
    REAL*8, DIMENSION(NM) :: &
        MUMU0, & ! mu*mu0
        SMUMU0   ! sqrt(1-mu*mu)*sqrt(1-mu0*mu0)
!===============================================================================
!
    MUMU0 = MU*MU0
    SMU0 = DSQRT(1.0D0 - MU0*MU0)
    SMUMU0 = DSQRT(1.0D0 - MU*MU)*SMU0
    CAZ = DCOS(AZI)
    DO IA = 1, NA; MUS(:, IA) = MUMU0 + SMUMU0*CAZ(IA); END DO
!
END SUBROUTINE CSCATANG
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 08Sep14 - Created and tested as part of SORD
!===============================================================================