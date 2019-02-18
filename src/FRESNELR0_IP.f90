SUBROUTINE FRESNELR0_IP(REFRE, REFIM, MUI, MUR, AZI, NA, RFR0)
!===============================================================================
! PURPOSE:
!   To compute the Fresnel reflection, including the effect of rotation, for the
!   case of reflection of unpolarized beam. Refer to the COMMENTS section for
!   definition of the system of coordinates.
!
! INPUT:
!   REFRE   D(1)    Real part of the refractive index
!   REFIM   D(1)    Imaginary part of the refractive index
!   MUI     D(1)    Cosine of incident zenith angle, MUI > 0
!   MUR     D(1)    Cosine of reflected zenith angle, MUR < 0
!   AZI     D(NA)   Relative azimuths in radians. AZI=[0:2PI]
!   NA      I(1)    Number of azimuths
!
! OUTPUT:
!   RFR0  D(NA, 3)   1st column of the Fresnel reflection matrix for all AZI
!
! TREE:
!   FRESNELR0
!           |
!           +-ROTATOR2 (To compute the second rotation)
!
! COMMENTS:
!   This subroutine is a simplified version of the FRESNELR_IP subroutine.
!   Refer to FRESNELR_IP.F90 for definition of the coordinate system.
!
!   ROTATOR2 computes rotation for a single azimuth. Consider vectorization.
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NA
    REAL*8, INTENT(IN) :: REFRE, REFIM, MUI, MUR
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NA), INTENT(IN) :: AZI
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NA, 3), INTENT(OUT) :: RFR0
!
! LOCAL VARIABLES
    INTEGER &
        IA ! Loop index over azimuth
    REAL*8 &
        MUI2, & ! MUI squared
        MUR2, & ! MUR squared
        NRI2, & ! Squared module of the refractive index
        SI,   & ! Sine of the angle of incidence
        SR      ! Sine of the angle of reflection
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NA) :: &
        C2,    & ! Cosine of 2*Angle_of_rotation
        CA,    & ! CA = DSQRT(NRI2 - 1.0D0 + C2INC)
        CAZ,   & ! cos(AZI)
        CB,    & ! CB = NRI2*CINC
        CINC,  & ! Cosine of angle of incidence
        CINC2, & ! Cosine of the doubled angle of incidence
        C2INC, & ! Squared cosine of the angle of incidence
        F1,    & ! [11] = [22] element of the Fresnel reflection matrix
        F2,    & ! [12] = [21] element of the Fresnel reflection matrix
        RL,    & ! Parallel component of reflection
        RR,    & ! Perpendicular component of reflection
        RL2,   & ! RL2 = RL*RL
        RR2,   & ! RR2 = RR*RR
        S2       ! Sine of 2*Angle_of_rotation
!===============================================================================
!
    CAZ = DCOS(AZI)
!
    MUI2 = MUI*MUI
    MUR2 = MUR*MUR
    SI = DSQRT(1.0D0 - MUI2)
    SR = DSQRT(1.0D0 - MUR2)
!
!   Scattering angle = 2* Angle of incidence
!   Angle between v and -v' is needed
    CINC2 = -MUI*MUR - SI*SR*CAZ
!
!   Angle of incidence = Angle of reflection
    C2INC = 0.5D0 + 0.5D0*CINC2
    CINC = DSQRT(C2INC)
!
!   Fresnel components
    NRI2 = REFRE*REFRE + REFIM*REFIM ! Following SHARM
    CA = DSQRT(NRI2 - 1.0D0 + C2INC)
    CB = NRI2*CINC
    RL = (CB - CA)/(CB + CA)
    RR = (CINC - CA)/(CINC + CA)
!   Components of the pure Fresnel matrix
    RL2 = RL*RL
    RR2 = RR*RR
    F1 = 0.5D0*(RL2 + RR2)
    F2 = 0.5D0*(RL2 - RR2)
!
!   Compute the rotation matrix (think about vectorization)
    DO IA = 1, NA; CALL ROTATOR2(MUI, MUR, AZI(IA), S2(IA), C2(IA)); END DO
!   Apply the rotator
    RFR0(:, 1) = F1
!   Non polarizing reflection: RFR0(:, 2:3) = 0
    RFR0(:, 2) = F2*C2
    RFR0(:, 3) = F2*S2
!
END SUBROUTINE FRESNELR0_IP
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 10Sep14 - Created from FRESNELR0 for a single point {SZA, VZA, AZA}. Tested
!           against the old version for SZA = VZA = 0:1:89, AZI=0:1:360,
!           REFRE = 1.33, REFIM = 0.005. Coincidence within REAL*8 precision.
!           This version is ~1.3 times faster than the old one.
!===============================================================================