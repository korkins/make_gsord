SUBROUTINE ROTATOR2(MUI, MUS, AZI, S2, C2)
!===============================================================================
! PURPOSE:
!   To compute SIN and COS of the doubled rotation angles for inverse rotation
!   in VRTE.
!
! INPUT:
!   MUI   D(1)   Cosine of 180-Zenith_of_incidence, MUI /= 0
!   MUS   D(1)   Cosine of 180-Zenith_of_scattering, MUS /= 0
!   AZI   D(1)   Relative azimuth (fi - fi'). AZI = [0:2PI] radians
!
! OUTPUT:
!   S2   D(1)   sin(2*X2)
!   C2   D(1)   cos(2*X2)
!
! TREE:
!   -
!
! COMMENTS:
!    For the vector raditive transfer applications: positive MUI and MUS are
!    defined in the direction of increasing of optical depth ("down"). Refer to
!    [1, p.70, Eq.(3.18)]: MUI = U = -cos(Zenith_of_incidence), and 
!    MUR = U' = -cos(Zenith_of_scattering). Zenith = 0 corresponds to "up",
!    Zenith = 180 corresponds to "down".
!
!    This subroutine is a simplified version of ROTATOR.f. It computes only the
!    second rotator, R(X2), from the following equation
!
!    B = R(X2) A R(X1)                                                       (1)
!
!    The coordinate system is defined in [1, p.66]. 'Theta' means zenith angle.
!
!    As defined in [1, p.11, Eq.(1.51)], the rotation matrix through an angle
!    X > 0
!
!        |1   0         0         0|   |1  0   0   0|
!        |0   cos(2X)   sin(2X)   0|   |0  C2  S2  0|
!    R = |0  -sin(2X)   cos(2X)   0| = |0 -S2  C2  0|                        (2)
!        |0   0         0         1|   |0  0   0   1|
!
!    In this subroutine, not the matrix R but the elements S2, C2 are computed.
!    The '-' in (1) is must be considered in the calling program.
!
!    Formulas are given in [1, pp.69-70, Eqs.(3.10)-(3.23)]. For the cases when
!    AZI = 0, PI, 2PI no rotation is performed. Two special cases must be
!    considered using L'Hopital rule when [1, Eqs.(3.15)-(3.16)] can not be used
!    [2].
!
!    Scattering in the horizontal directions is not considered. Refer to [1, 2]
!    for details.
!
!    a) For the case AZ = 0, PI, 2PI no rotation needed.
!    b) For the case when incident AND(!) scattered rays travel along the
!    zenith/nadir direction (0, PI) no rotation is needed
!    c) For the case when only incident ray travels along the zenith/nadir, the
!    SECOND rotator, R2, is not needed (scattering plane contains the nadir/
!    zenith direction and the direction of scattering).
!    d) For the case when only scattered ray travels along the zenith/nadir
!    direction, the FIRST rotator, R1, is not needed (scattering plane contains
!    the nadir/zenith and the direction of incidence)
!
!    All the cases (a)-(d) are included in this subroutine.
!
! REFERENCES:
!    1. Hovenier JW, van der Mee C, Domke H. Transfer of polarized light in
!       planetary atmospheres. Basic concepts and practical methods. Dordrecht:
!       Kluwer Academic Publishers, 2004.
!    2. Rozanov VV, Bremen University. Private communication.
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    REAL*8, PARAMETER :: &
        PI   = 3.1415926535897932D0, & ! The number PI
        PI2  = 6.2831853071795865D0, & ! The number 2PI
        TINY = 1.0D-8                  ! To compare doubles
!
! INPUT VARIABLES
    REAL*8, INTENT(IN) :: MUI, MUS, AZI
!
! OUTPUT VARIABLES
    REAL*8, INTENT(OUT) :: S2, C2
!
! LOCAL VARIABLES
    REAL*8 &
        CAZ,  & ! cos(AZI)
        CSCA, & ! Cosine of the scattering angle
        CX2,  & ! cos(X2)
        C2X2, & ! cos(X2)*cos(X2)
        SAZ,  & ! -sin(AZI), NOTE: '-'
        SNI,  & ! Sine of the zenith of incidence
        SNS,  & ! Sine of the zenith of scattering
        SSCA, & ! Sine of the scattering angle
        SX2     ! sin(X2)
!===============================================================================
!
!   No rotation by default and in particular cases
    C2 = 1.0D0
    S2 = 0.0D0 
    IF (1.0D0 - DABS(MUI) < TINY) RETURN
!
!   At least one rotation is needed in the following situations
    IF (AZI > TINY .AND. DABS(AZI - PI) > TINY &
                                                .AND. PI2 - AZI > TINY) THEN
        CAZ =  DCOS(AZI)
        SAZ = -DSIN(AZI)      ! for [1, p.70, Eq.(3.23)]: -sin(fi-fi')       (*)
        IF (1.0D0 - DABS(MUS) < TINY) THEN
            IF (1.0D0 + MUS < TINY) CX2 = -CAZ   ! '-' [2], -cos(Theta) = -1
            IF (1.0D0 - MUS < TINY) CX2 =  CAZ   ! '-' [2], -cos(Theta) = +1
            SX2 = SAZ                            ! [1, p.69, Eq.(3.13)] and  (*)
        ELSE
            SNI = DSQRT(1.0D0 - MUI*MUI)
            SNS = DSQRT(1.0D0 - MUS*MUS)
            CSCA = MUI*MUS + SNI*SNS*CAZ
            SSCA = DSQRT(1.0D0 - CSCA*CSCA)
            CX2 = (MUS*CSCA - MUI)/SNS/SSCA    ! [1, p.70, Eq.(3.21)]
            SX2 = SNI*SAZ/SSCA ! [1, p.70, Eq.(3.23)] and (*)
        END IF ! 1.0D0 - DABS(MUS) < TINY
 !
        C2X2 = CX2*CX2
        C2 = 2.0D0*C2X2 - 1.0D0
        S2 = 2.0D0*SX2*CX2
    END IF ! AZI > TINY .OR. DABS(AZI-PI) > TINY .OR. PI2 - AZI > TINY
!
END SUBROUTINE ROTATOR2
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 05Sep14 - New format of declaration of variables. Tested as part of SORD
!
! 25Aug14 - Cosines of angles, istead of angles, are now on input. Major changes
!           in the code. Tested against previous version of ROTATOR2 for
!           SZA = VZA = 0:PI, excluding PI/2, and AZI = 0:2PI. Step 1 degree.
!
! 06Jan14 - PI =  3.1415 9265 3589 7932D0, PI2 = 6.2831 8530 7179 5865D0
!
! 04Jan14 - THETAI, THETAS are noe in radians. Parameter D2R was removed.
!
! 26Dec13 - Unused variables CX1, C2X1, and SX1 were removed
!           AZA is now in radians, AZR is removed.
!           Tested vs the same scenario (4 solar, 5 azimuth angles) as a part
!           of RT code (atmosphere + ocean surface) before AZA was redefined.
!
! 27Aug13 - Tested as a part of the IQUVDOM: Rayleigh atmosphere over ocean.
!           SZA = 60, AZA = 45, UDA = [0:10:80]. Exact agreement (txt, 8 digits)
!
! 26Aug13 - Tested against ROTATOR.f for SZA = [0:1:80], and SZA = 180 - SZA,
!           AZA=[0:1:180], VZA = [0:1:180] excluding VZA = 90. Exact agreement.
!===============================================================================