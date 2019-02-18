SUBROUTINE ROTATOR(MUI, MUS, AZI, S1, C1, S2, C2)
!===============================================================================
! PURPOSE:
!   To compute SIN and COS of the doubled rotation angles for direct and inverse
!   rotations in VRTE
!
! INPUT:
!   MUI   D(1)   Cosine of 180-Zenith_of_incidence, MUI /= 0. Downward: MUI > 0
!   MUS   D(1)   Cosine of 180-Zenith_of_scattering, MUS /= 0. Upward: MUS < 0
!   AZI   D(1)   Relative azimuth (fi - fi'). AZI = [0:2PI] radians
!
! OUTPUT:
!   S1   D(1)   sin(2*X1)
!   C1   D(1)   cos(2*X1)
!   S2   D(1)   sin(2*X2)
!   C2   D(1)   cos(2*X2)
!
! TREE:
!   -
!
! COMMENTS:
!    Positive MU are measured opposed to positive Z-direction, i.e. MU > 0(down)
!    corresponds to zenith angle Theta = (pi/2:pi].
!
!    The coordinate system is defined in [1, p.66]. 'Theta' means zenith angle.
!
!    As defined in [1, p.11, Eq.(1.51)], the rotation matrix through an angle
!    X > 0
!
!        |1   0         0         0|   |1  0  0  0|
!        |0   cos(2X)   sin(2X)   0|   |0  C  S  0|
!    R = |0  -sin(2X)   cos(2X)   0| = |0 -S  C  0|                          (1)
!        |0   0         0         1|   |0  0  0  1|
!
!    In this subroutine, not the matrix R but the elements S, C are computed.
!    The '-' in (1) is must be considered in the calling program. S1, C1 and S2,
!    C2 are defined for X1 and X2, respectively. Both rotations are considered
!
!    B = R(X2) A R(X1)                                                       (2)
!
!    yet the first one, R(X1), is not necessary, for instance, in case of single
!    scattering of natural light.
!
!    Formulas are given in [1, pp.69-70, Eqs.(3.10)-(3.23)]. For the cases when
!    AZA = 0, PI, 2PI no rotation is performed. Two special cases must be
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
    REAL*8, PARAMETER :: PI =  3.1415926535897932D0 ! The number PI
    REAL*8, PARAMETER :: PI2 = 6.2831853071795865D0 ! The number 2PI
    REAL*8, PARAMETER :: TINY = 1.0D-8              ! A small angle, rad.
!
! INPUT VARIABLES
    REAL*8, INTENT(IN) :: MUI, MUS, AZI
!
! OUTPUT VARIABLES
    REAL*8, INTENT(OUT) :: S1, C1, S2, C2
!
! LOCAL VARIABLES
    REAL*8 &
        CAZ,  & ! cos(AZA)
        CSCA, & ! Cosine of the scattering angle
        CX1,  & ! cos(X1)
        CX2,  & ! cos(X2)
        C2X1, & ! cos(X1)*cos(X1)
        C2X2, & ! cos(X2)*cos(X2)
        SAZ,  & ! -sin(AZA), NOTE: '-'
        SNI,  & ! sin(THETAI)
        SNS,  & ! sin(THETAS)
        SSCA, & ! Sine of the scattering angle
        SX1,  & ! sin(X1)
        SX2     ! sin(X2)
!===============================================================================
!
!   No rotation by default
    C1 = 1.0D0
    S1 = 0.0D0
    C2 = 1.0D0
    S2 = 0.0D0
!
!   Particular cases: no rotation, R1 = R2 = diag[1 1 1 1]
!   Zenith_i = 0 (mu_i = -1) & Zenith_s = 0 (mu_s = -1)
    IF (1.0D0 + MUI < TINY .AND. 1.0D0 + MUS < TINY) RETURN
!   Zenith_i = 0 (mu_i = -1) & Zenith_s = pi (mu_s = +1)
    IF (1.0D0 + MUI < TINY .AND. 1.0D0 - MUS < TINY) RETURN
!   Zenith_i = pi (mu_i = +1) & Zenith_s = pi (mu_s = +1)
    IF (1.0D0 - MUI < TINY .AND. 1.0D0 - MUS < TINY) RETURN
!   Zenith_i = pi (mu_i = +1) & Zenith_s = 0 (mu_s = -1)
    IF (1.0D0 - MUI < TINY .AND. 1.0D0 + MUS < TINY) RETURN
!
!   At least one rotation is needed in the following situations
    IF (AZI > TINY .AND. DABS(AZI - PI) > TINY &
                                                .AND. PI2 - AZI > TINY) THEN
        CAZ =  DCOS(AZI)
        SAZ = -DSIN(AZI)      ! for [1, p.70, Eq.(3.23)]: -sin(fi-fi')       (*)
        IF (1.0D0 - DABS(MUI) < TINY) THEN
            IF (1.0D0 + MUI < TINY) CX1 = -CAZ       ! '-' [2], cos(Theta') = +1
            IF (1.0D0 - MUI < TINY) CX1 =  CAZ       ! '+' [2], cos(Theta') = -1
            SX1 =  SAZ                            ! [1, p.69, Eq.(3.13)] and (*)
            CX2 =  1.0D0
            SX2 =  0.0D0
        ELSEIF (1.0D0 - DABS(MUS) < TINY) THEN
            CX1 =  1.0D0
            SX1 =  0.0D0
            IF (1.0D0 + MUS < TINY) CX2 = -CAZ        ! '-' [2], cos(Theta) = +1
            IF (1.0D0 - MUS < TINY) CX2 =  CAZ        ! '-' [2], cos(Theta) = -1
            SX2 =  SAZ                            ! [1, p.69, Eq.(3.13)] and (*)
        ELSE
            SNI = DSQRT(1.0D0 - MUI*MUI)
            SNS = DSQRT(1.0D0 - MUS*MUS)
            CSCA = MUI*MUS + SNI*SNS*CAZ
            SSCA = DSQRT(1.0D0 - CSCA*CSCA)
            CX1 = (MUI*CSCA - MUS)/SNI/SSCA               ! [1, p.70, Eq.(3.20)]
            CX2 = (MUS*CSCA - MUI)/SNS/SSCA               ! [1, p.70, Eq.(3.21)]
            SX1 = SNS*SAZ/SSCA                    ! [1, p.70, Eq.(3.23)] and (*)
            SX2 = SNI*SAZ/SSCA                    ! [1, p.70, Eq.(3.23)] and (*)
        END IF ! 1.0D0 - DABS(MUI)) < TINY
 !
        C2X1 = CX1*CX1
        C2X2 = CX2*CX2
        C1 = 2.0D0*C2X1 - 1.0D0
        C2 = 2.0D0*C2X2 - 1.0D0
        S1 = 2.0D0*SX1*CX1
        S2 = 2.0D0*SX2*CX2
    END IF ! AZI > TINY .OR. DABS(AZI-PI) > TINY .OR. PI2 - AZI > TINY
!
END SUBROUTINE ROTATOR
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 11Sep14 - Cosines of zenith angles instead of angles are now used. New format
!           of declaration of variables. Tested against previous version of
!           ROTATOR for SZA=VZA=0:1:180 (excluding 90), AZI=0:1:360. Major
!           changes in the code. Exact agreement with the old version.
!           New version is faster: time_old/time_new = 1.3.
!
! 06Jan14 - PI =  3.1415 9265 3589 7932D0, PI2 = 6.2831 8530 7179 5865D0
!
! 04Jan14 - THETAI, THETAS are now in radians. Parameter D2R was removed.
!
! 26Dec13 - AZA is now in radians. AZR has been removed.
!           Tested vs the same scenario (4 solar, 5 azimuth angles) as a part
!           of RT code (atmosphere + ocean surface) before AZA was redefined.
!
! 07Dec12 - Used to compute the Single Scattering approximation and tested
!           against SCIATRAN. OK.
!
! 09Nov12 - The subroutine is used in FRESNELR subroutine that was tested
!           against SCIATRAN [2] for SZA = 45, VZA = 0:10:80, AZA = 0:45:180. OK
!           Additional test for SZA = VZA = 0, AZ = 0, 45, VZA = SZA = 60 and
!           AZ = 135, 180. OK.
!===============================================================================