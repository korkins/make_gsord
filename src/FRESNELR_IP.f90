SUBROUTINE FRESNELR_IP(REFRE, REFIM, MUI, MUR, AZI, NA, RFR)
!===============================================================================
! PURPOSE:
!   To compute the Fresnel reflection, including rotation. Refer to the COMMENTS
!   section for definition of the system of coordinates. Ellipticity is
!   neglected: size(Fresnel matrix) = 3-by-3.
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
!   RFR   D(NA, 9)   The Fresnel reflection matrix for all azimuths
!
! TREE:
!   FRESNELR
!          |
!          +-ROTATOR (To compute the elements C1, S1, C2, S2 of a rotator)
!
! COMMENTS:
!   The effect of roughness of the ocean surface and mutual shadowing are
!   NOT included in this subroutine. The ocean is assumed to be black (not
!   emitting) surface Horizontal traces are not included (refer to [1,
!   Eqs.(4.25)-(4.26)]) if needed.
!
!   The book by Hovenier et al.[1] is followed here for definition of the
!   geometry, the Stokes vector S = [I, Q, U, V], the direction of positive
!   rotation and so on. In particularly, Z-axis is pointed upwards (to the sky),
!   VZA = 0 is the zenith (looking UP to the sky IN the positive Z direction).
!   The relative azimuth is measured clockwise looking IN the positive Z
!   direction [1, p.64].
!
!   The following relations are used in this subroutine:
!   Fresnel matrix [2, Eq.(A.1)]:
!
!      |rp  rm  0   0 |
!      |rm  rp  0   0 |
!   F =|0   0   ro  0 |,                                                     (1)
!      |0   0   0   ro|
!
!    where
!
!    rp = (rl*rl + rr*rr); rm = (rl*rl - rr*rr); ro = Real(rl*rr) = rl*rr,
!    rl = (n2*cos(i) - S)/(n2*cos(i) + S); rr = (cos(i) - S)/(cos(i) + S)
!    S = sqrt(n2 - 1 + cos2(i)), cos2(i) = cos(i)*cos(i), n2 = n*n,
!
!    'n' is the index of refraction, and 'i' is the angle of incidence.
!    The rotation is described as following (note the signs!) [1, Eq.(3.7)]
!
!    r = R(pi - s2)*F*R(-s1),                                                (2)
!
!    where the rotation matrix through an angle 's' > 0 is [1, Eq.(1.51)]
!
!        |1   0         0         0|
!        |0   cos(2s)   sin(2s)   0|
!    F = |0  -sin(2s)   cos(2s)   0|.                                        (3)
!        |0   0         0         1|
!
!    Using the following notation (these elements are computed by ROTANG12)
!
!    C1 = cos(2s1), S1 = sin(2s1),                                           (4)
!    C2 = cos(2s2), S2 = sin(2s2),                                           (5)
!
!    the elements of 'r' are defined as follows [1, Eq.(3.9)] (Maple)
!
!        |rp     rm*C1                -rm*S1                 0 |
!        |rm*C2  rp*C1*C2 - ro*S1*S2  -rp*S1*C2 - ro*C1*S2   0 |
!    r = |rm*S2  rp*C1*S2 + ro*S1*C2  -rp*S1*S2 + ro*C1*C2   0 |,            (6)
!        |0      0                     0                     ro|
!
!    Given the geometry defined above the angle of incidence (the angle between
!    v and -v') is calculated as follows:
!
!    cos(2i) = -mu*mu' - smu*smu'*cos(fi - fi')                              (7)
!
!    Note that the angle of reflection, 2i, is the angle complementary to the
!    angle of scattering, S, defined as usual [1, Eq.(3.14)]
!
!    cos(S) = mu*mu' + smu*smu'*cos(fi - fi').                               (8)
!
! REFERENCES:
!    1. Hovenier JW, van der Mee C, Domke H. Transfer of polarized light in
!       planetary atmospheres. Basic concepts and practical methods. Dordrecht:
!       Kluwer Academic Publishers, 2004.
!    2. Y. Ota et al. JQSRT (2010), V.111, p.878-894
!    3. ocean.phase.f from http://www.giss.nasa.gov/staff/mmishchenko/brf/
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
    REAL*8, DIMENSION(NA, 9), INTENT(OUT) :: RFR
!
! LOCAL SCALARS
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
        C1,    & ! Cosine of doubled first rotation angle
        C2,    & ! Cosine of doubled second rotation angle
        CA,    & ! CA = DSQRT(NRI2 - 1.0D0 + C2INC)
        CAZ,   & ! cos(AZI)
        CB,    & ! CB = NRI2*CINC
        CINC,  & ! Cosine of angle of incidence
        CINC2, & ! Cosine of the doubled angle of incidence
        C2INC, & ! Squared cosine of the angle of incidence
        F1,    & ! [11] = [22] element of the Fresnel reflection matrix
        F2,    & ! [12] = [21] element of the Fresnel reflection matrix
        F3,    & ! [33] = [44] element of the Fresnel reflection matrix
        RL,    & ! Parallel component of reflection
        RR,    & ! Perpendicular component of reflection
        RL2,   & ! RL2 = RL*RL
        RR2,   & ! RR2 = RR*RR
        S1,    & ! Sine of doubled first rotation angle
        S2       ! Sine of doubled second rotation angle
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
    CINC2 = -MUI*MUR - SI*SR*CAZ ! Eq.(7)
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
    F3 = RL*RR
!
!   Compute the rotation matrix (think about vectorization)
    DO IA = 1, NA
        CALL ROTATOR(MUI, MUR, AZI(IA), S1(IA), C1(IA), S2(IA), C2(IA))
    END DO ! IA = 1, NA
!
!   Apply the rotator       Position of an element in the 3-by-3 Fresnel matrix:
    RFR(:, 1) =  F1                                                     ! [1, 1]
!   Non polarizing reflection: RFR0(:, 2:9) = 0
    RFR(:, 2) =  F2*C2                                                  ! [2, 1]
    RFR(:, 3) =  F2*S2                                                  ! [3, 1]
    RFR(:, 4) =  F2*C1                                                  ! [1, 2]
    RFR(:, 5) =  F1*C1*C2 - F3*S1*S2                                    ! [2, 2]
    RFR(:, 6) =  F1*C1*S2 + F3*S1*C2                                    ! [3, 2]
    RFR(:, 7) = -F2*S1                                                  ! [1, 3]
    RFR(:, 8) = -F1*S1*C2 - F3*C1*S2                                    ! [2, 3]
    RFR(:, 9) = -F1*S1*S2 + F3*C1*C2                                    ! [3, 3]
!
END SUBROUTINE FRESNELR_IP
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 11Sep14 - ROTATOR(cosine) is now used instead of ROTATOR(angle). Retested: Ok. 
!
! 10Sep14 - Created from FRESNELR for a single point {SZA, VZA, AZA}. Tested
!           against the old version for SZA = VZA = 0:1:89, AZI=0:1:360,
!           REFRE = 1.33, REFIM = 0.005. Coincidence within REAL*8 precision.
!===============================================================================