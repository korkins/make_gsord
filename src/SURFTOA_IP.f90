SUBROUTINE SURFTOA_IP(MU0, TAU0, MUR, NR, ISRF, PSRF, AZI, NA, IS0, QS0, US0)
!===============================================================================
! PURPOSE:
!   To compute contribution from reflection of the direct solar beam from
!   surface on TOA. The incident beam is assumed to be [2pi 0 0].
!
! INPUT:
!   MU0    D(1)      Cosine of solar zenith angle, mu0 > 0
!   TAU0   D(1)      Total optical thickness of atmosphere
!   MUR    D(NR)     Cosines of reflected view zeniths, mu_r < 0
!   NR     I(1)      Number of reflected zeniths
!   ISRF   I(1)      Index of the surface model, ISRF > 0
!   PSRF   D(NSPR)   Vector of parameters of a surface
!   AZI    D(NA)     View azimuths, [0:2PI], in radians
!   NA     I(1)      Number of view azimuths
!
! OUTPUT:
!   IS0, QS0, US0   D(NR, NA)   Surface contribution on TOA from each component
!
! TREE:
!   SURFTOA_IP
!            |
!            +-RPV (the Rahman-Pinty-Verstraete model for land reflection)
!            |
!            +-RTLS (the Ross-Thick Lee-Sparse model for land reflection)
!            |
!            +-FRESNELR_IP (Fresnel reflection matrix with rotation)
!            |
!            +-ROUGHSRF (Ocean waves)
!
! COMMENTS:
!   Note that output is normalized by 2pi to agree with the main code.
!
!   Contribution from each component is computed using a simple expression:
!
!   IS0(mu, mu0, azi) = mu0*exp(Tau/mu)*SURF_REF(mu, mu0, azi)*exp(-Tau0/mu0),
!
!   where Tau0 is the total optical thickness of atmosphere.
!
!   Current values of ISPR:
!       1 - Lambertian;
!       2 - RPV;
!       3 - mRPV;
!       4 - RTLS;
!       5 - modified RTLS;
!       6 - Azimuthally symmetric ocean reflection (Nakajima-Tanaka) model;
!       7 - Wind driven ocean reflection (Cox-Munk & Gramm-Charlie) model;
!
!       10 - linear mixing of 4 & 6 (RTLS & NT-ocean);
!       11 - linear mixing of 5 & 6 (mRTLS & NT-ocean);
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    INTEGER, PARAMETER :: &
        NSPR = 21              ! Max number of parameters for surface model
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: ISRF, NR, NA
    REAL*8, INTENT(IN) :: MU0, TAU0
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NSPR), INTENT(IN) :: PSRF
    REAL*8, DIMENSION(NR), INTENT(IN) :: MUR
    REAL*8, DIMENSION(NA), INTENT(IN) :: AZI
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NR, NA), INTENT(OUT) :: IS0, QS0, US0
!
! LOCAL VARIABLES
    INTEGER &
        MSRF, & ! MSRF=1: RPV, RTLS; MSRF=2: MRPV, MRTLS ...
        IA,   & ! Loop index over azimuth
        IR      ! Loop index over reflected beams
    REAL*8 &
        SW1, & ! Weight of the 1st surface in mixture
        SW2    ! Weight of the 2nd surface in mixture, 1.0-SW1
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NR) :: &
        ATTN ! mu0*exp(-Tau0/mu0)*exp(Tau0/mu)
    REAL*8, DIMENSION(NA) :: &
        BRDF, & ! Scalar surface reflection as a function of azimuth
        WSRF    ! Waves on surface
    REAL*8, DIMENSION(NA, 3) :: &
        BPDF ! Mueller matrix of surface as a function of azimuth
!===============================================================================
!
    MSRF = 1 ! By default: RPV, RTLS, NT, Vegetation, ...
    IF (PRODUCT((/3, 5, 7, 9, 11/)-ISRF) == 0) MSRF = 2 ! mRPV, mRTLS, ...
!
    IF (ISRF == 10 .OR. ISRF == 11) THEN
        SW1 = PSRF(NSPR)
        SW2 = 1.0D0 - SW1
    END IF ! ISRF = 10, 11
!
    SELECT CASE (ISRF)
        CASE (1) ! Lambert
            IS0 = PSRF(1)
        CASE (2:3) ! RPV or MRPV
            DO IR = 1, NR
                CALL RPV(MSRF, PSRF(1), PSRF(1), PSRF(2), PSRF(3), MU0, &
                                MUR(IR), AZI, NA, BRDF)
                IS0(IR, :) = BRDF  
            END DO
        CASE (4:5) ! RTLS or MRTLS
            DO IR = 1, NR
                CALL RTLS(MSRF, PSRF(1), PSRF(2), PSRF(3), MU0, &
                                MUR(IR), AZI, NA, BRDF)
                IS0(IR, :) = BRDF
            END DO ! IR = 1, NR
        CASE (6:7) ! NT or GC ocean
            DO IR = 1, NR
                CALL FRESNELR0_IP(PSRF(1), PSRF(2), MU0, &
                                MUR(IR), AZI, NA, BPDF)
                CALL ROUGHSRF(MSRF, PSRF(3), PSRF(4), MU0, &
                                MUR(IR), AZI, NA, BRDF)
                IS0(IR, :) = BRDF*BPDF(:, 1) ! I-component
                QS0(IR, :) = BRDF*BPDF(:, 2) ! Q-component
                US0(IR, :) = BRDF*BPDF(:, 3) ! U-component
            END DO ! IR = 1, NR
        CASE (8:9) ! Nadal-Breon
            ! CALL NADBRESRF - to be added
        CASE (10:11) ! Mixture of RTLS and NT ocean
            DO IR = 1, NR
                CALL FRESNELR0_IP(PSRF(11), PSRF(12), MU0, MUR(IR), AZI, &
                                  NA, BPDF)
                CALL ROUGHSRF(MSRF, PSRF(13), PSRF(14), MU0, MUR(IR), AZI, &
                                  NA, WSRF)
                CALL RTLS(MSRF, PSRF(1), PSRF(2), PSRF(3), MU0, MUR(IR), AZI, &
                                  NA, BRDF)
                IS0(IR, :) = SW1*BRDF + SW2*WSRF*BPDF(:, 1)
                QS0(IR, :) = SW2*WSRF*BPDF(:, 2)
                US0(IR, :) = SW2*WSRF*BPDF(:, 3)
            END DO ! IR = 1, NR
!
    END SELECT ! CASE (ISRF)
!
!   Attenuation factor
    ATTN = MU0*DEXP( TAU0*(MU0 - MUR)/(MU0*MUR) )
!
!   The factor of 2.0 comes from 1/pi(surface integral)*2pi(incident beam) = 2
    DO IA = 1, NA; IS0(:, IA) = IS0(:, IA)*ATTN; END DO
    IS0 = 2.0D0*IS0
    IF (ISRF > 5) THEN
        DO IA = 1, NA
            QS0(:, IA) = QS0(:, IA)*ATTN
            US0(:, IA) = US0(:, IA)*ATTN
        END DO ! IA = 1, NA
        QS0 = 2.0D0*QS0
        US0 = 2.0D0*US0
    END IF ! ISRF > 5
!
END SUBROUTINE SURFTOA_IP
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 02May15 - Mixture of RTLS&NT (ISRF=10) is modified. ISRF=11 is added.
!
! 16Nov14 - The RPV surface was added and tested in the Lambertian mode against
!           SORD+true Lambertian surface. Also tested RPV, R0=RC=0.219, K=0.463,
!           ASC=-0.287. See SURFACEM0_IP for details.
!
! 11Oct14 - Renamed: TYP -> MSRF, NPAR -> NSPR. Checked for ISRF = 4 against old
!           version.
!
! 22Sep14 - Created and tested as part of SORD for Lambert (ISRF = 1), RTLS
!           (ISRF = 4), NT Ocean (RTLS = 6), and NT+RTLS (ISRF=10). Perfect
!           agreement with the results reported in the SORD_IP subroutine.
!===============================================================================