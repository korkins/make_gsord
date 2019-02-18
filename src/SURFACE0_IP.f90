SUBROUTINE SURFACE0_IP(MU0, ISRF, PSRF, MUR, NR, AZI, NAZ, R01, R02, R03)
!===============================================================================
! PURPOSE:
!   Precomputes BRDF or 1st column of BPDF for a given direction of incidence  &
!   all directions of reflection.
!
! INPUT:
!   MU0    D(1)      mu0 = cos(SZA) > 0
!   ISRF   I(1)      Index of the surface model, ISRF > 1
!   PSRF   D(NSPR)   Vector of the surface parameters
!   MUR    I(NR)     Cosines of reflected zenith angles, MUR < 0
!   NR     I(1)      Number of reflected zeniths
!   AZI    D(NAZ)    Relative azimuth, radians
!   NAZ    I(1)      Number of azimuths
!
! OUTPUT:
!   R01, R02, R03   D(NAZ, NR)   1st column of the Mueller matrix of surface
!
! TREE:
!   SURFACE0_IP
!             |
!             +-RTLS
!             |
!             +-FRESNELR0_IP
!             |
!             +-ROUGHSRF
!
! COMMENTS:
!   Azimuth comes in 1st dimension to allow for faster Fourier expansion, i.e.
!   integration over azimuth, for any given angle of reflection (observation).
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    INTEGER, PARAMETER :: &
        NSPR = 21 ! Max number of parameters for surface model
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: ISRF, NR, NAZ
    REAL*8, INTENT(IN) :: MU0
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NSPR), INTENT(IN) :: PSRF
    REAL*8, DIMENSION(NR), INTENT(IN) :: MUR
    REAL*8, DIMENSION(NAZ), INTENT(IN) :: AZI
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NAZ, NR), INTENT(OUT) :: R01, R02, R03
!
! LOCAL VARIABLES
    INTEGER &
        IR, & ! Loop index over reflected zeniths
        MSRF  ! Modification of surface models
    REAL*8 &
        SW1, & ! Weight of the 1st surface in mixture
        SW2    ! Weight of the 2nd surface in mixture, 1.0-SW1
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NAZ) :: &
        BRDF, & ! Scalar reflection
        WSRF    ! Waves on surface
    REAL*8, DIMENSION(NAZ, 3) :: &
        BPDF ! Vector reflection (1st column only)
!===============================================================================
!
!   Original or modified surface model
    MSRF = 1
    IF (PRODUCT((/3, 5, 7, 9, 11/)-ISRF) == 0) MSRF = 2
!
    IF (ISRF == 10 .OR. ISRF == 11) THEN
        SW1 = PSRF(NSPR)
        SW2 = 1.0D0 - SW1
    END IF ! ISRF = 10, 11
!
    SELECT CASE (ISRF)
        CASE (4:5) ! RTLS & mRTLS
            DO IR = 1, NR
                CALL RTLS(MSRF, PSRF(1), PSRF(2), PSRF(3), MU0, MUR(IR), AZI, &
                          NAZ, R01(:, IR))
            END DO ! IR = 1, NR
        CASE (6:7) ! NT & GC Ocean
            DO IR = 1, NR
                CALL FRESNELR0_IP(PSRF(1), PSRF(2), MU0, MUR(IR), AZI, &
                                  NAZ, BPDF)
                CALL ROUGHSRF(MSRF, PSRF(3), PSRF(4), MU0, MUR(IR), AZI, &
                              NAZ, BRDF)
                R01(:, IR) = BRDF*BPDF(:, 1) ! I-component
                R02(:, IR) = BRDF*BPDF(:, 2) ! Q-component
                R03(:, IR) = BRDF*BPDF(:, 3) ! U-component
            END DO ! IR = 1, NR
        CASE (10:11) ! Mix RTLS & Ocean
            DO IR = 1, NR
                CALL FRESNELR0_IP(PSRF(11), PSRF(12), MU0, MUR(IR), AZI, &
                                  NAZ, BPDF)
                CALL ROUGHSRF(MSRF, PSRF(13), PSRF(14), MU0, MUR(IR), AZI, &
                                  NAZ, WSRF)
                CALL RTLS(MSRF, PSRF(1), PSRF(2), PSRF(3), MU0, MUR(IR), AZI, &
                                  NAZ, BRDF)
                R01(:, IR) = SW1*BRDF + SW2*WSRF*BPDF(:, 1)
                R02(:, IR) = SW2*WSRF*BPDF(:, 2)
                R03(:, IR) = SW2*WSRF*BPDF(:, 3)
            END DO ! IR = 1, NR
!  
    END SELECT ! CASE (ISRF)
!
END SUBROUTINE SURFACE0_IP
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 02May15 - ISRF = 10, 11 are added.
!
! 01May15 - Ocean is added and tested.
!
! 24Apr15 - Azimuth weights, WAZ, are now excluded from this subroutine.
!           The subroutine is renamed: WGASURF0_IP -> SURFACE0_IP
!
! 23Apr15 - Output is redefined: 3 2D arrays are now used. Tested.
!
! 14Nov14 - First created and tested for RTLS only.
!===============================================================================