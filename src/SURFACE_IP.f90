SUBROUTINE SURFACE_IP(ISRF, PSRF, MUI, NG1, MUR, NR, AZI, NAZ, &
                      R11, R12, R13, R21, R22, R23, R31, R32, R33)
!===============================================================================
! PURPOSE:
!   To compute the 3-by-3 Mueller matrix for surface reflection (pBRDF) for a
!   set of directions of incidence, MUI (Gauss nodes), observation, MUR, and
!   relative azimuths, AZI. ISRF specifies the surface.
!
! INPUT:
!   ISRF   I(1)      Index of the surface model, ISRF > 1
!   PSRF   D(NSPR)   Vector of parameters of the surface
!   MUI    D(NG1)    Cosines of incident zenith angles, MUI > 0
!   NG1    I(1)      Number of incident zeniths (Gauss nodes), MUI > 0
!   MUR    I(NR)     Cosines of reflected zenith angles, MUR < 0
!   NR     I(1)      Number of directions of reflection
!   AZI    D(NAZ)    Azimuths of observation, [0:2PI] rad
!   NAZ    I(1)      Number of relative azimuths, AZI
!
! OUTPUT:
!   Rij   D(NAZ, NG1, NR)   pBRDF of the surface for the given geometries
!
! TREE:
!   SURFACE_IP
!            |
!            +-RTLS (the Ross-Thick Lee-Sparse model for land reflection)
!            |
!            +-FRESNELR_IP (Fresnel reflection matrix with rotation)
!            |
!            +-ROUGHSRF (Ocean waves)
!
! COMMENTS:
!   *** To be added: RPV (ISRF=2) & modified RPV (ISRF=3) ***
!   *** To be added: Nadal-Breon (ISRF=8)                 ***
!
!   Current values of ISPR:
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
        NSPR = 21 ! Max number of parameters for surface model
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: ISRF, NG1, NR, NAZ
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NSPR), INTENT(IN) :: PSRF
    REAL*8, DIMENSION(NG1), INTENT(IN) :: MUI
    REAL*8, DIMENSION(NR), INTENT(IN) :: MUR
    REAL*8, DIMENSION(NAZ), INTENT(IN) :: AZI
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NAZ, NG1, NR), INTENT(OUT) :: R11, R12, R13, &
                                                    R21, R22, R23, &
                                                    R31, R32, R33
!
! LOCAL VARIABLES
    INTEGER &
        IG,   & ! Loop index over incident beams
        IR,   & ! Loop index over reflected beams
        MSRF, & ! Modified surfaces: mRPV, mRTLS, GC Ocean ...
        NRU     ! Number of user-defined directions of reflection, NR-NG1
    REAL*8 &
        SW1, & ! Weight of the 1st surface in mixture
        SW2    ! Weight of the 2nd surface in mixture, 1.0-SW1
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NAZ) :: &
        BRDF, & ! Scalar reflection
        WSRF    ! Waves on ocean
    REAL*8, DIMENSION(NAZ, 9) :: &
        BPDF  ! Vector reflection
    REAL*8, DIMENSION(NAZ, NG1, NR) :: &
        WSRFS ! Same as BRDF but stores values in RAM to apply symmetry
!===============================================================================
!
    MSRF = 1
    NRU = NR-NG1
    IF (PRODUCT((/3, 5, 7, 9, 11/)-ISRF) == 0) MSRF = 2
!
    IF (ISRF == 10 .OR. ISRF == 11) THEN
        SW1 = PSRF(NSPR)
        SW2 = 1.0D0 - SW1
    END IF ! ISRF = 10, 11
!
    SELECT CASE (ISRF)
        CASE (4:5) ! RTLS & mRTLS
!           No symmetry applied
            DO IR = 1, NRU ! if NRU = 0, the loop is skipped
                DO IG = 1, NG1
                    CALL RTLS(MSRF, PSRF(1), PSRF(2), PSRF(3), MUI(IG), &
                                MUR(IR), AZI, NAZ, R11(:, IG, IR))
                END DO ! IG = 1, NG1
            END DO ! IR = 1, NR
!           Symmetry applied
            DO IR = 1, NG1
                DO IG = IR, NG1 ! IG < IR: skip surface computation
                    CALL RTLS(MSRF, PSRF(1), PSRF(2), PSRF(3), MUI(IG), &
                              MUR(IR+NRU), AZI, NAZ, R11(:, IG, IR+NRU))
                END DO ! IG = 1, NG1
            END DO ! IMU = 1, NMU
!
        CASE (6:7) ! NT & GC Ocean
            DO IR = 1, NRU
                DO IG = 1, NG1
                    CALL FRESNELR_IP(PSRF(1), PSRF(2), MUI(IG), MUR(IR),    &
                                        AZI, NAZ, BPDF)
                    CALL ROUGHSRF(MSRF, PSRF(3), PSRF(4), MUI(IG), MUR(IR), &
                                        AZI, NAZ, WSRF)
                    R11(:, IG, IR) = WSRF*BPDF(:, 1)
                    R21(:, IG, IR) = WSRF*BPDF(:, 2)
                    R31(:, IG, IR) = WSRF*BPDF(:, 3)
                    R12(:, IG, IR) = WSRF*BPDF(:, 4)
                    R22(:, IG, IR) = WSRF*BPDF(:, 5)
                    R32(:, IG, IR) = WSRF*BPDF(:, 6)
                    R13(:, IG, IR) = WSRF*BPDF(:, 7)
                    R23(:, IG, IR) = WSRF*BPDF(:, 8)
                    R33(:, IG, IR) = WSRF*BPDF(:, 9)          
                END DO ! IG = 1, NG1
            END DO ! IR = 1, NR
!
            DO IR = 1, NG1
                DO IG = 1, NG1
                    CALL FRESNELR_IP(PSRF(1), PSRF(2), MUI(IG), MUR(IR+NRU), &
                                        AZI, NAZ, BPDF)
                    IF (IG < IR) THEN ! IG < IR: do not compute the diagonal
                        WSRFS(:, IG, IR) = WSRFS(:, IR, IG)
                    ELSE ! IG >= IR; Symmetry: diagonal elements ONLY                     
                        CALL ROUGHSRF(MSRF, PSRF(3), PSRF(4), MUI(IG), &
                                        MUR(IR+NRU), AZI, NAZ, WSRFS(:, IG, IR))
                        R11(:, IG, IR+NRU) = WSRFS(:, IG, IR)*BPDF(:, 1)
                        R22(:, IG, IR+NRU) = WSRFS(:, IG, IR)*BPDF(:, 5)
                        R33(:, IG, IR+NRU) = WSRFS(:, IG, IR)*BPDF(:, 9)
                    END IF ! IG < IR                 
!
                    R21(:, IG, IR+NRU) = WSRFS(:, IG, IR)*BPDF(:, 2)
                    R31(:, IG, IR+NRU) = WSRFS(:, IG, IR)*BPDF(:, 3)
                    R12(:, IG, IR+NRU) = WSRFS(:, IG, IR)*BPDF(:, 4)
                    R32(:, IG, IR+NRU) = WSRFS(:, IG, IR)*BPDF(:, 6)
                    R13(:, IG, IR+NRU) = WSRFS(:, IG, IR)*BPDF(:, 7)
                    R23(:, IG, IR+NRU) = WSRFS(:, IG, IR)*BPDF(:, 8)
                END DO ! IG = 1, NG1
            END DO ! IR = 1, NG1
!
        CASE (10:11) ! Mix RTLS & Ocean
            DO IR = 1, NRU
                DO IG = 1, NG1
                    CALL FRESNELR_IP(PSRF(11), PSRF(12), MUI(IG), MUR(IR),    &
                                        AZI, NAZ, BPDF)
                    CALL ROUGHSRF(MSRF, PSRF(13), PSRF(14), MUI(IG), MUR(IR), &
                                        AZI, NAZ, WSRF)
                    CALL RTLS(MSRF, PSRF(1), PSRF(2), PSRF(3), MUI(IG), &
                                        MUR(IR), AZI, NAZ, BRDF)
                    R11(:, IG, IR) = SW1*BRDF + &
                                     SW2*WSRF*BPDF(:, 1)
                    R21(:, IG, IR) = SW2*WSRF*BPDF(:, 2)
                    R31(:, IG, IR) = SW2*WSRF*BPDF(:, 3)
                    R12(:, IG, IR) = SW2*WSRF*BPDF(:, 4)
                    R22(:, IG, IR) = SW2*WSRF*BPDF(:, 5)
                    R32(:, IG, IR) = SW2*WSRF*BPDF(:, 6)
                    R13(:, IG, IR) = SW2*WSRF*BPDF(:, 7)
                    R23(:, IG, IR) = SW2*WSRF*BPDF(:, 8)
                    R33(:, IG, IR) = SW2*WSRF*BPDF(:, 9)          
                END DO ! IG = 1, NG1
            END DO ! IR = 1, NR
!
            DO IR = 1, NG1
                DO IG = 1, NG1
                    CALL FRESNELR_IP(PSRF(11), PSRF(12), MUI(IG), MUR(IR+NRU), &
                                        AZI, NAZ, BPDF)
                    IF (IG < IR) THEN ! IG < IR: do not compute the diagonal
                        WSRFS(:, IG, IR) = WSRFS(:, IR, IG)
                    ELSE ! IG >= IR; Symmetry: diagonal elements ONLY                     
                        CALL ROUGHSRF(MSRF, PSRF(13), PSRF(14), MUI(IG), &
                                        MUR(IR+NRU), AZI, NAZ, WSRFS(:, IG, IR))
                        CALL RTLS(MSRF, PSRF(1), PSRF(2), PSRF(3), MUI(IG), &
                              MUR(IR+NRU), AZI, NAZ, BRDF)
                        R11(:, IG, IR+NRU) = SW1*BRDF + &
                                          SW2*WSRFS(:, IG, IR)*BPDF(:, 1)        
                        R22(:, IG, IR+NRU) = SW2*WSRFS(:, IG, IR)*BPDF(:, 5)
                        R33(:, IG, IR+NRU) = SW2*WSRFS(:, IG, IR)*BPDF(:, 9)
                    END IF ! IG < IR                 
!
                    R21(:, IG, IR+NRU) = SW2*WSRFS(:, IG, IR)*BPDF(:, 2)
                    R31(:, IG, IR+NRU) = SW2*WSRFS(:, IG, IR)*BPDF(:, 3)
                    R12(:, IG, IR+NRU) = SW2*WSRFS(:, IG, IR)*BPDF(:, 4)
                    R32(:, IG, IR+NRU) = SW2*WSRFS(:, IG, IR)*BPDF(:, 6)
                    R13(:, IG, IR+NRU) = SW2*WSRFS(:, IG, IR)*BPDF(:, 7)
                    R23(:, IG, IR+NRU) = SW2*WSRFS(:, IG, IR)*BPDF(:, 8)
                END DO ! IG = 1, NG1
            END DO ! IR = 1, NG1
!
    END SELECT ! CASE (ISRF)
!
END SUBROUTINE SURFACE_IP
!===============================================================================
! 01May16 - IRU is replaced with IR+NRU
!
! 22Apr16 - Minor changes in comments
!
! 02May15 - ISRF = 10 & 11 are added.
!
! 01May15 - Ocean reflection, with symmetry for diagonal elements, is added &
!           tested.
!
! 17Apr14 - Created from SURFACEM_IP and tested for ISRF = 4 (RTLS). 
!===============================================================================