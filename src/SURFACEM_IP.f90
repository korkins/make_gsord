SUBROUTINE SURFACEM_IP(IM, ISRF, MUWG, WCMA, WSMA, NG1, NR, NGA, NM,           &
                                                                R11, R12, R13, &
                                                                R21, R22, R23, &
                                                                R31, R32, R33, &
                                                             RM11, RM12, RM13, &
                                                             RM21, RM22, RM23, &
                                                             RM31, RM32, RM33)
!===============================================================================
! PURPOSE:
!   To compute the m-th Fourier expansion coefficient for diffuse surface
!   reflection. Zenith Gauss weights for integration over incident beams are
!   already included in Rij.
!
! INPUT:
!   IM     I(1)              The Fourier moment, IM = 1, 2, ... (m = 0, 1, ...)
!   ISRF   I(1)              Index of the surface model, ISRF > 0
!   MUWG   D(NG1)            MUI*WG for integration over incident beams
!   WCMA   D(NGA, NM)        Gauss_weight*cos(m*az) for all m and az
!   WSMA   D(NGA, NM)        Same as WCMA except for WSMA
!   NG1    I(1)              Number of incident zeniths (Gauss nodes)
!   NR     I(1)              Number of reflected zeniths
!   NGA    I(1)              Number of Gauss nodes for integration over azimuth
!   NM     I(1)              Total number of the Fourier moments
!   Rij    D(NGA, NG1, NR)   Weighted surface reflection matrix.
!
! OUTPUT:
!   RMij   D(NG1, NR)   The m-th Fourier expansion moment of each element of the
!                       Mueller matrix of surface for all NG1 incident & all NR
!                       reflected beams
!
! TREE:
!   -
!
! COMMENTS:
!   -
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: IM, ISRF, NGA, NG1, NR, NM
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NGA, NM), INTENT(IN) :: WCMA, WSMA
    REAL*8, DIMENSION(NG1), INTENT(IN) :: MUWG
    REAL*8, DIMENSION(NGA, NG1, NR), INTENT(IN) :: R11, R12, R13, &
                                                   R21, R22, R23, &
                                                   R31, R32, R33
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NG1, NR), INTENT(OUT) :: RM11, RM12, RM13, &
                                               RM21, RM22, RM23, &
                                               RM31, RM32, RM33
!
! LOCAL VARIABLES
    INTEGER &
        IG,  & ! Loop index over incident beams
        IR,  & ! Loop index over reflected beams
        NRU    ! Number of user-defined directions of reflection, NR-NG1
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NG1) :: &
        MUWG2 ! 2.0*MUWG for convenience
!===============================================================================
!
    NRU = NR-NG1
    MUWG2 = 2.0D0*MUWG
!
    IF (ISRF < 6) THEN
        DO IR = 1, NRU ! if NRU = 0, the loop is skipped
            DO IG = 1, NG1
                RM11(IG, IR) = SUM(R11(:, IG, IR)*WCMA(:, IM))
            END DO ! IG = 1, NG1
        END DO ! IR = 1, NR
        DO IR = 1, NG1
            DO IG = 1, NG1
                IF (IG < IR) THEN
                    RM11(IG, IR+NRU) = RM11(IR, IG+NRU)
                ELSE ! IG >= IR
                    RM11(IG, IR+NRU) = SUM(R11(:, IG, IR+NRU)*WCMA(:, IM))
                END IF ! IG < IR
            END DO ! IG = 1, NG1
        END DO ! IR = 1, NR
!
!       RMij must be weighted AFTER symmetry is applied
        DO IR = 1, NR 
            RM11(:, IR) = MUWG2*RM11(:, IR)
        END DO ! IR = 1, NR
!
    ELSE ! ISRF > 5
!
        IF (IM == 1) THEN
            DO IR = 1, NRU
                DO IG = 1, NG1
                    RM11(IG, IR) = SUM(R11(:, IG, IR)*WCMA(:, IM))
                    RM12(IG, IR) = SUM(R12(:, IG, IR)*WCMA(:, IM))
                    RM21(IG, IR) = SUM(R21(:, IG, IR)*WCMA(:, IM))
                    RM22(IG, IR) = SUM(R22(:, IG, IR)*WCMA(:, IM))
                    RM33(IG, IR) = SUM(R33(:, IG, IR)*WCMA(:, IM))
                    RM31(IG, IR) = 0.0D0
                    RM32(IG, IR) = 0.0D0
                    RM13(IG, IR) = 0.0D0
                    RM23(IG, IR) = 0.0D0
                END DO ! IG = 1, NG1
            END DO ! IR = 1, NR
!      
            DO IR = 1, NG1
                DO IG = 1, NG1
                    IF (IG < IR) THEN
                        RM11(IG, IR+NRU) = RM11(IR, IG+NRU)
                        RM22(IG, IR+NRU) = RM22(IR, IG+NRU)
                        RM33(IG, IR+NRU) = RM33(IR, IG+NRU)
                    ELSE ! IG >= IR
                        RM11(IG, IR+NRU) = SUM(R11(:, IG, IR+NRU)*WCMA(:, IM))
                        RM22(IG, IR+NRU) = SUM(R22(:, IG, IR+NRU)*WCMA(:, IM)) 
                        RM33(IG, IR+NRU) = SUM(R33(:, IG, IR+NRU)*WCMA(:, IM))
                    END IF ! IG < IR
                    RM12(IG, IR+NRU) = SUM(R12(:, IG, IR+NRU)*WCMA(:, IM))
                    RM21(IG, IR+NRU) = SUM(R21(:, IG, IR+NRU)*WCMA(:, IM))
                    RM31(IG, IR+NRU) = 0.0D0
                    RM32(IG, IR+NRU) = 0.0D0
                    RM13(IG, IR+NRU) = 0.0D0
                    RM23(IG, IR+NRU) = 0.0D0
                END DO ! IG = 1, NG1
            END DO ! IR = 1, NR
!
            DO IR = 1, NR
                RM11(:, IR) = MUWG2*RM11(:, IR)
                RM12(:, IR) = MUWG2*RM12(:, IR)
                RM21(:, IR) = MUWG2*RM21(:, IR)
                RM22(:, IR) = MUWG2*RM22(:, IR)
                RM33(:, IR) = MUWG2*RM33(:, IR)
            END DO ! IR = 1, NR
!          
        ELSE ! IM > 1
!
            DO IR = 1, NRU
                DO IG = 1, NG1
                    RM11(IG, IR) =  SUM(R11(:, IG, IR)*WCMA(:, IM))
                    RM12(IG, IR) =  SUM(R12(:, IG, IR)*WCMA(:, IM))
                    RM21(IG, IR) =  SUM(R21(:, IG, IR)*WCMA(:, IM))
                    RM22(IG, IR) =  SUM(R22(:, IG, IR)*WCMA(:, IM))
                    RM33(IG, IR) =  SUM(R33(:, IG, IR)*WCMA(:, IM))
                    RM31(IG, IR) =  SUM(R31(:, IG, IR)*WSMA(:, IM))
                    RM32(IG, IR) =  SUM(R32(:, IG, IR)*WSMA(:, IM))
                    RM13(IG, IR) = -SUM(R13(:, IG, IR)*WSMA(:, IM))
                    RM23(IG, IR) = -SUM(R23(:, IG, IR)*WSMA(:, IM))
                END DO ! IG = 1, NG1
            END DO ! IR = 1, NR
            DO IR = 1, NG1
                DO IG = 1, NG1
                    IF (IG < IR) THEN
                        RM11(IG, IR+NRU) = RM11(IR, IG+NRU)
                        RM22(IG, IR+NRU) = RM22(IR, IG+NRU)
                        RM33(IG, IR+NRU) = RM33(IR, IG+NRU)
                    ELSE ! IG >= IR
                        RM11(IG, IR+NRU) = SUM(R11(:, IG, IR+NRU)*WCMA(:, IM))
                        RM22(IG, IR+NRU) = SUM(R22(:, IG, IR+NRU)*WCMA(:, IM)) 
                        RM33(IG, IR+NRU) = SUM(R33(:, IG, IR+NRU)*WCMA(:, IM))
                    END IF ! IG < IR
                    RM12(IG, IR+NRU) =  SUM(R12(:, IG, IR+NRU)*WCMA(:, IM))
                    RM21(IG, IR+NRU) =  SUM(R21(:, IG, IR+NRU)*WCMA(:, IM)) 
                    RM31(IG, IR+NRU) =  SUM(R31(:, IG, IR+NRU)*WSMA(:, IM))
                    RM32(IG, IR+NRU) =  SUM(R32(:, IG, IR+NRU)*WSMA(:, IM))
                    RM13(IG, IR+NRU) = -SUM(R13(:, IG, IR+NRU)*WSMA(:, IM))
                    RM23(IG, IR+NRU) = -SUM(R23(:, IG, IR+NRU)*WSMA(:, IM))
                END DO ! IG = 1, NG1
            END DO ! IR = 1, NR
!
            DO IR = 1, NR
                RM11(:, IR) = MUWG2*RM11(:, IR)
                RM12(:, IR) = MUWG2*RM12(:, IR)
                RM13(:, IR) = MUWG2*RM13(:, IR)
                RM21(:, IR) = MUWG2*RM21(:, IR)
                RM22(:, IR) = MUWG2*RM22(:, IR)
                RM23(:, IR) = MUWG2*RM23(:, IR)
                RM31(:, IR) = MUWG2*RM31(:, IR)
                RM32(:, IR) = MUWG2*RM32(:, IR)
                RM33(:, IR) = MUWG2*RM33(:, IR)
            END DO ! IR = 1, NR
!
        END IF ! IM == 1
!
    END IF ! ISRF < 6
!
END SUBROUTINE SURFACEM_IP
!===============================================================================
! 01May16 - IRU is repelced with IR+NRU
!
! 22Apr16 - Minor changes in comments
!
! 01May15 - Tested for RTLS & NT ocean.
!
! 24Apr15 - First created
!===============================================================================