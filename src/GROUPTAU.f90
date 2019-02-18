SUBROUTINE GROUPTAU(NZ, ZKM, NL, TAU, NX, XKM, NG, TAUG)
!===============================================================================
! PURPOSE:
!   To rearrange the input TAU(ZKM) on a new grid, XKM
!
! INPUT:
!   NZ    I(1)    Number of heights for the input TAU
!   ZKM   D(NZ)   Heights for TAU, TOA to BOA; ZKM(1)=XKM(1) & ZKM(NZ)=XKM(NX)
!   NL    I(1)    Number of input layers, NZ-1
!   TAU   D(NL)   Thickness for each layer
!   NX    I(1)    Number of heights for the output TAUG
!   XKM   D(NX)   Heights for TAUG, TOA to BOA; XKM(1)=ZKM(1) & XKM(NX)=ZKM(NZ)
!   NG    I(1)    Number of output layers, NX-1
!
! OUTPUT:
!   TAUG  D(NG)   TAU(ZKM) rearranged according to XKM
!
! TREE:
!   -
!
! COMMENTS:
!   ZKM and XKM must coincide on TOA & BOA, respectively.
!
!   The idea is to compute contribution (weight) of each input layer TAU(IL),
!   located in ZKM(IL):ZKM(IL+1), to each output layer TAUG(IG), located in
!   XKM(IG):XKM(IG+1). There are 4 options:
!
!   Case 0: TAU(IL) & TAUG(IG) do not overlap - zero weight, W = 0
!   Case 1: TAU(IL) completely fits into TAUG(IG). W = 1.
!           XKM(IG) >= ZKM(IL) > ZKM(IL+1) >= XKM(IG+1)
!   Case 2: Opposite to Case 1, TAU(IG) completely fits into TAU(IL)
!           ZKM(IL) > XKM(IG) > XKM(IG+1) > ZKM(IL+1)
!           W = (XKM(IG) - XKM(IG+1))/DZ
!   Case 3: TAU(IL) partially overlaps TAU(IG) from the top
!           ZKM(IL) > XKM(IG) > ZKM(IL+1) > XKM(IG+1)
!           W = (XKM(IG) - ZKM(IL+1))/DZ
!   Case 4: TAU(IL) partially overlaps TAU(IG) from the bottom
!           XKM(IG) > ZKM(IL) > XKM(IG+1) >= ZKM(IL+1)
!           W = (ZKM(IL) - XKM(IG+1))/DZ
!
!   In Cases 1-4, DZ = ZKM(IL) - ZKM(IL+1).
!
!   TOP & BOT indices for IX & IZ are introduced for convenience only. 
!
! REFERENCES:
!   -
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NZ, NL, NX, NG
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NZ), INTENT(IN) :: ZKM
    REAL*8, DIMENSION(NL), INTENT(IN) :: TAU
    REAL*8, DIMENSION(NX), INTENT(IN) :: XKM
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NG), INTENT(OUT) :: TAUG
!
! LOCAL VARIBLES
    INTEGER &
        IG,             & ! Loop index over TAUG(IG) ("steady" layer)
        IL,             & ! Loop index over TAU(L)   ("moving" layer)
        IX_TOP, IX_BOT, & ! Define position of TAUG(IG)
        IZ_TOP, IZ_BOT    ! Define position of TAU(IL)
!===============================================================================
!
    DO IG = 1, NG-1
        IX_TOP = IG
        IX_BOT = IG+1
        TAUG(IG) = 0.0D0 ! Initialize & account for non-overlapping layers
        DO IL = 1, NL
            IZ_TOP = IL
            IZ_BOT = IL+1
            IF (ZKM(IZ_TOP) <= XKM(IX_TOP) .AND. &
                ZKM(IZ_BOT) >= XKM(IX_BOT)) THEN     ! Case (1)
                TAUG(IG) = TAUG(IG) + TAU(IL)
            ELSE IF (ZKM(IZ_TOP) > XKM(IX_TOP) .AND. &
                     ZKM(IZ_BOT) < XKM(IX_BOT)) THEN ! Case (2)
                 TAUG(IG) = TAUG(IG) + TAU(IL)*(XKM(IX_TOP) - XKM(IX_BOT))/ &
                                               (ZKM(IZ_TOP) - ZKM(IZ_BOT))
            ELSE IF (ZKM(IZ_TOP) > XKM(IX_TOP) .AND. &
                     ZKM(IZ_BOT) < XKM(IX_TOP)) THEN ! Case (3)
                 TAUG(IG) = TAUG(IG) + TAU(IL)*(XKM(IX_TOP) - ZKM(IZ_BOT))/ &
                                               (ZKM(IZ_TOP) - ZKM(IZ_BOT))
            ELSE IF (ZKM(IZ_TOP) > XKM(IX_BOT) .AND. &
                     ZKM(IZ_BOT) < XKM(IX_BOT)) THEN ! Case (4)
                 TAUG(IG) = TAUG(IG) + TAU(IL)*(ZKM(IZ_TOP) - XKM(IX_BOT))/ &
                                               (ZKM(IZ_TOP) - ZKM(IZ_BOT))
            END IF ! Cases 1-4
        END DO ! IL = 1, NL
    END DO ! IG = 1, NG-1
!
    TAUG(NG) = SUM(TAU) - SUM(TAUG(1:NG-1))
!
END SUBROUTINE GROUPTAU
!===============================================================================
! 22Apr16 - Minor changes in comments. Tested as part of SORD_IP vs IPRT cases
!           for 30 layers (Emde C. et al, JQSRT, 2015).
!
! 15Mar16 - Created & tested:
!           - Test 1 -
!              ZKM  = [28 21 15 10 6 3 1 0]
!              TAU  = [1 2 3 4 5 6 7]
!              XKM  = [28 14 7 2 1.5 0]
!              TAUG = [3.6 5.4 9 1.5 8.5] - ok
!
!           - Test 2 -
!              ZKM  = [10 0]
!              TAU  = 1
!              XKM  = [10 8 6 4 2 0]
!              TAUG = [0.2 0.2 0.2 0.2 0.2] - ok
!
!           - Test 3 -
!              ZKM  = [10 0]
!              TAU  = 1
!              XKM  = [10 8 2 0]
!              TAUG = [0.2 0.6 0.2] - ok
!
!           - Test 4 -
!              ZKM  = [28 21 15 10 6 3 1 0]
!              TAU  = [1 2 3 4 5 6 7]
!              XKM  = [28 0]
!              TAUG = sum(TAU) = 28 - ok
!===============================================================================