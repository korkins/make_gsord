SUBROUTINE SPLITTAU(DTAU, NX, TAU, NL, ZKM, NZ, XKM)
!===============================================================================
! PURPOSE:
!   To compute heights, XKM, splitting the input Tau(z) into dTau-thick pieces
!
! INPUT:
!   DTAU  D(1)    Thickness between XKM(IX) & XKM(IX+1) (Increment)
!   NX    I(1)    Number of boundaries at the dTau-layers: sum(TAU)=dTau*(NX-1)
!   TAU   D(NL)   Thicknesses of each input layer
!   NL    I(1)    Number of input layers, NZ-1
!   ZKM   D(NZ)   Heights at Tau, from TOA to BOA: ZKM(I) > ZKM(I+1)
!   NZ    I(1)    Number of input heights, NL+1
!
! OUTPUT:
!   XKM   D(NX)   Heights, TOA to BOA, which split TAU with increment DTAU.
!                 Note, XKM(1) = ZKM(1) & XKM(NB) = ZKM(NZ)
!
! TREE:
!   -
!
! COMMENTS:
!   On input, the user must satisfy the following condition:
!
!       DTAU*(NX-1) = SUM(DTAU),
!
!   e.g. by doing DTAU=TAU0/CEILING(TAU0/DTAU). Here TAU0=SUM(TAU) and
!   CEILING(X) returns the least integer .ge. X [1].
!
!   There are only two possible cases. Case A: the current layer TAU(IL=IZ-1)
!   fits in the residual thickness entirely. Crop TAU(IL) from the residual, get
!   XKM(IX) and proceed to the next IX. Case B: there is at least one layer,
!   TAU(IL), between the known XKM(IX-1) and the next XKM(IX) which is yet to be
!   determined. Accumulate layer(s), TAU(IL), in this case. Note, the case
!   XKM(IX+1) = ZKM(IZ) falls under case A.
!
! REFERENCES:
!   1. https://gcc.gnu.org/onlinedocs/gfortran/CEILING.html#CEILING
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NX, NL, NZ
    REAL*8, INTENT(IN) :: DTAU
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NL), INTENT(IN) :: TAU
    REAL*8, DIMENSION(NZ), INTENT(IN) :: ZKM
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NX), INTENT(OUT) :: XKM
!
! LOCAL VARIABLES
    INTEGER &
        IX, & ! Loop index over XKM
        IZ    ! Loop index over ZKM
    REAL*8 &
        EXT,  & ! Extinction in the layer IL=IZ-1
        TACC, & ! Accumulated Tau
        TCOM, & ! Complimentary thickness, DTAU - TACC
        TRES    ! Residual thickness for the current layer IL=IZ-1
!===============================================================================
!
!   TOA
    XKM(1) = ZKM(1)
!
!   Internal XKM, if any
    IZ = 2
    EXT = TAU(1)/(ZKM(1) - ZKM(2))
    DO IX = 2, NX-1
        TRES = EXT*(XKM(IX-1) - ZKM(IZ))
        IF (TRES >= DTAU) THEN
!           Case A:
!               1. Residual of the layer IZ-1 is thick enough to hold dTau;
!               2. XKM(IX) belongs to the layer IZ-1, i.e. XKM(IX)>=ZKM(IZ);
!               3. Do not move ZKM(IZ) in this case;
!               4. Proceed to next XKM(IX) from the loop.
            XKM(IX) = XKM(IX-1) - DTAU/EXT
        ELSE ! TRES < DTAU
!           Case B:
!               1. Residual of the layer IZ-1 is thin;
!               2. Move ZKM(IZ) in this case;
!               3. Combine with the next layer(s) until thick enough;
!               4. XKM(IX) belongs to the layer that violates TACC < DTAU;
!               5. Complementary thickness, TCOM, is cropped from this layer;
!               6. TAU(IL=IZ-1) - TCOM becomes the residual thickness, TRES;
!               7. Extinction, EXT, is redefined:
!                                        EXT = TAU(IZ-1)/( ZKM(IZ-1) - ZKM(IZ) );
!               8. Proceed to next XKM(IX) from the loop.
            TACC = TRES
            DO WHILE (TACC < DTAU)
                IZ = IZ+1               ! Take the next layer ...
                TACC = TACC + TAU(IZ-1) ! ... and keep accumulating
            END DO ! WHILE (TACC < DTAU)
            TACC = TACC - TAU(IZ-1)  ! Overaccumulation: go one step back
            TCOM = DTAU - TACC       ! Out of TAU(IZ-1) only TCOM is needed
            EXT = TAU(IZ-1)/(ZKM(IZ-1) - ZKM(IZ)) ! New extinction for Case A
            XKM(IX) = ZKM(IZ-1) - TCOM/EXT
        END IF ! TRES >= DTAU
    END DO ! IX = 2, NX-1
!
!   BOA
    XKM(NX) = ZKM(NZ)
!
END SUBROUTINE SPLITTAU
!===============================================================================
! 22Apr16 - Minor changes in comments. Tested as part of SORD against IPRT tests
!           (Emde C et al, JQSRT 2015, Cases B2 & B3)
!
! 18Mar16 - Created & tested:
!
!           - Common input for Tests 1:4 -
!             DTAU = 1.0D0
!             TAU0 = SUM(TAU) = 10.0
!             NL = 10, NZ = NZ = 11
!             NX = CEILING(TAU0/DTAU)+1 = 11
!             ZKM = [10 9 8 7 6 5 4 3 2 1 0], NZ = 11
!
!           - Test 1 -
!             TAU = [0.1 0.2 0.3 0.4 1.0 >> 4.0 << 2.0 1.0 0.5 0.5]
!             XKM = [10.0 6.0 5.0 4.75 4.5 4.25 4.0 3.5 3.0 2.0 0.0] - ok
!
!           - Test 2 -
!             TAU = [4.0 << 0.1 0.2 0.3 0.4 1.0 2.0 1.0 0.5 0.5]
!             XKM = [10.0 9.75 9.5 9.25 9.0 5.0 4.0 3.5 3.0 2.0 0.0] - ok
!
!           - Test 3 -
!             TAU = [0.1 0.2 0.3 0.4 1.0 2.0 1.0 0.5 0.5 >> 4.0 <<]
!             XKM = [10.0 6.0 5.0 4.5 4.0 3.0 1.0 0.75 0.5 0.25 0.0] - ok
!
!           - Test 4 -
!             TAU = [0.1 0.2 0.3 >> 4.0 << 0.4 1.0 0.5 2.0 1.0 0.5]
!             XKM = [10.0 6.9 6.65 6.40 6.15 5.0 4.0 2.75 2.25 1.5 0.0] - ok
!
!           - Test 5 -
!             Similar to Test 1, except for
!             ZKM = [30 25 21 18 16 15                  10   6 3 1 0], NZ = 11
!             XKM = [30          16 15 13.75 12.5 11.25 10 8 6 3   0] - ok
!
!             See Test 1:
!             ZKM = [10 9 8 7 6 5               4     3 2 1 0], NZ = 11
!             XKM = [10       6 5 4.75 4.5 4.25 4 3.5 3 2   0] - ok
!
!           - Test 6 -
!             Similar to Test 4, except for
!             ZKM = [30 25 21 18 16 15 10 6 3 1 0], NZ = 11
!             XKM = [30             15 10             0] - ok (for shown values)
!
!             See Test 4:
!             ZKM = [10 9   8    7    6    5 4 3    2    1   0], NZ = 11
!             XKM = [10 6.9 6.65 6.40 6.15 5 4 2.75 2.25 1.5 0] - ok
!
!           For Tests 7-9, we start from the solution, XKM, dTau, & Tau. From
!           that we derives input ZKM, compute XKM' and compare XKM' and XKM,
!           which must exactly coincide.
!
!           - Test 7 -
!             Solution: NZ < NX
!             XKM = [10 9 8 7 6 5 4 3 2 1 0], dTau = 0.5
!             Tau = [1.5 1.0 2.5]
!             ZKM = [10 7 5 0] -> XKM = [10 9 8 7 6 5 4 3 2 1 0] - ok
!
!           - Test 8 -
!             Solution: NZ = NX
!             TAU = [dTau dTau dTau ... dTau]
!             ZKM = XKM - ok
!
!           - Test 9 -
!             Solution: NZ > NX
!             XKM = [10 9 8 7 6 5 4 3 2 1 0], dTau = 0.5
!             Tau = [dTau/2 dTau/2 ... dTau/2]
!             ZKM = [10: -0.5 : 0] -> XKM = [10 9 8 7 6 5 4 3 2 1 0] - ok
!
!           - Test 10 -
!             NL = 1, NZ = 2, Z = [10 7], Tau = 1
!             dTau = 1/3, NX = 4
!             XKM = [10 9 8 7] - ok
!
!           - Test 11 -
!             NL = 3, NZ = 4, Z = [10 9 8 7], Tau = [1 1 1]
!             dTau = 10 > sum(Tau), NX = 2, XKM = [10 7] - ok
!===============================================================================