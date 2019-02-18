SUBROUTINE SUMKA1B1(A1K, B1K, NK, MUS, NM, NA, A1, B1)
!===============================================================================
! PURPOSE:
!   To compute two elements of the scattering matrix, a1 and b1, using their
!   moments a1k and b1k.
!
! INPUT:
!   A1K   D(NK)       Moments of a1, including the (2k+1) factor
!   B1K   D(NK)       Moments of b1, including the (2k+1) factor
!   NK    I(1)        Total number of moments, k = 0, 1, 2, ... NK-1
!   MUS   D(NM, NA)   Cosine of scattering angle
!   NM    I(1)        Number of zeniths
!   NA    I(1)        Number of azimuths
!
! OUTPUT:
!   A1   D(NM, NA)   Element a1
!   B1   D(NM, NA)   Element b1
!
! TREE:
!   SUMKA1B1
!          |
!          +-QX (Q polynomial)
!
! COMMENTS:
!   This subroutine was created from the IPOL's SCATMATAB.F90 subroutine.
!
!   On input, A1K and B1K are scaled by a factor of (2k+1), k=0...K-1.
!
!   B1K(1) = B1K(2) == 0
!
!   The scattering matrix [1, p.49, Eq.(2.135)] and the matrix of moments
!   [similar bot not exactly equal to 1, p.91, Eq.(3.122)] are defined as:
!
!               |a1 b1  0  0 |
!               |b1 a2  0  0 |
!           X = |0  0   a3 b2|                                               (1)
!               |0  0  -b2 a4|
!
!                | a1k -b1k  0    0  |
!                |-b1k  a2k  0    0  |
!           Xk = | 0    0    a3k -b2k|                                       (2)
!                | 0    0    b2k  a4k|
!
!   Note that Xk, Eq.(2), contains the expansion moments of Eq.(1) over matrix
!   polynomials PI(mu) [1, p.91, Eq.(3.125)]. On input, Eq.(2) b1k is defined
!   WITHOUT "-".
!
!   Eqs.(2.152), (2.153), (2.156)-(2.159) [1, pp.53, 54] are used. Note that
!   the [(l-2)!/(l+2)!]^1/2 is included in the Q polynomial computed in the
!   QRTm subroutine. Making use of the following relation [1, Eq.(B.20),
!   p.209]
!
!           Pk20 = (-1)*Qk2                                                  (3)
! 
!   we write the following expressions that are used in this subroutine
!
!           a1(mus) =  sum(a1k*Qk0(mus), k = 0..K)                           (4)
!           b1(mus) = -sum(b1k*Qk2(mus), k = 2..K)                           (5)
!
! REFERENCES:
!   1. Hovenier JW, van der Mee C, Domke H. Transfer of polarized light in
!      planetary atmospheres. Basic concepts and practical methods. Kluwer, 2004
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NK, NM, NA
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NK), INTENT(IN) :: A1K, B1K
    REAL*8, DIMENSION(NM, NA), INTENT(IN) :: MUS
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NM, NA), INTENT(OUT) :: A1, B1
!
! LOCAL VARIABLES
    INTEGER &
        IA, & ! Loop index over azimuth
        IK    ! Loop index over expansion moments, k=0:NK-1 -> IK=1:NK
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NM, NK) :: &
        Q ! Q-polynomial
!===============================================================================
!
    DO IA = 1, NA
!
!       m = 0
        CALL QX(0, MUS(:, IA), NM, NK, Q)
        DO IK = 1, NK
            IF (IK == 1) THEN
                A1(:, IA) = A1K(IK)*Q(:, IK)    ! Initialize according to Eq.(4)
            ELSE ! IF IK == 1
                A1(:, IA) = A1(:, IA) + A1K(IK)*Q(:, IK) ! Accumulate A1, Eq.(4)
            END IF ! IK == 1
        END DO ! IK = 1, NK
!
!       m = 2
        CALL QX(2, MUS(:, IA), NM, NK, Q)
        DO IK = 3, NK
            IF (IK == 3) THEN
                B1(:, IA) = -B1K(IK)*Q(:, IK)   ! Initialize according to Eq.(5)
            ELSE ! IF IK == 3
                B1(:, IA) = B1(:, IA) - B1K(IK)*Q(:, IK) ! Accumulate B1, Eq.(5)
            END IF ! IK == 3
        END DO ! IK = 3, NK
    END DO ! IA = 1, NA
!
END SUBROUTINE SUMKA1B1
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 30Nov15 - Renamed: A1B1FROMXK -> SUMKA1B1
!
! 10Apr15 - QPOL is replaced with QX
!
! 08Sep14 - Computation of cosine of scattering angle was moved to a separate
!           subroutine, CSCATANG. Test from 07Sep14 was repeated
!
! 07Sep14 - Created and tested as part of SORD for the phase matrix from
!           Wauben et al, 1993: Astron. Astroph. 276, 589-602.
!===============================================================================