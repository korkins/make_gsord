SUBROUTINE ZAERM0_IP(M, MU0, A1K, B1K, NK, MU, NMU, ZP01, ZP02, ZP03)
!===============================================================================
! PURPOSE:
!   To compute the 1st column of the phase matrix for a given Fourier order M,
!   one incident and several scattering directions, MU0 and MU, respectively.
!
! INPUT:
!   M     I(1)     The Fourier order, M = 0, 1, 2 (3 values)
!   MU0   D(1)     cos(SZA) defines incident direction
!   A1K   D(NK)    Phase function expansion moments
!   B1K   D(NK)    Off-diagonal expansion moments
!   NK    I(1)     Total number of expansion moments
!   MU    D(NMU)   Directions of scattering (observation)
!   NMU   I(1)     Length of MU
!
! OUTPUT:
!   ZP01, ZP02, ZP03   D(NMU)   1st column of the phase matrix
!
! TREE:
!   -
!
! COMMENTS:
!   A1K(1) = 1. A1K & B1K include the factor of 2k+1.
!
!   The 1st column of the phase matrix is [1]
!
!               | Qkm(mu) a1k Qkm(mu0)|
!   ZP0M = sum( |-Rkm(mu) b1k Qkm(mu0)| ),                                   (1)
!               | Tkm(mu) b1k Qkm(mu0)|
!
!   where summation is taken over k = 0, 1, 2, ... NK-1
!
! REFERENCES:
!   1. HovenierJW et al. Transfer of polarized light, Kluwer 2004
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: M, NK, NMU
    REAL*8, INTENT(IN) :: MU0
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NK), INTENT(IN) :: A1K, B1K
    REAL*8, DIMENSION(NMU), INTENT(IN) :: MU
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NMU), INTENT(OUT) :: ZP01, ZP02, ZP03
!
! LOCAL VARIABLES
    INTEGER &
        IK, & ! Loop index over order of polynomial
        IMU   ! Loop index over MU
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NK) :: &
        Q0         ! Polynomials for all orders k and a single direction mu0
    REAL*8, DIMENSION(NMU, NK) :: &
        Q, R, T    ! Polynomials for all orders k and all directions mu
    REAL*8, DIMENSION(NK, NMU) :: &
        QT, RT, TT ! Transposed QRT for fast summation over k = 0:NK-1
!===============================================================================
!
!   Polynomials for all directions
    CALL QX0(M, MU0, NK, Q0)
    CALL QX(M, MU, NMU, NK, Q)
    CALL RTX(M, MU, NMU, NK, R, T)
!
!   Transposed to allow for fast summation over k (order of polynomial)
!   Symmetry can be applied !!
    DO IMU = 1, NMU
        QT(:, IMU) = Q(IMU, :) ! QT = TRANSPOSE(Q)
        RT(:, IMU) = R(IMU, :) ! RT = TRANSPOSE(R)
        TT(:, IMU) = T(IMU, :) ! TT = TRANSPOSE(T)
    END DO
!
    IF (M == 0) THEN ! Tkm == 0
        DO IMU = 1, NMU
!           IK = 1
            ZP01(IMU) =  QT(1, IMU)*A1K(1)*Q0(1)
            ZP02(IMU) = -RT(1, IMU)*B1K(1)*Q0(1)
            ZP03(IMU) =  0.0D0
            DO IK = 2, NK
                ZP01(IMU) = ZP01(IMU) + QT(IK, IMU)*A1K(IK)*Q0(IK)
                ZP02(IMU) = ZP02(IMU) - RT(IK, IMU)*B1K(IK)*Q0(IK)
            END DO ! IK = 1, NK
        END DO ! IMU = 1, NMU
    ELSE ! M > 0
        DO IMU = 1, NMU
!           IK = 1
            ZP01(IMU) =  QT(1, IMU)*A1K(1)*Q0(1)
            ZP02(IMU) = -RT(1, IMU)*B1K(1)*Q0(1)
            ZP03(IMU) =  TT(1, IMU)*B1K(1)*Q0(1)
            DO IK = 2, NK
                ZP01(IMU) = ZP01(IMU) + QT(IK, IMU)*A1K(IK)*Q0(IK)
                ZP02(IMU) = ZP02(IMU) - RT(IK, IMU)*B1K(IK)*Q0(IK)
                ZP03(IMU) = ZP03(IMU) + TT(IK, IMU)*B1K(IK)*Q0(IK)
            END DO ! IK = 1, NK
        END DO ! IMU = 1, NMU
    END IF ! M = 0
!
END SUBROUTINE ZAERM0_IP
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 17Jul15 - Renamed to ZAERM0_IP
!
! 16May15 - Tested for Test 5 (aerosol)
!
! 14May15 - First created
!===============================================================================