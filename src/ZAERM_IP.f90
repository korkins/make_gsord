SUBROUTINE ZAERM_IP(M, XK, NK, WG, NG2, MU, NMU, NRU, ZP11, ZP12, ZP13,   &
                                                      ZP21, ZP22, ZP23,   &
                                                      ZP31, ZP32, ZP33)
!===============================================================================
! PURPOSE:
!   To compute a phase matrix of a group of particles for a given Fourier order,
!   M. On output, the phase matrix is weighted with Gauss weights for
!   integration over zenith angle.
!
! INPUT:
!   M      I(1)       The Fourier order, M = 0, 1, 2 (3 values)
!   XK     D(NK, 4)   Expansion moments: a1k, a2k, a3k, b1k
!   NK     I(1)       Total number of expansion moments
!   WG     D(NG2)     Gauss weights
!   NG2    I(1)       Number of Gauss nodes per sphere
!   MU     D(NMU)     Directions of scattering (observation)
!   NMU    I(1)       Length of MU
!   NRU    I(1)       Number of reflected user defined (dummy) nodes
!
! OUTPUT:
!   ZPij   D(NG2, NMU)   Elements of the phase matrix
!
! TREE:
!   ZPARM_IP
!          |
!          +-QX (...)
!          |
!          +-RTX (...)
!
! COMMENTS:
!   XK(1, 1) = 1. XK(k, :) include the factor of 2k+1.
!
!   Theoretical background for the subroutine is given in [1]. Note that
!   polynomials used in [1] and in this subroutine differ by a factor of
!   sqrt((k+m)!/(k-m)!). Nevertheless, the final result for the phase matrix, W
!   [1, p.92, Eq.(3.128)], coincide exactly with the result provided by this
!   subroutine.
!
!   The phase matrix is (element-by-element, mu' correspond to Gauss nodes)
!
!   ZPM11 = sum{  Qkm(mu) a1k Qkm(mu') }                                     (1)
!   ZPM21 = sum{ -Rkm(mu) b1k Qkm(mu') }
!   ZPM31 = sum{  Tkm(mu) b1k Qkm(mu') }, Tk0 == 0 !!
!
!   ZPM12 = sum{ -Qkm(mu) b1k Rkm(mu') }
!   ZPM22 = sum{  Rkm(mu) a2k Rkm(mu') + Tkm(mu) a3k Tkm(mu') }, Tk0 == 0 !!
!   ZPM32 = sum{ -Tkm(mu) a2k Rkm(mu') - Rkm(mu) a3k Tkm(mu') }, Tk0 == 0 !!
!
!   ZPM13 = sum{  Qkm(mu) b1k Tkm(mu') }, Tk0 == 0 !!
!   ZPM23 = sum{ -Rkm(mu) a2k Tkm(mu') - Tkm(mu) a3k Rkm(mu') }, Tk0 == 0 !!
!   ZPM33 = sum{  Tkm(mu) a2k Tkm(mu') + Rkm(mu) a3k Rkm(mu') }, Tk0 == 0 !!
!
!   where summations are taken over k = 0, 1, ... NK-1.
!
! REFERENCES:
!   1. HovenierJW et al. Transfer of polarized light, Kluwer 2004
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: M, NK, NG2, NMU, NRU
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NK, 4), INTENT(IN) :: XK
    REAL*8, DIMENSION(NG2), INTENT(IN) :: WG
    REAL*8, DIMENSION(NMU), INTENT(IN) :: MU
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NG2, NMU), INTENT(OUT) :: ZP11, ZP12, ZP13, &
                                                ZP21, ZP22, ZP23, &
                                                ZP31, ZP32, ZP33
!
! LOCAL VARIABLES
    INTEGER &
        IMU, & ! Loop index over MU
        IG     ! Loop index over Gauss nodes
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NK) :: &
        A1Q, B1R, B1T, B1Q, & ! Product of a polynomial & expansion moment, e.g.
        A2R, A3T, A2T, A3R    ! A1Q = a1k*Qkm(mu) for all k, mu, and a given m
    REAL*8, DIMENSION(NMU, NK) :: &
        Q, R, T    ! Polynomials
    REAL*8, DIMENSION(NK, NMU) :: &
        QT, RT, TT ! Transposed QRT for fast summation over k = 0:NK-1
!===============================================================================
!
!   Polynomials for all directions
    CALL QX( M, MU, NMU, NK, Q )
    CALL RTX( M, MU, NMU, NK, R, T )
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
            DO IG = 1, NG2
                A1Q = QT(:, NRU+IG)*XK(:, 1)
                A2R = RT(:, NRU+IG)*XK(:, 2)
                A2T = TT(:, NRU+IG)*XK(:, 2)
                A3T = TT(:, NRU+IG)*XK(:, 3)
                A3R = RT(:, NRU+IG)*XK(:, 3)
                B1R = RT(:, NRU+IG)*XK(:, 4)
                B1T = TT(:, NRU+IG)*XK(:, 4)
                B1Q = QT(:, NRU+IG)*XK(:, 4)
!
                ZP11(IG, IMU) =  WG(IG)*SUM( QT(:, IMU)*A1Q )
                ZP12(IG, IMU) = -WG(IG)*SUM( QT(:, IMU)*B1R )
                ZP13(IG, IMU) =  0.0D0
                ZP21(IG, IMU) = -WG(IG)*SUM( RT(:, IMU)*B1Q )
                ZP22(IG, IMU) =  WG(IG)*SUM( RT(:, IMU)*A2R )
                ZP23(IG, IMU) =  0.0D0
                ZP31(IG, IMU) =  0.0D0
                ZP32(IG, IMU) =  0.0D0
                ZP33(IG, IMU) =  WG(IG)*SUM( RT(:, IMU)*A3R )
            END DO ! IG = 1, NG2
        END DO ! IMU = 1, NMU
    ELSE ! M > 0
        DO IMU = 1, NMU       
            DO IG = 1, NG2
                A1Q = QT(:, NRU+IG)*XK(:, 1)
                A2R = RT(:, NRU+IG)*XK(:, 2)
                A2T = TT(:, NRU+IG)*XK(:, 2)
                A3T = TT(:, NRU+IG)*XK(:, 3)
                A3R = RT(:, NRU+IG)*XK(:, 3)
                B1R = RT(:, NRU+IG)*XK(:, 4)
                B1T = TT(:, NRU+IG)*XK(:, 4)
                B1Q = QT(:, NRU+IG)*XK(:, 4)
!
                ZP11(IG, IMU) =  WG(IG)*SUM( QT(:, IMU)*A1Q )
                ZP12(IG, IMU) = -WG(IG)*SUM( QT(:, IMU)*B1R )
                ZP13(IG, IMU) =  WG(IG)*SUM( QT(:, IMU)*B1T )
                ZP21(IG, IMU) = -WG(IG)*SUM( RT(:, IMU)*B1Q )
                ZP22(IG, IMU) =  WG(IG)*( SUM( RT(:, IMU)*A2R ) +              &
                                          SUM( TT(:, IMU)*A3T ) )
                ZP23(IG, IMU) = -WG(IG)*( SUM( RT(:, IMU)*A2T ) +              &
                                          SUM( TT(:, IMU)*A3R ) )
                ZP31(IG, IMU) =  WG(IG)*SUM( TT(:, IMU)*B1Q )
                ZP32(IG, IMU) = -WG(IG)*( SUM( TT(:, IMU)*A2R ) +              &
                                          SUM( RT(:, IMU)*A3T ) )
                ZP33(IG, IMU) =  WG(IG)*( SUM( TT(:, IMU)*A2T ) +              &
                                          SUM( RT(:, IMU)*A3R ) )
            END DO ! IG = 1, NG2
        END DO ! IMU = 1, NMU
    END IF ! M = 0
!
END SUBROUTINE ZAERM_IP
!===============================================================================
! 09Nov16 - MUG removed from input as redundant
!
! 01May16 - JG is replaced with NRU+IG
!
! 22Apr16 - Minor changes in comments
!
! 27Sep15 - DDOT5(A, B, N) is replaced with SUM(A(:)*B(:))
!
! 17Jul15 - Renamed to ZAERM_IP
!
! 16May15 - Tested for Test 5 (aerosol)
!
! 14May15 - First created
!===============================================================================