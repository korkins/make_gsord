SUBROUTINE ZRAYM_IP(M, DEPF, MUG, WG, NG2, MU, NMU, ZR11, ZR12, ZR13,          &
                                                    ZR21, ZR22, ZR23,          &
                                                    ZR31, ZR32, ZR33)
!===============================================================================
! PURPOSE:
!   To compute the Rayleigh phase matrix for a given Fourier order, M = 0, 1, 2.
!   On output, the phase matrix is weighted with Gauss weights for integration
!   over cosine of incident zenith angle, mu'.
!
! INPUT:
!   M      I(1)     The Fourier order, M = 0, 1, 2 (3 values)
!   DEPF   D(1)     Depolarization factor, DEPF = [0 ... ~6/7]
!   MUG    D(NG2)   Gauss nodes (incident rays)
!   WG     D(NG2)   Gauss weights
!   NG2    I(1)     Number of Gauss nodes per sphere
!   MU     D(NMU)   Directions of observation (scattered rays)
!   NMU    I(1)     Length of MU
!
! OUTPUT:
!   ZRij   D(NG2, NMU)   Elements of the phase matrix
!
! TREE:
!   -
!
! COMMENTS:
!   The phase matrix has symmetry with respect to Gauss nodes. The symmetry
!   has not been implemented in this subroutine at this point.
!
!   Theoretical background for the subroutine is given in [1, p.93,Section
!   3.4.4]. One should note that polynomials used in [1] and in this subroutine
!   differ by a factor of sqrt((k+m)!/(k-m)!). Nevertheless, the final result
!   for the phase matrix, W [1, p.92, Eq.(3.128)], coincide exactly with the
!   result provided by this subroutine.
!
!   The phase matrix is (element-by-element, mu' correspond to Gauss nodes)
!
!   ZRM11 = sum{  Qkm(mu) a1k Qkm(mu') }                                     (1)
!   ZRM21 = sum{ -Rkm(mu) b1k Qkm(mu') }
!   ZRM31 = sum{  Tkm(mu) b1k Qkm(mu') }
!
!   ZRM12 = sum{ -Qkm(mu) b1k Rkm(mu') }
!   ZRM22 = sum{  Rkm(mu) a2k Rkm(mu') + Tkm(mu) a3k Tkm(mu') }, a3k == 0 !!
!   ZRM32 = sum{ -Tkm(mu) a2k Rkm(mu') - Rkm(mu) a3k Tkm(mu') }, a3k == 0 !!
!
!   ZRM13 = sum{  Qkm(mu) b1k Tkm(mu') }
!   ZRM23 = sum{ -Rkm(mu) a2k Tkm(mu') - Tkm(mu) a3k Rkm(mu') }, a3k == 0 !!
!   ZRM33 = sum{  Tkm(mu) a2k Tkm(mu') + Rkm(mu) a3k Rkm(mu') }, a3k == 0 !!
!
!   where summations are taken over k = 0, 1, 2 (3 values). The expansion
!   moments are [1, p.56, Section 2.9]
!
!                  a1k = [1;  0;  D/2]                                       (2)
!                  a2k = [0;  0;  3D]                                        (3)
!                  a3k = [0;  0;  0]                                         (4)
!                  b1k = [0;  0;  (sqrt(6)/2)D]                              (5)
!
!   where
!
!   D = (1 - d)/(1 + d/2),                                                   (6)
!
!   and d is the depolarization factor on input; d = 0 (D = 1) corresponds to
!   pure Rayleigh scattering without depolarization.
!
!   d = Il(SCAT_ANG = 90)/Ir(SCAT_ANG = 90)                                  (7)
!
!   Range of the depolarization factor, DEPF, is [0..6/7] [1, p.56]. DEPF is
!   often assumed spectrally independent [2, p.3494].
!
!   Explicit values for polynomials involved in computations are
!
!   m = 0:                                                                   (8)
!       Q00 = 1; Q10 = mu; Q20 = (3*mu*mu - 1)/2
!       R00 = 0; R10 = 0; R20 = sqrt(3/8)*(1 - mu*mu)
!       Tk0 = 0 for all k.
!   m = 1:                                                                   (9)
!       Q01 = 0; Q11 = sqrt(1 - mu*mu)/2; Q21 = 3*mu*sqrt((1-mu*mu)/6)
!       R01 = 0; R11 = 0; R21 = -mu*sqrt(1 - mu*mu)/2
!       T01 = 0; T11 = 0; T21 = -sqrt(1 - mu*mu)/2
!   m = 2:                                                                  (10)
!       Q02 = 0; Q12 = 0; Q22 = 3*(1-mu*mu)/(2*sqrt(6))
!       R02 = 0; R12 = 0; R22 = (1 + mu*mu)/4
!       T02 = 0; T12 = 0; T22 = mu/2
!
! REFERENCES:
!   1. HovenierJW et al. Transfer of polarized light, Kluwer 2004
!   2. Bates DR, 1984, Planet. Space Sci., V32,pp.785–790
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: M, NG2, NMU
    REAL*8, INTENT(IN) :: DEPF
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NG2), INTENT(IN) :: MUG, WG
    REAL*8, DIMENSION(NMU), INTENT(IN) :: MU
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NG2, NMU), INTENT(OUT) :: ZR11, ZR12, ZR13, &
                                                ZR21, ZR22, ZR23, &
                                                ZR31, ZR32, ZR33
!
! LOCAL VARIABLES
    INTEGER &
        IG,  & ! Loop index over Gauss nodes
        IMU    ! Loop index over MU
    REAL*8 &
        A, B, C, D, E, F, & ! Temporary variables to compute Zij
        DF,               & ! Depolarization constant, Eq.(6)
        MU2,              & ! MU(IMU)*MU(IMU)
        MUG2,             & ! MUG(IG)*MUG(IG)
        X,                & ! Temporary variable 1, function of MU(IMU)
        Y,                & ! Temporary variable 2, function of MU(IMU)
        Z                   ! Temporary variable 3, function of MU(IMU)
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NG2) :: &
        DWG ! DF*WG*C where C = 1/8, 3/4, and 3/16 for m = 0, 1, 2, respectively
!===============================================================================
!
!   Depolarization constant, Eq.(6)
    DF = (1.0D0 - DEPF)/(1.0D0 + 0.5D0*DEPF)
!
    IF (M == 0) THEN        ! [1, p.93, Eq.(3.134)]
        DWG = 0.125D0*DF*WG ! 0.125 = 1/8
        DO IMU = 1, NMU
            MU2 = MU(IMU)*MU(IMU)
            X   = 3.0D0*MU2 - 1.0D0
            Y   = MU2 - 1.0D0
            DO IG = 1, NG2
                MUG2 = MUG(IG)*MUG(IG)
                A = DWG(IG)*X
                B = 3.0D0*DWG(IG)*Y
                C = 3.0D0*MUG2 - 1.0D0
                D = 3.0D0*(MUG2 - 1.0D0)
!
                ZR11(IG, IMU) = WG(IG) + A*C
                ZR21(IG, IMU) = B*C
                ZR31(IG, IMU) = 0.0D0
                ZR12(IG, IMU) = A*D
                ZR22(IG, IMU) = B*D
                ZR32(IG, IMU) = 0.0D0
                ZR13(IG, IMU) = 0.0D0
                ZR23(IG, IMU) = 0.0D0
                ZR33(IG, IMU) = 0.0D0
            END DO ! IMU = 1, NMU
        END DO ! IG = 1, NG2
    ELSEIF (M == 1) THEN   ! [1, p.94, Eq.(3.138)]
        DWG = 0.75D0*DF*WG ! 0.75 = 3/4
        DO IMU = 1, NMU
            MU2 = MU(IMU)*MU(IMU)
            X = DSQRT(1.0D0 - MU2)
            DO IG = 1, NG2
                MUG2 = MUG(IG)*MUG(IG)
                A =  DWG(IG)*X*DSQRT(1.0D0 - MUG2)
                B =  A*MU(IMU)*MUG(IG)
                C = -A*MUG(IG)
                D = -A*MU(IMU)
!
                ZR11(IG, IMU) = B
                ZR21(IG, IMU) = B
                ZR31(IG, IMU) = C
                ZR12(IG, IMU) = B
                ZR22(IG, IMU) = B
                ZR32(IG, IMU) = C
                ZR13(IG, IMU) = D
                ZR23(IG, IMU) = D
                ZR33(IG, IMU) = A
            END DO ! IMU = 1, NMU
        END DO ! IG = 1, NG2
    ELSEIF (M == 2) THEN     ! [1, p.94, Eq.(3.140)]
        DWG = 0.1875D0*DF*WG ! 0.1875 = 3/16
        DO IMU = 1, NMU
            MU2 = MU(IMU)*MU(IMU)
            X = 1.0D0 - MU2
            Y = 1.0D0 + MU2
            Z = 2.0D0*MU(IMU)
            DO IG = 1, NG2
                MUG2 = MUG(IG)*MUG(IG)
!
                A =  DWG(IG)*X
                B = -DWG(IG)*Y
                C =  DWG(IG)*Z
                D =  1.0D0 - MUG2
                E = -(1.0D0 + MUG2)
                F =  2.0D0*MUG(IG)
!
                ZR11(IG, IMU) = A*D
                ZR21(IG, IMU) = B*D
                ZR31(IG, IMU) = C*D
                ZR12(IG, IMU) = A*E
                ZR22(IG, IMU) = B*E
                ZR32(IG, IMU) = C*E
                ZR13(IG, IMU) = A*F
                ZR23(IG, IMU) = B*F
                ZR33(IG, IMU) = C*F
            END DO ! IMU = 1, NMU
        END DO ! IG = 1, NG2
    END IF ! M = 0
!
END SUBROUTINE ZRAYM_IP
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 27Sep15 - MUG(IG) is used instead of MU(JG); JG & NRU (input) are not used.
!           Some minor improvements: f(MU) are moved out of the IG-loop; X, Y,
!           and Z temporary variables are now used.
!
! 15May15 - Compared numerically vs previous version of SORD for the Rayleigh
!           scenario from Kokhanovsky et al., JQSRT, 2010.
!
! 13May15 - First created and compared analytically against [1, p.94].
!===============================================================================