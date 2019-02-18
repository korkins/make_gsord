SUBROUTINE ZRAYM0_IP(M, DEPF, MU0, MU, NMU, ZR01, ZR02, ZR03)
!===============================================================================
! PURPOSE:
!   To compute the 1st column of the Rayleigh phase matrix for a given Fourier
!   order M = 0, 1, 2, one incident and several scattering directions, MU0 and
!   MU, respectively.
!
! INPUT:
!   M      I(1)     The Fourier order, M = 0, 1, 2 (3 values)
!   DEPF   D(1)     Depolarization factor, DEPF = [0 ... ~6/7]
!   MU0    D(1)     cos(SZA) defines incident direction
!   MU     D(NMU)   cos(VZA) defines directions of scattering (observation)
!   NMU    I(1)     Length of MU
!
! OUTPUT:
!   ZR01, ZR02, ZR03   D(NMU)   The 1st column of the Rayleigh phase matrix
!
! TREE:
!   -
!
! COMMENTS:
!   Theoretical background for this subroutine is described in [1, p.93, Section
!   3.4.4]. Note that polynomials used in [1] and in this subroutine differ by a
!   factor of sqrt((k+m)!/(k-m)!). Nevertheless, the final result for the phase
!   matrix, W [1, p.92, Eq.(3.128)], coincide exactly with the result of this
!   subroutine.
!
!   The 1st column of the phase matrix is
!
!               | Qkm(mu) a1k Qkm(mu0)|
!   ZR0M = sum( |-Rkm(mu) b1k Qkm(mu0)| ),                                   (1)
!               | Tkm(mu) b1k Qkm(mu0)|
!
!   where summation is taken over k = 0, 1, 2 (3 values). The expansion moments
!   are [1, p.56, Section 2.9]
!
!   a1k = [1  0   D/2],                                                      (2)
!   b1k = [0  0  (sqrt(6)/2)D],                                              (3)
!
!   where
!
!   D = (1 - d)/(1 + d/2),                                                   (4)
!
!   and d is the depolarization factor on input; d = 0 (D = 1) corresponds to
!   pure Rayleigh scattering without depolarization.
!
!   d = Il(SCAT_ANG = 90)/Ir(SCAT_ANG = 90)                                  (5)
!
!   Range of the depolarization factor, DEPF, is [0...6/7] [1, p.56]. DEPF is
!   often assumed spectrally independent [2, p.3494].
!
!   Explicit values for polynomilas involved in computations are
!
!   m = 0:                                                                   (6)
!       Q00 = 1; Q10 = mu; Q20 = (3*mu*mu - 1)/2
!       R00 = 0; R10 = 0; R20 = sqrt(3/8)*(1 - mu*mu)
!       Tk0 = 0 for all k.
!   m = 1:                                                                   (7)
!       Q01 = 0; Q11 = sqrt(1 - mu*mu)/2; Q21 = 3*mu*sqrt((1-mu*mu)/6)
!       R01 = 0; R11 = 0; R21 = -mu*sqrt(1 - mu*mu)/2
!       T01 = 0; T11 = 0; T21 = -sqrt(1 - mu*mu)/2
!   m = 2:                                                                   (8)
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
    INTEGER, INTENT(IN) :: M, NMU
    REAL*8, INTENT(IN) :: DEPF, MU0
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NMU), INTENT(IN) :: MU
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NMU), INTENT(OUT) :: ZR01, ZR02, ZR03
!
! LOCAL VARIABLES
    REAL*8 &
        DF ! Eq.(4)
!===============================================================================
!
!   Depolarization constant
    DF = (1.0D0 - DEPF)/(1.0D0 + 0.5D0*DEPF)
!
    IF (M == 0) THEN       ! [1, p.93, Eq.(3.134)]
        ZR01 =  1.0D0 + DF*(3.0D0*MU*MU - 1.0D0)*(3.0D0*MU0*MU0 - 1.0D0)/8.0D0
        ZR02 = -0.375D0*DF*(1.0D0 - MU*MU)*(3.0D0*MU0*MU0 - 1.0D0)!3/8=0.375
        ZR03 =  0.0D0
    ELSEIF (M == 1) THEN   ! [1, p.94, Eq.(3.138)]
        ZR03 =  0.75D0*DF*MU0*DSQRT((1.0D0 - MU*MU)*(1.0D0 - MU0*MU0)) !3/4=0.75
        ZR01 =  ZR03*MU
        ZR02 =  ZR01
        ZR03 = -ZR03
    ELSEIF (M == 2) THEN   ! [1, p.94, Eq.(3.140)]
        ZR03 =  0.1875D0*DF*(1.0D0 - MU0*MU0) !3/16=0.1875
        ZR01 =  ZR03*(1.0D0 - MU*MU)
        ZR02 = -ZR03*(1.0D0 + MU*MU)
        ZR03 =  2.0D0*MU*ZR03
    END IF ! M = 0
!
END SUBROUTINE ZRAYM0_IP
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 17Jul15 - Renamed to ZRAYM0_IP
!
! 15May15 - Compared numerically vs previous version of SORD for the Rayleigh
!           scenario from Kokhanovsky et al., JQSRT, 2010.
!
! 13May15 - First created and compared analyticaly against [1, p.94].
!===============================================================================