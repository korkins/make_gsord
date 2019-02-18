SUBROUTINE SSCAT1M_IP(M, MU0, TAU, ZM1, ZM2, ZM3, MU, NMU, NUP, NDN, &
                      I1M, Q1M, U1M)
!===============================================================================
! PURPOSE:
!   To compute reflected and transmitted path radiance in the single scattering
!   approximation for a given Fourier order M, single layer and one direction of
!   the Solar beam. The incident Stokes vector is assumed to be [2pi 0 0 0].
!
! INPUT:
!   M     I(1)     The Fourier order, M = 0 ... NM-1
!   MU0   D(1)     Cosine of the solar zenith angle, mu0 = cos(SZA) > 0
!   TAU   D(1)     Optical thickness of the layer
!   ZMi   D(NMU)   Fourier moment of the phase matrix (1st column) at all MU
!   MU    D(NMU)   Cos(VZA), VZA /= 90. On TOA: mu < 0. On BOA: mu > 0
!   NMU   I(1)     Number of view directions. Length of MU. NMU = NUP + NDN > 0
!   NUP   I(1)     Number of directions with mu < 0, if any
!   NDN   I(1)     Number of directions with mu > 0, if any
!
! OUTPUT:
!   I1M, Q1M, U1M   D(NMU)   Components of the Stokes vector, normalized by 2pi
!
! TREE:
!   SSCAT1M_IP
!            |
!            +-QX0 (Polynomials Qkm(mu0))
!            |
!            +-QX (Polynomials Qkm(mu))
!            |
!            +-RTX (Polynomials Rkm(mu), Tkm(mu))
!
! COMMENTS:
!   IMPORTANT 1: on input, A1K and B1K contain the factor of (2k+1)SSA/2, where
!   SSA is the single scattering albedo of the layer.
!
!   IMPORTANT 2: on output, I1M, Q1M, and U1M are normalized to incident beam
!   Io = [2pi 0 0 0].
!
!   IMPORTANT 3: mu < 0, if any, MUST come first, followed by mu > 0, if any.
!
!   The system of coordinates used in this subroutine is exactly the same as
!   defined in [1].
!
!   The following expressions are used in this subroutine:
!
!   I = 1/2 * SSA * Z(:, 1) * G                                              (1)
!
!   where
!
!   G = mu0/(mu0-mu)*[1 - exp{(mu0-mu)*Tau/(mu0*mu)}], mu < 0                (2)
!   G = mu0/(mu0-mu)*[exp(-Tau/mu0) - exp(-Tau/mu)],  mu > 0, mu /= mu0      (3)
!   G = Tau/mu * exp(-Tau/mu),  mu = mu0                                     (4)
!
!   The phase matrix and the expansion moments are (note the "-" signs!!)
!
!        |a1  b1   0   0|
!        |b1  a2   0   0|
!   F  = | 0   0  a3  b2|,                                                   (5)
!        | 0   0 -b2  a4|
!
!        | a1k -b1k  0    0  |
!        |-b1k  a2k  0    0  |
!   Fk = | 0    0    a3k -b2k|,                                              (6)
!        | 0    0    b2k  a4k|
!
!   and matrix polynomilas are [2]
!
!                    | Qkm(mu)  0        0        0      |
!                    | 0        Rkm(mu) -Tkm(mu)  0      |
!   Pki = Pkm(mui) = | 0       -Tkm(mu)  Rkm(mu)  0      |,                  (7)
!                    | 0        0        0        Qkm(mu)|
!
!   where QRT polynomials are defined following [1, 2, see QRTm.f for details].
!
!   Note that on input, b1k and b2k come WITHOUT the "-" sign!! For instance,
!   the Rayleigh scattering is described by (note the "+" signs!!)
!
!   k=0: |1 0 0 0| k=1: |0 0 0 0  | k=2: | 1/2       +sqrt(6)/2 0 0|
!        |0 0 0 0|      |0 0 0 0  |      |+sqrt(6)/2  3         0 0|
!        |0 0 0 0|      |0 0 0 0  |      | 0          0         0 0|         (8)
!        |0 0 0 0|      |0 0 0 3/2|      | 0          0         0 0|
!
!   Finally, the first column of Pkm(mu)*Fk*Pkm(mu0), omitting the 4th element
!   which equals zero, is
!
!   Z(:, 1) = sum{Zk(:, 1), k = 0..(KN-1)},                                  (9)
!
!   where (note the sign at the 2nd element)
!
!              |+Qkm(mu)*a1k*Qkm(mu0)|
!   Zk(:, 1) = |-Rkm(mu)*b1k*Qkm(mu0)|.                                     (10)
!              |+Tkm(mu)*b1k*Qkm(mu0)|
!
! REFERENCES:
!   1. Hovenier JW et al., 2004: Transfer of Polarized Light in Planetary
!      Atmosphere. Basic Concepts and Practical Methods, Dordrecht: Kluwer
!      Academic Publishers;
!   2. Siewert CE, JQSRT(2000), V.64, P.227.
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    REAL*8, PARAMETER :: TINY = 1.0D-8   ! A small number to compare doubles
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: M, NMU, NUP, NDN
    REAL*8, INTENT(IN) :: TAU, MU0
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NMU), INTENT(IN) :: MU, ZM1, ZM2, ZM3
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NMU), INTENT(OUT) :: I1M, Q1M, U1M
!
! LOCAL VARIABLES
    INTEGER &
        IMU ! Loop index for mu > 0: IMU = 1, NDN
    REAL*8 &
        DMU, & ! DMU = DMUDN(IMU) = MU0 - MUDN(IMU)
        TMU    ! Tau/mu in case of mu -> mu0
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NDN) :: &
        DMUDN, & ! mu0 - mu for mu > 0
        GDN, &   ! Geometric factor for mu > 0
        MUDN     ! mu > 0
    REAL*8, DIMENSION(NUP) :: &
        DMUUP, & ! mu0 - mu for mu < 0
        GUP,   & ! Geometric factor for mu < 0
        MUUP     ! mu < 0
!===============================================================================
!
!   Process upward directions, if any
    IF (NUP /= 0) THEN
        MUUP = MU(1:NUP)
        DMUUP = MU0 - MUUP
        GUP = MU0/DMUUP*( 1.0D0 - DEXP(TAU*DMUUP/(MU0*MUUP)) )
        I1M(:NUP) = GUP*ZM1(:NUP)
        Q1M(:NUP) = GUP*ZM2(:NUP)
        IF (M > 0) THEN
            U1M(:NUP) = GUP*ZM3(:NUP)
        ELSE ! M = 0
            U1M(:NUP) = 0.0D0
        END IF
    END IF ! NUP /= 0
!
!   Process downward directions, if any
    IF (NDN /= 0) THEN
        MUDN = MU(NUP+1:NMU)
        DMUDN = MU0 - MUDN
        DO IMU = 1, NDN
            DMU = DMUDN(IMU)
            IF (DABS(DMU) < TINY) THEN
!               Singular point: mu -> mu0
                TMU = TAU/MUDN(IMU)
                GDN(IMU) = TMU*DEXP(-TMU)
            ELSE ! DABS(DMU) < TINY
!               Regular point
                GDN(IMU) = MU0/DMU*( DEXP(-TAU/MU0) - DEXP(-TAU/MUDN(IMU)) )
            END IF ! DABS(DMU) < TINY
        END DO ! IMU = 1, NDN
        I1M(NUP+1:NMU) = GDN*ZM1(NUP+1:NMU)
        Q1M(NUP+1:NMU) = GDN*ZM2(NUP+1:NMU)
        IF (M > 0) THEN
            U1M(NUP+1:NMU) = GDN*ZM3(NUP+1:NMU)
        ELSE ! M = 0
            U1M(NUP+1:NMU) = 0.0D0
        END IF
    END IF ! NDN /= 0
!
END SUBROUTINE SSCAT1M_IP
!===============================================================================
! 01May16 - IDN is replaced with NUP+1
!
! 22Apr16 - Minor changes in comments
!
! 14May15 - Elements of the phase matrix are now on input.
!
! 11Apr15 - QRTm was replaced with QX & RTX. Tested.
!
! 07Apr15 - QPOL1 was replaced with QX0 (modification of QPOL1). Tested.
!
! 05Sep14 - New format of declaration of variables. Tested as part of SORD
!
! 25Aug14 - Renamed from SSCAT1M to SSCAT1M_IP, where 'IP' stands for Intensity
!           and (linear) Polarization. Tested with SORD_IPV.
!
! 14Aug14 - This subroutine was tested against the similar one, but using
!           WHERE .. END WHERE statements. They both were tested against each
!           other - perfect agreement (KN = 100, VZA = 0:1:180, M = 0, 1, 2).
!           Then, each subroutine was run NRUNS times with the following result:
!
!           NRUNS      time(this sub)/time(no WHERE)       time(this subroutine)
!             1000     1.05                                0.042 seconds
!            10000     1.07                                0.425 seconds
!           100000     1.06                                4.250 seconds
!
!           The WHERE .. END WHERE statement takes ~3-5% more time but makes the
!           subroutine shorter.
!
! 06Aug14 - Tests for Rayleigh scattering, no depolarization:
!
!           Test 1: SSA = 0.99999999, mu0 = cos(0), A1K = 0.5*SSA*[1 0 0.5],
!           B1K = 0.5*SSA*[0 0 +sqrt(6)/2], Tau=0.5. VZA = [0:10:180] (exluding
!           90), mu = -cos(VZA). Note "-" here. M = 0. Tested against SGLSCAT
!           from the code IPOL. Perfect agreement: 6 significant digits in txt.
!           Note that singular point VZA = SZA = 0 is included.
!
!           Test 2: the same as Test 1 except for SSA = 0.7, TAU = 1.5. Ok.
!
!           Test 3: the same as Test 2 excpet for TAU = 0.1. Ok.
!
!           Test 4: the same as Test 3, except for SZA = 70, TAU = 0.5, and
!           SSA = 0.9. Tested against single scattering in the IQUVDOM from the
!           code IPOL (not possible to use SGLSCAT beacuse it includes all
!           M = 0, 1, 2). Perfect agrrement of 6 digits printed in a txt file.
!           Note that singular point MU = MU0 = cos(70) is included.
!
!           Test 5: the same as in Test 4 but for M = 1. Unlike in Tests 1-4,
!           U1M is tested here. Perfect agreement.
!
!           Test 6: the same as in Test 5 but for M = 2. Perfect agreement.
!
!           Test 7: the same as in Test 6, but for different arrangement of MU.
!           only negative, only positive, negative first, positive first, and
!           mixed togeteher. Compared against SSCAT1M with MU arranged as in
!           Tests 1-6. Ok.
!===============================================================================