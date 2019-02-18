SUBROUTINE SSCAT1_IP(MU0, TAU, A1, B1, MU, AZI, NMU, NUP, NDN, NAZ, I1, Q1, U1)
!===============================================================================
! PURPOSE:
!   To compute reflected and transmitted path radiance in the single scattering
!   approximation for a single layer and one direction of the Solar beam.
!   The incident Stokes vector is assumed to be [2pi 0 0 0].
!
! INPUT:
!   MU0   D(1)          Cosine of the solar zenith angle, mu0 = cos(SZA) > 0
!   TAU   D(1)          Optical thickness of the layer
!   A1    D(NMU, NAZ)   The [11]-element of the phase matrix*SSA/2
!   B1    D(NMU, NAZ)   The [12]=[21]-element of the phase matrix*SSA/2
!   MU    D(NMU)        Cos(VZA), VZA /= 90. On TOA: mu < 0. On BOA: mu > 0
!   AZI   D(NAZ)        Relative azimuth, [0:2PI] rad., where A1 and B1 are given
!   NMU   I(1)          Number of view directions, NUP + NDN > 0. Length of MU
!   NUP   I(1)          Number of directions with mu < 0, if any
!   NDN   I(1)          Number of directions with mu > 0, if any
!   NAZ   I(1)          Number of relative azimuths where A1 and B1 are given
!
! OUTPUT:
!   I1, Q1, U1   D(NMU, NAZ)   Components of the Stokes vector, normalized by 2pi
!
! TREE:
!   SSCAT1_IP
!           |
!           +-ROTATOR2 (Rotate from scattering plane to meridian plane)
!
! COMMENTS:
!   IMPORTANT 1: on input, A1 and B1 contain the factor of SSA/2, where
!   SSA is the single scattering albedo of the layer.
!
!   IMPORTANT 2: on output, I1, Q1, and U1 are normalized to incident beam
!   Io = [2pi 0 0 0].
!
!   IMPORTANT 3: mu < 0, if any, MUST come first, followed by mu > 0, if any.
!
!   The system of coordinates used in this subroutine is exactly the same as
!   defined in [1].
!
!   The following expressions are used in this subroutine:
!
!   I = 1/2 * SSA * R(pi-x2) * F(:, 1) * G                                   (1)
!
!   where
!
!   G = mu0/(mu0-mu)*[1 - exp{(mu0-mu)*Tau/(mu0*mu)}], mu < 0                (2)
!   G = mu0/(mu0-mu)*[exp(-Tau/mu0) - exp(-Tau/mu)],  mu > 0, mu /= mu0      (3)
!   G = Tau/mu * exp(-Tau/mu),  mu = mu0                                     (4)
!
!   The phase matrix and the rotation matrix are (note the "-" signs!!)
!
!        |a1  b1   0   0|
!        |b1  a2   0   0|
!   F  = | 0   0  a3  b2|,                                                   (5)
!        | 0   0 -b2  a4|
!
!               |1  0   0   0|
!               |0  C2 -S2  0|
!   R(pi-x2)  = |0  S2  C2  0|,                                              (6)
!               |0  0   0   1|
!
!   where C2 = cos(x2) and S2 = sin(x2). Eqs. (5), (6) and unpolarized direct
!   solar beam give for the 3 components
!
!                        |1 |
!   R(pi-x2)*F*Io  = 2pi*|C2|.                                               (7)
!                        |S2|
!
! REFERENCES:
!   1. Hovenier JW et al., 2004: Transfer of Polarized Light in Planetary
!      Atmosphere. Basic Concepts and Practical Methods, Dordrecht: Kluwer
!      Academic Publishers.
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    REAL*8, PARAMETER :: TINY = 1.0D-8   ! A small number to compare doubles
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NMU, NUP, NDN, NAZ
    REAL*8, INTENT(IN) :: MU0, TAU
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NMU, NAZ), INTENT(IN) :: A1, B1
    REAL*8, DIMENSION(NAZ), INTENT(IN) :: AZI
    REAL*8, DIMENSION(NMU), INTENT(IN) :: MU
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NMU, NAZ), INTENT(OUT) :: I1, Q1, U1
!
! LOCAL VARIABLES
    INTEGER &
        IA,  & ! Loop index over azimuth angle
        IMU    ! Loop index for mu > 0: IMU = 1, NDN
    REAL*8 &
        C2,  & ! cos(2X) in the rotator
        S2,  & ! sin(2X) in the rotator
        DMU, & ! mu0-mu for a single value of mu
        TMU    ! Tau/mu in case of mu -> mu0
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NDN) :: &
        DMUDN, & ! mu0 - mu for mu > 0
        GDN,   & ! Geometric factor for mu > 0
        MUDN     ! mu > 0
    REAL*8, DIMENSION(NUP) :: &
        DMUUP, & ! mu0 - mu for mu < 0
        GUP, &   ! Geometric factor for mu < 0
        MUUP     ! mu < 0
!===============================================================================
!
    I1 = A1
    Q1 = B1
    U1 = 0.0D0
!
!   Process upward directions, if any
    IF (NUP /= 0) THEN
        MUUP = MU(:NUP)
        DMUUP = MU0 - MUUP
        GUP = MU0/DMUUP*( 1.0D0 - DEXP(TAU*DMUUP/(MU0*MUUP)) )
        DO IA = 1, NAZ
            I1(:NUP, IA) = GUP*I1(:NUP, IA)
            Q1(:NUP, IA) = GUP*Q1(:NUP, IA)
        END DO ! IA = 1, NAZ
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
        DO IA = 1, NAZ
            I1(NUP+1:NMU, IA) = GDN*I1(NUP+1:NMU, IA)
            Q1(NUP+1:NMU, IA) = GDN*Q1(NUP+1:NMU, IA)
        END DO ! IA = 1, NAZ
    END IF ! NDN /= 0
!
!   Rotation of the reference plane
    DO IA = 1, NAZ
        DO IMU = 1, NMU
            CALL ROTATOR2(MU0, MU(IMU), AZI(IA), S2, C2)
            U1(IMU, IA) = Q1(IMU, IA)*S2
            Q1(IMU, IA) = Q1(IMU, IA)*C2
        END DO ! IMU = 1, NMU
    END DO ! IA = 1, NAZ
!
END SUBROUTINE SSCAT1_IP
!===============================================================================
! 01May16 - IDN is replaced with NUP+1
!
! 22Apr16 - Minor changes in comments
!
! 05Sep14 - New format of declaration of variables. Tested as part of SORD
!
! 28Aug14 - First created and tested against SGLSCAT from IPOL for Rayleigh
!           scattering, depol=0.0, SSA=0.9, Tau=0.2, SZA=50, TOA_Flux=1.0.
!           VZA=[0:1:80]&[100:1:180], AZI=[0:45:360]. Exact coincidence of 6
!           significant digits (txt file).
!
!           Timing (Intel i7 2.2GHz, W7, 64bit, ifort 11.0.3453.2008, no
!           optimization): 0.013s. for the same case as above, except for
!           AZI = 0:1:360, NAZ=361.
!===============================================================================