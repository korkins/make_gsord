SUBROUTINE FORD_IP(MU0, TAU, MU, AZI, A1, B1, NL, NIB, NMU, NUP, NDN, NAZ,     &
                                                                     I1, Q1, U1)
!===============================================================================
! PURPOSE:
!   To compute reflected on TOA and transmitted on BOA path radiance in the
!   single scattering approximation for a set of layers and one direction of
!   the Solar beam. The incident Stokes vector is assumed to be [2pi 0 0 0].
!
! INPUT:
!   MU0   D(1)              Cosine of the solar zenith angle, mu0=cos(SZA)>0
!   TAU   D(NL)             Optical thickness of layers
!   MU    D(NMU)            Cos(VZA), VZA/=90. On TOA: mu<0. On BOA: mu>0
!   AZI   D(NAZ)            View (relative) azimuth, [0:2PI] rad
!   A1    D(NMU, NAZ, NL)   Phase function
!   B1    D(NMU, NAZ, NL)   12=21 element of the phase matrix
!   NL    I(1)              Number of layers
!   NIB   I(1)              Number of internal boundaries
!   NMU   I(1)              Number of view directions, NUP+NDN>0. Length of MU
!   NUP   I(1)              Number of directions with mu<0, if any
!   NDN   I(1)              Number of directions with mu>0, if any
!   NAZ   I(1)              Number of relative azimuths where A1 & B1 are given
!
! OUTPUT:
!   I1, Q1, U1   D(NMU, NAZ)   The Stokes vector normalized by 2pi
!
! TREE:
!   FORD_IP
!         |
!         +-SSCAT1_IP (First order of scattering by a single layer)
!
! COMMENTS:
!   On output, I1, Q1, and U1 are normalized to incident beam Io = [2pi 0 0 0].
!
!   Directions of reflection, mu < 0, if any, MUST come first, followed by
!   mu > 0, if any.
!
!   TODO: remove rotation from SSCAT1_IP and apply in the end of FORD_IP to
!   compute the U. Advantage: one does not have to accumulate U over boundaries.
!
! REFERENCES:
!   -
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NL, NIB, NMU, NUP, NDN, NAZ
    REAL*8, INTENT(IN) :: MU0
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NMU, NAZ, NL), INTENT(IN) :: A1, B1
    REAL*8, DIMENSION(NAZ), INTENT(IN) :: AZI
    REAL*8, DIMENSION(NMU), INTENT(IN) :: MU
    REAL*8, DIMENSION(NL), INTENT(IN) :: TAU
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NMU, NAZ), INTENT(OUT) :: I1, Q1, U1
!
! LOCAL VARIABLES
    INTEGER &
        IA,  & ! Loop index over azimuth angle
        IB,  & ! Loop index over internal boundaries
        IDN, & ! Initial index for downward directions, mu > 0
        IL     ! Loop index over layers
    REAL*8 &
        TAC    ! Tau accumulated from TOA
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NIB) :: &
        E0 ! Attenuation of the direct beam
    REAL*8, DIMENSION(NDN, NIB) :: &
        EDN ! Bouguer attenuation for mu > 0
    REAL*8, DIMENSION(NUP, NIB) :: &
        EUP ! Bouguer attenuation for mu < 0
    REAL*8, DIMENSION(NMU, NAZ) :: &
        I1L, & ! Intensity for a given layer IL
        Q1L, & ! Q-component for a given layer IL
        U1L    ! U-component for a given layer IL
!===============================================================================
!
!   1st layer
    IL = 1
    CALL SSCAT1_IP(MU0, TAU(IL), A1(:, :, IL), B1(:, :, IL), MU, AZI, NMU, &
                   NUP, NDN, NAZ, I1, Q1, U1)
!
    IF (NL > 1) THEN
!       Attenuation of the direct beam and the beams with mu < 0, if any
        TAC = TAU(1)
        E0(1) = DEXP(-TAC/MU0)
        IF (NUP /= 0) EUP(:, 1) = DEXP(TAU(1)/MU(:NUP))
        DO IB = 2, NIB
            TAC = TAC + TAU(IB)
            E0(IB) = DEXP(-TAC/MU0)
            IF (NUP /= 0) EUP(:, IB) = DEXP(TAC/MU(:NUP))
        END DO ! IB = 2, NIB
!
!       Attenuation of the beams with mu > 0, if any
        IF (NDN /= 0) THEN
            IDN = NUP+1
            TAC = SUM(TAU) - TAU(1)
            EDN(:, 1) = DEXP(-TAC/MU(IDN:))
            DO IA = 1, NAZ
                I1(IDN:, IA) = I1(IDN:, IA)*EDN(:, 1)
                Q1(IDN:, IA) = Q1(IDN:, IA)*EDN(:, 1)
                U1(IDN:, IA) = U1(IDN:, IA)*EDN(:, 1)
            END DO ! IA = 1, NA
            DO IB = 2, NIB
                TAC = TAC - TAU(IB)
                EDN(:, IB) = DEXP(-TAC/MU(IDN:))
            END DO ! IB = 2, NIB
        END IF ! NDN /= 0
!
        DO IL = 2, NL
            CALL SSCAT1_IP(MU0, TAU(IL), A1(:, :, IL), B1(:, :, IL), MU, AZI, &
                                             NMU, NUP, NDN, NAZ, I1L, Q1L, U1L)
!
!           Account for attenuation of the direct beam
            I1L = I1L*E0(IL-1)
            Q1L = Q1L*E0(IL-1)
            U1L = U1L*E0(IL-1)
!
!           Accumulate attenuated radiation for each azimuth from all layers
            IF (NUP/=0) THEN
!               mu < 0
                DO IA = 1, NAZ
                    I1(:NUP, IA) = I1(:NUP, IA) + I1L(:NUP, IA)*EUP(:, IL-1)
                    Q1(:NUP, IA) = Q1(:NUP, IA) + Q1L(:NUP, IA)*EUP(:, IL-1)
                    U1(:NUP, IA) = U1(:NUP, IA) + U1L(:NUP, IA)*EUP(:, IL-1)
                END DO ! IA = 1, NA
            END IF ! NUP/=0
!
            IF (NDN /= 0) THEN
!               mu > 0
                IF (IL == NL) THEN
!                   Bottom layer - no attenuation for mu > 0
                    DO IA = 1, NAZ
                        I1(IDN:, IA) = I1(IDN:, IA) + I1L(IDN:, IA)
                        Q1(IDN:, IA) = Q1(IDN:, IA) + Q1L(IDN:, IA)
                        U1(IDN:, IA) = U1(IDN:, IA) + U1L(IDN:, IA)
                    END DO ! IA = 1, NA
                ELSE ! IF (IL == NL)
!                   Other layers starting from IL=2
                    DO IA = 1, NAZ
                        I1(IDN:, IA) = I1(IDN:, IA) + I1L(IDN:, IA)*EDN(:, IL)
                        Q1(IDN:, IA) = Q1(IDN:, IA) + Q1L(IDN:, IA)*EDN(:, IL)
                        U1(IDN:, IA) = U1(IDN:, IA) + U1L(IDN:, IA)*EDN(:, IL)
                    END DO ! IA = 1, NA
                 END IF ! IL == NL
            END IF ! NDN /= 0
!
        END DO ! IL = 2, NL
    END IF ! NL > 1
!
END SUBROUTINE FORD_IP
!===============================================================================
! 01May16 - IL1=IL-1 is replaced with IL-1
!
! 22Apr16 - Minor changes in comments
!
! 11Oct15 - Mixing of components is moved to separate subroutine.
!           Input is changed
!
! 15May15 - Input is changed: Rayleigh and scattering component are now mixed
!           in this subroutine using ratio, RTO, as well as weighted by SSA/2.
!
! 05Sep14 - New format of declaration of variables. Tested as part of SORD
!
! 28Aug14 - First created and tested in the following scenario against SGLSCAT
!           from IPOL: Rayleigh, dep.factor = 0.0, SZA=50, AZA=0:45:360,
!           TAU=[0.1 0.2 0.4 0.3 0.2], SSA=[1.0 0.9 0.8 0.8 0.7], TOA Flux = 1,
!           VZA=[0:1:80]&[100:1:180]. Exact agreement of 6 digits (txt file).
!           Time = 0.002s (NAZ=9) and  0.0087 (NAZ=361).
!===============================================================================