SUBROUTINE SORDM_IP(IM, NOS, ISRF, EPSI, I1B, Q1B, U1B, ZM11, ZM12, ZM13,      &
                                                        ZM21, ZM22, ZM23,      &
                                                        ZM31, ZM32, ZM33,      &
                                                             RM11, RM12, RM13, &
                                                             RM21, RM22, RM23, &
                                                             RM31, RM32, RM33, &
                    EUP, EDN, NEWLR, NL, NB, NLR, NMU, NUP, NDN, NG1, NG2,     &
                    I2M, Q2M, U2M)
!===============================================================================
! PURPOSE:
!   To compute the m-th Fourier moment of solution of the vector radiative
!   transfer equation using the method of successive orders of scattering [1-4]
!
! INPUT:
!   IM      I(1)               Number of the Fourier harmonic, IM = 1, 2, 3
!   ISRF    I(1)               Index of the surface model, ISRF = 0, 1, 2 ...
!   EPSI    D(1)               Absolute accuracy for I, or manual mode
!   ZMij    D(NG2, NMU, NLR)   M-th Fourier moment of the phase matrix
!   RMij    D(NG1, NUP)        Fourier moments for surface reflection
!   EUP     D(NUP, NL)         Attenuation of ascending radiation at all levels
!   EDN     D(NDN)             Attenuation of descending radiation through dTau
!   NEWLR   I(NL)              Mask for new OPTICAL layer
!   NL      I(1)               Number of microlayers
!   NB      I(1)               Number of boundaries at microlayers, NL+1
!   NLR     I(1)               Number of optical layers, NLR <= NL
!   NMU     I(1)               Total number of zeniths, Gauss & user defined
!   NUP     I(1)               Total number of negative Gauss & dummy nodes
!   NDN     I(1)               Total number of positive Gauss & dummy nodes
!   NG1     I(1)               Number of Gaussian nodes in hemisphere
!   NG2     I(1)               Number of Gaussian nodes in sphere, 2*NG1
!
! INOUT:
!   NOS     I(1)               Number of orders of scattering (if EPSI < 0)
!                              or number of accumulated orders (if EPSI > 0)
!   I1B     D(NG2, NB)         m-th Fourier moment of the Stokes vector ..
!   Q1B     D(NG2, NB)         .. for the 1st order at Gauss nodes & all ..
!   U1B     D(NG2, NB)         .. boundaries. *** REDEFINED on output ***
!
! OUTPUT:
!   IM, QM, UM   D(NMU, 2)     Solution of the VRTE at all zenith angles,
!                              on TOA & BOA
!
! TREE:
!   -
!
! COMMENTS:
!   IM=1 is the azimuthally averaged value, i.e. corresponds to m=0 in
!   theoretical evaluations.
!
!   At least 2nd order of scattering is assumed necessary. Check must be
!   completed by a calling subroutine.
!
!   NEWLR = [1 1 1 2 2 3 4 4 4 4 4 4] stands for 12 microlayers and 4 layers
!   with different optical properties. Thickness of the 3rd one equals dTau.
!
!   RMij must contain the factor of 2.0 which comes from normalization
!   1/pi * 2pi = 2. Here 1/pi comes from the form of the bottom boundary
!   condition, and 2pi from normalization of the direct beam on TOA.
!
!   Dummy (user defined) nodes are optional.
!
! REFERENCES:
!   1. van de Hulst HC, 1980: Multiple light scattering. Tables, Formulas, and
!      Applications. Vol.1, Section 4.3, p. 46
!   2. Lenoble J (Ed.), 1985: Radiative transfer in scattering and absorbing
!      atmospheres: standard computational procedures. Part 1, Section 3.7 (Rus)
!   3. Lenoble J et al, 2007: JQSRT, V107, P479 (in [3], mu > 0 means "up")
!   4. Wauben WMF et al, 1993: Astron. Astrophys., V276, P589
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    INTEGER, PARAMETER :: NOS_MAX = 50  ! Maximum number of orders
    REAL*8, PARAMETER :: TINY = 1.0D-12 ! Defines small I2B relative to average
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: IM, ISRF, NL, NB, NLR, NMU, NUP, NDN, NG1, NG2
    REAL*8, INTENT(IN) :: EPSI
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NG2, NMU, NLR), INTENT(IN) :: ZM11, ZM12, ZM13,          &
                                                    ZM21, ZM22, ZM23,          &
                                                    ZM31, ZM32, ZM33
    REAL*8, DIMENSION(NG1, NUP), INTENT(IN) :: RM11, RM12, RM13,               &
                                               RM21, RM22, RM23,               &
                                               RM31, RM32, RM33
    INTEGER, DIMENSION(NL), INTENT(IN) :: NEWLR
    REAL*8, DIMENSION(NUP, NL), INTENT(IN) :: EUP
    REAL*8, DIMENSION(NDN), INTENT(IN) :: EDN
!
! DUAL INTENT VARIABLES
    INTEGER, INTENT(INOUT) :: NOS
    REAL*8, DIMENSION(NG2, NB), INTENT(INOUT) :: I1B, Q1B, U1B
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NMU, 2), INTENT(OUT) :: I2M, Q2M, U2M
!
! LOCAL VARIABLES
    LOGICAL &
        CONT   ! Continue or stop accumulation of orders
    INTEGER &
        IB,  & ! Loop index over boundaries of microlayers, 1:NB
        IG,  & ! Loop index over Gauss nodes, 1:NG2
        IG1, & ! mu(IG1) is the first Gauss node (negative, i.e. upward)
        IG2, & ! mu(IG2) is the last Gauss node (positive, i.e. downward)
        IL,  & ! Loop index over microlayers, 1:NL
        ILR, & ! Loop index over (macro)layers
        IMU, & ! Loop index over mu, IMU = 1:NMU
        I0,  & ! Initial value of IMU on BOA for geometric series
        IO     ! Loop index over orders of scattering
    REAL*8 &
        AVI    ! Average intensity for IM > 1, used in geometric series
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NG2, NL) :: &
        I1A, Q1A, U1A ! Average of two components at two adjacent boundaries
    REAL*8, DIMENSION(NUP) :: &
        I2R, Q2R, U2R ! Diffuse reflection from surface
    REAL*8, DIMENSION(NDN) :: &
        CEDN          ! Complementary EDN, 1 - EDN
    REAL*8, DIMENSION(NUP) :: &
        CEUP          ! Complementary EUP, 1 - EUP(:, 1)
    REAL*8, DIMENSION(NMU) :: &
        GSTOA, GSBOA  ! Common ratio for geometric series
    REAL*8, DIMENSION(NMU, NL) :: &
        ZI, ZQ, ZU    ! Convolution of Z and the Stokes vector in each level
    REAL*8, DIMENSION(NMU, NB) :: &
        I2B, Q2B, U2B ! Diffuse IQU at each boundary at Gauss & dummy nodes
    REAL*8, DIMENSION(NMU, 2) :: & ! NMU=126
        I2GS, Q2GS, U2GS ! Save previous on TOA & BOA order for geometric series
!
!!   For Aeronet-Linux: uncomment the line below
!    REAL*8 DDOT ! BLAS (See also OpenBLAS & BLIS - faster ??)
!===============================================================================
!
!   No U for m = 0 (IM = m+1 = 1)
    IF (IM == 1) U2M = 0.0D0
!
!   Initialize
    IO = 1        ! Count orders
    CONT = .TRUE. ! Compute at least 2nd order of scattering
!
!   Locate elements in arrays
    IG1 = NUP-NG1+1
    IG2 = IG1+NG2-1
!
!   1 - exp(-dTau/|mu|)
    CEDN = 1.0D0 - EDN
    CEUP = 1.0D0 - EUP(:, 1)
!
!   MAIN LOOP OVER ORDERS OF SCATTERING
!
    DO WHILE (CONT)
        IO = IO+1 ! IO = 2 is the first order to be processed
!        WRITE(*, *) 'IO = ', IO
!
        DO IL = 1, NL
            I1A(:, IL) = 0.5D0*(I1B(:, IL) + I1B(:, IL+1))
            Q1A(:, IL) = 0.5D0*(Q1B(:, IL) + Q1B(:, IL+1))
        END DO ! IL = 1, NL
!
        IF (IM > 1) THEN ! U-component
            DO IL = 1, NL
                U1A(:, IL) = 0.5D0*(U1B(:, IL) + U1B(:, IL+1))
            END DO ! IL = 1, NL
        END IF ! IM > 1
!
!       INTEGRATE Z*[I Q U] OVER mu' AT EACH BOUNDARY
!
        DO IL = 1, NL
            ILR = NEWLR(IL)
            DO IMU = 1, NMU
                ZI(IMU, IL) = SUM(ZM11(:, IMU, ILR)*I1A(:, IL)) +              &
                              SUM(ZM12(:, IMU, ILR)*Q1A(:, IL)) +              &
                              SUM(ZM13(:, IMU, ILR)*U1A(:, IL))
                ZQ(IMU, IL) = SUM(ZM21(:, IMU, ILR)*I1A(:, IL)) +              &
                              SUM(ZM22(:, IMU, ILR)*Q1A(:, IL)) +              &
                              SUM(ZM23(:, IMU, ILR)*U1A(:, IL))
!! For Aeronet-Linux
!                ZI(IMU, IL) = DDOT(NG2, ZM11(:, IMU, ILR), 1, I1A(:, IL), 1) + &
!                              DDOT(NG2, ZM12(:, IMU, ILR), 1, Q1A(:, IL), 1) + &
!                              DDOT(NG2, ZM13(:, IMU, ILR), 1, U1A(:, IL), 1)
!                ZQ(IMU, IL) = DDOT(NG2, ZM21(:, IMU, ILR), 1, I1A(:, IL), 1) + &
!                              DDOT(NG2, ZM22(:, IMU, ILR), 1, Q1A(:, IL), 1) + &
!                              DDOT(NG2, ZM23(:, IMU, ILR), 1, U1A(:, IL), 1)
!! End Aeronet-Linux
            END DO ! IMU = 1, NMU
        END DO ! IL = 1, NL
!
        IF (IM > 1) THEN ! U-component
            DO IL = 1, NL
                ILR = NEWLR(IL)
                DO IMU = 1, NMU
                    ZU(IMU, IL) = SUM(ZM31(:, IMU, ILR)*I1A(:, IL)) +          &
                                  SUM(ZM32(:, IMU, ILR)*Q1A(:, IL)) +          &
                                  SUM(ZM33(:, IMU, ILR)*U1A(:, IL))
!! For Aeronet-Linux
!!              ZU(IMU, IL) = DDOT(NG2, ZM31(:, IMU, ILR), 1, I1A(:, IL), 1) + &
!!                            DDOT(NG2, ZM32(:, IMU, ILR), 1, Q1A(:, IL), 1) + &
!!                            DDOT(NG2, ZM33(:, IMU, ILR), 1, U1A(:, IL), 1)
!! End Aeronet-Linux
                END DO ! IMU = 1, NMU
            END DO ! IL = 1, NL
        END IF ! IM > 1
!
!       INTEGRATE DOWNWARD RADIATION OVER TAU
!
!       Integral from TOA to TOA, mu > 0
        I2B(NUP+1:, 1) = 0.0D0
        Q2B(NUP+1:, 1) = 0.0D0
        IF (IM > 1) U2B(NUP+1:, 1) = 0.0D0
!       Integral from TOA to each boundary other than TOA:
!       .. 2nd boundary (accounts for zeros on the 1st boundary)
        I2B(NUP+1:, 2) = CEDN*ZI(NUP+1:, 1)
        Q2B(NUP+1:, 2) = CEDN*ZQ(NUP+1:, 1)
        IF (IM > 1) U2B(NUP+1:, 2) = CEDN*ZU(NUP+1:, 1)
!       .. all other boundaries, if any
        DO IB = 3, NB
            I2B(NUP+1:, IB) = EDN*I2B(NUP+1:, IB-1) + CEDN*ZI(NUP+1:, IB-1)
            Q2B(NUP+1:, IB) = EDN*Q2B(NUP+1:, IB-1) + CEDN*ZQ(NUP+1:, IB-1)
        END DO ! IB = 3, NB
!
        IF (IM > 1) THEN ! U-component
            DO IB = 3, NB
                U2B(NUP+1:, IB) = EDN*U2B(NUP+1:, IB-1) + CEDN*ZU(NUP+1:, IB-1)
            END DO ! IB = 3, NB
        END IF ! IM > 1     
!
!       INTEGRATE UPWARD RADIATION OVER TAU
!
!       Integral from BOA to BOA, mu < 0
        I2B(:NUP, NB) = 0.0D0
        Q2B(:NUP, NB) = 0.0D0
        IF (IM > 1) U2B(:NUP, NB) = 0.0D0
!       Integral from BOA to each boundary other than BOA:
!       .. (NB-1)-boundary (accounts for zeros on the last boundary)
        I2B(:NUP, NL) = CEUP*ZI(:NUP, NL)
        Q2B(:NUP, NL) = CEUP*ZQ(:NUP, NL)
        IF (IM > 1) U2B(:NUP, NL) = CEUP*ZU(:NUP, NL)
!       .. all other boundaries, if any
        DO IB = NL-1, 1, -1 ! NL-1 = NB-2
            I2B(:NUP, IB) = EUP(:, 1)*I2B(:NUP, IB+1) + CEUP*ZI(:NUP, IB)
            Q2B(:NUP, IB) = EUP(:, 1)*Q2B(:NUP, IB+1) + CEUP*ZQ(:NUP, IB)
        END DO ! IB = NL-1, 1, -1
!
        IF (IM > 1) THEN ! U-component
            DO IB = NL-1, 1, -1 ! NL-1 = NB-2
                U2B(:NUP, IB) = EUP(:, 1)*U2B(:NUP, IB+1) + CEUP*ZU(:NUP, IB)
            END DO ! IB = NL-1, 1, -1
        END IF ! IM > 1
!
!       SURFACE REFLECTION
!
        IF (ISRF > 0 .AND. ISRF < 6) THEN ! Scalar (depolarizing) reflection
!
            DO IMU = 1, NUP
                I2R(IMU) = SUM(RM11(:, IMU)*I1B(NG1+1:, NB))
!! For Aeronet-Linux
!!                I2R(IMU) = DDOT(NG1, RM11(:, IMU), 1, I1B(NG1+1:, NB), 1)
!! End Aeronet-Linux
            END DO ! IMU = 1, NMU
!
!           Add the attenuated reflected component to the path radiance
            I2B(:NUP, NB) = I2B(:NUP, NB) + I2R ! No attenuation on BOA
!
            DO IB = 1, NL ! NL = NB-1, attenuation at all boundaries
                I2B(:NUP, IB) = I2B(:NUP, IB) + I2R*EUP(:, NB-IB)
            END DO ! IB = 1, NL
!
        ELSE IF (ISRF > 5) THEN ! Polarized reflectance
!
            DO IMU = 1, NUP
                I2R(IMU) = SUM(RM11(:, IMU)*I1B(NG1+1:, NB)) +                   &
                           SUM(RM12(:, IMU)*Q1B(NG1+1:, NB)) +                   &
                           SUM(RM13(:, IMU)*U1B(NG1+1:, NB))
                Q2R(IMU) = SUM(RM21(:, IMU)*I1B(NG1+1:, NB)) +                   &
                           SUM(RM22(:, IMU)*Q1B(NG1+1:, NB)) +                   &
                           SUM(RM23(:, IMU)*U1B(NG1+1:, NB))
!! For Aeronet-Linux
!!                I2R(IMU) = DDOT(NG1, RM11(:, IMU), 1, I1B(NG1+1:, NB), 1) +      &
!!                           DDOT(NG1, RM12(:, IMU), 1, Q1B(NG1+1:, NB), 1) +      &
!!                           DDOT(NG1, RM13(:, IMU), 1, U1B(NG1+1:, NB), 1)
!!                Q2R(IMU) = DDOT(NG1, RM21(:, IMU), 1, I1B(NG1+1:, NB), 1) +      &
!!                           DDOT(NG1, RM22(:, IMU), 1, Q1B(NG1+1:, NB), 1) +      &
!!                           DDOT(NG1, RM23(:, IMU), 1, U1B(NG1+1:, NB), 1)
!! End Aeronet-Linux
            END DO ! IMU = 1, NMU
            I2B(:NUP, NB) = I2B(:NUP, NB) + I2R
            Q2B(:NUP, NB) = Q2B(:NUP, NB) + Q2R
!
            DO IB = 1, NL
                I2B(:NUP, IB) = I2B(:NUP, IB) + I2R*EUP(:, NB-IB)
                Q2B(:NUP, IB) = Q2B(:NUP, IB) + Q2R*EUP(:, NB-IB)
            END DO ! IB = 1, NL
!
            IF (IM > 1) THEN ! U-component
                DO IMU = 1, NUP
                    U2R(IMU) = SUM(RM31(:, IMU)*I1B(NG1+1:, NB)) +               &
                               SUM(RM32(:, IMU)*Q1B(NG1+1:, NB)) +               &
                               SUM(RM33(:, IMU)*U1B(NG1+1:, NB))
!! For Aeronet-Linux
!!                    U2R(IMU) = DDOT(NG1, RM31(:, IMU), 1, I1B(NG1+1:, NB), 1) + &
!!                               DDOT(NG1, RM32(:, IMU), 1, Q1B(NG1+1:, NB), 1) + &
!!                               DDOT(NG1, RM33(:, IMU), 1, U1B(NG1+1:, NB), 1)
!! End Aeronet-Linux
                END DO ! IMU = 1, NMU
                U2B(:NUP, NB) = U2B(:NUP, NB) + U2R
                DO IB = 1, NL
                    U2B(:NUP, IB) = U2B(:NUP, IB) + U2R*EUP(:, NB-IB)
                END DO ! IB = 1, NL
            END IF ! IM > 1
!
        END IF ! ISRF > 0 .AND. ISRF < 6
!
!       ACCUMULATE SCATTERING ORDERS
!
        IF (IO == 2) THEN
!           Initialize output
            I2M(:, 1) = I2B(:, 1)   ! TOA
            I2M(:, 2) = I2B(:, NB)  ! BOA
            Q2M(:, 1) = Q2B(:, 1)
            Q2M(:, 2) = Q2B(:, NB)
            IF (IM > 1) THEN
                U2M(:, 1) = U2B(:, 1)
                U2M(:, 2) = U2B(:, NB)
            END IF ! IM > 1
        ELSE ! IO > 2
!           Accumulate output
            I2M(:, 1) = I2M(:, 1) + I2B(:, 1)
            I2M(:, 2) = I2M(:, 2) + I2B(:, NB)
            Q2M(:, 1) = Q2M(:, 1) + Q2B(:, 1)
            Q2M(:, 2) = Q2M(:, 2) + Q2B(:, NB)
            IF (IM > 1) THEN
                U2M(:, 1) = U2M(:, 1) + U2B(:, 1)
                U2M(:, 2) = U2M(:, 2) + U2B(:, NB)
            END IF ! IM > 1
        END IF ! IO == 2
!
!       CHECK STOP CRITERION
!
        CONT = .FALSE. ! By default, do not compute next order
        IF (EPSI < 0.0D0) THEN ! Manual mode: use NOS
            IF (IO < NOS) CONT = .TRUE.
        ELSE ! EPSI > 0: Automatic mode, use EPSI or NOS_MAX
            NOS = IO ! Saved for output, no other purpose
            IF (IO < NOS_MAX) THEN ! Maximum number of orders not yet reached    
!               Check magnitude of the current order at all boundaries but
!               at Gauss nodes only (dummy nodes are not used for next order)
                DO IB = 1, NB
                    DO IG = IG1, IG2
                        IF (DABS(I2B(IG, IB)) > EPSI) THEN
                            CONT = .TRUE.
                            EXIT ! Exit the IG-loop
                        END IF ! II1B(IG, 1) > EPSI .OR. I1B(IG, NB) > EPSI
                    END DO ! IG = IG1, IG2
                    IF (CONT) EXIT ! Exit from the IB-loop
                END DO ! IB = 1, NB
            END IF ! IO < NOS_MAX
        END IF ! EPSI < 0.0D0
!
!       REDEFINE PREVIOUS SCATTERING ORDER TO GET THE NEXT ONE
!
        IF (CONT) THEN
!           Only Gauss nodes are used; I1B, Q1B, U1B are redefined
            I1B = I2B(IG1:IG2, :)
            Q1B = Q2B(IG1:IG2, :)
            IF (IM > 1) U1B = U2B(IG1:IG2, :)
!           Save current order for geometric series (automatic mode only)
            IF (IO > 1 .AND. EPSI > 0.0D0) THEN
                I2GS(:, 1) = I2B(:,  1) ! TOA
                I2GS(:, 2) = I2B(:, NB) ! BOA
                Q2GS(:, 1) = Q2B(:,  1) ! TOA
                Q2GS(:, 2) = Q2B(:, NB) ! BOA
                IF (IM > 1) THEN
                    U2GS(:, 1) = U2B(:,  1) ! TOA
                    U2GS(:, 2) = U2B(:, NB) ! BOA
                END IF ! IM > 1
            END IF ! (IO > 1 .AND. EPSI > 0.0D0)   
        END IF ! CONT
!
    END DO ! WHILE (CONT == .TRUE.)
!
!   GEOMETRIC SERIES: IN AUTOMATIC MODE ONLY
!
    IF (IO > 2 .AND. EPSI > 0.0D0) THEN
!       Compute infinite sum of the series starting from the last computed
!       order, I2B ( ***** note, same ratio is used for Q & U *****)
!               S = I2B/(1 - r).                                             (A)
!       Ratio is computed using the previous order, I2GS
!               r = I2B/I2GS                                                 (B)
!       Also, I2B has already been added to the series. Hence
!               S0 = S - I2B                                                 (C)
!       is to be computed. Combining (A), (B), and (C), one gets
!               S0 = I2B*I2B/(I2GS - I2B).                                   (D)
        IF (ISRF > 0) THEN
            I0 = 1     ! Process UP & DN directions on BOA
        ELSE
            I0 = NUP+1 ! UP = 0.0, process only DN directions
        END IF
!
        GSTOA = 0.0D0
        GSBOA = 0.0D0
!
        IF (IM == 1) THEN
!           Simply compute the ratio
            GSTOA(:NUP) = I2B(:NUP,  1)/(I2GS(:NUP, 1) - I2B(:NUP,  1))
            GSBOA(I0:NMU) = I2B(I0:NMU, NB)/(I2GS(I0:NMU, 2) - I2B(I0:NMU, NB))
        ELSE ! IM > 1
!           For higher Fourier moments, the ratio I{n+1}/I{n} is not defined
!           at the zenith on TOA and nadir on BOA. One has to find these
!           directions and avoid geometric series for them.
!
!           TOA:
            AVI = SUM( ABS(I2B(:NUP, 1)) )/NUP
            WHERE ( DABS(I2B(:NUP, 1))/AVI > TINY ) GSTOA(:NUP) = &
                I2B(:NUP,  1)/(I2GS(:NUP, 1) - I2B(:NUP,  1))
!           BOA:
            AVI = SUM( ABS(I2B(I0:NMU, NB)) )/(NMU+1-I0)
            WHERE ( DABS(I2B(I0:NMU, NB))/AVI > TINY ) GSBOA(I0:NMU) = &
                I2B(I0:NMU, NB)/(I2GS(I0:NMU, 2) - I2B(I0:NMU, NB))
!
        END IF ! IM == 1
!
!       I-component: see Eq.(D) above
        I2B(:,  1) = I2B(:,  1)*GSTOA
        I2B(:, NB) = I2B(:, NB)*GSBOA
        I2M(:, 1) = I2M(:, 1) + I2B(:, 1)
        I2M(:, 2) = I2M(:, 2) + I2B(:, NB)
!       Q-component: corrected using geometric series from the intensity
        Q2B(:,  1) = Q2B(:,  1)*GSTOA
        Q2B(:, NB) = Q2B(:, NB)*GSBOA
        Q2M(:, 1) = Q2M(:, 1) + Q2B(:, 1)
        Q2M(:, 2) = Q2M(:, 2) + Q2B(:, NB)
!       U-component: again, ratio from I is used
        U2B(:,  1) = U2B(:,  1)*GSTOA
        U2B(:, NB) = U2B(:, NB)*GSBOA
        U2M(:, 1) = U2M(:, 1) + U2B(:, 1)
        U2M(:, 2) = U2M(:, 2) + U2B(:, NB)
!
    END IF ! IO > 2 .AND. EPSI > 0.0D0
!
END SUBROUTINE SORDM_IP
!===============================================================================
! 01May16 - All indices "I1=I+1" are replaced with explicit "I+1"
!
! 22Apr16 - Minor changes in comments
!
! 11Nov15 - Geometric series have been reorganized and tested. Now it works for
!           SZA > 0 as well.
!
! 31Oct15 - Geometric progression from I is now applied to Q and U. Tested.
!           IM is now on input to process U-component only if IM > 1.
!           Significant changes in the text. Tested.
!
! 14Oct15 - Renamed to SORDM_IP. Tested.
!
! 19Sep15 - On output, SORD_IP now provides values at user & Gauss nodes
!           CRIT is not on input, EPSI < 0 is used if number of orders is
!           defined by a user. Geometric progression is implemented. Several
!           other modifications are made.
!
! 17Jul15 - 1st scattering order has been totally removed from this subroutine.
!           Fluxes are temporary removed. Input has been changed significantly.
!
! 07Jul15 - Integration over Tau has been slightly modified: instead of trapz,
!           exp(-dTau/|mu|) is now integrated analytically, while convolution
!           time is averaged over two neighbor boundaries. Tested vs previous
!           version.
!
! 02Jul15 - Integration over Tau was carefully validated, no errors found. But
!           Gauss seems to show better performance as compared to trapz. Some
!           cosmetic changes were made.
!
! 28May15 - Corrected lines 422-424
!
! 30Apr15 - DDOT is replaced with DDOT5. Tested for Rayleigh + RTLS
!
! 27Apr15 - The subroutine, including input parameters, has been significantly
!           modified and tested using Rayleigh atmosphere over RTLS
!
! 14Apr15 - Surface reflection of the direct beam is restructured. Tested using
!           RTLS.
!
! 11Apr15 - QRTm is replaced with QX & RTX. Tested.
!
! 02Apr15 - SORD_IP was reorganized: IQU are now considered separately, not as
!           elements of a single vector. Step 3 (components) is omitted.
!           Only BLAS_DDOT is now used (no external libraries needed). Tested
!           using test 5. Ok.
!
! 20Mar15 - SORD_IP now provides time on output.
!
! 03Mar15 - In IF (ISRF > 0), SVUP(1:NG3) is now used instead of just SVUP.
!           All test have been repeated. Ok.
!
! 19Feb15 - FLUXI is now included. Tested vs RT3 only the case of Lambertian
!           surface, Rayleigh scattering and three different Tau.
!
! 16Feb15 - NEWLR is now an INTEGER array. Tested vs previous version. Ok.
!
! 23Nov14 - a4k is not used. Only 4 elements are needed. Renamed: NX5 -> NX4,
!           NLR5 -> NLR4. Removed: ILA4. Tested against old version.
!
! 13Nov14 - Automatic criteria is now used as well as the manual mode. In the
!           manual mode, results are exactly as before for 10 scatterings.
!           For the same test case but in the automatic mode, results are as
!           follows:
!
!
!           AZ        SORD, NOS=10   SORD, EPSI=0.001   SORD, EPSI=0.0001
!             0 Iup   0.326743E+00    0.327824E+00       0.327979E+00
!               Idn   0.155616E+01    0.155655E+01       0.155678E+01
!               Qup  -0.463167E-01   -0.460766E-01      -0.461790E-01
!               Qdn   0.688248E-02    0.687937E-02       0.693961E-02
!            90 Iup   0.333687E+00    0.334714E+00       0.334890E+00
!               Idn   0.151043E+00    0.151616E+00       0.151708E+00
!               Qup   0.276127E-01    0.276466E-01       0.277042E-01
!               Qdn   0.271820E-01    0.272139E-01       0.272785E-01
!               Uup  -0.450346E-01   -0.448034E-01      -0.449556E-01
!               Udn  -0.438792E-01   -0.436226E-01      -0.437914E-01
!           180 Iup   0.442629E+00    0.443475E+00       0.443788E+00
!               Idn   0.137329E+00    0.137959E+00       0.138026E+00
!               Qup   0.638373E-02    0.639001E-02       0.644167E-02
!               Qdn  -0.446624E-01   -0.443911E-01      -0.445119E-01
!
!           Time, s:  4.3             1.69               1.96
!
!           M=0:20, NOS(0.001)  = 13 5 5 4 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2
!           M=0:20, NOS(0.0001) = 17 6 7 5 5 4 4 3 3 3 2 2 2 2 2 2 2 2 2 2 2
!                                    m=1: NOS = 6, m=2: NOS = 7
!
! 28Oct14 - NKS is now used in the scattering integral instead of NG2. Test from
!           20Sep14, NKS=min(NK=117, NG2=48)=48 was repeated. Ok. The same
!           numbers. Another test: Rayleigh, dep.factir=0.03, wavelength=0.46,
!           NKS=min(NK=3, NG2=48)=3. Compared vs IPOL, 100 microlayers.
!           TauR~=0.2.
!
!           AZ        IPOL           SORD           Err%
!             0 Iup   0.686642E-01   0.686744E-01  -0.01
!               Idn   0.116664E+00   0.116677E+00  -0.01
!               Qup  -0.453188E-01  -0.453162E-01   0.01
!               Qdn   0.447305E-02   0.447851E-02  -0.12
!            90 Iup   0.727651E-01   0.727763E-01  -0.02
!               Idn   0.716958E-01   0.717070E-01  -0.02
!               Qup   0.238105E-01   0.238154E-01  -0.02
!               Qdn   0.234436E-01   0.234484E-01  -0.02
!               Uup  -0.415167E-01  -0.415190E-01  -0.01
!               Udn  -0.408466E-01  -0.408489E-01  -0.01
!           180 Iup   0.118484E+00   0.118497E+00  -0.01
!               Idn   0.676483E-01   0.676585E-01  -0.02
!               Qup   0.450118E-02   0.450665E-02  -0.12
!               Qdn  -0.445428E-01  -0.445402E-01   0.01
!
! 27Oct14 - Renamed: NO->NOS. Tested against previous version.
!
! 22Sep14 - Linear combination of RTLS & Ocean was tested in case of RTLS -> 0 &
!           Ocean -> 0 against the results below. Ok. General case: 30% RTLS,
!           70% Ocean IPOL_SE was used to generate benchmark results:
!
!           AZ        IPOL           SORD1          Err1%   SORD2         Err2%
!             0 Iup   0.612458E+00   0.610992E+00    0.24   0.612502E+00   -0.01
!               Idn   0.194318E+01   0.194247E+01    0.04   0.194321E+01   -0.00
!               Qup  -0.466480E+00  -0.499403E+00   -7.06  -0.499700E+00   -7.12
!               Qdn  -0.632291E-03  -0.113180E-02  -79.00  -0.123846E-02  -95.87
!            90 Iup   0.100996E+00   0.998933E-01    1.09   0.100986E+00    0.01
!               Idn   0.294339E-01   0.289555E-01    1.63   0.294363E-01   -0.01
!               Qup   0.903498E-03  -0.143614E-05  100.16  -0.102982E-03  111.40
!               Qdn   0.695479E-03   0.586104E-03   15.73   0.556365E-03   20.00
!               Uup  -0.152627E-02  -0.178319E-02  -16.83  -0.181287E-02  -18.78
!               Udn  -0.127056E-02  -0.124370E-02    2.11  -0.124257E-02    2.20
!           180 Iup   0.137962E+00   0.136971E+00   -0.72   0.138025E+00   -0.05
!               Idn   0.130876E-01   0.126666E-01    3.22   0.130976E-01   -0.08
!               Qup  -0.257075E-03  -0.709042E-03 -175.85  -0.788879E-03 -206.87
!               Qdn  -0.145651E-02  -0.153194E-02   -5.18  -0.155310E-02   -6.63
!
!            ***************************************************************
!            * WARNING: CONVERGENCE IN THIS CASE IS NOT SATISFACTORY. WHY? *
!            * A PROBLEM WITH IPOL_SE IS POSSIBLE. RE-RUN THE TEST AGAIN.  *
!            ***************************************************************
!
! 21Sep14 - Ocean, ISRF = 6 (NT). Same atmosphere & geometry as on 20Sep14.
!           NH20=1.33, KH20=0.005, U=2.0 m/s.
!
!           AZ        IPOL           SORD1          Err1%   SORD2         Err2%
!             0 Iup   0.747522E+00   0.746458E+00    0.14   0.747537E+00  -0.00
!               Idn   0.194110E+01   0.194048E+01    0.03   0.194111E+01  -0.00
!               Qup  -0.713250E+00  -0.712812E+00    0.06  -0.713248E+00   0.00
!               Qdn  -0.185676E-02  -0.169662E-02    8.62  -0.185744E-02  -0.04
!            90 Iup   0.122534E-01   0.118579E-01    3.22   0.122555E-01  -0.02
!               Idn   0.262421E-01   0.259778E-01    1.01   0.262457E-01  -0.01
!               Qup  -0.245256E-03  -0.141585E-03   42.27  -0.243352E-03  -0.78
!               Qdn   0.459536E-03   0.500718E-03   -8.96   0.460438E-03  -0.20
!               Uup  -0.210307E-02  -0.205447E-02    2.31  -0.210202E-02   0.05
!               Udn  -0.120761E-02  -0.121201E-02   -0.36  -0.120777E-02  -0.01
!           180 Iup   0.129106E-01   0.126926E-01    1.69   0.129131E-01  -0.02
!               Idn   0.922837E-02   0.906824E-02    1.74   0.923016E-02  -0.02
!               Qup  -0.864547E-03  -0.816653E-03    5.54  -0.865840E-03  -0.15
!               Qdn  -0.165096E-02  -0.162906E-02    1.33  -0.165069E-02   0.02
!
! 20Sep14 - Surface reflection was added. Tested against IPOL (V was neglected),
!           for atmosphere from [4]: mu0 = 0.6,  mu = +/- mu0. SORD1 corresponds
!           to 5 scattering orders & 100 microlayers; SORD2 - twice as much
!           orders. dTau = 0.2/100. Note, SSA = 1.0. Err% < 0: SORD > IPOL
!
!           Lambert, ro = 0.0
!
!           AZ        IPOL           SORD1         Err1%   SORD2         Err2%
!             0 Iup   0.278921E-01   0.276617E-01   0.83   0.279095E-01  -0.06
!               Idn   0.193763E+01   0.193741E+01   0.01   0.193764E+01  -0.00
!               Qup  -0.129097E-02  -0.129259E-02  -0.12  -0.129089E-02   0.01
!               Qdn   0.121312E-03   0.118669E-03   2.18   0.121450E-03  -0.11
!            90 Iup   0.909198E-02   0.898380E-02   1.19   0.909721E-02  -0.06
!               Idn   0.253195E-01   0.252153E-01   0.41   0.253247E-01  -0.02
!               Qup   0.654675E-03   0.651131E-03   0.54   0.654872E-03  -0.03
!               Qdn   0.733838E-03   0.730223E-03   0.49   0.734050E-03  -0.03
!               Uup  -0.111103E-02  -0.110902E-02   0.18  -0.111118E-02   0.01
!               Udn  -0.126710E-02  -0.126482E-02   0.18  -0.126729E-02  -0.02
!           180 Iup   0.115451E-01   0.114803E-01   0.56   0.115480E-01  -0.03
!               Idn   0.840196E-02   0.833875E-02   0.75   0.840506E-02  -0.04
!               Qup  -0.367889E-04  -0.386375E-04  -5.02  -0.368490E-04  -0.16
!               Qdn  -0.128986E-02  -0.129037E-02  -0.04  -0.129003E-02  -0.01
!
!           Lambert, ro = 0.1
!
!           AZ        IPOL           SORD1         Err1%   SORD2         Err2%
!             0 Iup   0.844288E-01   0.839691E-01   0.54   0.844492E-01  -0.02
!               Idn   0.193947E+01   0.193916E+01   0.02   0.193949E+01  -0.00
!               Qup  -0.129518E-02  -0.129802E-02  -0.22  -0.129508E-02   0.01
!               Qdn   0.125523E-03   0.121578E-03   3.14   0.125678E-03  -0.12
!            90 Iup   0.656287E-01   0.652911E-01   0.51   0.656369E-01  -0.01
!               Idn   0.271626E-01   0.269625E-01   0.74   0.271699E-01  -0.03
!               Qup   0.650464E-03   0.645710E-03   0.73   0.650680E-03  -0.03
!               Qdn   0.738049E-03   0.733132E-03   0.67   0.738279E-03  -0.03
!               Uup  -0.111103E-02  -0.110902E-02   0.18  -0.111118E-02  -0.01
!               Udn  -0.126710E-02  -0.126482E-02   0.18  -0.126729E-02  -0.02
!           180 Iup   0.680818E-01   0.677876E-01   0.43   0.680877E-01  -0.01
!               Idn   0.102451E-01   0.100860E-01   1.55   0.102503E-01  -0.05
!               Qup  -0.409995E-04  -0.440587E-04  -7.46  -0.410412E-04  -0.10
!               Qdn  -0.128565E-02  -0.128746E-02  -0.14  -0.128580E-02  -0.01
!
!           Lambert, ro = 0.9
!
!           AZ        IPOL           SORD1          Err1%   SORD2         Err2%
!             0 Iup   0.556689E+00   0.550112E+00   1.18    0.556677E+00  0.00
!               Idn   0.195487E+01   0.195343E+01   0.07    0.195489E+01 -0.00
!               Qup  -0.133035E-02  -0.134344E-02  -0.98   -0.133022E-02  0.01
!               Qdn   0.160694E-03   0.144499E-03  10.08    0.160854E-03 -0.10
!            90 Iup   0.537889E+00   0.531434E+00   1.20    0.537865E+00  0.00
!               Idn   0.425591E-01   0.412339E-01   3.11    0.425700E-01 -0.03
!               Qup   0.615293E-03   0.600282E-03   2.44    0.615538E-03 -0.04
!               Qdn   0.773220E-03   0.756053E-03   2.22    0.773455E-03 -0.03
!               Uup  -0.111103E-02  -0.110902E-02   0.18   -0.111118E-02 -0.01
!               Udn  -0.126710E-02  -0.126482E-02   0.18   -0.126729E-02 -0.01
!           180 Iup   0.540342E+00   0.533931E+00   1.19    0.540316E+00  0.00
!               Idn   0.256416E-01   0.243573E-01   5.00    0.256504E-01 -0.03
!               Qup  -0.761709E-04  -0.894872E-04 -17.48   -0.761833E-04 -0.02
!               Qdn  -0.125048E-02  -0.126454E-02  -1.12   -0.125063E-02 -0.01
!
!           RTLS, [KL KV KG] = [0.539625D0, 0.330065D0, 0.045027]
!
!           AZ        IPOL           SORD1          Err1%   SORD2         Err2%
!             0 Iup   0.304249E+00   0.300276E+00   1.31    0.304241E+00  0.00
!               Idn   0.194830E+01   0.194721E+01   0.06    0.194832E+01 -0.00
!               Qup  -0.133066E-02  -0.134356E-02  -0.10   -0.133065E-02  0.00
!               Qdn   0.204034E-03   0.188038E-03   7.84    0.204092E-03 -0.03
!            90 Iup   0.313980E+00   0.310059E+00   1.24    0.313963E+00 -0.01
!               Idn   0.370013E-01   0.359614E-01   2.81    0.370078E-01 -0.02
!               Qup   0.685632E-03   0.669839E-03   2.30    0.685811E-03 -0.03
!               Qdn   0.823654E-03   0.806342E-03   2.10    0.823828E-03 -0.02
!               Uup  -0.116720E-02  -0.116428E-02   0.25   -0.116738E-02 -0.02
!               Udn  -0.131820E-02  -0.131528E-02   0.22   -0.131840E-02 -0.02
!           180 Iup   0.435179E+00   0.431166E+00   0.92    0.435165E+00  0.00
!               Idn   0.221697E-01   0.210877E-01   4.88    0.221778E-01 -0.04
!               Qup   0.331487E-06  -0.135360E-04 4183.4    0.252881E-06 23.71
!               Qdn  -0.126323E-02  -0.127667E-02  -1.06   -0.126344E-02 -0.02
!
! 04Sep14 - Several layers are now supported. Tested against single layer model
!           as on 25Aug14, AZ=90, for the following cases (dTau=0.002):
!
!           a) NLR = 2, Tau = [50*dTau  50*dTau]
!           b) NLR = 2, Tau = [ 1*dTau  99*dTau]
!           c) NLR = 2, Tau = [99*dTau   1*dTau]
!           d) NLR = 2, Tau = [30*dTau  70*dTau]
!           e) NLR = 2, Tau = [70*dTau  30*dTau]
!           
!           f) NLR = 3, Tau = [15*dTau 30*dTau 55*dTau]
!           g) NLR = 3, Tau = [55*dTau 30*dTau 15*dTau]
!
!           h) NLR = 10, Tau = [10*dTau 10*dTau ... 10*dTau]
!
!           i) NLR = NL = 100, Tau = [dTau dTau ... dTau]
!
!           In all mentioned cases, the same result as for single layer was
!           computed (see below).
!
!           Time: NLR = 2,3 t~3.6s; NLR = 10 t~4.6s; NLR = 100 t~11s.
!
! 25Aug14 - Created from SORD_IPV and tested for Higher orders of scattering,
!           computed from [4, Table 1, p.596] as Total-First for NL=100, NO=10.
!           The following results were obtained:
!
!               TOTAL-FIRST    SORD_IPV    SORD_IP (equal to SORD_IPV)
!            AZ   TOA
!           I  0  0.015018     0.015036    0.015036
!           Q  0 -0.000236    -0.000235   -0.000235
!                 BOA
!           I  0  0.121190     0.121207    0.121207
!           Q  0  0.000121     0.000121    0.000121
!                 TOA
!           I 90  0.004462     0.004467    0.004467
!           Q 90  0.000238     0.000239    0.000239
!           U 90 -0.000331    -0.000331   -0.000331
!                 BOA
!           I 90  0.009406     0.009411    0.009411
!           Q 90  0.000244     0.000244    0.000244
!           U 90 -0.000349    -0.000349   -0.000349
!                 TOA
!           I180  0.003881     0.003884    0.003884
!           Q180  0.000037    -0.000037   -0.000037
!                 BOA
!           I180  0.003616     0.003618    0.003618
!           Q180  0.000416    -0.000416   -0.000416
!
!           TIME, seconds:     4.6         4.0  for all IM=0:20, NM=21, NG1 = 24
!===============================================================================