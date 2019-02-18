SUBROUTINE SORD_IP(NM, NOS, NC, ISUR, EPSI, MU0, TAU0, DTAU, DEPF, TAU, SSA,   &
                   RTO, XK, F1, F2, R1, R2, PSRF, MU, MUG, WG, AZI, ZAZ, WCMA, &
                   WSMA, NEWLR, NK, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD, NRU, &
                   NG1, NG2, NGA, II, QQ, UU, FLXI)
!===============================================================================
! PURPOSE:
!   To solve the 1D monochromatic vector radiative transfer equation on TOA and
!   BOA only using the method of successive orders of scattering. V-component of
!   the Stokes vector is ignored.
!
! INPUT:
!   NM      I(1)          Total number of the Fourier harmonics, m = 0:NM-1
!   NOS     I(1)          Number of orders of scattering
!   NC      I(1)          Number of components: 1-Rayleigh, 2-R & Aerosol
!   ISUR    I(1)          Surface index: 0-black, 1-Lambertian, 2-...
!   EPSI    D(1)          Absolute accuracy for I (if CRIT = 2)
!   MU0     D(1)          cos(SZA) > 0
!   TAU0    D(1)          Total atmosphere optical thickness
!   DTAU    D(1)          Step for integration over Tau
!   DEPF    D(1)          Depolarization factor for Rayleigh scattering
!   TAU     D(NLR)        Optical thicknesses of layers from TOA to BOA
!   SSA     D(NLR)        The same as TAU but for single scattering albedo
!   RTO     D(NLR)        Aerosol-Rayleigh ratio (0 - no aerosol)
!   XK      D(NK, 4)      The phase matrix expansion moments*(2k+1)
!   F1      D(NUD, NAZ)   The [11]-element of the phase matrix
!   F2      D(NUD, NAZ)   The [12]-element of the phase matrix
!   R1      D(NUD, NAZ)   The [11]-element of the phase matrix
!   R2      D(NUD, NAZ)   The [12]-element of the phase matrix
!   PSRF    D(NSPR)       Set of surface parameters; NSPR is in PARAMETERS
!   MU      D(NMU)        Gaussian and dummy nodes, mu < 0 come first
!   MUG     D(NG2)        Gaussian nodes, mug < 0 come first
!   WG      D(NG2)        Gaussian weights corresponding to mug
!   AZI     D(NAZ)        View (relative) azimuth, radians
!   ZAZ     D(NGA)        Azimuthal Gauss nodes, radians
!   WCMA    D(NGA, NM)    Gauss_weight*cos(m*Azimuth_nodes)
!   WSMA    D(NGA, NM)    Similar to WCMA except for sin
!   NEWLR   I(NL)         Profile of microlayers in each optical layer
!   NK      I(1)          Number of k-moments in the scattering integral
!   NLR     I(1)          Number of layers with different optical properties
!   NL      I(1)          Number of microlayers for integration over Tau
!   NB      I(1)          Number of boundaries
!   NMU     I(1)          Total number of Gaussian and dummy nodes
!   NAZ     I(1)          Number of view azimuths
!   NUP     I(1)          Total number of negative Gaussian and dummy nodes
!   NDN     I(1)          Total number of positive Gaussian and dummy nodes
!   NUD     I(1)          Number of user defined (dummy) nodes
!   NRU     I(1)          Number of reflected user defined zeniths
!   NG1     I(1)          Number of Gaussian nodes per hemisphere
!   NG2     I(1)          Number of nodes per sphere, 2*NG1
!   NGA     I(1)          Order of integration over azimuth
!
! OUTPUT:
!   II, QQ, UU   D(NUD, NAZ)   Three components of the Stokes vector;
!                                   mu < 0, if any, come first
!   FLXI   D(4)   Fluxes for I: Up TOA, Up BOA, Down BOA, and Down+Direct BOA
!
! TREE:
!   SORD1WL_IP
!            |
!            +-FORD_IP (Computes first order of scattering on TOA and BOA)
!            |
!            +-SURFTOA_IP (Computes surface correction on TOA)
!            |
!            +-SURFREFM_IP (Computes matrix operator for surface reflection)
!            |
!            +-SURFACEM_IP (Computes m-th Fourier moment for BRDF or BPDF)
!            |
!            +- ... SURFACEM0_IP (Same as SURFACEM_IP but for the direct beam)
!            |
!            +-SORD_IP (Successive orders of scattering)
!
! COMMENTS:
!   The output result is normalized as [2pi 0 0 0].
!
!   This subroutine includes single scattering correction of path radiance
!   and surface correction on TOA.
!
!   The phase function is not truncated.
!
!   Nodes, MU(1:NMU), are organized as follows:
!
!         1 <negative dummy nodes> NRU          NR1 = NRU+1
!       NR1 <negative Gauss nodes> NUP          NT1 = NUP+1
!       NT1 <positive Gauss nodes> NT2          NT2 = NUP+NG1
!       NT3 <positive dummy nodes> NMU          NT3 = NT2+1
!
! REFERENCES:
!   -
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    INTEGER, PARAMETER :: &
        NSPR = 21, & ! Number of surface parameters
        NFL = 4      ! Number of fluxes
    REAL*8, PARAMETER :: &
        PI2 = 6.2831853071795865D0 ! 2pi for net fluxes
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NM, NK, NLR, NL, NB, ISUR, NMU, NUP, &
                           NDN, NUD, NRU, NG1, NG2, NAZ, NGA, NC
    REAL*8, INTENT(IN) :: EPSI, MU0, TAU0, DTAU, DEPF
!
! INPUT ARRAYS
    INTEGER, DIMENSION(NL), INTENT(IN) :: NEWLR
    REAL*8, DIMENSION(NSPR), INTENT(IN) :: PSRF
    REAL*8, DIMENSION(NMU), INTENT(IN) :: MU
    REAL*8, DIMENSION(NG2), INTENT(IN) :: MUG, WG
    REAL*8, DIMENSION(NAZ), INTENT(IN) :: AZI
    REAL*8, DIMENSION(NLR), INTENT(IN) :: TAU, SSA, RTO
    REAL*8, DIMENSION(NGA), INTENT(IN) :: ZAZ
    REAL*8, DIMENSION(NGA, NM), INTENT(IN) :: WCMA, WSMA
    REAL*8, DIMENSION(NK, 4), INTENT(IN) :: XK
    REAL*8, DIMENSION(NUD, NAZ), INTENT(IN) :: F1, F2, R1, R2
!
! DUAL INTENT VARIABLES
    INTEGER, INTENT(INOUT) :: NOS
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NFL), INTENT(OUT) :: FLXI
    REAL*8, DIMENSION(NUD, NAZ), INTENT(OUT) :: II, QQ, UU
!
! LOCAL VARIABLES
    LOGICAL   &
        CONT    ! Compute or skip 2+ orders for the given Fourier moment 
    INTEGER   &
        IA,   & ! Loop index over azimuths
        IB,   & ! Loop index over boundaries
        IG,   & ! Loop index over Gauss nodes
        IL,   & ! Loop index over microlayers
        IM,   & ! Loop index over Fourier orders
        IM_1, & ! IM-1
        IMU,  & ! Loop index over view angles
        ISRF, & ! Copy of ISUR: if (ISRF == 1) ISRF = 0
        NF,   & ! Number of actually computed F-moments, NF <= NM
        NTU,  & ! Number of user-defined transmitted directions
        NT2,  & ! NUP+NG1, last positive Gauss node
        NT3     ! NT2+1, 1st positive user-defined node 
    REAL*8 &
        CAZ,   & ! 2*cos(m*az)
        MU0E0, & ! mu0*exp(-Tau0/mu0)
        SAZ      ! 2*sin(m*az)
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NUD) :: &
        MUD ! Dummy nodes, MUD < 0 comes first (if any)
    REAL*8, DIMENSION(NG1) :: &
        EDNG,  & ! Bouguer attenuation for descending radiation through dTau
        MUWGU, & ! MUG*WG for all mu_gauss < 0
        MUWGD    ! MUG*WG for all mu_gauss > 0
    REAL*8, DIMENSION(NL) :: &
        E0 ! Bouguer attenuation of the solar beam by dTau, exp(-IL*dTau/mu0)
    REAL*8, DIMENSION(NG1, 3) :: &
        RM0 ! m-th Fourier moment of the direct beam reflected from surface
    REAL*8, DIMENSION(NG1, NL) :: &
        EUPG ! Bouguer attenuation for ascending radiation at each level 
    REAL*8, DIMENSION(NUP, NL) :: &
        EUP ! Bouguer attenuation for ascending radiation from each layer
    REAL*8, DIMENSION(NDN, NL) :: &
        EDN ! Bouguer attenuation for descending radiation from each layer
    REAL*8, DIMENSION(NRU, NAZ) :: &
        IS0, QS0, US0 ! Surface correction on TOA
    REAL*8, DIMENSION(NG1, NUP) :: &
        RM11, RM12, RM13, & ! Fourier moments of the Mueller matrix
        RM21, RM22, RM23, & !   for diffuse radiation
        RM31, RM32, RM33
    REAL*8, DIMENSION(NUD, NM) :: &
        I2M, Q2M, U2M ! IQU for higher orders of scattering & all F-moments
    REAL*8, DIMENSION(NGA, NG1) :: &
        R01, R02, R03 ! Direct solar beam reflected from surface
    REAL*8, DIMENSION(NUD, NAZ, NLR) :: &
        Z1, Z2        ! Mixed elements of the phase matrix in each layer
    REAL*8, DIMENSION(NGA, NG1, NUP) :: &
        R11, R12, R13, & ! Precomputed Mueller matrix of surface
        R21, R22, R23, &
        R31, R32, R33
    REAL*8, DIMENSION(NG2) :: &
        ZR01, ZR02, ZR03, & ! Rayleigh scattering of the direct solar beam
        ZA01, ZA02, ZA03    !
    REAL*8, DIMENSION(NMU, 2) :: &
        I2MGU, Q2MGU, U2MGU ! Output of SORD_IP
    REAL*8, DIMENSION(NG2, NLR) :: &
        ZM01, ZM02, ZM03
    REAL*8, DIMENSION(NG2, NMU) :: &
        ZR11, ZR12, ZR13, & !
        ZR21, ZR22, ZR23, & !
        ZR31, ZR32, ZR33, & !
        ZA11, ZA12, ZA13, & !
        ZA21, ZA22, ZA23, & !
        ZA31, ZA32, ZA33
    REAL*8, DIMENSION(NG2, NMU, NLR) :: &
        ZM11, ZM12, ZM13, & !
        ZM21, ZM22, ZM23, & !
        ZM31, ZM32, ZM33
    REAL*8, DIMENSION(NG2, NB) :: &
        I1B, Q1B, U1B
!===============================================================================
!
!   Number of F-moments for summation is not yet defined
    NF = -1
!
    ISRF = ISUR
    NTU = NUD-NRU
    NT2 = NUP+NG1
    NT3 = NT2+1
!
    MUWGU = MUG(:NG1)*WG(:NG1)
    MUWGD = MUG(NG1+1:)*WG(NG1+1:)
!
!   Compute the first order for path radiance and, optionally, surface correction
    IF (NRU > 0) MUD(:NRU) = MU(:NRU)
    IF (NTU > 0) MUD(NRU+1:) = MU(NT3:) ! NT3=NUP+NG1+1
!
    CALL Z0_IP(NC, SSA, RTO, NLR, F1, F2, R1, R2, NUD, NAZ, Z1, Z2)
    CALL FORD_IP(MU0, TAU, MUD, AZI, Z1, Z2, NLR, NLR-1, NUD, NRU, NTU, NAZ, &
                                                                  II, QQ, UU)
!           
!   Compute direct surface reflection on TOA if needed
    IF (NRU > 0 .AND. ISRF > 0) THEN
        CALL SURFTOA_IP(MU0, TAU0, MU(:NRU), NRU, ISRF, PSRF, AZI, NAZ, IS0,   &
                        QS0, US0)
        II(:NRU, :) = II(:NRU, :) + IS0
        IF (ISRF > 5) THEN
            QQ(:NRU, :) = QQ(:NRU, :) + QS0
            UU(:NRU, :) = UU(:NRU, :) + US0
        END IF ! ISRF > 5
    END IF ! NRU > 0
!
    IF ( EPSI > 0.0D0 .OR. (EPSI < 0.0D0 .AND. NOS > 1) ) THEN ! Higher orders
        MU0E0 = MU0*DEXP(-TAU0/MU0)
!
!       Exps for the Bouguer attenuation of the Sun beam & along view directions
        DO IL = 1, NL
            E0(IL) = 0.5D0*(DEXP(-(IL-1)*DTAU/MU0) + DEXP(-IL*DTAU/MU0))
        END DO
        EUP(:, 1) = DEXP(DTAU/MU(:NUP))
        EDN(:, 1) = DEXP(-DTAU/MU(NUP+1:))
        EUPG(:, 1) = DEXP(DTAU/MUG(:NG1))
        EDNG = DEXP(-DTAU/MUG(NG1+1:))
        DO IL = 2, NL
            EUP(:, IL) = EUP(:, IL-1)*EUP(:, 1)
            EDN(:, IL) = EDN(:, IL-1)*EDN(:, 1)
            EUPG(:, IL) = EUPG(:, IL-1)*EUPG(:, 1)
        END DO ! IL = 2, NL
!
!       SURFACE REFLECTION
!
        IF (ISRF == 1) THEN     ! Isotropic reflection
            RM0(:, 1) = 2.0D0*PSRF(1)
            RM11(:, 1) = RM0(1, 1)*MUWGD
            DO IMU = 2, NUP; RM11(:, IMU) = RM11(:, 1); END DO
        ELSE IF (ISRF > 1) THEN ! Anisotropic reflection
            CALL SURFACE0_IP(MU0, ISRF, PSRF, MUG(:NG1), NG1, ZAZ, NGA,        &
                                                         R01, R02, R03)
            CALL SURFACE_IP(ISRF, PSRF, MUG(NG1+1:), NG1, MU(:NUP),              &
                                          NUP, ZAZ, NGA, R11, R12, R13,        &
                                                         R21, R22, R23,        &
                                                         R31, R32, R33)
        END IF ! ISRF = 1
!
!       FOURIER SUMMATION
!
        DO IM = 1, NM
            IM_1 = IM-1
!
!           Rayleigh
            CALL ZRAYM0_IP(IM_1, DEPF, MU0, MUG, NG2, ZR01, ZR02, ZR03)
            CALL ZRAYM_IP(IM_1, DEPF, MUG, WG, NG2, MU, NMU, ZR11, ZR12, ZR13, &
                                                             ZR21, ZR22, ZR23, &
                                                             ZR31, ZR32, ZR33)
!
!           Aerosol
            IF (NC > 1) THEN
                CALL ZAERM0_IP(IM_1, MU0, XK(:, 1), XK(:, 4), NK, MUG, NG2,    &
                                                               ZA01, ZA02, ZA03)
                CALL ZAERM_IP(IM_1, XK, NK, WG, NG2, MU, NMU, NRU,             &
                                                             ZA11, ZA12, ZA13, &
                                                             ZA21, ZA22, ZA23, &
                                                             ZA31, ZA32, ZA33)
            END IF ! NC > 1
!
!           Mix Rayleigh & Aerosol in each optical layer, include SSA/2
            CALL ZM0_IP(IM_1, NC, SSA, RTO, NLR, NG2, ZR01, ZR02, ZR03,        &
                                                       ZA01, ZA02, ZA03,       &
                                                        ZM01, ZM02, ZM03)
            CALL ZM_IP(IM_1, NC, SSA, RTO, NLR, NG2, NMU, ZR11, ZR12, ZR13,    &
                                                          ZR21, ZR22, ZR23,    &
                                                          ZR31, ZR32, ZR33,    &
                                                           ZA11, ZA12, ZA13,   &
                                                           ZA21, ZA22, ZA23,   &
                                                           ZA31, ZA32, ZA33,   &
                                                             ZM11, ZM12, ZM13, &
                                                             ZM21, ZM22, ZM23, &
                                                             ZM31, ZM32, ZM33)
!
!           Surface (ISRF = 1 is processed separately; see above)
            IF (ISRF > 1) THEN
                CALL SURFACE0M_IP(IM, ISRF, R01, R02, R03, WCMA, WSMA, NGA,    &
                                      NG1, NM, RM0(:, 1), RM0(:, 2), RM0(:, 3))
                CALL SURFACEM_IP(IM, ISRF, MUWGD, WCMA, WSMA, NG1,             &
                                               NUP, NGA, NM,                   &
                                                   R11, R12, R13,              &
                                                   R21, R22, R23,              &
                                                   R31, R32, R33,              &
                                                             RM11, RM12, RM13, &
                                                             RM21, RM22, RM23, &
                                                             RM31, RM32, RM33)
            END IF ! ISRF > 1
!
!           Single scattering: path radiance
            CALL FORDM_IP(ZM01, ZM02, ZM03, NEWLR, E0, EUPG(:, 1), EDNG, NL,   &
                                               NB, NLR, NG1, NG2, I1B, Q1B, U1B)
!           Single scattering: surface contribution
            IF (ISRF > 0) THEN
                CALL SURFBNDM_IP(MU0E0, ISRF, RM0, I1B(:NG1, :), Q1B(:NG1, :), &
                                 U1B(:NG1, :), EUPG, NL, NB, NG1)
            END IF ! ISRF > 0
!           Single scattering: fluxes
            IF (IM == 1) THEN
                FLXI(1) = -SUM(I1B(:NG1,  1)*MUWGU) ! TOA up
                FLXI(2) = -SUM(I1B(:NG1, NB)*MUWGU) ! BOA up
                FLXI(3) =  SUM(I1B(NG1+1:, NB)*MUWGD) ! BOA dn
                FLXI(4) =  FLXI(3) + MU0E0          ! BOA dn 1st+direct
            END IF ! IM == 1
!
!           Check magnitude of the 1st order of the current Fourier moment
            IF (EPSI > 0.0D0) THEN ! Automatic mode
                CONT = .FALSE.     ! By default, 1st order is weak
                DO IB = 1, NB
                    DO IG = 1, NG2
                        IF (DABS(I1B(IG, IB)) > EPSI) THEN
                            CONT = .TRUE.
                            EXIT ! Exit the IG-loop
                        END IF ! DABS(I1B(IG, IB)) > EPSI
                    END DO ! IG = 1, NG2
                    IF (CONT) EXIT ! Exit from the IB-loop
                END DO ! IB = 1, NB
            ELSE ! EPSI < 0.0D0 - manual mode
                CONT = .TRUE. ! Compute as many orders as requested
            END IF ! EPSI > 0.0D0      
!
!           Multiple scattering
            IF (CONT) THEN
                NF = IM ! Update number of the actually used F-moments
                CALL SORDM_IP(IM, NOS, ISRF, EPSI, I1B, Q1B, U1B,              &
                                 ZM11, ZM12, ZM13,                             &
                                 ZM21, ZM22, ZM23,                             &
                                 ZM31, ZM32, ZM33,                             &
                                                   RM11, RM12, RM13,           &
                                                   RM21, RM22, RM23,           &
                                                   RM31, RM32, RM33,           &
                                  EUP, EDN(:, 1), NEWLR, NL, NB, NLR, NMU,     &
                                  NUP, NDN, NG1, NG2, I2MGU, Q2MGU, U2MGU)         
                IF (IM == 1) THEN
                    IF (ISRF == 1) ISRF = 0 ! No Lambert for m > 0 (IM > 1)
                    FLXI(1) = FLXI(1) - SUM( I2MGU(NRU+1:NUP, 1)*MUWGU ) ! TOA up
                    FLXI(2) = FLXI(2) - SUM( I2MGU(NRU+1:NUP, 2)*MUWGU ) ! BOA up
                    FLXI(4) = FLXI(4) + SUM( I2MGU(NUP+1:NT2, 2)*MUWGD ) ! total
                    FLXI(3) = FLXI(4) - MU0E0                          ! BOA dn
                END IF ! IM == 1
!
!               Extract user defined zeniths
                IF (NRU > 0) THEN           ! Up
                    I2M(:NRU, IM) = I2MGU(:NRU, 1)
                    Q2M(:NRU, IM) = Q2MGU(:NRU, 1)
                    U2M(:NRU, IM) = U2MGU(:NRU, 1)
                END IF ! NRU > 0
                IF (NTU > 0) THEN           ! Down
                    I2M(NRU+1:, IM) = I2MGU(NT3:, 2)
                    Q2M(NRU+1:, IM) = Q2MGU(NT3:, 2) ! NT3 = NUP+NG1+1
                    U2M(NRU+1:, IM) = U2MGU(NT3:, 2)
                END IF ! NTU > 0
            ELSE ! CONT = FALSE - do not compute 2+ orders for current IM
                EXIT ! Skip all higher IMs if 1st order is weak for the given IM
            END IF ! CONT = TRUE (compute next scattering order for the given IM)
!
!           WRITE(*, *) 'IM = ', IM
        END DO ! IM = 0, NM
!
!       FOURIER SUMMATION AT USER DEFINED NODES
!
!       IM = 1 (m = 0): average over azimuth
        DO IA = 1, NAZ
            II(:, IA) = II(:, IA) + I2M(:, 1)
            QQ(:, IA) = QQ(:, IA) + Q2M(:, 1)
        END DO ! IA = 1, NAZ
!
        DO IM = 2, NF ! NF < NM if the 1st order of NF+1 is negligible
            DO IA = 1, NAZ
!               (2.0 - Kronecker delta) = 2.0 for m > 0 (IM > 1)
                CAZ = 2.0D0*DCOS( (IM-1)*AZI(IA) )
                SAZ = 2.0D0*DSIN( (IM-1)*AZI(IA) )
                II(:, IA) = II(:, IA) + I2M(:, IM)*CAZ
                QQ(:, IA) = QQ(:, IA) + Q2M(:, IM)*CAZ
                UU(:, IA) = UU(:, IA) + U2M(:, IM)*SAZ               
            END DO ! IA = 1, NAZ
        END DO ! IM = 2, NM
! 
    END IF ! EPSI > 0.0D0 .OR. (EPSI < 0.0D0 .AND. NOS > 1)
!
!   THINKME: normalization of fluxes
!
!   Integration over azimuth, m=0, gives 2pi. For the direct beam,
!   2pi is assumed on TOA. Hence, the factor of 2pi must be included
!   in fluxes.
    FLXI = FLXI*PI2
!
END SUBROUTINE SORD_IP
!===============================================================================
! 09Nov16 - MUG is removed from ZAERM_IP as redundant 
!
! 01May16 - All INTEGER variables of the form NGX=NG1+1 are replaced with
!           corresponding explicit expression, NG1+1
!
! 30Apr16 - DCMA & DSMA are removed from input and replaced with CAZ & SAZ.
!           Checked under W+ifort & L+pgf90 - no difference for 53 tests.
!
! 22Apr16 - Minor changes in comments
!
! 20Jan15 - FOURIER SUMMATION AT USER DEFINED NODES was moved inside the
!           IF ( EPSI > 0.0D0 .OR. (EPSI < 0.0D0 .AND. NOS > 1) ) condition
!           (i.e., inside computation of high orders)
!
! 14Nov15 - Fourier accumulation is now in a separate loop. A new variable NF
!           is used; I2M, Q2M, U2M are now 2D arrays storing all IM = 1, NF 
!
! 31Oct15 - SORDM_IP: IM is now on input
!
! 15Oct15 - SURFBNDM_IP added & tested in 11 published tests.
!
! 14Oct15 - Renamed to SORD_IP. Tested.
!
! 30Sep15 - CRIT is added to compute\skip next scattering order
!
! 20Sep15 - Significantly modified: RTE is solved for one mu0; SORD_IP computes
!           IQU at Gauss & user defined nodes; this subroutine now computes
!           fluxes, and so on
!
! 22Jul15 - Several unused local variables & arrays are removed.
!
! 28Apr15 - The subroutine, including input parameters, has been significantly
!           modified and tested using Rayleigh atmosphere over RTLS.
!
! 14Apr15 - SURFACEM0_IP is removed, SRFM0 is redefined. Tested for RTLS.
!
! 02Apr15 - Timing and several minor modifications.
!
! 20Mar15 - SORD_IP now provides TIME on output.
!
! 10Mar15 - IF (NO > 1) -> IF (NO > 1 .OR. CRIT == 2)
!
! 06Mar15 - IF (ISRF == 1) ISRF = 0 has been moved from the body of the
!           DO IS = 1, NS loop.
!
! 03Mar15 - Lambertian surface, IF (ISRF == 1), the following lines are added:
!           SRFM0(2:NG3:3) = 0.0D0, SRFM0(3:NG3:3) = 0.0D0. All tests have been
!           repeated. Ok.
!
! 21Feb15 - An error (typo) was fixed: IF (NRU > 0 .AND. ISRF > 1) replaced with
!           IF (NRU > 0 .AND. ISRF > 0) - to include ISRF=1 - Lambert. Tested
!           vs RT3.
!
! 19Feb15 - FLUXI is added.
!
! 16Feb15 - NEWLR is now an INTEGER array. Tested vs previous version. Ok.
!
! 14Feb15 - cos(m*phi) & sin(m*phi) are now computed with trigonometric addition
!           formulas. The 2D arrays CMAZ, SMAZ and real variables IM8 and MAZ
!           are removed. Tested against previous version - ok.
!
! 11Feb15 - Compute direct surface reflection on TOA if needed: IF (NRU > 0) ->
!           IF (NRU > 0 .AND. ISRF > 1). Skip computation for a black surface.
!           
!           IF (NO > 1) THEN ... END IF was added to avoid computation of higher
!           orders when only 1st order is needed.
!
! 28Dec14 - MICROPROFILE was removed, NEWLR is now on input. Tested against
!           old version for the same input. Ok.
!
! 23Nov14 - a4k is not used. Renamed: NLR5 -> NLR4. Added: parameter NX = 4.
!           Tested against old version.
!
! 13Nov14 - Automatic criterion is now used as well as the manual mode. Tested
!           against previous results for 10 orders. See SORD_IP for details.
!
! 30Oct14 - TAU0 & DTAU are now on input. CALL MICROPROFILE(...) was modified.
!           Tested as part of SORD.
!
! 28Oct14 - NKS is now used as the input parameter in SORD_IP (see tests there)
!
! 11Oct14 - Renamed: NPAR -> NSPR. Checked for ISRF = 4 against old version.
!
! 09Oct14 - Order of input variables was changed. Tested against old format for
!           RTLS and 6 layers.
!
! 08Oct14 - Tests for 2, 3, 5, 6, and 10 layers with optical depths as below:
!               Tau2 = [0.1 0.1];
!               Tau3 = [0.05 0.1 0.05];
!               Tau5 = [0.05 0.02 0.06 0.02 0.05];
!               Tau6_asc = [1 2 3 5 4 5]/100;
!               Tau6_des = [5 4 5 3 2 1]/100;
!               Tau10 = [2 2 2 2 2 2 2 2 2 2]/100.
!           Compared against SORD, single layer with Tau0 = 0.2. Ok.
!
! 07Oct14 - Tests from SORD_IP for ISRF = 0, 1 (ro=0.1 & 0.9), 4, 6, 10 were
!           repeated. Same results observed.
!
! 06Oct14 - Created and tested against previous version SORD for 1 layer case
!           and RTLS surface. Ok.
!===============================================================================