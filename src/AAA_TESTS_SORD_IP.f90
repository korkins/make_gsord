PROGRAM AAA_TESTS_SORD_IP
!===============================================================================
! PURPOSE:
!   (a) To compare SORD_IP against published results.
!   (b) To provide examples of input for different scenarious.
!
! COMMENTS:
!   Refer to RT_TESTS.xlsx for a complete list of tests with brief description.
!
!   Meaning of all parameters are explained in each test.
!   'N/U' stands for a parameter which is used in a particular test.
!
!   When installed, update the PATH & LPATH variables below.
!===============================================================================
!
! DECLARATION OF VARIABLES
    IMPLICIT NONE
!
! PARAMETERS
    INTEGER, PARAMETER :: &
        LPATH = 21, & ! Lenght of the PATH variable
        LLOGF =  7, & ! Length of the log-file name & extension (log.txt = 7)
        LTXTF = 13, & ! Lenghth of 'TXXX/TXXX.txt'
        LLOGP = LPATH+LLOGF, & ! Length of the full path to the log file
        LTXTP = LPATH+LTXTF ! Length of the full path to benchmark results
    CHARACTER(LEN=LPATH), PARAMETER :: PATH = '/home/sk/gsord/TESTS/'
    CHARACTER(LEN=LLOGP), PARAMETER :: FNAME = PATH//'log.txt'
    CHARACTER(LEN=LLOGP), PARAMETER :: FSCRN = PATH//'scr.txt'
    INTEGER, PARAMETER :: &
        FNUM_XK = 11, & ! Phase function expansion moments file ID
        FNUM_BM = 12, & ! Benchmark data
        FNUM = 13,    & ! Report (log) file ID
        FSCR = 14,    & ! Print screen in a txt file
        NTST = 99       ! Maximum number of tests
    INTEGER, PARAMETER :: & ! Maximum number of ...
        NAZ_MAX = 181,  & ! ... view azimuths
        NGA_MAX = 360,  & ! ... Gauss azimuths for the surface Furier expansion
        NL_MAX  = 1000, & ! ... levels over heights (layers, boundaries, etc)
        NUD_MAX = 180,  & ! ... user-defined angles
        NG1_MAX = 1024, & ! ... Gauss nodes in HEMIsphere
        NG2_MAX = 2*NG1_MAX, & ! ... Gauss nodes in sphere
        NMU_MAX = NG2_MAX+NUD_MAX, & ! ... Gauss & user nodes
        NK_MAX = 2048,  &  ! ... phase function expansion moments, xk
        NX_MAX = 8,     &  ! ... columns in benchmark file (fixed)
        NY_MAX = NAZ_MAX*NUD_MAX, & ! ... lines in benchmark file
        NM_MAX = NG1_MAX ! ... Fourier components
    REAL*8, PARAMETER :: &
        PI  = DACOS(-1.0D0), & ! The number pi
        PI2 = 2.0D0*PI,      & ! The number 2pi
        D2R = PI/180.0D0,    & ! degrees-> radians, pi/180
        R2D = 180.0D0/PI,    & ! radians-> degrees, 180/pi
        TINY = 1.0D-12
!
! LOCAL SCALARS
    LOGICAL &
        ISFIRST, NEXT_IKM
    INTEGER &
        IA, IZ, ISRF, NKS, NG1, NG2, NGA, NUP, NDN, NUD, NMU, NL, &
        NLR, NML, NM, NO, NAZ, NRU, NTU, NK, ITEST, NX, NY, BEG, FIN,     &
        IBR, NC, IL, IL1, IL2, ILR, NLR0, NZ, IKM, NBR, IK, IY, NB, IRUN, &
        NKM
    REAL*4 &
        CPU_TIME1, CPU_TIME2, TIME
    REAL*8 &
        SZA, MU0, DEPF, TAU0, DTAU, EPSI, DELT, TAU_ACC, TAUG, TAU_DIF,  &
        DX, I_BMRK, Q_BMRK, U_BMRK, P_BMRK, I_SORD, Q_SORD, U_SORD, P_SORD,    &
        MAX_ERR_I, AVR_ERR_I, MAX_ERR_P, AVR_ERR_P, SSA_AER, TAU_SCA_AR
!
! LOCAL ARRAYS
    LOGICAL &
        TSTNUM(NTST)
    INTEGER &
        NEWLR(NL_MAX)
    REAL*8 & ! 1D arrays with fixed size
        PSRF(21), FLUXI(4)
    REAL*8 & ! 1D arrays
        ZAZ(NGA_MAX), WAZ(NGA_MAX), MUD(NUD_MAX), TAU(NL_MAX), SSA(NL_MAX),  &
        MUWG(NG1_MAX), Z(NG2_MAX), W(NG2_MAX), AZA(NAZ_MAX), AZI(NAZ_MAX),   &
        UDA(NUD_MAX), MU(NMU_MAX), RATIO(NL_MAX), ZKM(NL_MAX), XKM(NL_MAX), &
        TAU_RAY(NL_MAX), TAU_GAS(NL_MAX), TAU_ZKM(NL_MAX),                  &
        TAU_AER(NL_MAX), TAU_SCA(NL_MAX), TAU_SCA_RAY(NL_MAX),              &
        TAU_SCA_AER(NL_MAX)
    REAL*8 & ! 2D arrays
        WCMA(NGA_MAX, NM_MAX), WSMA(NGA_MAX, NM_MAX),   &
        FK(NK_MAX, 7), XK(NK_MAX, 4),                   &
        BMARK(NUD_MAX*NAZ_MAX, 8), OUT(NUD_MAX*NAZ_MAX, 15)
    REAL*8, DIMENSION(NUD_MAX, NAZ_MAX) :: & ! 3D arrays
        II, QQ, UU, F11, F12, R11, R12, MUS
!===============================================================================
!0
!   Select what tests to run: 0 - test OFF; 1 - test ON
    TSTNUM = 0
    TSTNUM(1:20) = 1 !1:25
!0
!   *** Skip these ***
    TSTNUM(29) = 0   ! Switch OFF: takes a long time to run, ~50min
    TSTNUM(47) = 0   ! Switch OFF: needs 1.5Gb RAM, takes ~15min, not IPRT bmrk
!0
!    TSTNUM(56) = 1
!0
!   Open file for screen printing
    OPEN(FSCR, FILE = FSCRN)
    WRITE(FSCR, 20) '|dI%|', '<|dI%|>', '|dP|', '<|dP|>', 'time, sec.'
    WRITE(FSCR, *)
!0
!   Open the log file. By default, an existing file is replaced
    OPEN(FNUM, FILE = FNAME)
    WRITE(FNUM, 10) 'SZA      mu0         Az       VZA      mu        ' //     &
                    'dI %    dP        Ib               Is               ' //  &
                    'Pb        Ps        Qb               Qs               ' //&
                    'Ub               Us'
!0
    WRITE(*, 20) '|dI%|', '<|dI%|>', '|dP|', '<|dP|>', 'time, sec.'
    WRITE(*, *)
!0
DO IRUN = 1, 1 ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
WRITE(FSCR, *) 'IRUN = ', IRUN
WRITE(*, *) 'IRUN = ', IRUN
!0
!*******************************************************************************
!
!   001. GOAL: To test single scattering approximation
!
!*******************************************************************************
!1
    ITEST = 1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!1
!       ACCURACY
!1
!       Number of Gauss nodes per hemisphere (N/U)
        NG1 = 2
!       Total number of the Fourier moments (N/U)
        NM = 1
!       Number of microlayers, defines dTau (N/U)
        NL = 1
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = -999.0D0
!       Number of scattering orders to be computed: single scattering 
        NO = 1
!1
!       GEOMETRY
!1
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.6D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 1       ! Number of reflected angles
        NTU = NRU     ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NUD) = (/SZA, 180.0D0 - SZA/) ! Reflection first
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!1
!       ATMOSPHERE
!1
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.2D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 117  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T001/Xk0117_0804.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!1
!       SURFACE
!1
!       Surface model (N/U)
        ISRF = 0
!       Set of parameters (N/U)
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion (N/U)
        NGA = 2
!1
!       COMPUTE AUXILIARY PARAMETERS
!1
!       Number of Gauss nodes per whole sphere (N/U)
        NG2 = NG1*2
!       Total number of nodes for ascending radiation (N/U)
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation (N/U)
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy (N/U)
        NMU = NG2+NUD
!       Number of moments in the scattering integral (N/U)
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau (N/U)
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!1
!       RUN RT CODE
!1
        CALL CPU_TIME(CPU_TIME1)
!1
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!1
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!1
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM),  &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!1
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!1
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!1
        CALL CPU_TIME(CPU_TIME2)
!1
        TIME = CPU_TIME2 - CPU_TIME1
!1
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T001/T001.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!1
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!1
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!1
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!1
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!1
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!1
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!1
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!1
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 001:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!1
        WRITE(*, 50) 'T001: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T001: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!1
    END IF ! TESTNUM(ITEST)
!1
!*******************************************************************************
!
!   002. GOAL: To test double scattering approximation
!
!*******************************************************************************
!2
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!2
!       ACCURACY
!2
!       Number of Gauss nodes per hemisphere
        NG1 = 24
!       Total number of the Fourier moments
        NM = 16
!       Number of microlayers, defines dTau
        NL = 10
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = -999.0D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!2
!       GEOMETRY
!2
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.6D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 1       ! Number of reflected angles
        NTU = NRU     ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NUD) = (/SZA, 180.0D0 - SZA/) ! Reflection first
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!2
!       ATMOSPHERE
!2
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.2D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 117  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T002/Xk0117_0804.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!2
!       SURFACE
!2
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!2
!       COMPUTE AUXILIARY PARAMETERS
!2
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!2
!       RUN RT CODE
!2
        CALL CPU_TIME(CPU_TIME1)
!2
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!2
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!2
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!2
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!2
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!2
        CALL CPU_TIME(CPU_TIME2)
!2
        TIME = CPU_TIME2 - CPU_TIME1
!2
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T002/T002.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!2
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!2
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!2
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!2
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!2
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!2
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!2
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!2
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 002:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!2
        WRITE(*, 50) 'T002: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T002: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!2
    END IF ! TESTNUM(ITEST)
!2
!*******************************************************************************
!
!   003. GOAL: To test double scattering for reflection only (NTU=0)
!
!*******************************************************************************
!3
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!3
!       ACCURACY
!3
!       Number of Gauss nodes per hemisphere
        NG1 = 24
!       Total number of the Fourier moments
        NM = 16
!       Number of microlayers, defines dTau
        NL = 10
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = -999.0D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!3
!       GEOMETRY
!3
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.6D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 1 ! Number of reflected angles
        NTU = 0 ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NUD) = SZA ! Reflection first
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!3
!       ATMOSPHERE
!3
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.2D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 117 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T003/Xk0117_0804.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!3
!       SURFACE
!3
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!3
!       COMPUTE AUXILIARY PARAMETERS
!3
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!3
!       RUN RT CODE
!3
        CALL CPU_TIME(CPU_TIME1)
!3
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!3
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!3
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!3
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!3
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!3
        CALL CPU_TIME(CPU_TIME2)
!3
        TIME = CPU_TIME2 - CPU_TIME1
!3
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T003/T003.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!3
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!3
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!3
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!3
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!3
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!3
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!3
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!3
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 003:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!3
        WRITE(*, 50) 'T003: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T003: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!3
    END IF ! TESTNUM(ITEST)
!3
!*******************************************************************************
!
!   004. GOAL: To test double scattering for transmission only (NRU=0)
!
!*******************************************************************************
!4
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!4
!       ACCURACY
!4
!       Number of Gauss nodes per hemisphere
        NG1 = 24
!       Total number of the Fourier moments
        NM = 16
!       Number of microlayers, defines dTau
        NL = 10
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = -999.0D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!4
!       GEOMETRY
!4
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.6D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 0 ! Number of reflected angles
        NTU = 1 ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NUD) = 180.0D0 - SZA ! Reflection first
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!4
!       ATMOSPHERE
!4
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.2D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 117  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T004/Xk0117_0804.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!4
!       SURFACE
!4
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!4
!       COMPUTE AUXILIARY PARAMETERS
!4
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau (N/U)
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!4
!       RUN RT CODE
!4
        CALL CPU_TIME(CPU_TIME1)
!4
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!4
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!4
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!4
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!4
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!4
        CALL CPU_TIME(CPU_TIME2)
!4
        TIME = CPU_TIME2 - CPU_TIME1
!4
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T004/T004.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!4
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!4
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!4
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!4
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!4
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!4
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!4
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!4
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 004:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!4
        WRITE(*, 50) 'T004: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T004: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!4
    END IF ! TESTNUM(ITEST)
!4
!*******************************************************************************
!
!   005. GOAL: To test aerosol scattering, high scattering orders and
!              automatic convergence.
!
!*******************************************************************************
!5
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!5
!       ACCURACY
!5
!       Number of Gauss nodes per hemisphere
        NG1 = 24
!       Total number of the Fourier moments
        NM = 16
!       Number of microlayers, defines dTau
        NL = 10
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!5
!       GEOMETRY
!5
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.6D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 1   ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NUD) = (/SZA, 180.0D0 - SZA/) ! Reflection first
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!5
!       ATMOSPHERE
!5
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.2D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 117  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T005/Xk0117_0804.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!5
!       SURFACE
!5
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!5
!       COMPUTE AUXILIARY PARAMETERS
!5
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!5
!       RUN RT CODE
!5
        CALL CPU_TIME(CPU_TIME1)
!5
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!5
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!5
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!5
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!5
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!5
        CALL CPU_TIME(CPU_TIME2)
!5
        TIME = CPU_TIME2 - CPU_TIME1
!5
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T005/T005.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!5
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!5
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!5
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!5
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!5
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!5
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!5
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!5
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 005:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!5
        WRITE(*, 50) 'T005: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T005: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!5
    END IF ! TESTNUM(ITEST)
!5
!*******************************************************************************
!
!   006. GOAL: To test Rayleigh scattering.
!
!*******************************************************************************
!6
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!6
!       ACCURACY
!6
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 20
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!6
!       GEOMETRY
!6
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 90 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(0.0D0 + (IK - 1.0D0), IK = 1, NRU)/) ! Reflection
        UDA(NRU+1:NUD) = 91.0D0 + UDA(1:NRU)                 ! Transmittance
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!6
!       ATMOSPHERE
!6
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.3262D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!6
!       SURFACE
!6
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!6
!       COMPUTE AUXILIARY PARAMETERS
!6
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!6
!       RUN RT CODE
!6
        CALL CPU_TIME(CPU_TIME1)
!6
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!6
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!6
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!6
        CALL CPU_TIME(CPU_TIME2)
!6
        TIME = CPU_TIME2 - CPU_TIME1
!6
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T006/T006.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!6
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!6
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!6
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!6
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!6
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!6
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!6
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!6
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 006:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!6
        WRITE(*, 50) 'T006: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T006: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!6
    END IF ! TESTNUM(ITEST)
!6
!*******************************************************************************
!
!   007. GOAL: To test fluxes.
!
!*******************************************************************************
!7
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!7
!       ACCURACY
!7
!       Number of Gauss nodes per hemisphere
        NG1 = 16
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!7
!       GEOMETRY
!7
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 30.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 16 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        CALL GAUSZW(0.0D0, 1.0D0, NRU, MUD(:NRU), UDA(:NRU)) ! UDA is redefined
        UDA(1:NRU) = DACOS(MUD(NRU:1:-1))*R2D          ! Reflectance
        UDA(NRU+1:NUD) = DACOS(-MUD(1:NRU))*R2D ! Transmittance
!       View azimuth angles, degrees
        NAZ = 5 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 45.0D0, 90.0D0, 135.0D0, 180.0D0/)
!7
!       ATMOSPHERE
!7
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.0D0
        SSA(1:NLR) = 0.9D0  
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!7
!       SURFACE
!7
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!7
!       COMPUTE AUXILIARY PARAMETERS
!7
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!7
!       RUN RT CODE
!7
        CALL CPU_TIME(CPU_TIME1)
!7
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!7
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!7
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!7
!       NORMALIZATION
        II = II/PI2/MU0
        QQ = QQ/PI2/MU0
        UU = UU/PI2/MU0
        FLUXI = FLUXI/PI2/MU0
!7
        CALL CPU_TIME(CPU_TIME2)
!7
        TIME = CPU_TIME2 - CPU_TIME1
!7
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T007/T007.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!7
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!7
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!7
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!7
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!7
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!7
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!7
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!7
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 007:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
!7
!       Compare fluxes
        WRITE(FNUM, *)
        WRITE(FNUM, 10) '           RT3             SORD            dF%'
        WRITE(FNUM, 60) 'F_up_dif', 0.294876D0, FLUXI(1), &
                        100.0D0*(1.0D0 - FLUXI(1)/0.294876D0)
        WRITE(FNUM, 60) 'F_dn_tot', 0.556938D0, FLUXI(4), &
                        100.0D0*(1.0D0 - FLUXI(4)/0.556938D0)
!7
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!7
        WRITE(*, 50) 'T007: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T007: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!7
    END IF ! TESTNUM(ITEST)
!7
!*******************************************************************************
!
!   008. GOAL: To test Lambertian surface reflection and fluxes.
!
!*******************************************************************************
!8
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!8
!       ACCURACY
!8
!       Number of Gauss nodes per hemisphere
        NG1 = 16
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 25
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!8
!       GEOMETRY
!8
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 45.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 16 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        CALL GAUSZW(0.0D0, 1.0D0, NRU, MUD(:NRU), UDA(:NRU)) ! UDA is redefined
        UDA(1:NRU) = DACOS(MUD(NRU:1:-1))*R2D          ! Reflectance
        UDA(NRU+1:NUD) = DACOS(-MUD(1:NRU))*R2D ! Transmittance
!       View azimuth angles, degrees
        NAZ = 5 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 45.0D0, 90.0D0, 135.0D0, 180.0D0/)
!8
!       ATMOSPHERE
!8
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.5D0
        SSA(1:NLR) = 0.9D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!8
!       SURFACE
!8
!       Surface model
        ISRF = 1
!       Set of parameters: Lambertian reflection, ro = PSRF(1)
        PSRF(1) = 0.7D0
        PSRF(2:) = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!8
!       COMPUTE AUXILIARY PARAMETERS
!8
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!8
!       RUN RT CODE
!8
        CALL CPU_TIME(CPU_TIME1)
!8
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!8
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!8
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!8
!       NORMALIZATION
        II = II/PI2/MU0
        QQ = QQ/PI2/MU0
        UU = UU/PI2/MU0
        FLUXI = FLUXI/PI2/MU0
!8
        CALL CPU_TIME(CPU_TIME2)
!8
        TIME = CPU_TIME2 - CPU_TIME1
!8
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T008/T008.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!8
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!8
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!8
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!8
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!8
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!8
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!8
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!8
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 008:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
!8
!       Compare fluxes
        WRITE(FNUM, *)
        WRITE(FNUM, 10) '           RT3             SORD            dF%'
        WRITE(FNUM, 60) 'F_up_dif', 0.609266D0, FLUXI(1), &
                        100.0D0*(1.0D0 - FLUXI(1)/0.609266D0)
        WRITE(FNUM, 60) 'F_dn_tot', 0.841381D0, FLUXI(4), &
                        100.0D0*(1.0D0 - FLUXI(4)/0.841381D0)
!8
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!8
        WRITE(*, 50) 'T008: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T008: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!8
    END IF ! TESTNUM(ITEST)
!8
!*******************************************************************************
!
!   009. GOAL: To test Lambertian surface reflection and fluxes.
!
!*******************************************************************************
!9
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!9
!       ACCURACY
!9
!       Number of Gauss nodes per hemisphere
        NG1 = 16
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 5
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!9
!       GEOMETRY
!9
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 75.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 16 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        CALL GAUSZW(0.0D0, 1.0D0, NRU, MUD(:NRU), UDA(:NRU)) ! UDA is redefined
        UDA(1:NRU) = DACOS(MUD(NRU:1:-1))*R2D          ! Reflectance
        UDA(NRU+1:NUD) = DACOS(-MUD(1:NRU))*R2D ! Transmittance
!       View azimuth angles, degrees
        NAZ = 5 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 45.0D0, 90.0D0, 135.0D0, 180.0D0/)
!9
!       ATMOSPHERE
!9
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.1D0
        SSA(1:NLR) = 0.9D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!9
!       SURFACE
!9
!       Surface model
        ISRF = 1
!       Set of parameters: Lambertian reflection, ro = PSRF(1)
        PSRF(1) = 0.7D0
        PSRF(2:) = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!9
!       COMPUTE AUXILIARY PARAMETERS
!9
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!9
!       RUN RT CODE
!9
        CALL CPU_TIME(CPU_TIME1)
!9
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!9
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!9
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!9
!       NORMALIZATION
        II = II/PI2/MU0
        QQ = QQ/PI2/MU0
        UU = UU/PI2/MU0
        FLUXI = FLUXI/PI2/MU0
!9
        CALL CPU_TIME(CPU_TIME2)
!9
        TIME = CPU_TIME2 - CPU_TIME1
!9
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T009/T009.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!9
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!9
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!9
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!9
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!9
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!9
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!9
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!9
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 009:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
!9
!       Compare fluxes
        WRITE(FNUM, *)
        WRITE(FNUM, 10) '           RT3             SORD            dF%'
        WRITE(FNUM, 60) 'F_up_dif', 0.690810D0, FLUXI(1), &
                        100.0D0*(1.0D0 - FLUXI(1)/0.690810D0)
        WRITE(FNUM, 60) 'F_dn_tot', 0.864513D0, FLUXI(4), &
                        100.0D0*(1.0D0 - FLUXI(4)/0.864513D0)
!9
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!9
        WRITE(*, 50) 'T009: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T009: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!9
    END IF ! TESTNUM(ITEST)
!9
!*******************************************************************************
!
!   010. GOAL: To test ocean (NT model) reflection.
!
!*******************************************************************************
!10
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!10
!       ACCURACY
!10
!       Number of Gauss nodes per hemisphere
        NG1 = 10
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 10
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!10
!       GEOMETRY
!10
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 9 ! Number of reflected angles
        NTU = 0 ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(0.0D0 + 10.0D0*(IK-1), IK = 1, NRU)/) ! Reflection
!       View azimuth angles, degrees
        NAZ = 2 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 45.0D0/)
!10
!       ATMOSPHERE
!10
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.1D0
        SSA(1:NLR) = 0.9D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!10
!       SURFACE
!10
!       Surface model
        ISRF = 6
!       Set of parameters
        PSRF(1) = 1.34D0 ! Refractive index of water, real part
        PSRF(2) = 0.00D0 ! Refractive index of water, imaginary part
        PSRF(3) = 1.00D0 ! Wind speed, m/s, at 10m above the surface
        PSRF(4:) = 0.0D0 ! Other surface parameters are not used   
!       Number of nodes for the Fourier expansion
        NGA = 180
!10
!       COMPUTE AUXILIARY PARAMETERS
!10
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!10
!       RUN RT CODE
!10
        CALL CPU_TIME(CPU_TIME1)
!10
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!10
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!10
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!10
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!10
        CALL CPU_TIME(CPU_TIME2)
!10
        TIME = CPU_TIME2 - CPU_TIME1
!10
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T010/T010.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!10
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!10
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!10
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!10
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!10
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!10
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!10
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!10
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 010:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!10
        WRITE(*, 50) 'T010: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T010: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!10
    END IF ! TESTNUM(ITEST)
!10
!*******************************************************************************
!
!   011. GOAL: To test ocean (NT model) reflection.
!
!*******************************************************************************
!11
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!11
!       ACCURACY
!11
!       Number of Gauss nodes per hemisphere
        NG1 = 10
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 20
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!11
!       GEOMETRY
!11
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 9 ! Number of reflected angles
        NTU = 0 ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(0.0D0 + 10.0D0*(IK-1), IK = 1, NRU)/) ! Reflection
!       View azimuth angles, degrees
        NAZ = 2 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 45.0D0/)
!11
!       ATMOSPHERE
!11
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.5D0
        SSA(1:NLR) = 0.9D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!11
!       SURFACE
!11
!       Surface model
        ISRF = 6
!       Set of parameters
        PSRF(1) = 1.34D0 ! Refractive index of water, real part
        PSRF(2) = 0.00D0 ! Refractive index of water, imaginary part
        PSRF(3) = 1.00D0 ! Wind speed, m/s, at 10m above the surface
        PSRF(4:) = 0.0D0 ! Other surface parameters are not used   
!       Number of nodes for the Fourier expansion
        NGA = 180
!11
!       COMPUTE AUXILIARY PARAMETERS
!11
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!11
!       RUN RT CODE
!11
        CALL CPU_TIME(CPU_TIME1)
!11
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!11
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!11
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!11
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!11
        CALL CPU_TIME(CPU_TIME2)
!11
        TIME = CPU_TIME2 - CPU_TIME1
!11
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T011/T011.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!11
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!11
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!11
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!11
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!11
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!11
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!11
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!11
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 011:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!11
        WRITE(*, 50) 'T011: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T011: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!11
    END IF ! TESTNUM(ITEST)
!11
!*******************************************************************************
!
!   012. GOAL: To test 2 optical layers, mixing of Rayleigh and aerosol, and the
!              depolarization factor.
!
!*******************************************************************************
!12
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!12
!       ACCURACY
!12
!       Number of Gauss nodes per hemisphere
        NG1 = 24
!       Total number of the Fourier moments
        NM = 16
!       Number of microlayers, defines dTau
        NL = 30
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!12
!       GEOMETRY
!12
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.5D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 4   ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU)     = DACOS(-(/-1.0D0, -0.8D0, -0.5D0, -0.2D0/))*R2D ! Upward
        UDA(NRU+1:NUD) = DACOS(-(/0.2D0, 0.5D0, 0.8D0, 1.0D0/))*R2D     ! Downward
!       View azimuth angles, degrees
        NAZ = 1 ! Number of azimuths
        AZA(1:NAZ) = 0.0D0
!12
!       ATMOSPHERE
!12
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0279D0
!       Number of optical layers
        NLR = 2
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = (/0.1D0, 0.5D0/)
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 117  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T012/Xk0117_0804.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!12
!       SURFACE
!12
!       Surface model
        ISRF = 1
!       Set of parameters
        PSRF(1) = 0.1D0
        PSRF(2:) = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!12
!       COMPUTE AUXILIARY PARAMETERS
!12
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!12
!       Mask microlayers in the optical layer No.2
        NEWLR(1:NL) = 2
!       Mask microlayers in the optical layer No.1
        NEWLR(1:INT(TAU(1)/DTAU)) = 1
!12
!       RUN RT CODE
!12
        CALL CPU_TIME(CPU_TIME1)
!12
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!12
!       No Aerosol in the top layer
        RATIO(1) = 0.0D0
!       Tau_a = 0.4; Tau_r = 0.1; Tau_abs = 0
        RATIO(2) = 0.4D0/(0.4D0 + 0.1D0)
!12
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!12
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!12
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!12
        CALL CPU_TIME(CPU_TIME2)
!12
        TIME = CPU_TIME2 - CPU_TIME1
!12
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T012/T012.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!12
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!12
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!12
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!12
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!12
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!12
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!12
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!12
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 012:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!12
        WRITE(*, 50) 'T012: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T012: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!12
    END IF ! TESTNUM(ITEST)
!12
!*******************************************************************************
!
!   013. GOAL: To test view directions close to horizon, mu = +/-0.1 (VZA ~= 84)
!
!*******************************************************************************
!13
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!13
!       ACCURACY
!13
!       Number of Gauss nodes per hemisphere
        NG1 = 24
!       Total number of the Fourier moments
        NM = 16
!       Number of microlayers, defines dTau
        NL = 100
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!13
!       GEOMETRY
!13
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.5D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 3   ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS((/1.0D0, 0.5D0, 0.1D0/))*R2D      ! Up
        UDA(NRU+1:NUD) = DACOS(-(/0.1D0, 0.5D0, 1.0D0/))*R2D ! Down
!       View azimuth angles, degrees
        NAZ = 2 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 30.0D0/)
!13
!       ATMOSPHERE
!13
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.0D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 117  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T013/Xk0117_0804.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!13
!       SURFACE
!13
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!13
!       COMPUTE AUXILIARY PARAMETERS
!13
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!13
!       RUN RT CODE
!13
        CALL CPU_TIME(CPU_TIME1)
!13
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!13
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!13
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!13
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!13
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!13
        CALL CPU_TIME(CPU_TIME2)
!13
        TIME = CPU_TIME2 - CPU_TIME1
!13
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T013/T013.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!13
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!13
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!13
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!13
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!13
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!13
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!13
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!13
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 013:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!13
        WRITE(*, 50) 'T013: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T013: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!13
    END IF ! TESTNUM(ITEST)
!13
!*******************************************************************************
!
!   014. GOAL: To test simulation for the low Sun, SZA ~= 84.26.
!
!*******************************************************************************
!14
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!14
!       ACCURACY
!14
!       Number of Gauss nodes per hemisphere
        NG1 = 24
!       Total number of the Fourier moments
        NM = 16
!       Number of microlayers, defines dTau
        NL = 100
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!14
!       GEOMETRY
!14
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.1D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 3   ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS((/1.0D0, 0.5D0, 0.1D0/))*R2D      ! Up
        UDA(NRU+1:NUD) = DACOS(-(/0.1D0, 0.5D0, 1.0D0/))*R2D ! Down
!       View azimuth angles, degrees
        NAZ = 2 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 30.0D0/)
!14
!       ATMOSPHERE
!14
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.0D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 117  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T014/Xk0117_0804.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!14
!       SURFACE
!14
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!14
!       COMPUTE AUXILIARY PARAMETERS
!14
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!14
!       RUN RT CODE
!14
        CALL CPU_TIME(CPU_TIME1)
!14
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!14
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!14
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!14
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!14
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!14
        CALL CPU_TIME(CPU_TIME2)
!14
        TIME = CPU_TIME2 - CPU_TIME1
!14
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T014/T014.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!14
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!14
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!14
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!14
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!14
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!14
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!14
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!14
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 014:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!14
        WRITE(*, 50) 'T014: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T014: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!14
    END IF ! TESTNUM(ITEST)
!14
!*******************************************************************************
!
!   015. GOAL: To test 2 optical layers and the U-component (missing in 012)
!
!*******************************************************************************
!15
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!15
!       ACCURACY
!15
!       Number of Gauss nodes per hemisphere
        NG1 = 24
!       Total number of the Fourier moments
        NM = 16
!       Number of microlayers, defines dTau
        NL = 30
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!15
!       GEOMETRY
!15
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.5D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 3   ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS((/1.0D0, 0.5D0, 0.1D0/))*R2D      ! Up
        UDA(NRU+1:NUD) = DACOS(-(/0.1D0, 0.5D0, 1.0D0/))*R2D ! Down
!       View azimuth angles, degrees
        NAZ = 2 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 30.0D0/)
!15
!       ATMOSPHERE
!15
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0279D0
!       Number of optical layers
        NLR = 2
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = (/0.1D0, 0.5D0/)
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 117  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T015/Xk0117_0804.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!15
!       SURFACE
!15
!       Surface model
        ISRF = 1
!       Set of parameters
        PSRF(1) = 0.1D0
        PSRF(2:) = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!15
!       COMPUTE AUXILIARY PARAMETERS
!15
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!15
!       Mask microlayers in the optical layer No.2
        NEWLR(1:NL) = 2
!       Mask microlayers in the optical layer No.1
        NEWLR(1:INT(TAU(1)/DTAU)) = 1
!15
!       RUN RT CODE
!15
        CALL CPU_TIME(CPU_TIME1)
!15
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!15
!       No Aerosol in the top layer
        RATIO(1) = 0.0D0
!       Tau_a = 0.4; Tau_r = 0.1; Tau_abs = 0
        RATIO(2) = 0.4D0/(0.4D0 + 0.1D0)
!15
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!15
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!15
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!15
        CALL CPU_TIME(CPU_TIME2)
!15
        TIME = CPU_TIME2 - CPU_TIME1
!15
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T015/T015.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!15
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!15
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!15
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!15
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!15
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!15
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!15
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!15
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 015:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!15
        WRITE(*, 50) 'T015: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T015: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!15
    END IF ! TESTNUM(ITEST)
!15
!*******************************************************************************
!
!   016. GOAL: To test 2 optical layers and the low Sun, SZA~=84.26.
!
!*******************************************************************************
!16
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!16
!       ACCURACY
!16
!       Number of Gauss nodes per hemisphere
        NG1 = 24
!       Total number of the Fourier moments
        NM = 16
!       Number of microlayers, defines dTau
        NL = 60
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!16
!       GEOMETRY
!16
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.1D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 3   ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS((/1.0D0, 0.5D0, 0.1D0/))*R2D      ! Up
        UDA(NRU+1:NUD) = DACOS(-(/0.1D0, 0.5D0, 1.0D0/))*R2D ! Down
!       View azimuth angles, degrees
        NAZ = 2 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 30.0D0/)
!16
!       ATMOSPHERE
!16
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0279D0
!       Number of optical layers
        NLR = 2
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = (/0.1D0, 0.5D0/)
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 117  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T016/Xk0117_0804.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!16
!       SURFACE
!16
!       Surface model
        ISRF = 1
!       Set of parameters
        PSRF(1) = 0.1D0
        PSRF(2:) = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!16
!       COMPUTE AUXILIARY PARAMETERS
!16
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!16
!       Mask microlayers in the optical layer No.2
        NEWLR(1:NL) = 2
!       Mask microlayers in the optical layer No.1
        NEWLR(1:INT(TAU(1)/DTAU)) = 1
!16
!       RUN RT CODE
!16
        CALL CPU_TIME(CPU_TIME1)
!16
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!16
!       No Aerosol in the top layer
        RATIO(1) = 0.0D0
!       Tau_a = 0.4; Tau_r = 0.1; Tau_abs = 0
        RATIO(2) = 0.4D0/(0.4D0 + 0.1D0)
!16
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!16
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!16
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!16
        CALL CPU_TIME(CPU_TIME2)
!16
        TIME = CPU_TIME2 - CPU_TIME1
!16
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T016/T016.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!16
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!16
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!16
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!16
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!16
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!16
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!16
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!16
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 016:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!16
        WRITE(*, 50) 'T016: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T016: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!16
    END IF ! TESTNUM(ITEST)
!16
!*******************************************************************************
!
!   017. GOAL: To test RT simulation in case of a thick layer. In particular,
!              approximation of scattering orders using geometric progression
!              is tested.  
!
!*******************************************************************************
!17
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!17
!       ACCURACY
!17
!       Number of Gauss nodes per hemisphere
        NG1 = 32
!       Total number of the Fourier moments
        NM = 1
!       Number of microlayers, defines dTau
        NL = 510
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!17
!       GEOMETRY
!17
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 1.0D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 31   ! Number of reflected angles
        NTU = NRU    ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS( (/(1.0D0 - 0.02D0*(IK-1), IK = 1, NRU)/) )*R2D ! Up
        UDA(NRU+1:NUD) = 180.0D0 - UDA(NRU:1:-1)                           ! Dn
!       View azimuth angles, degrees
        NAZ = 1 ! Number of azimuths
        AZA(1:NAZ) = 0.0D0
!17
!       ATMOSPHERE
!17
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0279D0
!       Number of optical layers
        NLR = 2
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = (/0.1D0, 10.1D0/)
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 117  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T017/Xk0117_0804.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!17
!       SURFACE
!17
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!17
!       COMPUTE AUXILIARY PARAMETERS
!17
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!17
!       Mask microlayers in the optical layer No.2
        NEWLR(1:NL) = 2
!       Mask microlayers in the optical layer No.1
        NEWLR(1:INT(TAU(1)/DTAU)) = 1
!17
!       RUN RT CODE
!17
        CALL CPU_TIME(CPU_TIME1)
!17
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!17
!       No Aerosol in the top layer
        RATIO(1) = 0.0D0
!       Tau_a = 10.0; Tau_r = 0.1; Tau_abs = 0
        RATIO(2) = 10.0D0/(10.0D0 + 0.1D0)
!17
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!17
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!17
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!17
        CALL CPU_TIME(CPU_TIME2)
!17
        TIME = CPU_TIME2 - CPU_TIME1
!17
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T017/T017.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!17
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!17
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!17
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!17
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!17
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!17
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!17
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!17
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 017:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!17
        WRITE(*, 50) 'T017: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T017: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!17
    END IF ! TESTNUM(ITEST)
!17
!*******************************************************************************
!
!   018. GOAL: To test RT simulation in case of a thick layer and SZA > 0.
!
!*******************************************************************************
!18
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!18
!       ACCURACY
!18
!       Number of Gauss nodes per hemisphere
        NG1 = 32
!       Total number of the Fourier moments
        NM = 24
!       Number of microlayers, defines dTau
        NL = 510
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!18
!       GEOMETRY
!18
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 45.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 31   ! Number of reflected angles
        NTU = NRU    ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS( (/(1.0D0 - 0.02D0*(IK-1), IK = 1, NRU)/) )*R2D ! Up
        UDA(NRU+1:NUD) = 180.0D0 - UDA(NRU:1:-1)                            ! Dn
!       View azimuth angles, degrees
        NAZ = 13 ! Number of azimuths
        AZA(1:NAZ) = (/(15.0D0*(IK-1), IK = 1, NAZ)/)
!18
!       ATMOSPHERE
!18
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0279D0
!       Number of optical layers
        NLR = 2
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = (/0.1D0, 10.1D0/)
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 117  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T018/Xk0117_0804.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!18
!       SURFACE
!18
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!18
!       COMPUTE AUXILIARY PARAMETERS
!18
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!18
!       Mask microlayers in the optical layer No.2
        NEWLR(1:NL) = 2
!       Mask microlayers in the optical layer No.1
        NEWLR(1:INT(TAU(1)/DTAU)) = 1
!18
!       RUN RT CODE
!18
        CALL CPU_TIME(CPU_TIME1)
!18
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!18
!       No Aerosol in the top layer
        RATIO(1) = 0.0D0
!       Tau_a = 10.0; Tau_r = 0.1; Tau_abs = 0
        RATIO(2) = 10.0D0/(10.0D0 + 0.1D0)
!18
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!18
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!18
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!18
        CALL CPU_TIME(CPU_TIME2)
!18
        TIME = CPU_TIME2 - CPU_TIME1
!18
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T018/T018.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!18
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!18
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!18
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!18
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!18
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!18
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!18
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!18
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 018:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!18
        WRITE(*, 50) 'T018: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T018: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!18
    END IF ! TESTNUM(ITEST)
!18
!*******************************************************************************
!
!   019. GOAL: To test pure Rayleigh scattering.
!
!*******************************************************************************
!19
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!19
!       ACCURACY
!19
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 1
!       Number of microlayers, defines dTau
        NL = 25
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!19
!       GEOMETRY
!19
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 0.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 17 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 100.0D0 + UDA(1:NRU)        ! Down
!       View azimuth angles, degrees
        NAZ = 1 ! Number of azimuths
        AZA(1:NAZ) = 0.0D0
!19
!       ATMOSPHERE
!19
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.5D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!19
!       SURFACE
!19
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!19
!       COMPUTE AUXILIARY PARAMETERS
!19
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!19
!       RUN RT CODE
!19
        CALL CPU_TIME(CPU_TIME1)
!19
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!19
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!19
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!19
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!19
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!19
        CALL CPU_TIME(CPU_TIME2)
!19
        TIME = CPU_TIME2 - CPU_TIME1
!19
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T019/T019.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!19
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!19
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!19
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!19
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!19
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!19
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!19
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!19
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 019:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!19
        WRITE(*, 50) 'T019: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T019: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!19
    END IF ! TESTNUM(ITEST)
!19
!*******************************************************************************
!
!   020. GOAL: To test Rayleigh scattering and depolarization factor.
!
!*******************************************************************************
!20
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!20
!       ACCURACY
!20
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 25
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!20
!       GEOMETRY
!20
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 30.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 17 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 100.0D0 + UDA(1:NRU)        ! Down
!       View azimuth angles, degrees
        NAZ = 73 ! Number of azimuths
        AZA(1:NAZ) = (/(5.0D0*(IK-1), IK = 1, NAZ)/)
!20
!       ATMOSPHERE
!20
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.5D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!20
!       SURFACE
!20
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!20
!       COMPUTE AUXILIARY PARAMETERS
!20
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!20
!       RUN RT CODE
!20
        CALL CPU_TIME(CPU_TIME1)
!20
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!20
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!20
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!20
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!20
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!20
        CALL CPU_TIME(CPU_TIME2)
!20
        TIME = CPU_TIME2 - CPU_TIME1
!20
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T020/T020.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!20
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!20
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!20
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!20
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!20
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!20
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!20
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!20
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 020:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!20
        WRITE(*, 50) 'T020: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T020: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!20
    END IF ! TESTNUM(ITEST)
!20
!*******************************************************************************
!
!   021. GOAL: To test solar azimuth, saz = 65 (0.0 is usually assumed in RT).
!
!*******************************************************************************
!21
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!21
!       ACCURACY
!21
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 25
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!21
!       GEOMETRY
!21
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 30.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 17 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 100.0D0 + UDA(1:NRU)        ! Down
!       View azimuth angles, degrees
        NAZ = 73 ! Number of azimuths
        AZA(1:13) = (/( 295.0D0 + 5.0D0*(IK-1), IK = 1, 13 )/)
        AZA(14:NAZ) = (/( 5.0D0*(IK-14), IK = 14, NAZ )/)    
!21
!       ATMOSPHERE
!21
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.1D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.5D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!21
!       SURFACE
!21
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!21
!       COMPUTE AUXILIARY PARAMETERS
!21
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!21
!       RUN RT CODE
!21
        CALL CPU_TIME(CPU_TIME1)
!21
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!21
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!21
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!21
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!21
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!21
        CALL CPU_TIME(CPU_TIME2)
!21
        TIME = CPU_TIME2 - CPU_TIME1
!21
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T021/T021.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!21
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!21
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!21
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!21
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!21
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!21
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!21
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!21
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 021:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!21
        WRITE(*, 50) 'T021: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T021: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!21
    END IF ! TESTNUM(ITEST)
!21
!*******************************************************************************
!
!   022. GOAL: To test Lambertian surface and Rayleigh with depolarization.
!
!*******************************************************************************
!22
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!22
!       ACCURACY
!22
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 10
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!22
!       GEOMETRY
!22
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 50.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 17 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 100.0D0 + UDA(1:NRU)        ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/(5.0D0*(IK-1), IK = 1, NAZ)/)
!22
!       ATMOSPHERE
!22
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.1D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!22
!       SURFACE
!22
!       Surface model: Lambertian
        ISRF = 1
!       Set of parameters
        PSRF(1) = 0.3D0  ! Surface albedo
        PSRF(2:) = 0.0D0 ! Other surface parameters
!       Number of nodes for the Fourier expansion
        NGA = 2
!22
!       COMPUTE AUXILIARY PARAMETERS
!22
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!22
!       RUN RT CODE
!22
        CALL CPU_TIME(CPU_TIME1)
!22
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!22
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!22
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!22
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!22
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!22
        CALL CPU_TIME(CPU_TIME2)
!22
        TIME = CPU_TIME2 - CPU_TIME1
!22
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T022/T022.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!22
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!22
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!22
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!22
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!22
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!22
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!22
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!22
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 022:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!22
        WRITE(*, 50) 'T022: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T022: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!22
    END IF ! TESTNUM(ITEST)
!22
!*******************************************************************************
!
!   023. GOAL: To test ocean surface for upward and downward view directions.
!
!*******************************************************************************
!23
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!23
!       ACCURACY
!23
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 10
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!23
!       GEOMETRY
!23
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 45.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 17 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 100.0D0 + UDA(1:NRU)        ! Down
!       View azimuth angles, degrees
        NAZ = 73 ! Number of azimuths
        AZA(1:NAZ) = (/(5.0D0*(IK-1), IK = 1, NAZ)/)
!23
!       ATMOSPHERE
!23
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.1D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!23
!       SURFACE
!23
!       Surface model: Nakajima-Tanaka ocen surface (no wind direction)
        ISRF = 6
!       Set of parameters
        PSRF(1) = 1.33D0  ! Refractive image of water: real part
        PSRF(2) = 0.00D0  ! Refractive image of water: imaginary part
        PSRF(3) = 2.00D0  ! Wind speed, m/s, at 10m above the surface
        PSRF(4:) = 0.0D0  ! Other surface parameters are not used
!       Number of nodes for the Fourier expansion
        NGA = 180
!23
!       COMPUTE AUXILIARY PARAMETERS
!23
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!23
!       RUN RT CODE
!23
        CALL CPU_TIME(CPU_TIME1)
!23
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!23
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!23
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)                           
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!23
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!23
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!23
        CALL CPU_TIME(CPU_TIME2)
!23
        TIME = CPU_TIME2 - CPU_TIME1
!23
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T023/T023.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!23
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!23
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!23
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!23
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!23
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!23
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!23
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!23
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 023:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!23
        WRITE(*, 50) 'T023: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T023: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!23
    END IF ! TESTNUM(ITEST)
!23
!*******************************************************************************
!
!   024. GOAL: To test light scattering by fine aerosol particles.
!
!*******************************************************************************
!24
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!24
!       ACCURACY
!24
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 6
!       Number of microlayers, defines dTau
        NL = 20
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!24
!       GEOMETRY
!24
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 40.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 17 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 100.0D0 + UDA(1:NRU)        ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/(5.0D0*(IK-1), IK = 1, NAZ)/)
!24
!       ATMOSPHERE
!24
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.2D0
        SSA(1:NLR) = 0.975683D0
!       Phase function
        NK = 10  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T024/Xk0010_0550.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!24
!       SURFACE
!24
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!24
!       COMPUTE AUXILIARY PARAMETERS
!24
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!24
!       RUN RT CODE
!24
        CALL CPU_TIME(CPU_TIME1)
!24
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!24
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!24
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!24
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!24
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!24
        CALL CPU_TIME(CPU_TIME2)
!24
        TIME = CPU_TIME2 - CPU_TIME1
!24
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T024/T024.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!24
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!24
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!24
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!24
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!24
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!24
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!24
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!24
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 024:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!24
        WRITE(*, 50) 'T024: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T024: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!24
    END IF ! TESTNUM(ITEST)
!24
!*******************************************************************************
!
!   025. GOAL: To test light scattering by dust (spheroids) particles.
!
!*******************************************************************************
!25
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!25
!       ACCURACY
!25
!       Number of Gauss nodes per hemisphere
        NG1 = 96
!       Total number of the Fourier moments
        NM = 64
!       Number of microlayers, defines dTau
        NL = 20
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!25
!       GEOMETRY
!25
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 40.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 17 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 100.0D0 + UDA(1:NRU)        ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/(5.0D0*(IK-1), IK = 1, NAZ)/)
!25
!       ATMOSPHERE
!25
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.2D0
        SSA(1:NLR) = 0.787581D0
!       Phase function
        NK = 1000 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T025/Xk1000_0838.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!25
!       SURFACE
!25
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!25
!       COMPUTE AUXILIARY PARAMETERS
!25
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!25
!       RUN RT CODE
!25
        CALL CPU_TIME(CPU_TIME1)
!25
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!25
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!25
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!25
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!25
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!25
        CALL CPU_TIME(CPU_TIME2)
!25
        TIME = CPU_TIME2 - CPU_TIME1
!25
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T025/T025.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!25
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!25
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!25
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!25
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!25
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!25
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!25
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!25
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 025:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!25
        WRITE(*, 50) 'T025: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T025: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!25
    END IF ! TESTNUM(ITEST)
!25
!*******************************************************************************
!
!   026. GOAL: To test cloud (thick) scenario & almucantar (AP) observation.
!
!*******************************************************************************
!26
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!26
!       ACCURACY
!26
!       Number of Gauss nodes per hemisphere
        NG1 = 128
!       Total number of the Fourier moments
        NM = 128
!       Number of microlayers, defines dTau
        NL = 250
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!26
!       GEOMETRY
!26
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 50.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 1 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = SZA               ! Up
        UDA(NRU+1:NUD) = 180.0D0 - SZA ! Down
!       View azimuth angles, degrees
        NAZ = 181 ! Number of azimuths
        AZA(1:NAZ) = (/(1.0D0*(IK-1), IK = 1, NAZ)/)
!26
!       ATMOSPHERE
!26
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 5.0D0
        SSA(1:NLR) = 0.999979D0
!       Phase function
        NK = 510 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T026/Xk0510_0859.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!26
!       SURFACE
!26
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!26
!       COMPUTE AUXILIARY PARAMETERS
!26
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!26
!       RUN RT CODE
!26
        CALL CPU_TIME(CPU_TIME1)
!26
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!26
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!26
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!26
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!26
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!26
        CALL CPU_TIME(CPU_TIME2)
!26
        TIME = CPU_TIME2 - CPU_TIME1
!26
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T026/T026.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!26
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!26
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!26
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!26
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!26
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!26
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!26
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!26
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 026:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!26
        WRITE(*, 50) 'T026: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T026: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!26
    END IF ! TESTNUM(ITEST)
!26
!*******************************************************************************
!
!   027. GOAL: To test cloud (thick) scenario & principal plane (PP) observation
!
!*******************************************************************************
!27
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!27
!       ACCURACY
!27
!       Number of Gauss nodes per hemisphere
        NG1 = 128
!       Total number of the Fourier moments
        NM = 128
!       Number of microlayers, defines dTau
        NL = 250
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!27
!       GEOMETRY
!27
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 50.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 81 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(1.0D0*(IK-1), IK = 1, NRU)/)               ! Up
        UDA(NRU+1:NUD) = (/(100.0D0 + 1.0D0*(IK-1), IK = 1, NTU)/) ! Down
!       View azimuth angles, degrees
        NAZ = 2 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 180.0D0/)
!27
!       ATMOSPHERE
!27
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 5.0D0
        SSA(1:NLR) = 0.999979D0
!       Phase function
        NK = 510 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T027/Xk0510_0859.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!27
!       SURFACE
!27
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!27
!       COMPUTE AUXILIARY PARAMETERS
!27
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!27
!       RUN RT CODE
!27
        CALL CPU_TIME(CPU_TIME1)
!27
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!27
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!27
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!27
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!27
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!27
        CALL CPU_TIME(CPU_TIME2)
!27
        TIME = CPU_TIME2 - CPU_TIME1
!27
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T027/T027.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!27
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!27
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!27
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!27
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!27
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!27
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!27
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!27
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 027:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!27
        WRITE(*, 50) 'T027: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T027: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!27
    END IF ! TESTNUM(ITEST)
!27
!*******************************************************************************
!
!   028. GOAL: To test scattering by coarse fraction of aerosol (Mie).
!
!*******************************************************************************
!28
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!28
!       ACCURACY
!28
!       Number of Gauss nodes per hemisphere
        NG1 = 128
!       Total number of the Fourier moments
        NM = 500 ! NM > 2*NG1=NG2. Why it works ??????
!       Number of microlayers, defines dTau
        NL = 20
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!28
!       GEOMETRY
!28
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 90 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(0.0D0 + (IK - 1.0D0), IK = 1, NRU)/) ! Reflection
        UDA(NRU+1:NUD) = 91.0D0 + UDA(1:NRU)                 ! Transmittance
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!28
!       ATMOSPHERE
!28
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.3262D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 919 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T028/Xk0919_0793.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!28
!       SURFACE
!28
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!28
!       COMPUTE AUXILIARY PARAMETERS
!28
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!28
!       RUN RT CODE
!28
        CALL CPU_TIME(CPU_TIME1)
!28
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!28
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!28
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!28
        CALL CPU_TIME(CPU_TIME2)
!28
        TIME = CPU_TIME2 - CPU_TIME1
!28
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T028/T028.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!28
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!28
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!28
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!28
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!28
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!28
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!28
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!28
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 0028:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!28
        WRITE(*, 50) 'T028: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T028: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!28
    END IF ! TESTNUM(ITEST)
!28
!*******************************************************************************
!
!   029. GOAL: To test scattering by thick cloud of water droplets (Mie).
!
!*******************************************************************************
!29
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!29
!       ACCURACY
!29
!       Number of Gauss nodes per hemisphere
        NG1 = 256
!       Total number of the Fourier moments
        NM = 256
!       Number of microlayers, defines dTau
        NL = 500
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!29
!       GEOMETRY
!29
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 90 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(0.0D0 + (IK - 1.0D0), IK = 1, NRU)/) ! Reflection
        UDA(NRU+1:NUD) = 91.0D0 + UDA(1:NRU)                 ! Transmittance
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!29
!       ATMOSPHERE
!29
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 5.0D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 1670 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T029/Xk1670_0861.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!29
!       SURFACE
!29
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!29
!       COMPUTE AUXILIARY PARAMETERS
!29
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!29
!       RUN RT CODE
!29
        CALL CPU_TIME(CPU_TIME1)
!29
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!29
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!29
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!29
        CALL CPU_TIME(CPU_TIME2)
!29
        TIME = CPU_TIME2 - CPU_TIME1
!29
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T029/T029.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!29
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!29
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!29
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!29
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!29
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!29
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!29
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!29
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 029:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!29
        WRITE(*, 50) 'T029: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T029: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!29
    END IF ! TESTNUM(ITEST)
!29
!*******************************************************************************
!
!   030. GOAL: To test light scattering by fine aerosol particles.
!
!*******************************************************************************
!30
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!30
!       ACCURACY
!30
!       Number of Gauss nodes per hemisphere
        NG1 = 6
!       Total number of the Fourier moments
        NM = 1
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!30
!       GEOMETRY
!30
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 1.0D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 10  ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS( (/(1.0D0 - 0.1D0*(IK-1), IK = 1, NRU)/) )*R2D ! Up
        UDA(NRU+1:NUD) = 180.0D0 - UDA(NRU:1:-1)                           ! Dn
!       View azimuth angles, degrees
        NAZ = 1 ! Number of azimuths
        AZA(1:NAZ) = 0.0D0
!30
!       ATMOSPHERE
!30
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.00D0
        SSA(1:NLR) = 0.99D0
!       Phase function
        NK = 12  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T030/Xk0012_0485.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!30
!       SURFACE
!30
!       Surface model
        ISRF = 1
!       Set of parameters
        PSRF(1) = 0.1D0  ! Surface albedo
        PSRF(2:) = 0.0D0 ! Other parameters
!       Number of nodes for the Fourier expansion
        NGA = 2
!30
!       COMPUTE AUXILIARY PARAMETERS
!30
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!30
!       RUN RT CODE
!30
        CALL CPU_TIME(CPU_TIME1)
!30
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!30
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!30
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!30
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!30
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!30
        CALL CPU_TIME(CPU_TIME2)
!30
        TIME = CPU_TIME2 - CPU_TIME1
!30
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T030/T030.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!30
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!30
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!30
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!30
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!30
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!30
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!30
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!30
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 030:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!30
        WRITE(*, 50) 'T030: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T030: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!30
    END IF ! TESTNUM(ITEST)
!30
!*******************************************************************************
!
!   031. GOAL: To test light scattering by fine aerosol particles.
!
!*******************************************************************************
!31
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!31
!       ACCURACY
!31
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 5
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!31
!       GEOMETRY
!31
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.6D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 10  ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS( (/(1.0D0 - 0.1D0*(IK-1), IK = 1, NRU)/) )*R2D ! Up
        UDA(NRU+1:NUD) = 180.0D0 - UDA(NRU:1:-1)                           ! Dn
!       View azimuth angles, degrees
        NAZ = 1 ! Number of azimuths
        AZA(1:NAZ) = 0.0D0
!31
!       ATMOSPHERE
!31
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.0D0
        SSA(1:NLR) = 0.973527D0
!       Phase function
        NK = 12  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T031/Xk0012_0701.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!31
!       SURFACE
!31
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!31
!       COMPUTE AUXILIARY PARAMETERS
!31
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!31
!       RUN RT CODE
!31
        CALL CPU_TIME(CPU_TIME1)
!31
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!31
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!31
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!31
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!31
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!31
        CALL CPU_TIME(CPU_TIME2)
!31
        TIME = CPU_TIME2 - CPU_TIME1
!31
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T031/T031.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!31
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!31
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!31
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!31
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!31
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!31
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!31
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!31
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 031:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!31
        WRITE(*, 50) 'T031: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T031: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!31
    END IF ! TESTNUM(ITEST)
!31
!*******************************************************************************
!
!   032. GOAL: To test low Sun, mu0=0.2 (SZA~=78.5).
!
!*******************************************************************************
!32
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!32
!       ACCURACY
!32
!       Number of Gauss nodes per hemisphere
        NG1 = 6
!       Total number of the Fourier moments
        NM = 6
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!32
!       GEOMETRY
!32
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.2D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 10  ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS( (/(1.0D0 - 0.1D0*(IK-1), IK = 1, NRU)/) )*R2D ! Up
        UDA(NRU+1:NUD) = 180.0D0 - UDA(NRU:1:-1)                           ! Dn
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!32
!       ATMOSPHERE
!32
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.00D0
        SSA(1:NLR) = 0.99D0
!       Phase function
        NK = 12  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T032/Xk0012_0485.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!32
!       SURFACE
!32
!       Surface model
        ISRF = 1
!       Set of parameters
        PSRF(1) = 0.1D0  ! Surface albedo
        PSRF(2:) = 0.0D0 ! Other parameters
!       Number of nodes for the Fourier expansion
        NGA = 2
!32
!       COMPUTE AUXILIARY PARAMETERS
!32
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!32
!       RUN RT CODE
!32
        CALL CPU_TIME(CPU_TIME1)
!32
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!32
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!32
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!32
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!32
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!32
        CALL CPU_TIME(CPU_TIME2)
!32
        TIME = CPU_TIME2 - CPU_TIME1
!32
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T032/T032.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!32
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!32
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!32
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!32
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!32
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!32
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!32
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!32
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 032:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!32
        WRITE(*, 50) 'T032: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T032: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!32
    END IF ! TESTNUM(ITEST)
!32
!*******************************************************************************
!
!   033. GOAL: To test thick atmopshere, Tau0=10.0.
!
!*******************************************************************************
!33
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!33
!       ACCURACY
!33
!       Number of Gauss nodes per hemisphere
        NG1 = 6
!       Total number of the Fourier moments
        NM = 1
!       Number of microlayers, defines dTau
        NL = 500
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!33
!       GEOMETRY
!33
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 1.0D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 10  ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS( (/(1.0D0 - 0.1D0*(IK-1), IK = 1, NRU)/) )*R2D ! Up
        UDA(NRU+1:NUD) = 180.0D0 - UDA(NRU:1:-1)                          ! Dn
!       View azimuth angles, degrees
        NAZ = 1 ! Number of azimuths
        AZA(1:NAZ) = 0.0D0
!33
!       ATMOSPHERE
!33
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 10.0D0
        SSA(1:NLR) = 0.99D0
!       Phase function
        NK = 12  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T033/Xk0012_0485.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!33
!       SURFACE
!33
!       Surface model
        ISRF = 1
!       Set of parameters
        PSRF(1) = 0.1D0  ! Surface albedo
        PSRF(2:) = 0.0D0 ! Other parameters
!       Number of nodes for the Fourier expansion
        NGA = 2
!33
!       COMPUTE AUXILIARY PARAMETERS
!33
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!33
!       RUN RT CODE
!33
        CALL CPU_TIME(CPU_TIME1)
!33
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!33
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!33
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!33
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!33
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!33
        CALL CPU_TIME(CPU_TIME2)
!33
        TIME = CPU_TIME2 - CPU_TIME1
!33
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T033/T033.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!33
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!33
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!33
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!33
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!33
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!33
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!33
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!33
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 033:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!33
        WRITE(*, 50) 'T033: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T033: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!33
    END IF ! TESTNUM(ITEST)
!33
!*******************************************************************************
!
!   034. GOAL: To test light scattering by mild aerosol particles.
!
!*******************************************************************************
!34
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!34
!       ACCURACY
!34
!       Number of Gauss nodes per hemisphere
        NG1 = 16
!       Total number of the Fourier moments
        NM = 1
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!34
!       GEOMETRY
!34
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 1.0D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 10  ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS( (/(1.0D0 - 0.1D0*(IK-1), IK = 1, NRU)/) )*R2D ! Up
        UDA(NRU+1:NUD) = 180.0D0 - UDA(NRU:1:-1)                          ! Dn
!       View azimuth angles, degrees
        NAZ = 1 ! Number of azimuths
        AZA(1:NAZ) = 0.0D0
!34
!       ATMOSPHERE
!34
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.00D0
        SSA(1:NLR) = 0.99D0
!       Phase function
        NK = 54  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T034/Xk0054_0677.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!34
!       SURFACE
!34
!       Surface model
        ISRF = 1
!       Set of parameters
        PSRF(1) = 0.1D0  ! Surface albedo
        PSRF(2:) = 0.0D0 ! Other parameters
!       Number of nodes for the Fourier expansion
        NGA = 2
!34
!       COMPUTE AUXILIARY PARAMETERS
!34
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!34
!       RUN RT CODE
!34
        CALL CPU_TIME(CPU_TIME1)
!34
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!34
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!34
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!34
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!34
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!34
        CALL CPU_TIME(CPU_TIME2)
!34
        TIME = CPU_TIME2 - CPU_TIME1
!34
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T034/T034.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!34
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!34
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!34
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!34
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!34
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!34
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!34
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!34
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 034:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!34
        WRITE(*, 50) 'T034: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T034: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!34
    END IF ! TESTNUM(ITEST)
!34
!*******************************************************************************
!
!   035. GOAL: To test low Sun, mu0=0.2 (SZA~=78.5)
!
!*******************************************************************************
!35
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!35
!       ACCURACY
!35
!       Number of Gauss nodes per hemisphere
        NG1 = 16
!       Total number of the Fourier moments
        NM = 24
!       Number of microlayers, defines dTau
        NL = 100
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!35
!       GEOMETRY
!35
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.2D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 10  ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS( (/(1.0D0 - 0.1D0*(IK-1), IK = 1, NRU)/) )*R2D ! Up
        UDA(NRU+1:NUD) = 180.0D0 - UDA(NRU:1:-1)                          ! Dn
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!35
!       ATMOSPHERE
!35
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.00D0
        SSA(1:NLR) = 0.99D0
!       Phase function
        NK = 54  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T035/Xk0054_0677.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!35
!       SURFACE
!35
!       Surface model
        ISRF = 1
!       Set of parameters
        PSRF(1) = 0.1D0  ! Surface albedo
        PSRF(2:) = 0.0D0 ! Other parameters
!       Number of nodes for the Fourier expansion
        NGA = 2
!35
!       COMPUTE AUXILIARY PARAMETERS
!35
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!35
!       RUN RT CODE
!35
        CALL CPU_TIME(CPU_TIME1)
!35
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!35
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!35
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!35
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!35
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!35
        CALL CPU_TIME(CPU_TIME2)
!35
        TIME = CPU_TIME2 - CPU_TIME1
!35
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T035/T035.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!35
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!35
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!35
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!35
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!35
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!35
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!35
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!35
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 035:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!35
        WRITE(*, 50) 'T035: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T035: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!35
    END IF ! TESTNUM(ITEST)
!35
!*******************************************************************************
!
!   036. GOAL: To test scattering by mild aerosol particles and thick atmopshere
!
!*******************************************************************************
!36
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!36
!       ACCURACY
!36
!       Number of Gauss nodes per hemisphere
        NG1 = 16
!       Total number of the Fourier moments
        NM = 1
!       Number of microlayers, defines dTau
        NL = 500
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!36
!       GEOMETRY
!36
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 1.0D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 10  ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS( (/(1.0D0 - 0.1D0*(IK-1), IK = 1, NRU)/) )*R2D ! Up
        UDA(NRU+1:NUD) = 180.0D0 - UDA(NRU:1:-1)                          ! Dn
!       View azimuth angles, degrees
        NAZ = 1 ! Number of azimuths
        AZA(1:NAZ) = 0.0D0
!36
!       ATMOSPHERE
!36
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 10.0D0
        SSA(1:NLR) = 0.99D0
!       Phase function
        NK = 54  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T036/Xk0054_0677.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!36
!       SURFACE
!36
!       Surface model
        ISRF = 1
!       Set of parameters
        PSRF(1) = 0.1D0  ! Surface albedo
        PSRF(2:) = 0.0D0 ! Other parameters
!       Number of nodes for the Fourier expansion
        NGA = 2
!36
!       COMPUTE AUXILIARY PARAMETERS
!36
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!36
!       RUN RT CODE
!36
        CALL CPU_TIME(CPU_TIME1)
!36
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!36
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!36
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!36
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!36
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!36
        CALL CPU_TIME(CPU_TIME2)
!36
        TIME = CPU_TIME2 - CPU_TIME1
!36
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T036/T036.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!36
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!36
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!36
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!36
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!36
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!36
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!36
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!36
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 036:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!36
        WRITE(*, 50) 'T036: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T036: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!36
    END IF ! TESTNUM(ITEST)
!36
!*******************************************************************************
!
!   037. GOAL: To test arbitrary number of optical layers, NLR=NL
!
!*******************************************************************************
!37
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!37
!       ACCURACY
!37
!       Number of Gauss nodes per hemisphere (NG1=20, see p.171 in the paper)
        NG1 = 20
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!37
!       GEOMETRY
!37
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.9D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees)
        NRU = 11 ! Number of reflected angles
        NTU = 0  ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NUD) = DACOS((/1.0D0, 0.9D0, 0.8D0, 0.7D0, 0.6D0, 0.5D0, & ! Up
                       0.4D0, 0.3D0, 0.2D0, 0.1D0, 0.05D0/))*R2D
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!37
!       ATMOSPHERE
!37
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = NL
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.0D0/NLR
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!37
!       SURFACE
!37
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!37
!       COMPUTE AUXILIARY PARAMETERS
!37
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers are different
        NEWLR = (/( IL, IL = 1, NL )/)
!37
!       RUN RT CODE
!37
        CALL CPU_TIME(CPU_TIME1)
!37
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!37
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!37
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!37
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!37
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!37
        CALL CPU_TIME(CPU_TIME2)
!37
        TIME = CPU_TIME2 - CPU_TIME1
!37
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T037/T037.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!37
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!37
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!37
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!37
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!37
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!37
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!37
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!37
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 037:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!37
        WRITE(*, 50) 'T037: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T037: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!37
    END IF ! TESTNUM(ITEST)
!37
!*******************************************************************************
!
!   038. GOAL: To test arbitrary number of optical layers, NLR
!
!*******************************************************************************
!38
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!38
!       ACCURACY
!38
!       Number of Gauss nodes per hemisphere (NG1=20, see p.171 in the paper)
        NG1 = 20
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!38
!       GEOMETRY
!38
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.9D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees)
        NRU = 11 ! Number of reflected angles
        NTU = 0  ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NUD) = DACOS((/1.0D0, 0.9D0, 0.8D0, 0.7D0, 0.6D0, 0.5D0, & ! Up
                       0.4D0, 0.3D0, 0.2D0, 0.1D0, 0.05D0/))*R2D
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!38
!       ATMOSPHERE
!38
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 4
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.0D0/NLR
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!38
!       SURFACE
!38
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!38
!       COMPUTE AUXILIARY PARAMETERS
!38
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!38
        IL1 = 1
        DO ILR = 1, NLR-1
            NLR0 = NINT(TAU(ILR)/DTAU)
            IL2 = IL1+NLR0-1
            DELT = TAU(ILR) - DTAU*NLR0
            TAU(ILR) = DTAU*NLR0
            TAU(ILR+1) = TAU(ILR+1) + DELT
            NEWLR(IL1:IL2) = ILR
            IL1 = IL2+1
        END DO ! ILR = 1, NLR-1
        NEWLR(IL1:NL) = NLR
!38
        TAU_ACC = 0.0D0 ! Accumulate Tau top to bottom
        DO ILR = 1, NLR
            SSA(ILR) = DEXP(-0.01D0*TAU_ACC) ! g = 0.01
            TAU_ACC = TAU_ACC + TAU(ILR)
            SSA(ILR) = (SSA(ILR) + DEXP(-0.01D0*TAU_ACC))/2.0D0
        END DO
!38
!       RUN RT CODE
!38
        CALL CPU_TIME(CPU_TIME1)
!38
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!38
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!38
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!38
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!38
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!38
        CALL CPU_TIME(CPU_TIME2)
!38
        TIME = CPU_TIME2 - CPU_TIME1
!38
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T038/T038.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!38
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!38
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!38
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!38
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!38
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!38
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!38
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!38
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 038:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!38
        WRITE(*, 50) 'T038: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T038: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!38
    END IF ! TESTNUM(ITEST)
!38
!*******************************************************************************
!
!   039. GOAL: To test arbitrary number of optical layers, NLR
!
!*******************************************************************************
!39
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!39
!       ACCURACY
!39
!       Number of Gauss nodes per hemisphere (NG1=20, see p.171 in the paper)
        NG1 = 20
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!39
!       GEOMETRY
!39
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.9D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees)
        NRU = 11 ! Number of reflected angles
        NTU = 0  ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NUD) = DACOS((/1.0D0, 0.9D0, 0.8D0, 0.7D0, 0.6D0, 0.5D0, & ! Up
                              0.4D0, 0.3D0, 0.2D0, 0.1D0, 0.05D0/))*R2D
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!39
!       ATMOSPHERE
!39
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 16
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.0D0/NLR
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!39
!       SURFACE
!39
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!39
!       COMPUTE AUXILIARY PARAMETERS
!39
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!39
        IL1 = 1
        DO ILR = 1, NLR-1
            NLR0 = NINT(TAU(ILR)/DTAU)
            IL2 = IL1+NLR0-1
            DELT = TAU(ILR) - DTAU*NLR0
            TAU(ILR) = DTAU*NLR0
            TAU(ILR+1) = TAU(ILR+1) + DELT
            NEWLR(IL1:IL2) = ILR
            IL1 = IL2+1
        END DO ! ILR = 1, NLR-1
        NEWLR(IL1:NL) = NLR
!39
        TAU_ACC = 0.0D0 ! Accumulate Tau top to bottom
        DO ILR = 1, NLR
            SSA(ILR) = DEXP(-0.1D0*TAU_ACC) ! g = 0.1
            TAU_ACC = TAU_ACC + TAU(ILR)
            SSA(ILR) = (SSA(ILR) + DEXP(-0.1D0*TAU_ACC))/2.0D0
        END DO
!39
!       RUN RT CODE
!39
        CALL CPU_TIME(CPU_TIME1)
!39
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!39
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!39
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!39
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!39
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!39
        CALL CPU_TIME(CPU_TIME2)
!39
        TIME = CPU_TIME2 - CPU_TIME1
!39
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T039/T039.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!39
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!39
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!39
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!39
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!39
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!39
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!39
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!39
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 039:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!39
        WRITE(*, 50) 'T039: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T039: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!39
    END IF ! TESTNUM(ITEST)
!39
!*******************************************************************************
!
!   040. GOAL: To test arbitrary number of optical layers, NLR
!
!*******************************************************************************
!40
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!40
!       ACCURACY
!40
!       Number of Gauss nodes per hemisphere (NG1=20, see p.171 in the paper)
        NG1 = 20
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!40
!       GEOMETRY
!40
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.9D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees)
        NRU = 11 ! Number of reflected angles
        NTU = 0  ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NUD) = DACOS((/1.0D0, 0.9D0, 0.8D0, 0.7D0, 0.6D0, 0.5D0, & ! Up
                              0.4D0, 0.3D0, 0.2D0, 0.1D0, 0.05D0/))*R2D
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/0.0D0, 90.0D0, 180.0D0/)
!40
!       ATMOSPHERE
!40
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = NL
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.0D0/NLR
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!40
!       SURFACE
!40
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!40
!       COMPUTE AUXILIARY PARAMETERS
!40
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!40
        IL1 = 1
        DO ILR = 1, NLR-1
            NLR0 = NINT(TAU(ILR)/DTAU)
            IL2 = IL1+NLR0-1
            DELT = TAU(ILR) - DTAU*NLR0
            TAU(ILR) = DTAU*NLR0
            TAU(ILR+1) = TAU(ILR+1) + DELT
            NEWLR(IL1:IL2) = ILR
            IL1 = IL2+1
        END DO ! ILR = 1, NLR-1
        NEWLR(IL1:NL) = NLR
!40
        TAU_ACC = 0.0D0 ! Accumulate Tau top to bottom
        DO ILR = 1, NLR
            SSA(ILR) = DEXP(-TAU_ACC)    ! g = 1.0
            TAU_ACC = TAU_ACC + TAU(ILR)
            SSA(ILR) = (SSA(ILR) + DEXP(-TAU_ACC))/2.0D0
        END DO
!40
!       RUN RT CODE
!40
        CALL CPU_TIME(CPU_TIME1)
!40
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!40
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!40
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!40
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!40
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!40
        CALL CPU_TIME(CPU_TIME2)
!40
        TIME = CPU_TIME2 - CPU_TIME1
!40
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T040/T040.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!40
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!40
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!40
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!40
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!40
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!40
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!40
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!40
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 040:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!40
        WRITE(*, 50) 'T040: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T040: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!40
    END IF ! TESTNUM(ITEST)
!40
!*******************************************************************************
!
!   041. GOAL: To test light scattering by fine prolate spheroids.
!
!*******************************************************************************
!41
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!41
!       ACCURACY
!41
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 11
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!41
!       GEOMETRY
!41
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.6D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 10  ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS( (/(1.0D0 - 0.1D0*(IK-1), IK = 1, NRU)/) )*R2D ! Up
        UDA(NRU+1:NUD) = 180.0D0 - UDA(NRU:1:-1)                          ! Dn
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/ 0.0D0, 90.0D0, 180.0D0 /)
!41
!       ATMOSPHERE
!41
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.0D0
        SSA(1:NLR) = 0.954174D0
!       Phase function
        NK = 26  ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T041/Xk0026_0806.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!41
!       SURFACE
!41
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!41
!       COMPUTE AUXILIARY PARAMETERS
!41
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!41
!       RUN RT CODE
!41
        CALL CPU_TIME(CPU_TIME1)
!41
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!41
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!41
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!41
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!41
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!41
        CALL CPU_TIME(CPU_TIME2)
!41
        TIME = CPU_TIME2 - CPU_TIME1
!41
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T041/T041.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!41
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!41
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!41
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!41
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!41
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!41
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!41
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!41
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 041:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!41
        WRITE(*, 50) 'T041: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T041: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!41
    END IF ! TESTNUM(ITEST)
!41
!*******************************************************************************
!
!   042. GOAL: To test light scattering by fine oblate spheroids
!
!*******************************************************************************
!42
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!42
!       ACCURACY
!42
!       Number of Gauss nodes per hemisphere
        NG1 = 6
!       Total number of the Fourier moments
        NM = 5
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!42
!       GEOMETRY
!42
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.6D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 10  ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS( (/(1.0D0 - 0.1D0*(IK-1), IK = 1, NRU)/) )*R2D ! Up
        UDA(NRU+1:NUD) = 180.0D0 - UDA(NRU:1:-1)                          ! Dn
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/ 0.0D0, 90.0D0, 180.0D0 /)
!42
!       ATMOSPHERE
!42
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.0D0
        SSA(1:NLR) = 0.973527D0
!       Phase function
        NK = 12 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T042/Xk0012_0701.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!42
!       SURFACE
!42
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!42
!       COMPUTE AUXILIARY PARAMETERS
!42
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!42
!       RUN RT CODE
!42
        CALL CPU_TIME(CPU_TIME1)
!42
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!42
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!42
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!42
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!42
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!42
        CALL CPU_TIME(CPU_TIME2)
!42
        TIME = CPU_TIME2 - CPU_TIME1
!42
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T042/T042.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!42
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!42
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!42
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!42
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!42
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!42
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!42
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!42
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 042:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!42
        WRITE(*, 50) 'T042: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T042: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!42
    END IF ! TESTNUM(ITEST)
!42
!*******************************************************************************
!
!   043. GOAL: To test light scattering by fine prolate spheroids.
!
!*******************************************************************************
!43
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!43
!       ACCURACY
!43
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 10
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!43
!       GEOMETRY
!43
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        MU0 = 0.6D0
        SZA = DACOS(MU0)*R2D
!       User-defined zenith angles, UDA (degrees), and/or cosine
        NRU = 10  ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = DACOS( (/(1.0D0 - 0.1D0*(IK-1), IK = 1, NRU)/) )*R2D ! Up
        UDA(NRU+1:NUD) = 180.0D0 - UDA(NRU:1:-1)                          ! Dn
!       View azimuth angles, degrees
        NAZ = 3 ! Number of azimuths
        AZA(1:NAZ) = (/ 0.0D0, 90.0D0, 180.0D0 /)
!43
!       ATMOSPHERE
!43
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 1.0D0
        SSA(1:NLR) = 0.951151D0
!       Phase function
        NK = 36 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T043/Xk0036_0695.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!43
!       SURFACE
!43
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!43
!       COMPUTE AUXILIARY PARAMETERS
!43
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!43
!       RUN RT CODE
!43
        CALL CPU_TIME(CPU_TIME1)
!43
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!43
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!43
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!43
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!43
!       The benchmark results are normalized by PI: divide output by 2.0
        II = 0.5D0*II; QQ = 0.5D0*QQ; UU = 0.5D0*UU
!43
        CALL CPU_TIME(CPU_TIME2)
!43
        TIME = CPU_TIME2 - CPU_TIME1
!43
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T043/T043.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!43
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!43
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!43
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!43
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!43
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!43
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!43
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!43
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 043:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!43
        WRITE(*, 50) 'T043: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T043: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!43
    END IF ! TESTNUM(ITEST)
!43
!*******************************************************************************
!
!   044. GOAL: To test the case of many optical layers, NLR=NL.
!
!*******************************************************************************
!44
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!44
!       ACCURACY
!44
!       Number of Gauss nodes per hemisphere
        NG1 = 6
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 20
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!44
!       GEOMETRY
!44
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 18 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 95.0D0 + UDA(1:NRU)         ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/( 5.0D0*(IK-1), IK = 1, NAZ )/)
!44
!       ATMOSPHERE
!44
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Number of optical layers
        NLR = NL
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.21880749390D0/NLR
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!44
!       SURFACE
!44
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!44
!       COMPUTE AUXILIARY PARAMETERS
!44
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!44
!       RUN RT CODE
!44
        CALL CPU_TIME(CPU_TIME1)
!44
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!44
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!44
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!44
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!44
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!44
        CALL CPU_TIME(CPU_TIME2)
!44
        TIME = CPU_TIME2 - CPU_TIME1
!44
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T044/T044.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!44
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!44
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!44
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!44
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!44
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!44
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!44
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!44
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 044:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!44
        WRITE(*, 50) 'T044: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T044: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!44
    END IF ! TESTNUM(ITEST)
!44
!*******************************************************************************
!
!   045. GOAL: To test the cloud (thick) over ocean scenario
!
!*******************************************************************************
!45
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!45
!       ACCURACY
!45
!       Number of Gauss nodes per hemisphere
        NG1 = 128
!       Total number of the Fourier moments
        NM = 128
!       Number of microlayers, defines dTau
        NL = 250
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!45
!       GEOMETRY
!45
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 18 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 95.0D0 + UDA(1:NRU)        ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/(5.0D0*(IK-1), IK = 1, NAZ)/)
!45
!       ATMOSPHERE
!45
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.0D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 5.0D0
        SSA(1:NLR) = 0.999979D0
!       Phase function
        NK = 510 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T045/Xk0510_0859.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!45
!       SURFACE
!45
!       Surface model: Nakajima-Tanaka ocen surface (no wind direction)
        ISRF = 6
!       Set of parameters
        PSRF(1) = 1.33D0  ! Refractive image of water: real part
        PSRF(2) = 0.00D0  ! Refractive image of water: imaginary part
        PSRF(3) = 2.00D0  ! Wind speed, m/s, at 10m above the surface
        PSRF(4:) = 0.0D0  ! Other surface parameters are not used
!       Number of nodes for the Fourier expansion
        NGA = 180
!45
!       COMPUTE AUXILIARY PARAMETERS
!45
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!45
!       RUN RT CODE
!45
        CALL CPU_TIME(CPU_TIME1)
!45
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!45
!       No Rayleigh
        RATIO(1:NLR) = 1.0D0
        R11(1:NUD, 1:NAZ) = 0.0D0
        R12(1:NUD, 1:NAZ) = 0.0D0
!45
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!45
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!45
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!45
        CALL CPU_TIME(CPU_TIME2)
!45
        TIME = CPU_TIME2 - CPU_TIME1
!45
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T045/T045.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!45
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!45
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!45
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!45
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!45
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!45
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!45
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!45
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 045:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!45
        WRITE(*, 50) 'T045: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T045: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!45
    END IF ! TESTNUM(ITEST)
!45
!*******************************************************************************
!
!   046. GOAL: To test the cloud (thick) mixed with Rayleigh over ocean scenario
!
!*******************************************************************************
!46
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!46
!       ACCURACY
!46
!       Number of Gauss nodes per hemisphere
        NG1 = 128
!       Total number of the Fourier moments
        NM = 128
!       Number of microlayers, defines dTau
        NL = 252
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!46
!       GEOMETRY
!46
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 18 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 95.0D0 + UDA(1:NRU)        ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/(5.0D0*(IK-1), IK = 1, NAZ)/)
!46
!       ATMOSPHERE
!46
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 5.0D0 + 0.02098185D0 ! Tau_Cloud + Tau_Rayleigh
        SSA(1:NLR) = (0.999979D0*5.0D0 + 0.02098185D0)/(5.0D0 + 0.02098185D0)
!       Phase function
        NK = 510 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T046/Xk0510_0859.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!46
!       SURFACE
!46
!       Surface model: Nakajima-Tanaka ocen surface (no wind direction)
        ISRF = 6
!       Set of parameters
        PSRF(1) = 1.33D0  ! Refractive image of water: real part
        PSRF(2) = 0.00D0  ! Refractive image of water: imaginary part
        PSRF(3) = 2.00D0  ! Wind speed, m/s, at 10m above the surface
        PSRF(4:) = 0.0D0  ! Other surface parameters are not used
!       Number of nodes for the Fourier expansion
        NGA = 180
!46
!       COMPUTE AUXILIARY PARAMETERS
!46
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!46
!       RUN RT CODE
!46
        CALL CPU_TIME(CPU_TIME1)
!46
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!46
!       Mix Cloud & Rayleigh
        RATIO(1:NLR) = 0.999979D0*5.0D0/(0.999979D0*5.0D0 + 0.02098185D0)
!46
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!46
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!46
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!46
        CALL CPU_TIME(CPU_TIME2)
!46
        TIME = CPU_TIME2 - CPU_TIME1
!46
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T046/T046.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!46
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!46
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!46
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!46
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!46
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!46
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!46
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!46
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 046:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!46
        WRITE(*, 50) 'T046: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T046: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!46
    END IF ! TESTNUM(ITEST)
!46
!*******************************************************************************
!
!   047. GOAL: To test the cloud (thick) mixed with Rayleigh over ocean scenario
!              and estimate influence of the Rayleigh profile.
!
!*******************************************************************************
!47
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!47
!       ACCURACY
!47
!       Number of Gauss nodes per hemisphere
        NG1 = 256
!       Total number of the Fourier moments
        NM = 128
!       Number of microlayers, defines dTau
        NL = 252
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 2
!47
!       GEOMETRY
!47
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 18 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 95.0D0 + UDA(1:NRU)        ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/(5.0D0*(IK-1), IK = 1, NAZ)/)
!47
!       ATMOSPHERE
!47
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 5.0D0 + 0.02098185D0 ! Tau_Cloud + Tau_Rayleigh
        SSA(1:NLR) = (0.999979D0*5.0D0 + 0.02098185D0)/(5.0D0 + 0.02098185D0)
!       Phase function
        NK = 510 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T047/Xk0510_0859.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!47
!       SURFACE
!47
!       Surface model: Nakajima-Tanaka ocen surface (no wind direction)
        ISRF = 6
!       Set of parameters
        PSRF(1) = 1.33D0  ! Refractive image of water: real part
        PSRF(2) = 0.00D0  ! Refractive image of water: imaginary part
        PSRF(3) = 2.00D0  ! Wind speed, m/s, at 10m above the surface
        PSRF(4:) = 0.0D0  ! Other surface parameters are not used
!       Number of nodes for the Fourier expansion
        NGA = 180
!47
!       COMPUTE AUXILIARY PARAMETERS
!47
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!47
!       RUN RT CODE
!47
        CALL CPU_TIME(CPU_TIME1)
!47
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
!47
!       Mix Cloud & Rayleigh
        RATIO(1:NLR) = 0.999979D0*5.0D0/(0.999979D0*5.0D0 + 0.02098185D0)
!47
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi
!47
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!47
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!47
        CALL CPU_TIME(CPU_TIME2)
!47
        TIME = CPU_TIME2 - CPU_TIME1
!47
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T047/T047.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!47
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!47
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!47
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!47
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!47
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!47
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
!47
!               ***************************************************
!               *** VZA = 80, 85, 95, 100 (4 total) are ignored ***
!               *** Corresponding IKs are 17, 18, 19, and 20    ***
!               ***************************************************
                IF (PRODUCT( (/17, 18, 19, 20/)-IK ) == 0) THEN
                    OUT(IY, 6) = 0.0D0 ! Skip dI%
                    OUT(IY, 7) = 0.0D0 ! Skip dP
                END IF
!47
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!47
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/(NY-4*NAZ) ! 4 ignored UDA ...
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/(NY-4*NAZ) ! ... per azimuth
!47
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 047:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!47
        WRITE(*, 50) 'T047: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T047: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!47
    END IF ! TESTNUM(ITEST)
!47
!*******************************************************************************
!
!   048. GOAL: To test the case with absrobtion profile. Single scattering only
!
!*******************************************************************************
!48
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!48
!       ACCURACY
!48
!       Number of Gauss nodes per hemisphere
        NG1 = 6
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 100
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = -999.0D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 1
!48
!       GEOMETRY
!48
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 18 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 95.0D0 + UDA(1:NRU)         ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/( 5.0D0*(IK-1), IK = 1, NAZ )/)
!48
!       ATMOSPHERE
!48
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Number of optical layers
        NLR = NL
!       Number of heights in km
        NKM = 31
!       Heights
        ZKM(1:NKM) = (/(IKM - 1.0D0, IKM = NKM, 1, -1)/)
!       Thicknesses: Rayleigh, gas absorbtion, and total
        TAU_RAY(1:NKM-1) = (/1.6635450D-03, 1.9419680D-03, 2.2679770D-03, &
                             2.6508105D-03, 3.0994260D-03, 3.6273125D-03, &
                             4.2496675D-03, 4.9817975D-03, 5.8451965D-03, &
                             6.8637535D-03, 8.0437120D-03, 9.4089160D-03, &
                             1.1006795D-02, 1.2874230D-02, 1.5060170D-02, &
                             1.7621605D-02, 2.0618885D-02, 2.4125770D-02, &
                             2.8215390D-02, 3.2443020D-02, 3.6710550D-02, &
                             4.1395475D-02, 4.6516100D-02, 5.2115720D-02, &
                             5.8220870D-02, 6.4854970D-02, 7.2055790D-02, &
                             7.9857485D-02, 8.8296165D-02, 9.7407360D-02/)
        TAU_GAS(1:NKM-1) = (/4.0070706D-03, 4.4480361D-03, 4.9292911D-03, &
                             5.4750735D-03, 6.0704857D-03, 6.5661138D-03, &
                             6.9236554D-03, 7.1938558D-03, 7.1913797D-03, &
                             7.0795558D-03, 6.8044136D-03, 6.2428236D-03, &
                             5.5932687D-03, 4.8454180D-03, 4.1839113D-03, &
                             3.7072194D-03, 3.3327532D-03, 3.0502997D-03, &
                             2.6909507D-03, 2.0556686D-03, 1.5176692D-03, &
                             1.1755900D-03, 9.7371357D-04, 9.1875740D-04, &
                             8.9990520D-04, 9.2304693D-04, 9.7556167D-04, &
                             1.0724266D-03, 1.1324552D-03, 1.1508426D-03/)
        TAU_ZKM(1:NKM-1) = TAU_RAY(1:NKM-1) + TAU_GAS(1:NKM-1)
!48
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!48
!       SURFACE
!48
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!48
!       COMPUTE AUXILIARY PARAMETERS
!48
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU_ZKM(1:NKM-1))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!48
!       REDEFINE ATMOSPHERIC PROFILE:
!           1.Split the total Tau = TauR+TauG with increment dTau, get XKM
        CALL SPLITTAU(DTAU, NB, TAU_ZKM(1:NKM-1), NKM-1, &
                      ZKM(1:NKM), NKM, XKM(1:NB))
!           2.Regroup Tau scattering (Rayleigh) from ZKM to XKM                   
        CALL GROUPTAU(NKM, ZKM(1:NKM), NKM-1, TAU_RAY(1:NKM-1), &
                      NB, XKM(1:NB), NL, TAU_SCA(1:NL))
!48
!       Compute SSA. DTAU includes Rayleigh scattering and Gas absorbtion
        SSA(1:NL) = TAU_SCA(1:NL)/DTAU
!48
!       All layers have the same thickness
        TAU(1:NL) = DTAU
!48
!       RUN RT CODE
!48
        CALL CPU_TIME(CPU_TIME1)
!48
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!48
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
        END DO ! IA = 1, NAZ
!48
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!48
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!48
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!48
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!48
        CALL CPU_TIME(CPU_TIME2)
!48
        TIME = CPU_TIME2 - CPU_TIME1
!48
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T048/T048.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!48
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!48
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!48
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!48
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!48
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!48
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!48
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!48
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 048:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!48
        WRITE(*, 50) 'T048: Completed', MAX_ERR_I, AVR_ERR_I,    &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T048: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!48
    END IF ! TESTNUM(ITEST)
!48
!*******************************************************************************
!
!   049. GOAL: To test the case with absrobtion and Rayleigh profiles.
!
!*******************************************************************************
!49
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!49
!       ACCURACY
!49
!       Number of Gauss nodes per hemisphere
        NG1 = 6
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!49
!       GEOMETRY
!49
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 18 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 95.0D0 + UDA(1:NRU)         ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/( 5.0D0*(IK-1), IK = 1, NAZ )/)
!49
!       ATMOSPHERE
!49
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Number of optical layers
        NLR = NL
!       Number of heights in km
        NKM = 31
!       Heights
        ZKM(1:NKM) = (/(IKM - 1.0D0, IKM = NKM, 1, -1)/)
!       Thicknesses: Rayleigh, gas absorbtion, and total
        TAU_RAY(1:NKM-1) = (/1.6635450D-03, 1.9419680D-03, 2.2679770D-03, &
                             2.6508105D-03, 3.0994260D-03, 3.6273125D-03, &
                             4.2496675D-03, 4.9817975D-03, 5.8451965D-03, &
                             6.8637535D-03, 8.0437120D-03, 9.4089160D-03, &
                             1.1006795D-02, 1.2874230D-02, 1.5060170D-02, &
                             1.7621605D-02, 2.0618885D-02, 2.4125770D-02, &
                             2.8215390D-02, 3.2443020D-02, 3.6710550D-02, &
                             4.1395475D-02, 4.6516100D-02, 5.2115720D-02, &
                             5.8220870D-02, 6.4854970D-02, 7.2055790D-02, &
                             7.9857485D-02, 8.8296165D-02, 9.7407360D-02/)
        TAU_GAS(1:NKM-1) = (/4.0070706D-03, 4.4480361D-03, 4.9292911D-03, &
                             5.4750735D-03, 6.0704857D-03, 6.5661138D-03, &
                             6.9236554D-03, 7.1938558D-03, 7.1913797D-03, &
                             7.0795558D-03, 6.8044136D-03, 6.2428236D-03, &
                             5.5932687D-03, 4.8454180D-03, 4.1839113D-03, &
                             3.7072194D-03, 3.3327532D-03, 3.0502997D-03, &
                             2.6909507D-03, 2.0556686D-03, 1.5176692D-03, &
                             1.1755900D-03, 9.7371357D-04, 9.1875740D-04, &
                             8.9990520D-04, 9.2304693D-04, 9.7556167D-04, &
                             1.0724266D-03, 1.1324552D-03, 1.1508426D-03/)
        TAU_ZKM(1:NKM-1) = TAU_RAY(1:NKM-1) + TAU_GAS(1:NKM-1)
!49
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!49
!       SURFACE
!49
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!49
!       COMPUTE AUXILIARY PARAMETERS
!49
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU_ZKM(1:NKM-1))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to different optical layers
        DO IL = 1, NL; NEWLR(IL) = IL; END DO
!49
!       REDEFINE ATMOSPHERIC PROFILE:
!           1.Split the total Tau = TauR+TauG with increment dTau, get XKM
        CALL SPLITTAU(DTAU, NB, TAU_ZKM(1:NKM-1), NKM-1, &
                      ZKM(1:NKM), NKM, XKM(1:NB))
!           2.Regroup Tau scattering (Rayleigh) from ZKM to XKM                   
        CALL GROUPTAU(NKM, ZKM(1:NKM), NKM-1, TAU_RAY(1:NKM-1), &
                      NB, XKM(1:NB), NL, TAU_SCA(1:NL))
!49
!       Compute SSA. DTAU includes Rayleigh scattering and Gas absorbtion
        SSA(1:NL) = TAU_SCA(1:NL)/DTAU
!49
!       All layers have the same thickness
        TAU(1:NL) = DTAU
!49
!       RUN RT CODE
!49
        CALL CPU_TIME(CPU_TIME1)
!49
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!49
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
        END DO ! IA = 1, NAZ
!49
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!49
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!49
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!49
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!49
        CALL CPU_TIME(CPU_TIME2)
!49
        TIME = CPU_TIME2 - CPU_TIME1
!49
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T049/T049.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!49
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!49
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!49
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!49
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!49
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!49
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!49
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!49
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 049:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!49
        WRITE(*, 50) 'T049: Completed', MAX_ERR_I, AVR_ERR_I,    &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T049: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!49
    END IF ! TESTNUM(ITEST)
!49
!*******************************************************************************
!
!   050. GOAL: To test aerosol (dust) + Rayleigh + absorbtion profile
!
!*******************************************************************************
!50
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!50
!       ACCURACY
!50
!       Number of Gauss nodes per hemisphere
        NG1 = 48
!       Total number of the Fourier moments
        NM = 24
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!50
!       GEOMETRY
!50
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 30.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 18 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 95.0D0 + UDA(1:NRU)         ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/( 5.0D0*(IK-1), IK = 1, NAZ )/)
!50
!       ATMOSPHERE
!50
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Aerosol SSA
        SSA_AER = 0.787581D0
!       Number of optical layers
        NLR = NL
!       Number of heights in km
        NKM = 31
!       Heights
        ZKM(1:NKM) = (/(IKM - 1.0D0, IKM = NKM, 1, -1)/)
!       Thicknesses: Rayleigh, gas absorbtion, and total
        TAU_RAY(1:NKM-1) = (/1.223920D-03, 1.436202D-03, 1.673156D-03, &
                             1.962176D-03, 2.278537D-03, 2.643477D-03, &
                             3.096986D-03, 3.630470D-03, 4.259550D-03, & 
                             5.000539D-03, 5.859522D-03, 6.853999D-03, &
                             8.017975D-03, 9.378394D-03, 1.097063D-02, &
                             1.283650D-02, 1.501988D-02, 1.757449D-02, &
                             2.055893D-02, 2.366180D-02, 2.677598D-02, &
                             3.019535D-02, 3.393684D-02, 3.802398D-02, &
                             4.247561D-02, 4.731790D-02, 5.257392D-02, &
                             5.826870D-02, 6.442852D-02, 7.106688D-02/)
        TAU_GAS(1:NKM-1) = (/1.171501D-04, 1.240940D-04, 1.312583D-04, &
                             1.389360D-04, 1.466986D-04, 1.483235D-04, &
                             1.448478D-04, 1.425411D-04, 1.374489D-04, &
                             1.291086D-04, 1.180679D-04, 1.040696D-04, &
                             8.707805D-05, 6.793898D-05, 4.908056D-05, &
                             3.260113D-05, 2.233359D-05, 1.841397D-05, &
                             1.635498D-05, 1.590512D-05, 1.650243D-05, &
                             1.723491D-05, 1.831206D-05, 2.030483D-05, &
                             2.259870D-05, 2.540327D-05, 2.870375D-05, &
                             3.280318D-05, 3.665020D-05, 4.005650D-05/)
        TAU_AER(1:NKM-1) = (/4.46335D-05,   6.4946D-05,    8.52585D-05,   &
                             0.000105571D0, 0.000125883D0, 0.000170829D0, &
                             0.000258579D0, 0.00037385D0,  0.000478735D0, &
                             0.000566486D0, 0.000607506D0, 0.000571678D0, &
                             0.000490678D0, 0.000419023D0, 0.000403446D0, &
                             0.0004346D0,   0.000497947D0, 0.000601275D0, &
                             0.000747699D0, 0.000995261D0, 0.00150682D0,  &
                             0.00264126D0,  0.0048852D0,   0.00709377D0,  &
                             0.00866111D0,  0.0108659D0,   0.0158173D0,   &
                             0.0265917D0,   0.0438819D0,   0.0700112D0/)
        TAU_ZKM(1:NKM-1) = TAU_RAY(1:NKM-1) + TAU_GAS(1:NKM-1) + TAU_AER(1:NKM-1)
!50
!       Phase function
        NK = 1000 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T050/Xk1000_0838.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!50
!       SURFACE
!50
!       Surface model
        ISRF = 0
!       Set of parameters
        PSRF = 0.0D0
!       Number of nodes for the Fourier expansion
        NGA = 2
!50
!       COMPUTE AUXILIARY PARAMETERS
!50
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU_ZKM(1:NKM-1))
!       Integration step over Tau
        DTAU = TAU0/NL
!50
!       REDEFINE ATMOSPHERIC PROFILE:
!           1.Split the total Tau = TauR+TauG with increment dTau, get XKM
        CALL SPLITTAU(DTAU, NB, TAU_ZKM(1:NKM-1), NKM-1, &
                      ZKM(1:NKM), NKM, XKM(1:NB))
!           2.Regroup Tau Rayleigh scattering from ZKM to XKM                 
        CALL GROUPTAU(NKM, ZKM(1:NKM), NKM-1, TAU_RAY(1:NKM-1), &
                      NB, XKM(1:NB), NL, TAU_SCA_RAY(1:NL))
!           3.Regroup Tau aerosol scattering from ZKM to XKM
        TAU_SCA(1:NKM-1) = SSA_AER*TAU_AER(1:NKM-1)                  
        CALL GROUPTAU(NKM, ZKM(1:NKM), NKM-1, TAU_SCA(1:NKM-1), &
                      NB, XKM(1:NB), NL, TAU_SCA_AER(1:NL))
!50
!       All microlayers belong to different optical layers;
!       Gas absorbtion is incliuded in DTAU;
!       RATIO = 0 - no aerosol.
        DO IL = 1, NL
            NEWLR(IL) = IL
            TAU_SCA_AR = TAU_SCA_AER(IL) + TAU_SCA_RAY(IL)
            SSA(IL) = TAU_SCA_AR/DTAU
            RATIO(IL) = TAU_SCA_AER(IL)/TAU_SCA_AR
        END DO ! IL = 1, NL
!50
!       All layers have the same total thickness
        TAU(1:NL) = DTAU
!50
!       RUN RT CODE
!50
        CALL CPU_TIME(CPU_TIME1)
!50
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
        END DO ! IA = 1, NAZ
!50
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!50
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!50
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!50
        CALL CPU_TIME(CPU_TIME2)
!50
        TIME = CPU_TIME2 - CPU_TIME1
!50
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T050/T050.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!50
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!50
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!50
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!50
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!50
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!50
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!50
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!50
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 050:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!50
        WRITE(*, 50) 'T050: Completed', MAX_ERR_I, AVR_ERR_I,    &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T050: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!50
    END IF ! TESTNUM(ITEST)
!50
!*******************************************************************************
!
!   051. GOAL: To test RTLS surface under Rayleigh atmopshere: first scattering
!
!*******************************************************************************
!51
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!51
!       ACCURACY
!51
!       Number of Gauss nodes per hemisphere
        NG1 = 6
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 20
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = -999.0D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 1
!51
!       GEOMETRY
!51
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 18 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 95.0D0 + UDA(1:NRU)         ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/( 5.0D0*(IK-1), IK = 1, NAZ )/)
!51
!       ATMOSPHERE
!51
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.21880749390D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!51
!       SURFACE
!51
!       Surface model
        ISRF = 4
!       Set of parameters
        PSRF(1:3) = (/0.330D0, 0.053D0, 0.066D0/)
!       Number of nodes for the Fourier expansion
        NGA = 180
!51
!       COMPUTE AUXILIARY PARAMETERS
!51
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!51
!       RUN RT CODE
!51
        CALL CPU_TIME(CPU_TIME1)
!51
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!51
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!51
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!51
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!51
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!51
        CALL CPU_TIME(CPU_TIME2)
!51
        TIME = CPU_TIME2 - CPU_TIME1
!51
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T051/T051.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!51
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!51
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!51
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!51
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!51
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!51
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!51
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!51
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 051:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!51
        WRITE(*, 50) 'T051: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T051: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!51
    END IF ! TESTNUM(ITEST)
!51
!*******************************************************************************
!
!   052. GOAL: To test RTLS surface under Rayleigh atmopshere
!
!*******************************************************************************
!52
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!52
!       ACCURACY
!52
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 20
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!52
!       GEOMETRY
!52
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 60.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 18 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 95.0D0 + UDA(1:NRU)         ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/( 5.0D0*(IK-1), IK = 1, NAZ )/)
!52
!       ATMOSPHERE
!52
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.21880749390D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!52
!       SURFACE
!52
!       Surface model
        ISRF = 4
!       Set of parameters
        PSRF(1:3) = (/0.330D0, 0.053D0, 0.066D0/)
!       Number of nodes for the Fourier expansion
        NGA = 180
!52
!       COMPUTE AUXILIARY PARAMETERS
!52
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!52
!       RUN RT CODE
!52
        CALL CPU_TIME(CPU_TIME1)
!52
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!52
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!52
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!52
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!52
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!52
        CALL CPU_TIME(CPU_TIME2)
!52
        TIME = CPU_TIME2 - CPU_TIME1
!52
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T052/T052.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!52
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!52
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!52
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!52
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!52
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!52
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!52
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!52
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 052:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!52
        WRITE(*, 50) 'T052: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T052: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!52
    END IF ! TESTNUM(ITEST)
!52
!*******************************************************************************
!
!   053. GOAL: To test realistic atmopshere (A+R+G) over RTLS surface
!
!*******************************************************************************
!53
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!53
!       ACCURACY
!53
!       Number of Gauss nodes per hemisphere
        NG1 = 48
!       Total number of the Fourier moments
        NM = 24
!       Number of microlayers, defines dTau
        NL = 50
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.000001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!53
!       GEOMETRY
!53
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 30.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 18 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 95.0D0 + UDA(1:NRU)         ! Down
!       View azimuth angles, degrees
        NAZ = 37 ! Number of azimuths
        AZA(1:NAZ) = (/( 5.0D0*(IK-1), IK = 1, NAZ )/)
!53
!       ATMOSPHERE
!53
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 2
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Aerosol SSA
        SSA_AER = 0.787581D0
!       Number of optical layers
        NLR = NL
!       Number of heights in km
        NKM = 31
!       Heights
        ZKM(1:NKM) = (/(IKM - 1.0D0, IKM = NKM, 1, -1)/)
!       Thicknesses: Rayleigh, gas absorbtion, and total
        TAU_RAY(1:NKM-1) = (/1.223920D-03, 1.436202D-03, 1.673156D-03, &
                             1.962176D-03, 2.278537D-03, 2.643477D-03, &
                             3.096986D-03, 3.630470D-03, 4.259550D-03, & 
                             5.000539D-03, 5.859522D-03, 6.853999D-03, &
                             8.017975D-03, 9.378394D-03, 1.097063D-02, &
                             1.283650D-02, 1.501988D-02, 1.757449D-02, &
                             2.055893D-02, 2.366180D-02, 2.677598D-02, &
                             3.019535D-02, 3.393684D-02, 3.802398D-02, &
                             4.247561D-02, 4.731790D-02, 5.257392D-02, &
                             5.826870D-02, 6.442852D-02, 7.106688D-02/)
        TAU_GAS(1:NKM-1) = (/1.171501D-04, 1.240940D-04, 1.312583D-04, &
                             1.389360D-04, 1.466986D-04, 1.483235D-04, &
                             1.448478D-04, 1.425411D-04, 1.374489D-04, &
                             1.291086D-04, 1.180679D-04, 1.040696D-04, &
                             8.707805D-05, 6.793898D-05, 4.908056D-05, &
                             3.260113D-05, 2.233359D-05, 1.841397D-05, &
                             1.635498D-05, 1.590512D-05, 1.650243D-05, &
                             1.723491D-05, 1.831206D-05, 2.030483D-05, &
                             2.259870D-05, 2.540327D-05, 2.870375D-05, &
                             3.280318D-05, 3.665020D-05, 4.005650D-05/)
        TAU_AER(1:NKM-1) = (/4.46335D-05,   6.4946D-05,    8.52585D-05,   &
                             0.000105571D0, 0.000125883D0, 0.000170829D0, &
                             0.000258579D0, 0.00037385D0,  0.000478735D0, &
                             0.000566486D0, 0.000607506D0, 0.000571678D0, &
                             0.000490678D0, 0.000419023D0, 0.000403446D0, &
                             0.0004346D0,   0.000497947D0, 0.000601275D0, &
                             0.000747699D0, 0.000995261D0, 0.00150682D0,  &
                             0.00264126D0,  0.0048852D0,   0.00709377D0,  &
                             0.00866111D0,  0.0108659D0,   0.0158173D0,   &
                             0.0265917D0,   0.0438819D0,   0.0700112D0/)
        TAU_ZKM(1:NKM-1) = TAU_RAY(1:NKM-1) + TAU_GAS(1:NKM-1) + TAU_AER(1:NKM-1)
!53
!       Phase function
        NK = 1000 ! Number of expansion moments in the phase function
        OPEN(FNUM_XK, FILE = PATH//'T053/Xk1000_0838.txt', ACTION = 'READ')
        DO IK = 1, NK; READ(FNUM_XK, *) FK(IK, 1:7); END DO
        CLOSE(FNUM_XK)
        XK(:, 1:3) = FK(:, 2:4) ! a1k, a2k, a3k
        XK(:, 4)   = FK(:, 6)   ! b1k
!53
!       SURFACE
!53
!       Surface model
        ISRF = 4
!       Set of parameters
        PSRF(1:3) = (/0.330D0, 0.053D0, 0.066D0/)
!       Number of nodes for the Fourier expansion
        NGA = 180
!53
!       COMPUTE AUXILIARY PARAMETERS
!53
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU_ZKM(1:NKM-1))
!       Integration step over Tau
        DTAU = TAU0/NL
!53
!       REDEFINE ATMOSPHERIC PROFILE:
!           1.Split the total Tau = TauR+TauG with increment dTau, get XKM
        CALL SPLITTAU(DTAU, NB, TAU_ZKM(1:NKM-1), NKM-1, &
                      ZKM(1:NKM), NKM, XKM(1:NB))
!           2.Regroup Tau Rayleigh scattering from ZKM to XKM                 
        CALL GROUPTAU(NKM, ZKM(1:NKM), NKM-1, TAU_RAY(1:NKM-1), &
                      NB, XKM(1:NB), NL, TAU_SCA_RAY(1:NL))
!           3.Regroup Tau aerosol scattering from ZKM to XKM
        TAU_SCA(1:NKM-1) = SSA_AER*TAU_AER(1:NKM-1)                  
        CALL GROUPTAU(NKM, ZKM(1:NKM), NKM-1, TAU_SCA(1:NKM-1), &
                      NB, XKM(1:NB), NL, TAU_SCA_AER(1:NL))
!53
!       All microlayers belong to different optical layers;
!       Gas absorbtion is incliuded in DTAU;
!       RATIO = 0 - no aerosol.
        DO IL = 1, NL
            NEWLR(IL) = IL
            TAU_SCA_AR = TAU_SCA_AER(IL) + TAU_SCA_RAY(IL)
            SSA(IL) = TAU_SCA_AR/DTAU
            RATIO(IL) = TAU_SCA_AER(IL)/TAU_SCA_AR
        END DO ! IL = 1, NL
!53
!       All layers have the same total thickness
        TAU(1:NL) = DTAU
!53
!       RUN RT CODE
!53
        CALL CPU_TIME(CPU_TIME1)
!53
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!       Get the F11 and F12 elements from the moments:
        CALL SUMKA1B1(XK(1:NK, 1), XK(1:NK, 4), NK, MUS(1:NUD, 1:NAZ), NUD,    &
                      NAZ, F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ))
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
        END DO ! IA = 1, NAZ
!53
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!53
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!53
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!53
        CALL CPU_TIME(CPU_TIME2)
!53
        TIME = CPU_TIME2 - CPU_TIME1
!53
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T053/T053.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!53
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!53
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!53
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!53
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!53
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!53
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!53
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!53
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 053:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!53
        WRITE(*, 50) 'T053: Completed', MAX_ERR_I, AVR_ERR_I,    &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T053: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!53
    END IF ! TESTNUM(ITEST)
!53
!*******************************************************************************
!
!   054. GOAL: To test the Nadal-Breon vegetation model.
!
!*******************************************************************************
!54
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!54
!       ACCURACY
!54
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 10
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!54
!       GEOMETRY
!54
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 45.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 17 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 100.0D0 + UDA(1:NRU)        ! Down
!       View azimuth angles, degrees
        NAZ = 73 ! Number of azimuths
        AZA(1:NAZ) = (/(5.0D0*(IK-1), IK = 1, NAZ)/)
!54
!       ATMOSPHERE
!54
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.1D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!54
!       SURFACE
!54
!       Surface model: Nadal-Breon vegetation model, ~ 1/(mu_i + mu_r)
        ISRF = 8
!       Set of parameters
        PSRF(1) = 1.50D0  ! Refractive index: real part
        PSRF(2:) = 0.0D0  ! Other surface parameters are not used
!       Number of nodes for the Fourier expansion
        NGA = 180
!54
!       COMPUTE AUXILIARY PARAMETERS
!54
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!54
!       RUN RT CODE
!54
        CALL CPU_TIME(CPU_TIME1)
!54
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!54
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!54
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)                           
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!54
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!54
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!54
        CALL CPU_TIME(CPU_TIME2)
!54
        TIME = CPU_TIME2 - CPU_TIME1
!54
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T054/T054.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!54
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!54
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!54
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!54
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!54
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!54
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!54
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!54
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 054:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!54
        WRITE(*, 50) 'T054: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T054: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!54
    END IF ! TESTNUM(ITEST)
!54
!*******************************************************************************
!
!   055. GOAL: To test the Nadal-Breon vegetation model.
!
!*******************************************************************************
!55
    ITEST = ITEST+1
    IF (TSTNUM(ITEST)) THEN ! Number of orders are selected manually
!55
!       ACCURACY
!55
!       Number of Gauss nodes per hemisphere
        NG1 = 8
!       Total number of the Fourier moments
        NM = 3
!       Number of microlayers, defines dTau
        NL = 10
!       Absolute error for automatic convergence EPSI < 0 => use NO
        EPSI = 0.00001D0
!       Number of scattering orders to be computed: only if EPSI < 0
        NO = 20
!55
!       GEOMETRY
!55
!       Solar angle, SZA (degrees), and/or mu0=cos(SZA)
        SZA = 45.0D0
        MU0 = DCOS(SZA*D2R)
!       User-defined zenith angles, UDA (degrees)
        NRU = 17 ! Number of reflected angles
        NTU = NRU ! Number of transmitted angles
        NUD = NRU+NTU ! Number of user-defiend angles
        UDA(1:NRU) = (/(5.0D0*(IK-1), IK = 1, NRU)/) ! Up
        UDA(NRU+1:NUD) = 100.0D0 + UDA(1:NRU)        ! Down
!       View azimuth angles, degrees
        NAZ = 73 ! Number of azimuths
        AZA(1:NAZ) = (/(5.0D0*(IK-1), IK = 1, NAZ)/)
!55
!       ATMOSPHERE
!55
!       Number of components: 1 - Rayleigh only, 2 - R&A
        NC = 1
!       Rayleigh depolarization factor (N/U)
        DEPF = 0.03D0
!       Number of optical layers
        NLR = 1
!       Optical thickness & SSA of EACH optical layer
        TAU(1:NLR) = 0.1D0
        SSA(1:NLR) = 1.0D0
!       Phase function
        NK = 3  ! Number of expansion moments in the phase function
!55
!       SURFACE
!55
!       Surface model: Nadal-Breon vegetation model, ~ 1/(mu_i + mu_r)
        ISRF = 9
!       Set of parameters
        PSRF(1) = 1.50D0  ! Refractive index: real part
        PSRF(2:) = 0.0D0  ! Other surface parameters are not used
!       Number of nodes for the Fourier expansion
        NGA = 180
!55
!       COMPUTE AUXILIARY PARAMETERS
!55
!       Number of Gauss nodes per whole sphere
        NG2 = NG1*2
!       Total number of nodes for ascending radiation
        NUP = NG1+NRU
!       Same as NUP, but for descending radiation
        NDN = NG1+NTU
!       Total number of nodes, Gauss & dummy
        NMU = NG2+NUD
!       Number of moments in the scattering integral
        NKS = MIN(NK, NG2)
!       Number of boundaries at microlayers dTau
        NB = NL+1
!       Total optical thickness of atmosphere
        TAU0 = SUM(TAU(1:NLR))
!       Integration step over Tau
        DTAU = TAU0/NL
!       All microlayers belong to the same optical layer No.1
        NEWLR(1:NL) = 1
!55
!       RUN RT CODE
!55
        CALL CPU_TIME(CPU_TIME1)
!55
!       Get zenith and azimuth nodes
        CALL NODES(NM, ISRF, UDA(1:NUD), AZA(1:NAZ), NUD, NRU, NTU, NG1, NG2,  &
                   NMU, NAZ, NGA, MU(1:NMU), MUD(1:NUD), Z(1:NG2), W(1:NG2),   &
                   AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM), WSMA(1:NGA, 1:NM))
!       Compute cosine of scattering angle
        CALL CSCATANG(MU0, MUD(1:NUD), NUD, AZI(1:NAZ), NAZ, MUS(1:NUD, 1:NAZ))
!
        DO IA = 1, NAZ
            CALL RAYLEIGH12(MUS(1:NUD, IA), NUD, DEPF, R11(1:NUD, IA),         &
                                                       R12(1:NUD, IA))
                            
        END DO ! IA = 1, NAZ
!55
!       No aerosol
        RATIO(1:NLR)      = 0.0D0
        XK(1:NKS, 1:4)    = 0.0D0
        F11(1:NUD, 1:NAZ) = 0.0D0
        F12(1:NUD, 1:NAZ) = 0.0D0
!55
!       Solve the VRTE. On output, the Stokes vector is normalized by 2PI
        CALL SORD_IP(NM, NO, NC, ISRF, EPSI, MU0, TAU0, DTAU, DEPF,            &
                     TAU(1:NLR), SSA(1:NLR), RATIO(1:NLR), XK(1:NKS, 1:4),     &
                     F11(1:NUD, 1:NAZ), F12(1:NUD, 1:NAZ), R11(1:NUD, 1:NAZ),  &
                     R12(1:NUD, 1:NAZ), PSRF, MU(1:NMU), Z(1:NG2), W(1:NG2),   &
                     AZI(1:NAZ), ZAZ(1:NGA), WCMA(1:NGA, 1:NM),                &
                     WSMA(1:NGA, 1:NM), &
                     NEWLR(1:NL), NKS, NLR, NL, NB, NMU, NAZ, NUP, NDN, NUD,   &
                     NRU, NG1, NG2, NGA, II(1:NUD, 1:NAZ), QQ(1:NUD, 1:NAZ),   &
                     UU(1:NUD, 1:NAZ), FLUXI)                           
!       II, QQ, UU are normalized by 2pi - coincide with pi/mu0 in the benchmark
!55
!       NORMALIZATION CONSTANT (FLUX ON TOA)
!55
!       The benchmark results are normalized to 1.0: devide by 2pi=6.28...
        II = II/PI2; QQ = QQ/PI2; UU = UU/PI2
!55
        CALL CPU_TIME(CPU_TIME2)
!55
        TIME = CPU_TIME2 - CPU_TIME1
!55
        NY = NUD*NAZ
        OPEN(FNUM_BM, FILE = PATH//'T055/T055.txt', ACTION = 'READ')
        DO IY = 1, NY; READ(FNUM_BM, *) BMARK(IY, :); END DO
        CLOSE(FNUM_BM)
!55
        IY = 0
        DO IA = 1, NAZ
            DO IK = 1, NUD
                IY = IY+1
                OUT(IY, 1:5) = BMARK(IY, 1:5) ! Geometry
!55
                I_BMRK = BMARK(IY, 6)
                Q_BMRK = BMARK(IY, 7)
                U_BMRK = BMARK(IY, 8)
                P_BMRK = DSQRT(Q_BMRK*Q_BMRK + U_BMRK*U_BMRK)/I_BMRK
!55
                I_SORD = II(IK, IA)
                Q_SORD = QQ(IK, IA)
                U_SORD = UU(IK, IA)
!55
                IF ( DABS(Q_SORD/I_SORD) < TINY ) Q_SORD = 0.0D0
                IF ( DABS(U_SORD/I_SORD) < TINY ) U_SORD = 0.0D0
!55
                P_SORD = DSQRT(Q_SORD*Q_SORD + U_SORD*U_SORD)/I_SORD
!55
                OUT(IY,  6) = 100.0D0*(1.0D0 - I_SORD/I_BMRK)
                OUT(IY,  7) = P_BMRK - P_SORD
                OUT(IY,  8) = I_BMRK
                OUT(IY,  9) = I_SORD
                OUT(IY, 10) = P_BMRK
                OUT(IY, 11) = P_SORD
                OUT(IY, 12) = Q_BMRK
                OUT(IY, 13) = Q_SORD
                OUT(IY, 14) = U_BMRK
                OUT(IY, 15) = U_SORD
            END DO ! IK = 1, NUD
        END DO ! IA = 1, NAZ
!55
        MAX_ERR_I = MAXVAL(DABS(OUT(1:NY, 6)))
        AVR_ERR_I = SUM(DABS(OUT(1:NY, 6)))/NY
        MAX_ERR_P = MAXVAL(DABS(OUT(1:NY, 7)))
        AVR_ERR_P = SUM(DABS(OUT(1:NY, 7)))/NY
!55
        WRITE(FNUM, *)
        WRITE(FNUM, *) 'TEST 055:'
        WRITE(FNUM, *)
        DO IY = 1, NY; WRITE(FNUM, 30) OUT(IY, :); END DO
        WRITE(FNUM, *)
        WRITE(FNUM, 40) 'time = ', TIME, ' sec.'
!55
        WRITE(*, 50) 'T055: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
        WRITE(FSCR, 50) 'T055: Completed', MAX_ERR_I, AVR_ERR_I, &
                                        MAX_ERR_P, AVR_ERR_P, TIME
!55
    END IF ! TESTNUM(ITEST)
!55
WRITE(FSCR, *)
WRITE(*, *)
END DO ! IRUN = 1, NRUN
!0
    CLOSE(FNUM) ! Close the log file
    CLOSE(FSCR) ! Close the print screen file
!0
    WRITE(*, *)
    WRITE(*, *) 'Done:'
    WRITE(*, *) '    '//PATH//'log.txt'
    WRITE(*, *) '    '//PATH//'scr.txt'
!    READ(*, *)
!0
10 FORMAT(1X, 256A)
20 FORMAT(19X, A, 4(3X, A))
30 FORMAT(1F6.2, 1F10.4, 2F10.2, 1F10.4, 1F8.2, 1F10.4, 2ES17.7, 2F10.4, 4ES17.7)
40 FORMAT(1X, A, F10.3, A)
50 FORMAT(1X, A, 2F8.2, 2F9.4, F12.3)
60 FORMAT(1X, A, 2ES16.7, 1F9.4)
!0
END PROGRAM AAA_TESTS_SORD_IP
