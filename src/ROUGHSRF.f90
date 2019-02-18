SUBROUTINE ROUGHSRF(MSRF, U, WAZ, MUI, MUR, AZI, NA, WSRF)
!===============================================================================
! PURPOSE:
!   To compute the scalar coefficient, WSRF, that describes the waved surface
!   as a function of geometry of irradiation and observation, and the wind
!   speed [m/s] at 10 meters above the surface. WAVED_SURFACE = WSRF*R2*F*R1,
!   where R2*F*R1 is the rotated Fresnel matrix (subroutine FRESNELR).
!
! INPUT:
!   MSRF   I(1)    1 for azimuthally independent (Nakajima-Tanaka) ocean
!                  surface, 2 for azimuthally dependent (Gram-Charlier)
!                  ocean surface
!   U      D(1)    Wind speed [m/s] 10 meters above the water surface
!   WAZ    D(1)    Relative azimuth, RAD, of the wind direction (for MSRF = 2)
!   MUI    D(1)    Cosine of incident zenith angle, MUI > 0
!   MUR    D(1)    Cosine of reflected zenith angle, MUR < 0
!   AZI    D(NA)   Relative azimuth of observation, radians. AZI = [0:2PI]
!   NA     I(1)    Number of azimuth angles
!
! OUTPUT:
!   WSRF   D(NA)   Coefficient that simulate waves as a function of azimuth
!
! TREE:
!   -
!
! COMMENTS:
!   This code was build using codes SHARM (Dr.Alexei Lyapustin, NASA GSFC) and
!   ocean.phase.f (Dr.Michael Mishchenko, NASA GISS, available online at
!   http://www.giss.nasa.gov/staff/mmishchenko/brf/).
!
!   The complementary error function, REAL*8 DERFC, is computed using intrinsic
!   function. If not available, derfc.f function should be downloaded from
!   www.netlib.org with dependencies). ERFC is defined as
!
!   erfc(x) = 1 - erf(x) = 2/sqrt(pi)*integral(exp(-t*t)dt, x, +inf)         (1)
!
!   The bidirectional shadowing function is
!
!   S = 1/(1 + F(mu) + F(mu')),                                              (2)
!
!   where
!
!   F = 1/2*{ exp(-v*v)/(sqrt(pi)*v) - erfc(v) },                            (3)
!
!   erfc(x) is defined in Eq.(1) and v = mu/(S*sqrt(1 - mu2)), S is the slope
!   distribution parameter usually defined as 'sigma'. Eq.(3) is taken from
!   [1, Eq.(18)] with a corrected misprinting (compare with [2, Eq.(6)] or with
!   [3, Eq.(88)]).
!
!   NOTE: F(mu) = 0 if no shadowing is considered.
!
!   NOTE: In this subroutine, S2 = sigma^2 is defined following code SHARM and
!   differs from that defined in ocean.phase.f by the factor of 2!
!
!   NOTE: Eq.(2) defined in [1, Eq.(17)] and in [3, Eq.(87)] differs from that
!   defined in  [2, Eq.(7)].
!
!   The probability distribution function is given in [2, Eq.(4)]
!
!   P = a^2/(pi*mu*S^2)*exp((1 - 2a)/S^2),                                   (4)
!   a = (1 + cos(2h))/(mu + mu')^2,                                          (5)
!   cos(2h) = mu*mu' - smu*smu'*cos(fi-fi'),                                 (6)
!
!   where 2h is the angle complementary to the scattering angle (h is the angle
!   of reflection = angle of incidence).
!
!   NOTE: the signs in Eqs.(4)-(6) are given for MUI > 0, and MUR > 0 [2].
!
!   The BRDF for the rough sea surface is defined as [2, Eq.(3)]
!
!   WSRF = P * S.                                                            (7)
!
!   In order to bring the boundary condition to the same form as over land, the
!   BRDF is defined following SHARM [4, Appendix, the last equation of section
!   'C.Ocean', unnumbered]
!
!   WSRF = P * S * pi/mu.                                                    (8)
!
!   NOTE: 1/pi in Eq.(4) and *pi in Eq.(8) are cancelled.
!
!   NOTE: the typical expression for 'sigma' is [1, Eq.(15); 2, p.4250]
!
!   S^2 = 0.00534*U,                                                         (9)
!
!   another definition was used in [3, Eq.18]
!
!   2*S^2 = 0.00512*U + 0.003,                                              (10)
!
!   and in SCIATRAN code. Definitions Eq.(9) and (10) are equivalent, i.e. give
!   the same result on output if S^2 = 2*S^2.
!
!   For the case of azimuth dependence, the Gram-Charlier distribution is used
!   [5, Eq.(18)]. It has a long coefficient, and is not typed in here.
!   Eq.(4) in case of Gram Charlier is also slightly changed [4, Eq.(A.24)].
!   All equations are used in the
!
!       IF (I == 0) THEN
!               NT model
!       ELSE
!           ->> GC model <<-
!       END IF
!
!   statement below explicitly. Computational relations are taken from SHARM [4]
!   Note that Eq.(18) from [5] is valid for ZU1(eta)=ZC1(ksi)=2.5 (next line,
!   immediately after Eq.(18) in [5]). This subroutine DOES NOT CHECK if the
!   Gram-Charlier representation is adequate for the given input parameters.
!
! REFERENCES:
!   1. Nakajima T, Tanaka M, JQSRT(1983), V.29, N.6, pp.521-537.
!   2. Gordon HR, Wang M, Appl.Opt.(1992), V.31, N.21, pp.4247-4260.
!   3. Mishchenko MI, Travis LD, JGR(1997), V.102, D14, pp.16989-17013.
!   4. Lyapustin AI, Appl.Opt.(2005), V.44, N.36, pp.7764-7772.
!   5. Cox C, Munk W, JOSA (1954), V.44, N.11, pp.838-850.
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    REAL*8, PARAMETER :: &
        C40 = 0.40D0,                 & ! Gram-Charlier
        C22 = 0.12D0,                 & !      expansion
        C04 = 0.23D0,                 & !          moments
        K = 1.5D0,                    & ! K = SGMU/SGMC, [4, Ap.D]
        NADIR = 1.0D-8,               & ! Near-nadir direction
        PI = 3.1415926535897932D0,    & ! The number PI
        SQRPI = 1.7724538509055160D0, & ! Square root of PI
        SR2 = 1.4142135623730950D0,   & ! Square root of 2.0
        S1 = 5.12D-3,                 & ! SGM2 = S1*U + S2
        S2 = 3.00D-3                    ! SGM2 = S1*U + S2
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: MSRF, NA
    REAL*8, INTENT(IN)  :: U, WAZ, MUI, MUR
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NA), INTENT(IN) :: AZI
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NA), INTENT(OUT) :: WSRF
!
! LOCAL VARIABLES
    REAL*8 &
        CWAZ,  & ! Cosine of azimuth of wind direction
        C21,   & ! Gram Charlier
        C03,   & !   expansion parameters
        FI,    & ! Shadowing factor for the incident direction
        FR,    & ! Shadowing factor for the reflected direction
        K2,    & ! K2 = K*K = SGMU2/SGMC2
        K3,    & ! K3 = 1.0D0/(1.0D0 + K2)
        MUP,   & ! -MUR; MUR < 0, MUP > 0
        SWAZ,  & ! Sine of azimuth of wind direction
        SMUI,  & ! Sine of the incident angle
        SMUR,  & ! Cosine of the incident angle
        SGM,   & ! Slope distribution parameter
        SGM2,  & ! Square of the slope distribution parameter
        SGMC,  & ! Slope distribution parameter: cross-wind
        SGMC2, & ! Square of the slope distribution parameter: cross-wind
        SGMU,  & ! Slope distribution parameter: up-wind
        SGMU2, & ! Square of the slope distribution parameter: up-wind
        X,     & ! Temporary variable for different purposes
        X2       ! Square of X
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NA) :: &
        A,     & ! Slope magnitude, Eq.(5)
        CAZ,   & ! Cosine of azimuth of observation
        C2INC, & ! Cosine of the doubled angle of incidence, Eq.(6)
        GC,    & ! Gram-Charlier series
        P,     & ! Probability density function, Eq.(4)
        SAZ,   & ! Sine of azimuth of observation
        ZC,    & ! Cross-wind
        ZC1,   & !   component
        ZC2,   & ! ZC2 = ZC1*ZC1
        ZU,    & ! Up-wind
        ZU1,   & !   component
        ZU2      ! ZU2 = ZU1*ZU1
!
! FUNCTION (uncomment for pgf90)
!    REAL*8 DERFC ! erfc(x) = 1 - erf(x)
!===============================================================================
!
    MUP = -MUR ! In equations below, MUR is assumed positive
    CAZ = DCOS(AZI)
!
    SGM2 = S1*U + S2
    SGM = DSQRT(SGM2)
!
! Shadowing factor for the incident direction
    SMUI = DSQRT(1.0D0 - MUI*MUI)
    IF (SMUI < NADIR) THEN
        FI = 0.0D0 ! No shadowing
    ELSE
        X =  MUI/SMUI/SGM
        X2 = X*X
        FI = DEXP(-X2)/(SQRPI*X)
        FI = 0.5D0*(FI - DERFC(X)) ! DERFC IS NOT ALWAYS AN INTRINSIC FUNCTION
    END IF ! SMUI < NADIR
!
! Shadowing factor for the reflected direction
    SMUR = DSQRT(1.0D0 - MUP*MUP)
    IF (SMUR < NADIR) THEN
        FR = 0.0D0 ! No shadowing
    ELSE
        X = MUP/SMUR/SGM
        X2 = X*X
        FR = DEXP(-X2)/(SQRPI*X)
        FR = 0.5D0*(FR - DERFC(X)) ! DERFC IS NOT ALWAYS AN INTRINSIC FUNCTION
    END IF ! SMUR < NADIR
!
!****************************UNCOMMENT FOR SHADOWS OFF**************************
!    FI = 0.0D0
!    FR = 0.0D0
!*******************************************************************************
!
!   Scattering angle (between v and -v') = 2*Angle_of_incidence
    C2INC = MUI*MUP - SMUI*SMUR*CAZ
    X = MUI + MUP
    X2 = X*X
    A = (1.0D0 + C2INC)/X2
!   The probability density function of the slope distribution
    IF (MSRF == 1) THEN
!       Azimuthally symmetric (NT) model
!       P == Y = 
!       = 1/(4 mu Sigma**2 mu_n**4) exp[ -(1 - mu_n**2)/ (Sigma**2 mu_n**2) ] =
!       = pi/(4 mu mu_n) P(mu_n) where P(mu_n) is defined in
!       [4, Eq. after (A21)]
        P = A*A*DEXP((1.0D0 - 2.0D0*A)/SGM2) / (MUP*SGM2)
    ELSE ! MSRF == 1
!       Observation
        SAZ = DSIN(AZI)
!       Wind direction
        CWAZ = DCOS(WAZ)
        SWAZ = DSIN(WAZ)
!
        ZC = SMUR*SAZ/X ! Cross-wind component. Refer to [4, Eq.(A22)]
        ZU = -(SMUR*CAZ - SMUI)/X ! Up-wind component Refer to [4, Eq.(A22)].
!       NOTE: Vladimir Rozanov reported different sigh at ZU on 03Jul12
!
!       Azimuthally dependent (GC) model
        K2 = K*K
        K3 = 1.0D0/(1.0D0 + K2)
        SGMC2 = SGM2/K3
        SGMU2 = SGM2*K2/K3
        SGMU = DSQRT(SGMU2)
        SGMC = DSQRT(SGMC2)
!
        ZU1 = ( CWAZ*ZU + SWAZ*ZC)/SGMU ! [4, Eq.(A24), eta]
        ZC1 = (-SWAZ*ZU + CWAZ*ZC)/SGMC ! [4, Eq.(A24), ksi]
!
        ZC2 = ZC1*ZC1
        ZU2 = ZU1*ZU1
!
        P = A*A*DEXP(-(ZC2+ZU2)/2.0D0) / (2.0D0*MUP*SGMU*SGMC) !/PI;[4,Eq.(A24)]
!
!       GC expansion coefficients
        C21 = 0.01D0 - 0.0086D0*U
        C03 = 0.04D0 - 0.0330D0*U
        GC = 1.0D0 - 0.5D0*C21*(ZC2-1.0D0)*ZU1 - C03/6.0D0*(ZU2-3.0D0)*ZU1 + &
			    C40/24.0D0*(ZC2*ZC2-6.0D0*ZC2 + 3.0D0) + &
			        C22/4.0D0*(ZC2-1.0D0)*(ZU2-1.0D0) + &
			            C04/24.0D0*(ZU2*ZU2 - 6.0D0*ZU2 + 3.0D0)
	    P = P*GC
    END IF ! MSRF == 1
!   BRDF, Eq.(8)
    WSRF = ( P/(1.0D0 + FI + FR) )/MUI
!   The factor of PI in P and WSRF is cancelled following code SHARM,
!   NT subroutine, BCondition.cpp.
!
END SUBROUTINE ROUGHSRF
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 11Oct14 - Renamed: TYP -> MSRF. Checked against old version.
!
! 10Sep14 - Created from ROUGHSRF for a single point {SZA, VZA, AZA}. Tested
!           against it for SZA=0:1:89, VZA=0:1:89, AZA=0:1:360(2PI) for
!           U=2.5, WAZ=75.0*PI/180, TYP = 1 & 2. Exact agreement in REAL*8. Ok.
!
!           The new version of ROUGHSRF is ~4 times faster.
!===============================================================================