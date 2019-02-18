SUBROUTINE ZM0_IP(M, NC, SSA, RTO, NLR, NG2, ZR01, ZR02, ZR03,                 &
                                                   ZA01, ZA02, ZA03,           &
                                                         ZM01, ZM02, ZM03)
!===============================================================================
! PURPOSE:
!   To compute a mixture of the Rayleigh (R) and Aerosol (A) phase matrices for
!   a given Fourier order, M, and optical layer, ILR, to be used for single
!   scattering of unpolarized light.
!
! INPUT:
!   M       I(1)     Fourier order, M = 0, 1, 2, 3 ...
!   NC      I(1)     Number of components: 1 - R; 2 - R & A
!   SSA     D(NLR)   Single scattering albedo in optical layers
!   RTO     D(NLR)   Aerosol-Rayleigh ratio in optical layers
!   NLR     I(1)     Number of optical layers
!   NG2     I(1)     Number of Gauss nodes per whole sphere
!   ZR0i    D(NG2)   Mth Fourier moment for the Rayleigh phase matrix (1st col.)
!   ZA0i    D(NG2)   Mth Fourier moment for the aerosol phase matrix (1st col.)
!
! OUTPUT:
!   ZM0i    D(NG2, NLR)   Elements of mixed phase matrix (1st column only)
!
! TREE:
!   -
!
! COMMENTS:
!   On output, ZM0i contain the factor of SSA/2.
!
!   The Fourier order, M, is used only to avoid mixing components for M > 2 not
!   existing in Rayleigh phase matrix.
!
!   RTO(ILR) = 0.0 means no Aerosol in the layer ILR (pure Rayleigh scattering)
!
! REFERENCES:
!   -
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: M, NC, NLR, NG2
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NLR), INTENT(IN) :: SSA, RTO
    REAL*8, DIMENSION(NG2), INTENT(IN) :: ZR01, ZR02, ZR03, &
                                          ZA01, ZA02, ZA03
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NG2, NLR), INTENT(OUT) :: ZM01, ZM02, ZM03
!
! LOCAL VARIABLES
    INTEGER &
        ILR      ! Loop index over optical layers
    REAL*8 &
        SSA05, & ! SSA/2
        WA,    & ! Weight of aerosol in mixture including SSA/2
        WR       ! Same as WA, except for Rayleigh
!===============================================================================
!
    IF (NC == 1) THEN ! Rayleigh only
        DO ILR = 1, NLR
            SSA05 = 0.5D0*SSA(ILR)
            ZM01(:, ILR) = SSA05*ZR01
            ZM02(:, ILR) = SSA05*ZR02
            ZM03(:, ILR) = SSA05*ZR03
        END DO ! ILR = 1, NLR
    ELSE ! NC > 1: Rayleigh & Aerosol
        IF (M < 3) THEN ! Rayleigh: M = 0, 1, 2
            DO ILR = 1, NLR
                SSA05 = 0.5D0*SSA(ILR)
                WA = SSA05*RTO(ILR)
                WR = SSA05 - WA
                ZM01(:, ILR) = WR*ZR01 + WA*ZA01
                ZM02(:, ILR) = WR*ZR02 + WA*ZA02
                ZM03(:, ILR) = WR*ZR03 + WA*ZA03
            END DO ! IL = 1, NL
        ELSE ! M > 2: no Rayleigh
            DO ILR = 1, NLR
                SSA05 = 0.5D0*SSA(ILR)
                WA = SSA05*RTO(ILR)
                ZM01(:, ILR) = WA*ZA01
                ZM02(:, ILR) = WA*ZA02
                ZM03(:, ILR) = WA*ZA03
            END DO ! ILR = 1, NLR
        END IF ! M < 3
    END IF ! NC = 1
!
END SUBROUTINE ZM0_IP
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 03Oct15 - NEWLR, NL are removed from input (no need to use them!)
!
! 07Jul15 - First created from ZM_IP
!===============================================================================