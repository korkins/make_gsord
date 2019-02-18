SUBROUTINE ZM_IP(M, NC, SSA, RTO, NLR, NG2, NMU, ZR11, ZR12, ZR13,             &
                                                 ZR21, ZR22, ZR23,             &
                                                 ZR31, ZR32, ZR33,             &
                                                       ZA11, ZA12, ZA13,       &
                                                       ZA21, ZA22, ZA23,       &
                                                       ZA31, ZA32, ZA33,       &
                                                             ZM11, ZM12, ZM13, &
                                                             ZM21, ZM22, ZM23, &
                                                             ZM31, ZM32, ZM33)
!===============================================================================
! PURPOSE:
!   To compute mixture of the Rayleigh (R) and Aerosol (A) phase matrices for a
!   given Fourier order, M, and for each optical layer, ILR.
!
! INPUT:
!   M       I(1)          Fourier order, M = 0, 1, 2, 3 ...
!   NC      I(1)          Number of components: 1 - R; 2 - R & A
!   SSA     D(NLR)        Single scattering albedo in optical layers
!   RTO     D(NLR)        Aerosol-Rayleigh ratio in optical layers
!   NLR     I(1)          Number of optical layers
!   NG2     I(1)          Number of Gauss nodes per whole sphere
!   NMU     I(1)          Number of Gauss and user defined view zeniths
!   ZRij    D(NG2, NMU)   Mth Fourier moment for the Rayleigh phase matrix
!   ZAij    D(NG2, NMU)   Same as ZRij except for Aerosol
!
! OUTPUT:
!   ZMij    D(NG2, NMU, NLR)   Elements of the mixture in each optical layer
!
! TREE:
!   -
!
! COMMENTS:
!   On output, Zij contains the factor of SSA/2.
!
!   The Fourier order, M, is used only to avoid mixing components for M > 2:
!   no Rayleigh in this case.
!
!   RTO(ILR) = 0.0 means no aerosol in the layer ILR (pure Rayleigh scattering).
!
! REFERENCES:
!   -
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: M, NC, NLR, NG2, NMU
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NLR), INTENT(IN) :: SSA, RTO
    REAL*8, DIMENSION(NG2, NMU), INTENT(IN) :: ZR11, ZR12, ZR13, &
                                               ZR21, ZR22, ZR23, &
                                               ZR31, ZR32, ZR33, &
                                               ZA11, ZA12, ZA13, &
                                               ZA21, ZA22, ZA23, &
                                               ZA31, ZA32, ZA33
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NG2, NMU, NLR), INTENT(OUT) :: ZM11, ZM12, ZM13, &
                                                     ZM21, ZM22, ZM23, &
                                                     ZM31, ZM32, ZM33
!
! LOCAL VARIABLES
    INTEGER &
        ILR      ! Loop index over optical layers
    REAL*8 &
        SSA05, & ! SSA/2
        WA,    & ! Weight of Aerosol in mixture including SSA/2
        WR       ! Same as WA except for Rayleigh
!===============================================================================
!
    IF (NC == 1) THEN ! Rayleigh only
        DO ILR = 1, NLR
            SSA05 = 0.5D0*SSA(ILR)
            ZM11(:, :, ILR) = SSA05*ZR11
            ZM12(:, :, ILR) = SSA05*ZR12
            ZM13(:, :, ILR) = SSA05*ZR13
            ZM21(:, :, ILR) = SSA05*ZR21
            ZM22(:, :, ILR) = SSA05*ZR22
            ZM23(:, :, ILR) = SSA05*ZR23
            ZM31(:, :, ILR) = SSA05*ZR31
            ZM32(:, :, ILR) = SSA05*ZR32
            ZM33(:, :, ILR) = SSA05*ZR33
        END DO ! ILR = 1, NLR
    ELSE ! NC > 1: Rayleigh & Aerosol
        IF (M < 3) THEN ! Rayleigh: m = 0, 1, 2
            DO ILR = 1, NLR
                SSA05 = 0.5D0*SSA(ILR)
                WA = SSA05*RTO(ILR)
                WR = SSA05 - WA
                ZM11(:, :, ILR) = WR*ZR11 + WA*ZA11
                ZM12(:, :, ILR) = WR*ZR12 + WA*ZA12
                ZM13(:, :, ILR) = WR*ZR13 + WA*ZA13
                ZM21(:, :, ILR) = WR*ZR21 + WA*ZA21
                ZM22(:, :, ILR) = WR*ZR22 + WA*ZA22
                ZM23(:, :, ILR) = WR*ZR23 + WA*ZA23
                ZM31(:, :, ILR) = WR*ZR31 + WA*ZA31
                ZM32(:, :, ILR) = WR*ZR32 + WA*ZA32
                ZM33(:, :, ILR) = WR*ZR33 + WA*ZA33
            END DO ! ILR = 1, NLR
        ELSE ! M > 2: no Rayleigh
            DO ILR = 1, NLR
                SSA05 = 0.5D0*SSA(ILR)
                WA = SSA05*RTO(ILR)
                ZM11(:, :, ILR) = WA*ZA11
                ZM12(:, :, ILR) = WA*ZA12
                ZM13(:, :, ILR) = WA*ZA13
                ZM21(:, :, ILR) = WA*ZA21
                ZM22(:, :, ILR) = WA*ZA22
                ZM23(:, :, ILR) = WA*ZA23
                ZM31(:, :, ILR) = WA*ZA31
                ZM32(:, :, ILR) = WA*ZA32
                ZM33(:, :, ILR) = WA*ZA33
            END DO ! ILR = 1, NLR
        END IF ! M < 3
    END IF ! NC = 1
!
END SUBROUTINE ZM_IP
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 03Oct15 - NEWLR, NL are removed from input (no need to use them!)
!
! 07Jul15 - First created from corresponding (totally modified) piece of SORD_IP
!===============================================================================