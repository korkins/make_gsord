SUBROUTINE Z0_IP(NC, SSA, RTO, NLR, A11, A12, R11, R12, NMU, NAZ, Z11, Z12)
!===============================================================================
! PURPOSE:
!   To compute a mixture of the Rayleigh (R) and Aerosol (A) phase matrices in
!   all optical layers to be used for single scattering of unpolarized light.
!
! INPUT:
!   NC    I(1)          Number of components: 1 - R; 2 - R & A
!   SSA   D(NLR)        Single scattering albedo in optical layers
!   RTO   D(NLR)        Aerosol-Rayleigh ratio in optical layers
!   NLR   I(1)          Number of optical layers
!   A11   D(NMU, NAZ)   The [11]-element of the phase matrix
!   A12   D(NMU, NAZ)   The [12]=[21]-element of the phase matrix
!   R11   D(NMU, NAZ)   The [11]-element of the Rayleigh phase matrix
!   R12   D(NMU, NAZ)   The [12]=[21]-element of the Rayleigh phase matrix
!   NMU   I(1)          Number of view zeniths
!   NAZ   I(1)          Number of azimuths
!
! OUTPUT:
!   Z11, Z12   D(NMU, NAZ, NLR)   Elements of mixture in each optical layer
!
! TREE:
!   -
!
! COMMENTS:
!   On output, Z11 & Z12, contain the factor of SSA/2
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
    INTEGER, INTENT(IN) :: NC, NLR, NMU, NAZ
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NMU, NAZ), INTENT(IN) :: A11, A12, R11, R12
    REAL*8, DIMENSION(NLR), INTENT(IN) :: SSA, RTO
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NMU, NAZ, NLR), INTENT(OUT) :: Z11, Z12
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
            Z11(:, :, ILR) = SSA05*R11
            Z12(:, :, ILR) = SSA05*R12
        END DO ! ILR = 1, NLR
    ELSE ! NC > 1: Rayleigh & Aerosol
        DO ILR = 1, NLR
            SSA05 = 0.5D0*SSA(ILR)
            WA = SSA05*RTO(ILR)
            WR = SSA05 - WA
            Z11(:, :, ILR) = WR*R11 + WA*A11
            Z12(:, :, ILR) = WR*R12 + WA*A12
        END DO ! IL = 1, NL
    END IF ! NC = 1
!
END SUBROUTINE Z0_IP
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 03Oct15 - First created using ZM0_IP
!===============================================================================