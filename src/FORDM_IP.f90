SUBROUTINE FORDM_IP(ZM01, ZM02, ZM03, NEWLR, E0, EUP, EDN, NL, NB, NLR, NG1,   &
                                                             NG2, I1B, Q1B, U1B)
!===============================================================================
! PURPOSE:
!   To compute the m-th Fourier moment of the 1st scattering order at Gauss
!   nodes and all levels (boundaries, IB = 1:NB). Path radiance is computed
!   here, surface contribution is computed elsewhere.
!
! INPUT:
!   ZM0i    D(NG2, NLR)   1st column of Fourier moment of the phase matrix
!   NEWLR   I(NL)         Profile of microlayers in each optical layer
!   E0      D(NL)         Average attenuation of the direct beam at two levels
!   EUP     D(NG1)        Bouguer attenuation for ascending beams at all levels
!   EDN     D(NG1)        Bouguer attenuation for descending beams through dTau
!   NL      I(1)          Number of microlayers, NL >= NLR
!   NB      I(1)          Number of boundaries at microlayers, NL+1
!   NLR     I(1)          Number of optical layers, NLR <= NL
!   NG1     I(1)          Number of Gauss nodes per hemisphere, NG2/2
!   NG2     I(1)          Number of Gauss nodes per sphere, NG1*2
!
! OUTPUT:
!   I1B, Q1B, U1B   D(NG2, NB)   Stokes vector at Gauss nodes at all levels
!
! TREE:
!   -
!
! COMMENTS:
!   The 1st scattering order is a source for higher scattering orders computed
!   in the SORD_IP subroutine.
!
!   Try using SSCAT1M (exact SS), splitted into UP & DN, instead of the French
!   technique [1] implemented below. For this, declare
!
!       REAL*8, DIMENSION(NG1) :: &
!       I1U, Q1U, U1U ! Single scattering by a single microlayer at Gauss nodes
!
!   and same for DN.
!
!   Following [1, p.497, Eq.(67)], the averaged intensity is
!
!       E0 = [exp(-i*dTau/mu0) + exp(-(i+1)*dTau/mu0)]/2
!
!   Note that exp(-dTau/cos(vza)) is missing from the equation by mistake (typo)
!
! REFERENCES:
!   1. Lenoble J et al, 2007: JQSRT, V107, P479
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NL, NB, NLR, NG1, NG2
!
! INPUT ARRAYS
    INTEGER, DIMENSION(NL), INTENT(IN) :: NEWLR
    REAL*8, DIMENSION(NL), INTENT(IN) :: E0
    REAL*8, DIMENSION(NG1), INTENT(IN) :: EUP, EDN
    REAL*8, DIMENSION(NG2, NLR), INTENT(IN):: ZM01, ZM02, ZM03
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NG2, NB), INTENT(OUT) :: I1B, Q1B, U1B
!
! LOCAL VARIABLES
    INTEGER &
        IB, & ! Loop index over boundaries at microlayers, 1:NB
        ILR   ! Loop index over (macro)layers
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NG1) :: &
        CEDN, & ! Complementary to EDN, 1 - EDN
        CEUP    ! Complementary to EUP, 1 - EUP
!===============================================================================
!
!   Precompute
    CEUP = 1.0D0 - EUP
    CEDN = 1.0D0 - EDN
!
!   Approximate 1st scattering at all boundaries and gauss nodes: *** UP ***
    I1B(:NG1, NB) = 0.0D0
    Q1B(:NG1, NB) = 0.0D0
    U1B(:NG1, NB) = 0.0D0
    DO IB = NB-1, 1, -1
        ILR = NEWLR(IB) ! Get index for optical layer
        I1B(:NG1, IB) = I1B(:NG1, IB+1)*EUP + ZM01(:NG1, ILR)*CEUP*E0(IB)
        Q1B(:NG1, IB) = Q1B(:NG1, IB+1)*EUP + ZM02(:NG1, ILR)*CEUP*E0(IB)
        U1B(:NG1, IB) = U1B(:NG1, IB+1)*EUP + ZM03(:NG1, ILR)*CEUP*E0(IB)
    END DO ! I = NB-1, 1, -1 
!
!   Approximate 1st scattering at all boundaries and gauss nodes: *** DOWN ***
    I1B(NG1+1:, 1) = 0.0D0
    Q1B(NG1+1:, 1) = 0.0D0
    U1B(NG1+1:, 1) = 0.0D0
    DO IB = 2, NB
        ILR = NEWLR(IB-1) ! Get index for optical layer
        I1B(NG1+1:, IB) = I1B(NG1+1:, IB-1)*EDN + &
                                                ZM01(NG1+1:, ILR)*CEDN*E0(IB-1)
        Q1B(NG1+1:, IB) = Q1B(NG1+1:, IB-1)*EDN + &
                                                ZM02(NG1+1:, ILR)*CEDN*E0(IB-1)
        U1B(NG1+1:, IB) = U1B(NG1+1:, IB-1)*EDN + &
                                                ZM03(NG1+1:, ILR)*CEDN*E0(IB-1)
    END DO ! I = 2, NB
!
END SUBROUTINE FORDM_IP
!===============================================================================
! 01May16 - NGX & IB1 are replaced with NG1+1 & IB+/-1, respectively
!
! 22Apr16 - Minor changes in comments
!
! 15Oct15 - Surface contribution is removed to a separate subroutine. Set of
!           input/output parameters is rearranged. Tested.
!
! 16Jul15 - First created from SORD_IP. Tested vs analytical SS for SZA=0, m=0,
!           32 single Gauss nodes, NL = 17 (NB = 18). Rayleigh over black. Ok.
!===============================================================================