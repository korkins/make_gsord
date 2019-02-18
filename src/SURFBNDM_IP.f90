SUBROUTINE SURFBNDM_IP(MU0E0, ISRF, RM0, I1UP, Q1UP, U1UP, EUP, NL, NB, NG1)
!===============================================================================
! PURPOSE:
!   To account for surface contribution to the m-th Fourier moment of the 1st
!   order of path radiance at Gauss nodes and all boundaries.
!
! INPUT:
!   MU0E0   D(1)          Attenuation of direct beam at surface scaled by mu0
!   ISRF    I(1)          Index of the surface model, ISRF = 1, 2, ...
!   RM0     D(NG1, 3)     1st column of Fourier moment of the surface matrix
!   EUP     D(NG1, NL)    Bouguer attenuation for ascending beams at all levels
!   NL      I(1)          Number of microlayers, NL >= NLR
!   NB      I(1)          Number of boundaries at microlayers, NL+1
!   NG1     I(1)          Number of Gauss nodes per hemisphere, NG2/2
!
! INOUT:
!   I1UP,   D(NG1, NB)    Ascending path radiance at Gauss nodes at all levels
!   Q1UP,
!   U1UP
!
! TREE:
!   -
!
! COMMENTS:
!   -
!
! REFERENCES:
!   -
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: ISRF, NL, NB, NG1
    REAL*8, INTENT(IN) :: MU0E0
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NG1, NL), INTENT(IN) :: EUP
    REAL*8, DIMENSION(NG1, 3), INTENT(IN) :: RM0
!
! DUAL INTENT ARRAYS
    REAL*8, DIMENSION(NG1, NB), INTENT(INOUT) :: I1UP, Q1UP, U1UP
!
! LOCAL VARIABLES
    INTEGER &
        IB ! Loop index over boundaries at microlayers, 1:NB
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NG1) :: &
        I1R, Q1R, U1R ! Exp-attenuated surface reflection of the direct beam
!===============================================================================
!
!   Surface reflection of the direct beam
    I1R = MU0E0*RM0(:, 1)
    IF (ISRF > 5) THEN ! Polarized reflectance
        Q1R = MU0E0*RM0(:, 2)
        U1R = MU0E0*RM0(:, 3)
    END IF ! ISRF > 5
!
!   Surface contribution at all boundaries
    I1UP(:, NB) = I1R
    DO IB = 1, NL ! NL=NB-1
        I1UP(:, IB) = I1UP(:, IB) + I1R*EUP(:, NB-IB)
    END DO ! IB = 1, NL
!   Polarized components of the vector surface
    IF (ISRF > 5) THEN ! Polarized reflectance
        Q1UP(:, NB) = Q1R
        U1UP(:, NB) = U1R
        DO IB = 1, NL
            Q1UP(:, IB) = Q1UP(:, IB) + Q1R*EUP(:, NB-IB)
            U1UP(:, IB) = U1UP(:, IB) + U1R*EUP(:, NB-IB)
        END DO ! IB = 1, NL
    END IF ! ISRF > 5
!
END SUBROUTINE SURFBNDM_IP
!===============================================================================
! 01May16 - ID is rpelaced with NB-IB
!
! 22Apr16 - Minor changes in comments
!
! 15Oct15 - First created and tested.
!===============================================================================