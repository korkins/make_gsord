SUBROUTINE RAYLEIGH12(MU, N, DF, R11, R12)
!===============================================================================
! PURPOSE:
!   To compute elements [11] and [12] for the Rayleigh scattering matrix.
!   The subroutine accounts for the Rayleigh depolarization factor, DF.
!
! INPUT:
!   MU   D(N)   Cosine of the scattering angle
!   N    I(1)   Length of MU
!   DF   D(1)   Depolarization factor
!
! OUTPUT:
!   R11, R12   D(N)   Elements of the scattering matrix
!
! TREE:
!   -
!
! COMMENTS:
!   The Rayleigh scattering matrix for isotropic spherical particles is
!   [1, Eq.(2.14)]
!
!           | 1 + MU2  -1 + MU2    0    0|
!           |-1 + MU2   1 + MU2    0    0|
!   R = 3/4 |       0         0  2MU    0|                                   (1)
!           |       0         0    0  2MU|
!
!   where MU2 = MU*MU, 2MU = 2*MU. The elements R11 = R(1, 1), R12 = R(1, 2),
!   R33 = R(3, 3) are computed. The phase function, R11, is normalized
!
!   1/2 integral(R11, -1...+1) = 1                                           (2)
!
!   The anisotropy of realistic randomly oriented molecules is accounted for as
!   follows [1, Eq.(2.15)]
!
!              | 1 + MU2  -1 + MU2    0     0  |
!              |-1 + MU2   1 + MU2    0     0  |
!   Ra = D 3/4 |       0         0  2MU     0  | +
!              |       0         0    0  D' 2MU|
!
!            |1 0 0 0|
!            |0 0 0 0|
!   +(1 - D) |0 0 0 0|,                                                      (3)
!            |0 0 0 0|
!
!   where [1, Eq.(2.16)]
!
!   D = (1 - d)/(1 + d/2), D' =  (1 - 2d)/(1 - d),                           (4)
!
!   and d is the depolarization factor, DF. For isotropic molecules, d = 0.
!   Some values of d are [1, p.542]: d(H2) = 0.02, d(N2) = 0.03, d(Air) = d(N2),
!   d(O2) = 0.06, d(CO2) = 0.09. See the reference to Penndorf(1957) in [1] or
!   [3].
!
!   In [4] Eq.(12):
!
!   Ra = 3/(4(1 + 2G))[(1 + 3G) + (1 - G)MU2]                                (5)
!
!   where
!
!   G = d/(2 - d)                                                            (6)
!
!   Eq.(5) is simply another (but equal!) form for R11 from Eq.(3).
!
!   In [5, p.56]: the depolarization factor, d = DF, belongs to [0..6/7] and
!   the typical values are [5, p.43]: d(N2) = 0.020, d(O2) = 0.058,
!   d(dry air) = 0.028, d(CO2) = 0.079. Depolarization factor is [5, Eq.(2.94)]
!
!   d = Il(MU=0)/Ir(MU=0).                                                   (7)
!
!   Note that:
!   1. The range of depolarization factor, DF, is [0...0.5] [1, p.542]
!      or [0..7/6 = 0.86] [5, p.56]
!   2. If DF > 0 then R33 <> R44 and R11 <> R22
!   3. DF is often assumed spectrally independent [2, p.3494]
!
! REFERENCES:
!   1. Hansen JE, Travis LD, 1974, Space.Sci.Rev., V16, pp.527-610
!   2. Bates DR, 1984, Planet. Space Sci., V32,pp.785–790
!   3. Caldas M, Semiano V, 2001, JOSA, V18, N4, pp.831-838
!   4. Bucholtz A, 1995, V34, N15, pp.2765-2773
!   5. HovenierJW et al. Transfer of polarized light, Kluwer 2004
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: N
    REAL*8, INTENT(IN) :: DF
!
! INPUT ARRAYS
    REAL*8, DIMENSION(N), INTENT(IN) :: MU
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(N), INTENT(OUT) :: R11, R12
!
! LOCAL VARIABLES
    REAL*8 &
        D,  & ! D  from Eq.(4)
        D0, & ! D*D' = (1 - 2d)/(1 + d/2)
        D1    ! D1 = 1 + d/2
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(N) :: &
        MU2 ! MU*MU
!===============================================================================
!
    MU2 = MU*MU
!
!   Eq.(1)
    R11 =  1.0D0 + MU2
    R12 = -1.0D0 + MU2
!
!   Eq.(2)
    R11 = 0.75D0*R11
    R12 = 0.75D0*R12
!
!   Anisotropic molecules
    IF (DF > 0.0D0) THEN
!       Depolarization constants
        D1 = 1.0D0 + 0.5D0*DF
        D = (1.0D0 - DF) / D1
        D0 = (1.0D0 - 2.0D0*DF) / D1
!       Recompute the elements
        R11 = D*(R11 - 1.0D0) + 1.0D0
        R12 = D*R12
    END IF ! DF > 0.0D0
!
END SUBROUTINE RAYLEIGH12
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 02Oct14 - Declaration of local scalars has slightly changed according to the
!           template from code SORD. Parameter O75 = 0.75D0 is removed.
!
! 14Oct13 - Created from the RAYLEIGH subroutine.
!===============================================================================