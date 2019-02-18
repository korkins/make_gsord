SUBROUTINE SURFACE0M_IP(IM, ISRF, R01, R02, R03, WCMA, WSMA, NGA, NR, NM, &
                        RM01, RM02, RM03)
!===============================================================================
! PURPOSE:
!   To compute the m-th Fourier expansion moment for a surface directly
!   irradiated by the Sun (single mu0).
!
! INPUT:
!   IM     I(1)         Fourier moment, m = 0, 1, 2 ...
!   ISRF   I(1)         Surface index, ISRF > 1
!   R01    D(NGA, NR)   BRDF at Gauss azimuth nodes
!   R02    D(NGA, NR)   Similar to R01 but for the 21 element (ISRF > 5)
!   R03    D(NGA, NR)   Similar to R01 but for the 31 element (ISRF > 5)
!   WCMA   D(NGA, NM)   Weighted cosines of Gauss azimuth nodes for all M
!   WCMA   D(NGA, NM)   Similar to WCMA except for sin-terms (ISRF > 5)
!   NGA    I(1)         Number of Gauss nodes
!   NR     I(1)         Number of reflected zeniths
!   NM     I(1)         Maximum Fourier order, m = 0, 1, 2 ... NM
!
! OUTPUT:
!   RM01   D(NR)   IM-moment of the BRDF at all zeniths of observation
!   RM02   D(NR)   Similar to RM01 but for the 21 element (ISRF > 5)
!   RM03   D(NR)   Similar to RM01 but for the 31 element (ISRF > 5)
!
! TREE:
!   -
!
! COMMENTS:
!
!   1/pi(surf.int.)*2pi(SORD norm.) = 2 - explanation of the factor of 2.0 below
!
!   The azimuth expansion of scalar surface is
!
!   SRF  = SUM[ (2 - d0m)*SRFM*cos(m*AZA), m = 0...M ].                      (1)
!
!   Eq.(1) determines the Fourier expansion for each element of SRFM as
!
!   SRFM = 1/(2*PI)*int(SRF(x)*cos(m*x)dx, 0..2*PI) ~=
!        ~= 1/(2*PI)*SUM[wj*SRF(2PI*xj)*cos(m*2PI*xj), j=1..NGA].            (2)
!
!   NOTE: in the code SHARM [1], file BCondition.cpp, 1/2PI constant is
!   omitted in Eq.(2). Here 1/2PI is included via Gauss knots as follows
!
!   1/2pi*int(f(fi)cos(m*fi)dfi, fi = 0..2pi) =
!               = int(f(2pi*z)cos(m*2pi*z)dz, z = 0..1)                      (3)
!
!   GAUSZW(0.0D0, 1.0D0, NAZ, Z, W) implements the Double Gauss scheme.
!   Numerical integration over azimuth is
!
!   I = int(F(x)dx, 0..2PI) ~= sum(wj*F(2PI*xj), j=1..N),                    (4)
!
!   where xj and wj are obtained using the Double Gauss scheme, xj = [0..1] and
!   2PI*xj = [0..2PI].
!
! REFERENCES:
!   1. Lyapustin A.I., Appl.Opt(2005), 44, N36, pp.7764-7772.
!===============================================================================
!
    IMPLICIT NONE
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: IM, ISRF, NGA, NR, NM
!
! INPUT ARRAYS
    REAL*8, DIMENSION(NGA, NR), INTENT(IN) :: R01, R02, R03
    REAL*8, DIMENSION(NGA, NM), INTENT(IN) :: WCMA, WSMA
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NR), INTENT(OUT) :: RM01, RM02, RM03
!
! LOCAL VARIABLES
    INTEGER &
        IR ! Loop index over reflected zeniths
!===============================================================================
!
    IF (ISRF < 6) THEN ! Reflection of I only
        DO IR = 1, NR
            RM01(IR) = 2.0D0*SUM(R01(:, IR)*WCMA(:, IM))
        END DO ! IR = 1, NR
    ELSE ! ISRF > 5: reflection of I & P
        IF (IM == 1) THEN ! No U-reflection
            DO IR = 1, NR
                RM01(IR) = 2.0D0*SUM(R01(:, IR)*WCMA(:, IM))
                RM02(IR) = 2.0D0*SUM(R02(:, IR)*WCMA(:, IM))
                RM03(IR) = 0.0D0
            END DO ! IR = 1, NR
        ELSE ! IM > 1 (m > 0)
            DO IR = 1, NR
                RM01(IR) = 2.0D0*SUM(R01(:, IR)*WCMA(:, IM))
                RM02(IR) = 2.0D0*SUM(R02(:, IR)*WCMA(:, IM))
                RM03(IR) = 2.0D0*SUM(R03(:, IR)*WSMA(:, IM))
            END DO ! IR = 1, NR
        END IF ! IM = 1
    END IF ! ISRF < 6
!
END SUBROUTINE SURFACE0M_IP
!===============================================================================
! 22Apr16 - Minor changes in comments
!
! 01May15 - Vector surface, ISRF > 5, is added and tested for NT ocean.
!
! 24Apr15 - First created.
!===============================================================================