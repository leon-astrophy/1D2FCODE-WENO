!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the WENO module that contains all necessary tool     !
! for WENO reconstruction of the numerical flux at the left    !
! and the right hand side of the boundary cell, asuming a flux !
! splitting procedure had been performed		       !			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE WENO_MODULE        
USE DEFINITION
IMPLICIT NONE

! The small parameter in WENO scheme that avoid the coefficient to be dividing by zero !
! (Do NOT change this unless you know what you are doing) !
REAL (DP), PARAMETER :: smallpara = 1.0E-40_DP

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do reconstruction using WENO interpolation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RECONSTRUCT_WENO
USE RIEMANN_MODULE
USE FLAME_MODULE
USE NUCLEAR_MODULE
USE DEFINITION  
IMPLICIT NONE

! integer parameter !
INTEGER :: i, j, in

! Dummy !
REAL (DP) :: left, right, dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do reconstruction steps for DM !
IF(RUNDM_flag == 1) THEN
	
	! Reconstruct hydrodynamics variables at cell boundaries for DM !
	DO j = -1, length_step_part_1 + 1
		CALL WENO (j, p1(j-2), p1(j-1), p1(j), p1(j+1), p1(j+2), p1R(j-1), p1L(j))
		CALL WENO (j, rho1(j-2), rho1(j-1), rho1(j), rho1(j+1), rho1(j+2), rho1R(j-1), rho1L(j))
		CALL WENO (j, vel1(j-2), vel1(j-1), vel1(j), vel1(j+1), vel1(j+2), vel1R(j-1), vel1L(j))
	END DO
	
	! For DM scalars !
	If(iminsca1 > 0) THEN
		DO i = iminsca1, imaxsca1
			DO j = -1, length_step_part_1 + 1
				CALL WENO (j, sca1(i,j-2), sca1(i,j-1), sca1(i,j), sca1(i,j+1), sca1(i,j+2), sca1R(i,j-1), sca1L(i,j))
			END DO
		END DO
	END IF

	! Get epsilon at boundary !
	DO j = -1, length_step_part_1 + 1
		CALL EOSEPSILON(rho1R(j), p1R(j), eps1R(j), 1)
		CALL EOSEPSILON(rho1L(j), p1L(j), eps1L(j), 1)
	END DO

	! Get speed of sound at boundary !
	DO j = -1, length_step_part_1 + 1
		CALL EOSSOUNDSPEED(p1R(j), rho1R(j), eps1R(j), cs1R(j), 1)
		CALL EOSSOUNDSPEED(p1L(j), rho1L(j), eps1L(j), cs1L(j), 1)
	END DO

	If(movinggriddm_flag == 1) THEN
		DO j = -1, length_step_part_1 + 1
			vf1R(j) = vel1_max*r1F(j)/radius1
			vf1L(j) = vel1_max*r1F(j)/radius1
		END DO
	END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reconstruct hydrodynamics variables at cell boundaries for NM !
DO j = -1, length_step_part_2 + 1
	CALL WENO (j, p2(j-2), p2(j-1), p2(j), p2(j+1), p2(j+2), p2R(j-1), p2L(j))
	CALL WENO (j, rho2(j-2), rho2(j-1), rho2(j), rho2(j+1), rho2(j+2), rho2R(j-1), rho2L(j))
	CALL WENO (j, vel2(j-2), vel2(j-1), vel2(j), vel2(j+1), vel2(j+2), vel2R(j-1), vel2L(j))
END DO

! Do extra reconstuctions for dual energy !
IF (dual_energy == 1) THEN
	DO j = -1, length_step_part_2 + 1
		CALL WENO (j, rhoe2(j-2), rhoe2(j-1), rhoe2(j), rhoe2(j+1), rhoe2(j+2), rhoe2R(j-1), rhoe2L(j))
	END DO
END IF

! For NM scalars !
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO j = -1, length_step_part_2 + 1
			CALL WENO (j, sca2(i,j-2), sca2(i,j-1), sca2(i,j), sca2(i,j+1), sca2(i,j+2), sca2R(i,j-1), sca2L(i,j))
		END DO
	END DO
END IF

! Get epsilon at boundary !
IF (dual_energy == 1) THEN
	DO j = -1, length_step_part_2 + 1
		eps2R(j) = rhoe2R(j)/rho2R(j)
		eps2L(j) = rhoe2L(j)/rho2L(j)
	END DO
ELSE
	DO j = -1, length_step_part_2 + 1
		CALL EOSEPSILON(rho2R(j), p2R(j), eps2R(j), 2)
		CALL EOSEPSILON(rho2L(j), p2L(j), eps2L(j), 2)
	END DO
END IF

! Get speed of sound at boundary !
IF(helmeos_flag == 1) THEN
	DO j = -1, length_step_part_2 + 1
		CALL WENO (j, cs2(j-2), cs2(j-1), cs2(j), cs2(j+1), cs2(j+2), cs2R(j-1), cs2L(j))
	END DO
ELSE
	DO j = -1, length_step_part_2 + 1
		CALL EOSSOUNDSPEED(p2R(j), rho2R(j), eps2R(j), cs2R(j), 2)
		CALL EOSSOUNDSPEED(p2L(j), rho2L(j), eps2L(j), cs2L(j), 2)
	END DO
END IF

! No reconstruction for frame velocity since we assume analytic continous form !
If(movinggridnm_flag == 1) THEN
	DO j = -1, length_step_part_2 + 1
		vf2R(j) = vel2_max*r2F(j)/radius2
		vf2L(j) = vel2_max*r2F(j)/radius2
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Coefficient crj that appears nautrally in the newtonian !
! interpolation using 3 point stencil to get a numerical  !
! flux at either the left or the right hand side of the	  !
! boundary cell. Must not change this			  !			 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CONSTC (c)
USE DEFINITION
IMPLICIT NONE

REAL (8), INTENT (OUT), DIMENSION (-1 : 2, 0 : 2) :: c

c (-1, 0) = 1.1E1_DP / 6.0E0_DP
c (-1, 1) = - 7.0E0_DP / 6.0E0_DP
c (-1, 2) = 1.0E0_DP / 3.0E0_DP

c (0, 0) = 1.0E0_DP / 3.0E0_DP
c (0, 1) = 5.0E0_DP / 6.0E0_DP
c (0, 2) = - 1.0E0_DP / 6.0E0_DP

c (1, 0) = - 1.0E0_DP / 6.0E0_DP
c (1, 1) = 5.0E0_DP / 6.0E0_DP
c (1, 2) = 1.0E0_DP / 3.0E0_DP

c (2, 0) = 1.0E0_DP / 3.0E0_DP
c (2, 1) = - 7.0E0_DP / 6.0E0_DP
c (2, 2) = 1.1E1_DP / 6.0E0_DP

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the coefficient dr that used to determine the convex !
! combination coefficient omega. This should not be change     !
! in case you provide another better WENO schem	  	       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CONSTD (d)
USE DEFINITION
IMPLICIT NONE

REAL (8), INTENT (OUT), DIMENSION (0 : 2) :: d

d (0) = 3.0E0_DP / 1.0E1_DP
d (1) = 3.0E0_DP / 5.0E0_DP
d (2) = 1.0E0_DP / 1.0E1_DP

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the WENO scheme for reconstructing the numerical flux at both the !
! left and right hand side located at the boundary cell. In this version, I !
! provide different WENO scheme that differ by their smoothness indicator   !
! I also include WENO scheme that use combination of high and low order     !
! polynominal as building block. Nonetheless, a monotonicity preserving     !
! limter option is provided so to make the solution to be MPW               !
! For details, please refer to the textbook with ISBN 3-540-65893-9         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WENO (i, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! Input integer !
INTEGER, INTENT(IN) :: i

! The input into the subroutine, including conservative variable and input flux function !
REAL (DP), INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL (DP), INTENT (OUT) :: vm_out, vp_out

! The coefficient c for langrange interpolation !
REAL (DP), DIMENSION (-1 : 2, 0 : 2) :: c

! For weight reconstruction !
REAL (DP), DIMENSION (0 : 2) :: d, td

! Temporal arrays !
REAL (DP), DIMENSION (0 : 2) :: vrhs, vlhs

! For assigning weights !
REAL (DP), DIMENSION (0 : 2) :: alpha, talpha, omega, tomega, beta

! Temporal arrays !
REAL (DP), DIMENSION (i - 2 : i + 2) :: v

! Integer !
INTEGER :: j, r, s

! Tempeorary parameter !
REAL (DP) :: tau5, temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We get the coefficient c !
CALL CONSTC (c)

! We get the coefficient d !
CALL CONSTD (d)

! Relating d with d-dash that used to calculate u at left and right boundary !
DO r = 0, 2
	td (r) = d (2 - r)
END DO

! Assign temporal arrays !
v(i - 2) = vm2
v(i - 1) = vm1
v(i) = vc
v(i + 1) = vp1
v(i + 2) = vp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! We calculate the value of u at each grid by the following loop !
! Do the right cell boundary !
DO r = 0, 2
	vrhs (r) = 0.0E0_DP
		
	! We calculate the value of u at right boundary !
	DO j = 0, 2
		vrhs (r) = vrhs (r) + c (r, j) * v (i - r + j)
	END DO
END DO

! Do the left cell boundary !
DO r = 0, 2
	vlhs (r) = 0.0E0_DP

	! Do the same for left boundary !
	DO j = 0, 2
		vlhs (r) = vlhs (r) + c (r - 1, j) * v (i - r + j)
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! These are essential parameter for further construction of u !
beta (0) = (1.3E1_DP / 1.2E1_DP) * (v (i) - 2 * v (i + 1) + v (i + 2)) ** 2 &
		+ (1.0E0_DP / 4.0E0_DP) * (3 * v (i) - 4 * v (i + 1) + v (i + 2)) ** 2
beta (1) = (1.3E1_DP / 1.2E1_DP) * (v (i - 1) - 2 * v (i) + v (i + 1)) ** 2 &
		+ (1.0E0_DP / 4.0E0_DP) * (v (i - 1) - v (i + 1)) ** 2
beta (2) = (1.3E1_DP / 1.2E1_DP) * (v (i - 2) - 2 * v (i - 1) + v (i)) ** 2 &
		+ (1.0E0_DP / 4.0E0_DP) * (v (i - 2) - 4 * v (i - 1) + 3 * v (i)) ** 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assigning tau for the WENO-Z corrections !
tau5 = abs (beta(0) -beta(2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do WENO-Z weight reconstructions !
DO r = 0, 2
	alpha (r) = d (r) * (1.0D0 + (tau5/(beta(r) + smallpara))**2)
END DO

temp = 0.0E0_DP
	
! The denominator in finding omega, a coefficient for the last step of reconstruction  !
DO s = 0, 2
	temp = temp + alpha (s)
END DO
	
! Find the omega !
DO r = 0, 2
	omega (r) = alpha (r) / temp
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the alpha, omega for the value of u at left grid boundary... !
DO r = 0, 2
	talpha (r) = td (r) * (1.0D0 + (tau5/(beta(r) + smallpara))**2)
END DO

temp = 0.0E0_DP

DO s = 0, 2
	temp = temp + talpha (s)
END DO

DO r = 0, 2
	tomega (r) = talpha (r) / temp
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

vp_out = 0.0E0_DP
	
! u at the left boundary !
DO r = 0, 2

	! Original WENO !
	vp_out = vp_out + omega (r) * vrhs (r)

END DO

vm_out = 0.0E0_DP
	
! u at the right boundary !
DO r = 0, 2

	! Original WENO !
	vm_out = vm_out + tomega (r) * vlhs (r)	

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

end module