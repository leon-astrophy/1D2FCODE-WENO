!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine did TVD reconstructions on edge state !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RECONSTRUCT_TVD
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
		CALL TVD (p1(j-1), p1(j), p1(j+1), p1R(j-1), p1L(j))
		CALL TVD (rho1(j-1), rho1(j), rho1(j+1), rho1R(j-1), rho1L(j))
		CALL TVD (vel1(j-1), vel1(j), vel1(j+1), vel1R(j-1), vel1L(j))
	END DO
	
	! For DM scalars !
	If(iminsca1 > 0) THEN
		DO i = iminsca1, imaxsca1
			DO j = -1, length_step_part_1 + 1
				CALL TVD (sca1(i,j-1), sca1(i,j), sca1(i,j+1), sca1R(i,j-1), sca1L(i,j))
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
	CALL TVD (p2(j-1), p2(j), p2(j+1), p2R(j-1), p2L(j))
	CALL TVD (rho2(j-1), rho2(j), rho2(j+1), rho2R(j-1), rho2L(j))
	CALL TVD (vel2(j-1), vel2(j), vel2(j+1), vel2R(j-1), vel2L(j))
END DO

! Do extra reconstuctions for dual energy !
IF (dual_energy == 1) THEN
	DO j = -1, length_step_part_2 + 1
		CALL TVD (rhoe2(j-1), rhoe2(j), rhoe2(j+1), rhoe2R(j-1), rhoe2L(j))
	END DO
END IF

! For NM scalars !
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO j = -1, length_step_part_2 + 1
			CALL TVD (sca2(i,j-1), sca2(i,j), sca2(i,j+1), sca2R(i,j-1), sca2L(i,j))
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
		CALL TVD (j, cs2(j-2), cs2(j-1), cs2(j), cs2(j+1), cs2(j+2), cs2R(j-1), cs2L(j))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine did TVD limiter on edge state !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TVD (vm1, v, vp1, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The numerical value of u at each grid center !
REAL (DP), INTENT (IN) :: vm1, v, vp1

! The numerical value of u at each grid boundary !
REAL (DP), INTENT (OUT) :: vm_out, vp_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Limited Slope !
REAL (DP) :: slope

! Get the limited slope !
slope = limiter(vp1 - v, v - vm1)

! interpolate !
vm_out = v - 0.5D0*slope
vp_out = v + 0.5D0*slope

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	real(DP) function limiter(alpha, beta)
	implicit none
	real(DP) :: alpha, beta, limier
	
	! Van Leer (van Leer 1974) !
	!limiter = (sign(1.0D0,alpha) + sign(1.0D0,beta))*(alpha*beta)/(abs(alpha) + abs(beta))
	
	! Monotonized Central (MC) (van Leer 1977) !
	limiter = 0.5D0*(sign(1.0D0,alpha) + sign(1.0D0,beta))*min(0.5D0*abs(alpha + beta), 2.0D0*abs(alpha), 2.0D0*abs(beta))

	! Superbee (Roe 1986) !
	!limiter = 0.5D0*(sign(1.0D0,alpha) + sign(1.0D0,beta))*max*(min(2.0D0*abs(alpha),abs(beta)),min(abs(alpha),2.0D0*abs(beta)))
	
	! Minmod (Roe 1986) !
	!limiter = 0.5D0*(sign(1.0D0,alpha) + sign(1.0D0,beta))*min(abs(alpha),abs(beta))

	end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE