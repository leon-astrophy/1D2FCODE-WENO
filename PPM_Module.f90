!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the PPM (Piecewise Parabolic Method) Module that reconstruct   !
! interface values of primitive variables. We can choose either extremum !
! preserving limiter or the original limiter. We use shock flattening	 !
! and contact steepening algorithm. For details, please refer to the	 !
! original PPM paper and the method paper by FLASH code. Also see do	 !
! refer to the PPM reconstruction procedure by GR1D code		 !     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE PPM_MODULE
USE DEFINITION
IMPLICIT NONE

! Flattening coefficients ! 
REAL (DP), DIMENSION(-4 : length_step_1 + 5) :: flat_1
REAL (DP), DIMENSION(-4 : length_step_2 + 5) :: flat_2

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine did the reconsturction steps of PPM !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RECONSTRUCT_PPM 
USE RIEMANN_MODULE
USE FLAME_MODULE
USE NUCLEAR_MODULE
USE DEFINITION 
IMPLICIT NONE

! integer parameter !
INTEGER :: i, j

! Dummy !
REAL (DP) :: left, right

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the same for DM !
IF(RUNDM_flag == 1) THEN
	DO j = -1, length_step_part_1 + 1
		CALL PPM_6THINTERPOLATE (j, rho1(j-3), rho1(j-2), rho1(j-1), rho1(j), rho1(j+1), rho1(j+2), rho1(j+3), rho1R(j-1), rho1L(j))
	END DO	
END IF

! We first interpolate density for NM !
DO j = -1, length_step_part_2 + 1
	CALL PPM_6THINTERPOLATE (j, rho2(j-3), rho2(j-2), rho2(j-1), rho2(j), rho2(j+1), rho2(j+2), rho2(j+3), rho2R(j-1), rho2L(j))
END DO	

! Do steepening on density profiles !
CALL PPM_STEEPEN

! Get flattening coefficients !
CALL GET_FLATTEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the remaining interpolation steps for DM !
IF(RUNDM_flag == 1) THEN
	DO j = -1, length_step_part_1 + 1
		CALL PPM_6THINTERPOLATE (j, p1(j-3), p1(j-2), p1(j-1), p1(j), p1(j+1), p1(j+2), p1(j+3), p1R(j-1), p1L(j))
		CALL PPM_6THINTERPOLATE (j, vel1(j-3), vel1(j-2), vel1(j-1), vel1(j), vel1(j+1), vel1(j+2), vel1(j+3), vel1R(j-1), vel1L(j))
	END DO
	If(iminsca1 > 0) THEN
		DO i = iminsca1, imaxsca1
			DO j = -1, length_step_part_1 + 1
				CALL PPM_6THINTERPOLATE (j, sca1(i,j-3), sca1(i,j-2), sca1(i,j-1), sca1(i,j), sca1(i,j+1), sca1(i,j+2), sca1(i,j+3), sca1R(i,j-1), sca1L(i,j))
			END DO
		END DO
	END IF
END IF

! Reconstruct remaining hydrodynamics variables at cell boundaries for NM !
DO j = -1, length_step_part_2 + 1
	CALL PPM_6THINTERPOLATE (j, p2(j-3), p2(j-2), p2(j-1), p2(j), p2(j+1), p2(j+2), p2(j+3), p2R(j-1), p2L(j))
	CALL PPM_6THINTERPOLATE (j, vel2(j-3), vel2(j-2), vel2(j-1), vel2(j), vel2(j+1), vel2(j+2), vel2(j+3), vel2R(j-1), vel2L(j))
END DO	

! Do extra reconstuctions for dual energy !
IF (dual_energy == 1) THEN
	DO j = -1, length_step_part_2 + 1
		CALL PPM_6THINTERPOLATE (j, rhoe2(j-3), rhoe2(j-2), rhoe2(j-1), rhoe2(j), rhoe2(j+1), rhoe2(j+2), rhoe2(j+3), rhoe2R(j-1), rhoe2L(j))
	END DO
END IF

! For NM scalars !
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO j = -1, length_step_part_2 + 1
			CALL PPM_6THINTERPOLATE (j, sca2(i,j-3), sca2(i,j-2), sca2(i,j-1), sca2(i,j), sca2(i,j+1), sca2(i,j+2), sca2(i,j+3), sca2R(i,j-1), sca2L(i,j))
		END DO
	END DO
END IF

! Speed of sound !
IF(helmeos_flag == 1) THEN
	DO j = -1, length_step_part_2 + 1
		CALL PPM_6THINTERPOLATE (j, cs2(j-3), cs2(j-2), cs2(j-1), cs2(j), cs2(j+1), cs2(j+2), cs2(j+3), cs2R(j-1), cs2L(j))
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Flatten the profile for DM !
IF(RUNDM_flag == 1) THEN
	DO j = -1, length_step_part_1 + 1
		rho1R(j-1) = flat_1(j)*rho1(j) + (1.0D0 - flat_1(j))*rho1R(j-1)
		rho1L(j) = flat_1(j)*rho1(j) + (1.0D0 - flat_1(j))*rho1L(j)
		p1R(j-1) = flat_1(j)*p1(j) + (1.0D0 - flat_1(j))*p1R(j-1)
		p1L(j) = flat_1(j)*p1(j) + (1.0D0 - flat_1(j))*p1L(j)
		vel1R(j-1) = flat_1(j)*vel1(j) + (1.0D0 - flat_1(j))*vel1R(j-1)
		vel1L(j) = flat_1(j)*vel1(j) + (1.0D0 - flat_1(j))*vel1L(j)
	END DO
	If(iminsca1 > 0) THEN
		DO i = iminsca1, imaxsca1
			DO j = -1, length_step_part_1 + 1
				sca1R(i,j-1) = flat_1(j)*sca1(i,j) + (1.0D0 - flat_1(j))*sca1R(i,j-1)
				sca1L(i,j) = flat_1(j)*sca1(i,j) + (1.0D0 - flat_1(j))*sca1L(i,j)
			END DO
		END DO
	END IF
END IF

! Do the same for NM !
DO j = -1, length_step_part_2 + 1
	rho2R(j-1) = flat_2(j)*rho2(j) + (1.0D0 - flat_2(j))*rho2R(j-1)
	rho2L(j) = flat_2(j)*rho2(j) + (1.0D0 - flat_2(j))*rho2L(j)
	p2R(j-1) = flat_2(j)*p2(j) + (1.0D0 - flat_2(j))*p2R(j-1)
	p2L(j) = flat_2(j)*p2(j) + (1.0D0 - flat_2(j))*p2L(j)
	vel2R(j-1) = flat_2(j)*vel2(j) + (1.0D0 - flat_2(j))*vel2R(j-1)
	vel2L(j) = flat_2(j)*vel2(j) + (1.0D0 - flat_2(j))*vel2L(j)
END DO
IF (dual_energy == 1) THEN
	DO j = -1, length_step_part_2 + 1
		rhoe2R(j-1) = flat_2(j)*rhoe2(j) + (1.0D0 - flat_2(j))*rhoe2R(j-1)
		rhoe2L(j) = flat_2(j)*rhoe2(j) + (1.0D0 - flat_2(j))*rhoe2L(j)
	END DO
END IF
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO j = -1, length_step_part_2 + 1
			sca2R(i,j-1) = flat_2(j)*sca2(i,j) + (1.0D0 - flat_2(j))*sca2R(i,j-1)
			sca2L(i,j) = flat_2(j)*sca2(i,j) + (1.0D0 - flat_2(j))*sca2L(i,j)
		END DO
	END DO
END IF
IF(helmeos_flag == 1) THEN
	DO j = -1, length_step_part_2 + 1
		cs2R(j-1) = flat_2(j)*cs2(j) + (1.0D0 - flat_2(j))*cs2R(j-1)
		cs2L(j) = flat_2(j)*cs2(j) + (1.0D0 - flat_2(j))*cs2L(j)
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Apply monotonicty constraints !
If(RUNDM_flag == 1) THEN
	DO j = -1, length_step_part_1 + 1
		CALL XPPM_MONOTONIZE (p1(j-2), p1(j-1), p1(j), p1(j+1), p1(j+2), p1R(j-1), p1L(j))
		CALL XPPM_MONOTONIZE (rho1(j-2), rho1(j-1), rho1(j), rho1(j+1), rho1(j+2), rho1R(j-1), rho1L(j))
		CALL XPPM_MONOTONIZE (vel1(j-2), vel1(j-1), vel1(j), vel1(j+1), vel1(j+2), vel1R(j-1), vel1L(j))
	END DO
	If(iminsca1 > 0) THEN
		DO i = iminsca1, imaxsca1
			DO j = -1, length_step_part_1 + 1
				CALL XPPM_MONOTONIZE (sca1(i,j-2), sca1(i,j-1), sca1(i,j), sca1(i,j+1), sca1(i,j+2), sca1R(i,j-1), sca1L(i,j))
			END DO
		END DO
	END IF
END IF

! Apply monotonicty constraints !
DO j = -1, length_step_part_2 + 1
	CALL XPPM_MONOTONIZE (p2(j-2), p2(j-1), p2(j), p2(j+1), p2(j+2), p2R(j-1), p2L(j))
	CALL XPPM_MONOTONIZE (rho2(j-2), rho2(j-1), rho2(j), rho2(j+1), rho2(j+2), rho2R(j-1), rho2L(j))
	CALL XPPM_MONOTONIZE (vel2(j-2), vel2(j-1), vel2(j), vel2(j+1), vel2(j+2), vel2R(j-1), vel2L(j))
END DO

! Do extra for dual energy !
IF (dual_energy == 1) THEN
	DO j = -1, length_step_part_2 + 1
		CALL XPPM_MONOTONIZE (rhoe2(j-2), rhoe2(j-1), rhoe2(j), rhoe2(j+1), rhoe2(j+2), rhoe2R(j-1), rhoe2L(j))
	END DO
END IF

! For NM scalars !
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO j = -1, length_step_part_2 + 1
			CALL XPPM_MONOTONIZE (sca2(i,j-2), sca2(i,j-1), sca2(i,j), sca2(i,j+1), sca2(i,j+2), sca2R(i,j-1), sca2L(i,j))
		END DO
	END DO
END IF

! Speed of sound !
IF(helmeos_flag == 1) THEN
	DO j = -1, length_step_part_2 + 1
		CALL XPPM_MONOTONIZE (cs2(j-2), cs2(j-1), cs2(j), cs2(j+1), cs2(j+2), cs2R(j-1), cs2L(j))
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the remaining !
IF(RUNDM_flag == 1) THEN

	! Get the epsilon at boundary !
	DO j = -1, length_step_part_1 + 1
		CALL EOSEPSILON (rho1R(j), p1R(j), eps1R(j), 1)
		CALL EOSEPSILON (rho1L(j), p1L(j), eps1L(j), 1)
	END DO

	! Get speed of sound at boundary !
	DO j = -1, length_step_part_1 + 1
		CALL EOSSOUNDSPEED (p1R(j), rho1R(j), eps1R(j), cs1R(j), 1)
		CALL EOSSOUNDSPEED (p1L(j), rho1L(j), eps1L(j), cs1L(j), 1)
	END DO

	If(movinggriddm_flag == 1) THEN
		DO j = -1, length_step_part_1 + 1
			vf1R(j) = vel1_max*r1F(j)/radius1
			vf1L(j) = vel1_max*r1F(j)/radius1
		END DO
	END IF

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
IF(helmeos_flag /= 1) THEN
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine did the 4th order interpolation to get the interface values !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_INTERPOLATE (i, vm2, vm1, v, vp1, vp2, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER, INTENT(IN) :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The input into the subroutine, including conservative variable and input flux function !
REAL (DP), INTENT (IN) :: vm2, vm1, v, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL (DP), INTENT (OUT) :: vp_out, vm_out

! Interpolation coefficients ! 
REAL (DP) :: wm1, wc, wp1, wp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! PPM coefficients !
wm1 = -1.0D0/1.2D1
wc = 7.0D0/1.2D1
wp1 = 7.0D0/1.2D1
wp2 = -1.0D0/1.2D1

! Interpolate and get the right and left states !
vp_out = wm1*vm1 + wc*v + wp1*vp1 + wp2*vp2
vm_out = wm1*vm2 + wc*vm1 + wp1*v + wp2*vp1

END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine did the 6th order interpolation to get the interface values !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_6THINTERPOLATE (i, vm3, vm2, vm1, v, vp1, vp2, vp3, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER, INTENT(IN) :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The input into the subroutine, including conservative variable and input flux function !
REAL (DP), INTENT (IN) :: vm3, vm2, vm1, v, vp1, vp2, vp3

! The output of the subroutine, the flux at cell boundary !
REAL (DP), INTENT (OUT) :: vp_out, vm_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Interpolate and get the right and left states !
vp_out = (3.7D1/6.0D1)*(vp1 + v) - (8.0D0/6.0D1)*(vp2 + vm1) + (1.0D0/6.0D1)*(vp3 + vm2)
vm_out = (3.7D1/6.0D1)*(v + vm1) - (8.0D0/6.0D1)*(vp1 + vm2) + (1.0D0/6.0D1)*(vp2 + vm3)

END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine steepen the density profile near contact discontinuity !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_STEEPEN
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! integer parameter !
INTEGER :: j

! Adiabatic Index !
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: gam_1
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: gam_2

! Altered face values of density !
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: rhop_1, rhom_1
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: rhop_2, rhom_2

! Numerical Second Derivative !
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: deltam_1, deltarho_1, delta2rho_1
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: deltam_2, deltarho_2, delta2rho_2

! Numerical Third Derivative !
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: eta_tilde_1, eta_out_1
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: eta_tilde_2, eta_out_2

! Dummy variables for criteria !
REAL (DP) :: criteria_1, criteria_2, criteria_3, criteria_4

! Do For NM first !
! Assign first derivative !
DO j = -3, length_step_part_2 + 4
	deltarho_2 (j) = 0.5D0*(rho2(j+1) - rho2(j-1))
	If((rho2(j+1) - rho2(j))*(rho2(j) - rho2(j-1)) > 0.0D0) THEN
		deltam_2 (j) = min(abs(deltarho_2 (j)),2.0D0*abs(rho2(j) - rho2(j-1)),2.0D0*abs(rho2(j) - rho2(j+1)))*sign(1.0D0, deltarho_2 (j))
	ELSE
		deltam_2 (j) = 0.0D0
	END IF
END DO

! Assign second derivative !
DO j = -3, length_step_part_2 + 4
	delta2rho_2 (j) = (rho2(j+1) - 2.0D0*rho2(j) + rho2(j-1))/(6.0D0*dx2**2)
END DO

! Get the eta !
DO j = -2, length_step_part_2 + 3
	criteria_1 = abs(rho2(j+1) - rho2(j-1))/min(rho2(j+1),rho2(j-1))
	criteria_2 = delta2rho_2 (j-1)*delta2rho_2 (j+1)
	If(criteria_1 < 0.01D0 .OR. criteria_2 > 0.0D0) THEN
		eta_tilde_2 (j) = 0.0D0
	ELSE
		eta_tilde_2 (j) = -(delta2rho_2 (j+1) - delta2rho_2 (j-1))/(rho2(j+1) - rho2(j-1))*dx2**2
	ENDIF
	eta_out_2 (j) = max(0.0D0,min(2.0D1*(eta_tilde_2 (j) - 0.05D0), 1.0D0))
END DO

! Shock detection !
DO j = -2, length_step_part_2 + 3
	CALL EOSGAMMA(p2(j), rho2(j), cs2(j), gam_2 (j))
	criteria_3 = abs(p2(j+1) - p2(j-1))/(min(p2(j+1),p2(j-1)))
	criteria_4 = 0.1D0*gam_2(j)*abs(rho2(j+1) - rho2(j-1))/(min(rho2(j+1),rho2(j-1)))
	If(criteria_3 - criteria_4 > 0.0D0) THEN
		eta_out_2 (j) = 0.0D0
	END IF
END DO

! Compute altered face values !
DO j = -2, length_step_part_2 + 3
	rhom_2 (j) = rho2(j-1) + 0.5D0*deltam_2 (j-1)
	rhop_2 (j) = rho2(j+1) - 0.5D0*deltam_2 (j+1)
END DO
DO j = -2, length_step_part_2 + 3
	rho2R(j-1) = rho2R(j-1)*(1.0D0 - eta_out_2 (j)) + rhom_2(j)*eta_out_2 (j)
	rho2L(j) = rho2L(j)*(1.0D0 - eta_out_2 (j)) + rhop_2(j)*eta_out_2 (j)
END DO

! Do all the same for DM !
IF(RUNDM_flag == 1) THEN
	DO j = -3, length_step_part_2 + 4
		deltarho_1 (j) = 0.5D0*(rho1(j+1) - rho1(j-1))
		If((rho1(j+1) - rho1(j))*(rho1(j) - rho1(j-1)) > 0.0D0) THEN
			deltam_1 (j) = min(abs(deltarho_1(j)),2.0D0*abs(rho1(j) - rho1(j-1)),2.0D0*abs(rho1(j) - rho1(j+1)))*sign(1.0D0, deltarho_1 (j))
		ELSE
			deltam_1 (j) = 0.0D0
		END IF	
	END DO
	DO j = -3, length_step_part_1 + 4
		delta2rho_1 (j) = (rho1(j+1) - 2.0D0*rho1(j) + rho1(j-1))/(6.0D0*dx1**2)
	END DO
	DO j = -2, length_step_part_1 + 3
		criteria_1 = abs(rho1(j+1) - rho1(j-1))/min(rho1(j+1),rho1(j-1))
		criteria_2 = delta2rho_1 (j-1)*delta2rho_1 (j+1)
		If(criteria_1 < 0.01D0 .OR. criteria_2 > 0.0D0) THEN
			eta_tilde_1 (j) = 0.0D0
		ELSE
			eta_tilde_1 (j) = -(delta2rho_1 (j+1) - delta2rho_1 (j-1))/(rho1(j+1) - rho1(j-1))*dx1**2
		ENDIF
		eta_out_1 (j) = max(0.0D0,min(2.0D1*(eta_tilde_1 (j) - 0.05D0), 1.0D0))
	END DO
	DO j = -2, length_step_part_1 + 3
		CALL EOSGAMMA(p1(j), rho1(j), cs1(j), gam_1 (j))
		criteria_3 = abs(p1(j+1) - p1(j-1))/(min(p1(j+1),p1(j-1)))
		criteria_4 = 0.1D0*gam_1(j)*abs(rho1(j+1) - rho1(j-1))/(min(rho1(j+1),rho1(j-1)))
		If(criteria_3 - criteria_4 > 0.0D0) THEN
			eta_out_1 (j) = 0.0D0
		END IF
	END DO
	DO j = -2, length_step_part_1 + 3
		rhom_1 (j) = rho1(j-1) + 0.5D0*deltam_1 (j-1)
		rhop_1 (j) = rho1(j+1) - 0.5D0*deltam_1 (j+1)
	END DO
	DO j = -2, length_step_part_1 + 3
		rho1R(j-1) = rho1R(j-1)*(1.0D0 - eta_out_1 (j)) + rhom_1 (j)*eta_out_1 (j)
		rho1L(j) = rho1L(j)*(1.0D0 - eta_out_1 (j)) + rhop_1 (j)*eta_out_1 (j)
	END DO
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculate the flattening coefficient !
! For profile flattening near strong shocks 	       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_FLATTEN
USE RIEMANN_MODULE
USE DEFINITION
IMPLICIT NONE

! integer parameter !
INTEGER :: j

! Criteria !
REAL (DP) :: criteria

! Velocity divergence !
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: div_v1
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: div_v2

! pressure graident !
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: p_grad1
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: p_grad2

! intermediate flattening coefficient !
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: f_tilde1
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: f_tilde2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do For NM First !
! Get the velocity divergence !
If(sp_dim == 0) THEN
	DO j = -3, length_step_part_2 + 4
		div_v2(j) = vel2(j+1) - vel2(j-1)
	END DO
ELSEIf(sp_dim == 1) THEN
	DO j = -3, length_step_part_2 + 4
		div_v2(j) = r2(j+1)*vel2(j+1) - r2(j-1)*vel2(j-1)
	END DO
ELSEIf(sp_dim == 2) THEN
	DO j = -3, length_step_part_2 + 4
		div_v2(j) = r2(j+1)**2*vel2(j+1) - r2(j-1)**2*vel2(j-1)
	END DO
END IF

! Get pressure gradient !
DO j = -2, length_step_part_2 + 3
	p_grad2(j) = (p2(j+1) - p2(j-1))/(p2(j+2) - p2(j-2))
END DO

! Set the criteria !
DO j = -2, length_step_part_2 + 3
	criteria = abs(p2(j+1) - p2(j-1))/min(p2(j+1),p2(j-1)) - 1.0D0/3.0D0
	If(criteria < 0.0D0 .OR. div_v2(j) > 0.0D0) THEN
		f_tilde2(j) = 0.0D0
	ELSE
		f_tilde2(j) = max(0.0D0,min(1.0D0,1.0D1*(p_grad2(j) - 0.75D0)))
	END IF
END DO

! Find the flattening coefficient !
DO j = -1, length_step_part_2 + 2
	If(p2(j+1) - p2(j-1) < 0.0D0) THEN
		flat_2 (j) = max(f_tilde2(j), f_tilde2(j+1))
	ELSE
		flat_2 (j) = max(f_tilde2(j), f_tilde2(j-1))
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! DO the same for DM !
IF(RUNDM_flag == 1) THEN
	If(sp_dim == 0) THEN
		DO j = -3, length_step_part_1 + 4
			div_v1(j) = vel1(j+1) - vel1(j-1)
		END DO
	ELSEIf(sp_dim == 1) THEN
		DO j = -3, length_step_part_1 + 4
			div_v1(j) = r1(j+1)*vel1(j+1) - r1(j-1)*vel1(j-1)
		END DO
	ELSEIf(sp_dim == 2) THEN
		DO j = -3, length_step_part_1 + 4
			div_v1(j) = r1(j+1)**2*vel1(j+1) - r1(j-1)**2*vel1(j-1)
		END DO
	END IF
	DO j = -2, length_step_part_1 + 3
		p_grad1(j) = (p1(j+1) - p1(j-1))/(p1(j+2) - p1(j-2))
	END DO
	DO j = -2, length_step_part_1 + 3
		criteria = abs(p1(j+1) - p1(j-1))/min(p1(j+1),p1(j-1)) - 1.0D0/3.0D0
		If(criteria < 0.0D0 .OR. div_v1(j) > 0.0D0) THEN
			f_tilde1(j) = 0.0D0
		ELSE
			f_tilde1(j) = max(0.0D0,min(1.0D0,1.0D1*(p_grad1(j) - 0.75D0)))
		END IF
	END DO
	DO j = -1, length_step_part_1 + 2
		If(p1(j+1) - p1(j-1) < 0.0D0) THEN
			flat_1 (j) = max(f_tilde1(j), f_tilde1(j+1))
		ELSE
			flat_1 (j) = max(f_tilde1(j), f_tilde1(j-1))
		END IF
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The Original PPM method by Colella (1984) named as:		     !
! The Piecewise Parabolic Method (PPM) for Gas-Dynamical Simulations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PPM_MONOTONIZE (vm2, vm1, v, vp1, vp2, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! integer parameter !
INTEGER :: j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The input into the subroutine, including conservative variable and input flux function !
REAL (DP), INTENT (IN) :: vm2, vm1, v, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL (DP), INTENT (INOUT) :: vm_out, vp_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Parameter for limiters !
REAL (DP), PARAMETER :: nppm = 2.0D0

! Local left and right state variables !
REAL (DP) :: vl, vr

! Some variables for differences !
REAL (DP) :: dv, dvp, dvm
REAL (DP) :: d1, dL, dR
REAL (DP) :: d1p, dLp, dRp
REAL (DP) :: d1m, dLm, dRm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Interpolate and get the right and left states !
vr = vp_out
vl = vm_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set the differences !
d1 = vp1 - vm1
dL = v - vm1
dR = vp1 - v
d1p = vp2 - v
dLp = vp1 - v
dRp = vp2 - vp1
d1m = v - vm2 
dLm = vm1 - vm2
dRm = v - vm1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set the Van-Leer limiter !
If(dL*dR > 0.0D0) THEN
	dv = sign(1.0D0,d1)*min(0.5D0*abs(d1),2.0D0*abs(dL),2.0D0*abs(dR))
ELSE
	dv = 0.0D0
ENDIF
If(dLp*dRp > 0.0D0) THEN
	dvp = sign(1.0D0,d1p)*min(0.5D0*abs(d1p),2.0D0*abs(dLp),2.0D0*abs(dRp))
ELSE
	dvp = 0.0D0
ENDIF
If(dLm*dRm > 0.0D0) THEN
	dvm = sign(1.0D0,d1m)*min(0.5D0*abs(d1m),2.0D0*abs(dLm),2.0D0*abs(dRm))
ELSE
	dvm = 0.0D0
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check monotonicity constraints !
If((v - vr)*(vp1 - vr) > 0.0D0) THEN
	vr = 0.5D0*(v + vp1) - (dvp - dv)/6.0D0
END IF

! For left state !
If((vm1 - vl)*(v - vl) > 0.0D0) THEN
	vl = 0.5D0*(vm1 + v) - (dv - dvm)/6.0D0
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check for extremum conditions !
If((vr - v)*(v - vl) <= 0.0D0 .OR. (vm1 - v)*(v - vp1) <= 0.0D0) THEN
	vr = v
	vl = v
ELSE

	! Limiter away from limiters !
	If(abs(vr - v)>=nppm*abs(v - vl)) THEN
		vr = v + nppm*(v - vl)
	ELSE
		CONTINUE
	END IF
	If(abs(vl - v)>=nppm*abs(v - vr)) THEN
		vl = v + nppm*(v - vr)
	ELSE
		CONTINUE
	END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign output !
vp_out = vr
vm_out = vl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The Modified PPM method by Colella to preserve extremum !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE XPPM_MONOTONIZE (vm2, vm1, v, vp1, vp2, vm_out, vp_out)
USE DEFINITION
IMPLICIT NONE

! integer parameter !
INTEGER :: j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The input into the subroutine, including conservative variable and input flux function !
REAL (DP), INTENT (IN) :: vm2, vm1, v, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL (DP), INTENT (INOUT) :: vp_out, vm_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Parameter for constants !
REAL (DP), PARAMETER :: const = 1.25D0

! Parameter for limiters !
REAL (DP), PARAMETER :: nppm = 2.0D0

! Local left and right state variables !
REAL (DP) :: vl, vr, alphal, alphar

! Temporal variables !
REAL (DP) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
REAL (DP) :: tmp7, tmp8, tmp9, tmp10, signs

! Some variables for differences !
REAL (DP) :: d2vp, d2vpL, d2vpR
REAL (DP) :: d2vm, d2vmL, d2vmR
REAL (DP) :: d2vplim, d2vmlim, d2vlim
REAL (DP) :: d2v, d2vC, d2vL, d2vR

! For alternate limiters !
REAL (DP) :: dfext, dv

! Right and left state !
REAL (DP) :: left, right

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Interpolate and get the right and left states !
vr = vp_out
vl = vm_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set the differences !
d2vp = 3.0D0*(v - 2.0D0*vr + vp1)
d2vpL = vm1 - 2.0D0*v + vp1
d2vpR = v - 2.0D0*vp1 + vp2
d2vm = 3.0D0*(vm1 - 2.0D0*vl + v)
d2vmL = vm2 - 2.0D0*vm1 + v
d2vmR = vm1 - 2.0D0*v + vp1
d2v = 6.0D0*(vl - 2.0D0*v + vr)
d2vC = vm1 - 2.0D0*v + vp1
d2vL = vm2 - 2.0D0*vm1 + v
d2vR = v - 2.0D0*vp1 + vp2

! Set the sign !
tmp1 = sign(1.0D0,d2vp)
tmp2 = sign(1.0D0,d2vpL)
tmp3 = sign(1.0D0,d2vpR)
tmp4 = sign(1.0D0,d2vm)
tmp5 = sign(1.0D0,d2vmL)
tmp6 = sign(1.0D0,d2vmR)
tmp7 = sign(1.0D0,d2v)
tmp8 = sign(1.0D0,d2vC)
tmp9 = sign(1.0D0,d2vL)
tmp10 = sign(1.0D0,d2vR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set the Van Leer limiter !
If(tmp1 == tmp2 .AND. tmp2 == tmp3) THEN
	d2vplim = tmp1*min(const*abs(d2vpL),const*abs(d2vpR),abs(d2vp))
ELSE
	d2vplim = 0.0D0
END IF
If(tmp4 == tmp5 .AND. tmp5 == tmp6) THEN
	d2vmlim = tmp4*min(const*abs(d2vmL),const*abs(d2vmR),abs(d2vm))
ELSE
	d2vmlim = 0.0D0
END IF
If(tmp7 == tmp8 .AND. tmp8 == tmp9 .AND. tmp9 == tmp10) THEN
	d2vlim = tmp7*min(const*abs(d2vL),const*abs(d2vR),const*abs(d2vC),abs(d2v))
ELSE
	d2vlim = 0.0D0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check monotonicity constraints !
If((v - vr)*(vp1 - vr) > 0.0D0) THEN
	vr = 0.5D0*(v + vp1) - d2vplim/6.0D0
END IF

! For left state !
If((vm1 - vl)*(v - vl) > 0.0D0) THEN
	vl = 0.5D0*(vm1 + v) - d2vmlim/6.0D0
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign alpha !
alphar = vr - v
alphal = vl - v

! Check local extremum constraints !
If((vr - v)*(v - vl) <= 0.0D0 .OR. (vm1 - v)*(v - vp1) <= 0.0D0) THEN
	If(d2v == 0.0D0) THEN
		vr = v
		vl = v
	ELSE
		vr = v + (vr - v)*d2vlim/d2v
		vl = v + (vl - v)*d2vlim/d2v
	END IF
ELSE
	
	! We can choose different limiter away from extremum !
	If(abs(vr - v) >= nppm*abs(v - vl)) THEN
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Mutted this because may be unstable !
		!dfext = -0.25D0*alphar/(alphar + alphal)
		!dv = vm1 - v
		!signs = sign(1.0D0,dv)
		!If(dfext >= signs*dv) THEN
		!	vr = v - 2.0D0*(dv + signs*sqrt(dv**2 - dv*alphal))
		!ELSE
		!	CONTINUE
		!END IF
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		vr = v + nppm*(v - vl)
	ELSE
		CONTINUE
	END IF
	If(abs(vl - v) >= nppm*abs(v - vr)) THEN
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Mutted this because may be unstable !
		!dfext = -0.25D0*alphal/(alphar + alphal)
		!dv = vp1 - v
		!signs = sign(1.0D0,dv)
		!If(dfext >= signs*dv) THEN
		!	vl = v - 2.0D0*(dv + signs*sqrt(dv**2 - dv*alphar))
		!ELSE
		!	CONTINUE
		!END IF
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		vl = v + nppm*(v - vr)
	ELSE
		CONTINUE
	END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!S

! Assign output !
vm_out = vl
vp_out = vr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

END MODULE