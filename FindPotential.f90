!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine find the total gravitational potential energy of NM and DM !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDPOTENTIAL
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
INTEGER :: j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We need to find the total gravitating mass !
!CALL findgravrho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We find the gravitational potential only if you set the precense of gravity !
! First, we guess the form of the potential !
IF (w_gravity_i == 1) THEN

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! These are old scheme where we used iterative poisson solvers, mutted !
	! We find the trial gravitational potential first !
	!temp (1) = (4.0E0_DP / 3.0E0_DP) * pi_old * r2(1) ** 3 * rho (1)
	!temp (2) = 9.0E0_DP * temp (1) + (4.0E0_DP / 3.0E0_DP) * pi_old * dx * r2(2) ** 2 * rho (2)
	!DO j = 3, length_step
	!	temp (j) = temp (j - 2) + 4.0E0_DP * pi_old * (dx / 3.0E0_DP) &
	!			* (rho (j - 2) * r2(j-2) ** 2 &
	!			+ 4.0E0_DP * rho (j - 1) * r2(j-1) ** 2 &
	!			+ rho (j) * r2(j) ** 2)
	!END DO
	!Do j = 1, length_step
	!	phip (j) = temp (j) / r2(j) ** 2
	!END DO
	! We copy the value to ghost shell !
	!CALL BOUNDARY1D (phip, odd)
	! Find the mass !
	!CALL FINDMASS
	! Fix the boundary condition, Dirichlet boundary condition !
	! At the edge of simulation box, the potential goes as GM/r !
	!If (DM_flag == 1) THEN
	!	phi (length_step) = - (mass1 + mass2) / r2(j)
	!ELSE	
	!	phi (length_step) = - mass2 / r2(j)
	!END IF 
	!DO j = length_step - 1, 1, -1
	!	phi (j) = phi (j + 1) - (phip (j) + phip (j + 1)) * dx / 2.0E0_DP
	!END DO
	! We copy the value to ghost shell !
	!CALL BOUNDARY1D (phi, even)
	! We call the relaxation iteration method to find the actual potential !
	!CALL POTENTIALRELAX
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! New method, use multipole expansion !
	IF(DM_flag == 1) THEN
		CALL MULTIPOLE_DM
	END IF
	CALL MULTIPOLE_NM
	CALL INTERPOLATION

	! To see if we need to use GR potential !
	IF(gr_flag == 1) THEN
		CALL GRPOTENTIAL
		phinm = phigr
	END IF
ELSE
	
	! No gravitational potential, set them to zero !
	DO j = -4, length_step_2 + 5
		phinm (j) = 0.0D0
		phinm_1 (j) = 0.0E0_DP	
		phinm_2 (j) = 0.0E0_DP	
		phip_nm (j) = 0.0E0_DP
	END DO
	
	! For running DM !
	IF(DM_flag == 1) THEN
		DO j = -4, length_step_1 + 5
			phidm (j) = 0.0D0
			phidm_1 (j) = 0.0E0_DP	
			phidm_2 (j) = 0.0E0_DP	
			phip_dm (j) = 0.0E0_DP
		END DO
	END IF

END IF

CALL GRPOTENTIAL

END SUBROUTINE