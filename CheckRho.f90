!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine check whether the density of NM or DM      !
! reached atmospheric density, and alter primitive variables !
! accordingly 						     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CHECKRHO
USE DEFINITION
USE NUCLEAR_MODULE
IMPLICIT NONE

! Integer variables !
INTEGER :: j, grid1, grid2

! Real Variables !
REAL(DP) :: rho_min_nm, rho_min_dm

! We do the NM first !
! Set NM Threshold Density !
rho_min_nm = 1.1E0_DP * rho2_a

! Assign the appropriate values of hydro variables !
DO j = -4, length_step_2 + 5
	IF (rho2 (j) <= rho_min_nm) THEN
		rho2 (j) = rho2_a
		epsilon2 (j) = epsilon2_a
		vel2 (j) = vel2_a
		ye2 (j) = ye2_a
		temp2 (j) = temp2_a
	END IF
END DO

! For Dual energy !
IF(dual_energy == 1) THEN
	DO j = -4, length_step_2 + 5
		IF (rho2 (j) <= rho_min_nm) THEN
			rhoe2(j) = rho2_a*epsilon2_a
		END IF
	END DO
	CALL BOUNDARY1D_NM (rhoe2, even) 
END IF

! Copy the variables to ghost shell !
CALL BOUNDARY1D_NM (rho2, even)
CALL BOUNDARY1D_NM (vel2, odd)
CALL BOUNDARY1D_NM (epsilon2, even)
CALL BOUNDARY1D_NM (ye2, even)  
CALL BOUNDARY1D_NM (temp2, even) 

if(xisotran_flag == 1) then
	DO j = -4, length_step_2 + 5
		IF (rho2 (j) <= rho_min_nm) THEN
			xiso(:, j) = xiso_a
		END IF
	END DO

	! Check whether the sum of mass fraction equals one !
	call checkxisotope
	
	! Copy the chemical composition to ghost cell !   
	CALL BOUNDARY2D_X  
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the same thing for DM only if DM is precense !
IF (DM_flag == 1) THEN

	! Set DM Threshold Density !
	rho_min_dm = 1.1E0_DP * rho1_a

	! Assign the appropriate values of hydro variables !
	DO j = -4, length_step_1 + 5
		IF (rho1 (j) <= rho_min_dm) THEN
			rho1 (j) = rho1_a
			epsilon1 (j) = epsilon1_a
		END IF
	END DO

	! Copy the variables to ghost shell !
	CALL BOUNDARY1D_DM (rho1, even)
	CALL BOUNDARY1D_DM (epsilon1, even)

	! assign value for velocity only if DM is movable !
	IF (RUNDM_flag == 1) THEN
		DO j = -4, length_step_1 + 5
			IF (rho1 (j) <= rho_min_dm) THEN
				vel1 (j) = vel1_a
			END IF
		END DO
		CALL BOUNDARY1D_DM (vel1, odd)
	END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we start to look for the minimum size of the box !
! which can contain the whole star, but minimize the calculation !
IF(checkstepnm_flag == 1) THEN
	
	! Obtain length_step_part_2 !
	! Initialize grid2 !
	grid2 = 0
	
	! Find the outermost grid which is the star surface !
	! Do for NM first !
	DO j = 1, length_step_2 - 1
		if(rho2 (j) > rho_min_nm .and. rho2 (j+1) == rho2_a) then
			if(j > grid2) then
				grid2 = j
			endif
		END IF
	END DO

	! Set the effective length_step_part_2 !
	IF(grid2 == 0) then
		length_step_part_2 = length_step_2
	ELSEIF(grid2 /= 0 .and. grid2 < length_step_2 - 10) then
		length_step_part_2 = grid2 + 10
	ELSE
		length_step_part_2 = length_step_2
	ENDIF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now do the same if DM is movable ! 
IF (checkstepdm_flag == 1) THEN

	! Obtain length_step_part_1 !
	! Initialize grid1 !
	grid1 = 0
	
	! Find the outermost grid which is the star surface !
	DO j = 1, length_step_1 - 1
		if(rho1 (j) > rho_min_dm .and. rho1 (j+1) == rho1_a) then
			if(j > grid1) then
				grid1 = j
			endif
		END IF
	END DO

	! Set the effective length_step_part_1 !
	IF(grid1 == 0) then
		length_step_part_1 = length_step_1
	ELSEIF(grid1 /= 0 .and. grid1 < length_step_1 - 10) then
		length_step_part_1 = grid1 + 10
	ELSE
		length_step_part_1 = length_step_1
	ENDIF

END IF	

END SUBROUTINE