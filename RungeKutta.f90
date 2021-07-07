!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the 4 step rungekutta time evolution of the hydro - equation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RUNGEKUTTA (n)
USE DEFINITION
USE WENO_MODULE
USE NUCLEAR_MODULE
USE FLAME_MODULE
IMPLICIT NONE

! Integer parameter !
INTEGER, INTENT (IN) :: n
INTEGER :: i, j

! Dummy !
REAL (DP) :: rhoaold, dummy

! For movinggrid !
REAL (DP) :: dx1_old, dx1_2, dx1_3, a1_3
REAL (DP) :: dx2_old, dx2_2, dx2_3, a2_3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the boundaries !
IF(movinggridnm_flag == 1) THEN
	CALL FINDBOUNDARY_NM
END IF
IF(movinggriddm_flag == 1) THEN
	CALL FINDBOUNDARY_DM
END IF

! assign the grid size !
IF(movinggridnm_flag == 1) THEN
	dx2_old = dx2
END IF
IF(movinggriddm_flag == 1) THEN
	dx1_old = dx1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First step in SSPRK(5,4) !
CALL BACKUPCONS (u_old1, u_temp1, u_old2, u_temp2)

! First calculate evaulate the flux difference !
CALL SPATIAL (u_old1, u_old2)

! Now do the first iteration step !
IF (RUNDM_flag == 1) THEN
	DO i = imin1, imax1
		DO j = 1, length_step_part_1
			u_new1 (i, j) = u_old1 (i, j) + 0.391752226571890D0 * dt * l1 (i, j)
		END DO
	END DO
END IF
DO i = imin2, imax2
	DO j = 1, length_step_part_2
		u_new2 (i, j) = u_old2 (i, j) + 0.391752226571890D0 * dt * l2 (i, j)
	END DO
END DO

! Update dx !
IF(movinggriddm_flag == 1) then
	dx1 = dx1_old + 0.391752226571890D0 * dt * vel1_max * dx1 / radius1
	call getgrid_DM
ENDIF
IF(movinggridnm_flag == 1) then
	dx2 = dx2_old + 0.391752226571890D0 * dt * vel2_max * dx2 / radius2
	call getgrid_NM
ENDIF

! We copy the value to ghost shell !
CALL BOUNDARY2D_DM (u_new1)
CALL BOUNDARY2D_NM (u_new2)

! Convert conservative variables to natural variables !
CALL FROMUTORVE (u_new1, u_new2)

! We update abar and zbar !
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF

! Update everything !
CALL UPDATE

! Do conversion again !
CALL FROMRVETOU (u_new1, u_new2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Second step in SSPRK(5,4) !
CALL BACKUPCONS (u_new1, u_temp1, u_new2, u_temp2)
CALL SPATIAL (u_new1, u_new2)

IF (RUNDM_flag == 1) THEN
	DO i = imin1, imax1
		DO j = 1, length_step_part_1
			u_new1 (i, j) = 0.444370493651235D0 * u_old1 (i, j) + 0.555629506348765D0 * u_new1 (i, j) + 0.368410593050371D0 * dt * l1 (i, j)
			u2_dm (i, j) = u_new1 (i, j) 	
		END DO
	END DO
END IF
DO i = imin2, imax2
	DO j = 1, length_step_part_2
		u_new2 (i, j) = 0.444370493651235D0 * u_old2 (i, j) + 0.555629506348765D0 * u_new2 (i, j) + 0.368410593050371D0 * dt * l2 (i, j)
		u2_nm (i, j) = u_new2 (i, j) 
	END DO
END DO

! Update dx !
if(movinggriddm_flag == 1) then
	dx1 = 0.444370493651235D0 * dx1_old + 0.555629506348765D0 * dx1 +  0.368410593050371D0 * dt * vel1_max * dx1 / radius1
	dx1_2 = dx1
	call getgrid_DM
endif
if(movinggridnm_flag == 1) then
	dx2 = 0.444370493651235D0 * dx2_old + 0.555629506348765D0 * dx2 +  0.368410593050371D0 * dt * vel2_max * dx2 / radius2
	dx2_2 = dx2
	call getgrid_NM
endif

CALL BOUNDARY2D_DM (u_new1)
CALL BOUNDARY2D_NM (u_new2)
CALL FROMUTORVE (u_new1,u_new2)
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF
CALL UPDATE
CALL FROMRVETOU (u_new1,u_new2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Third step in SSPRK(5,4) !
CALL BACKUPCONS (u_new1, u_temp1, u_new2, u_temp2)
CALL SPATIAL (u_new1, u_new2)

IF (RUNDM_flag == 1) THEN
	DO i = imin1, imax1
		DO j = 1, length_step_part_1
			u_new1 (i, j) = 0.620101851488403D0 * u_old1 (i, j) + 0.379898148511597D0 * u_new1 (i, j) + 0.251891774271694D0 * dt * l1 (i, j)
			u3_dm (i, j) = u_new1 (i, j)
			l3_dm (i, j) = l1 (i, j)
		END DO
	END DO
END IF
DO i = imin2, imax2
	DO j = 1, length_step_part_2
		u_new2 (i, j) = 0.620101851488403D0 * u_old2 (i, j) + 0.379898148511597D0 * u_new2 (i, j) + 0.251891774271694D0 * dt * l2 (i, j)
		u3_nm (i, j) = u_new2 (i, j)
		l3_nm (i, j) = l2 (i, j)
	END DO
END DO

! Update dx !
if(movinggriddm_flag == 1) then
	dx1 = 0.620101851488403D0 * dx1_old + 0.379898148511597D0 * dx1 +  0.251891774271694D0 * dt * vel1_max * dx1 / radius1
	dx1_3 = dx1
	a1_3 = vel1_max * dx1 / radius1
	call getgrid_DM
endif
if(movinggridnm_flag == 1) then
	dx2 = 0.620101851488403D0 * dx2_old + 0.379898148511597D0 * dx2 +  0.251891774271694D0 * dt * vel2_max * dx2 / radius2
	dx2_3 = dx2
	a2_3 = vel2_max * dx2 / radius2
	call getgrid_NM
endif

CALL BOUNDARY2D_DM (u_new1)
CALL BOUNDARY2D_NM (u_new2)
CALL FROMUTORVE (u_new1,u_new2)
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF
CALL UPDATE
CALL FROMRVETOU (u_new1,u_new2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Fourth step in SSPRK(5,4) !
CALL BACKUPCONS (u_new1, u_temp1, u_new2, u_temp2)
CALL SPATIAL (u_new1, u_new2)

IF (RUNDM_flag == 1) THEN
	DO i = imin1, imax1
		DO j = 1, length_step_part_1
			u_new1 (i, j) = 0.178079954393132D0 * u_old1 (i, j) + 0.821920045606868D0 * u_new1 (i, j) + 0.544974750228521D0 * dt * l1 (i, j)
		END DO
	END DO
END IF
DO i = imin2, imax2
	DO j = 1, length_step_part_2
		u_new2 (i, j) = 0.178079954393132D0 * u_old2 (i, j) + 0.821920045606868D0 * u_new2 (i, j) + 0.544974750228521D0 * dt * l2 (i, j)
	END DO
END DO

! Update dx !
if(movinggriddm_flag == 1) then
	dx1 = 0.178079954393132D0 * dx1_old + 0.821920045606868D0 * dx1 +  0.544974750228521D0 * dt * vel1_max * dx1 / radius1
	call getgrid_DM
endif
if(movinggridnm_flag == 1) then
	dx2 = 0.178079954393132D0 * dx2_old + 0.821920045606868D0 * dx2 +  0.544974750228521D0 * dt * vel2_max * dx2 / radius2
	call getgrid_NM
endif

CALL BOUNDARY2D_DM (u_new1)
CALL BOUNDARY2D_NM (u_new2)
CALL FROMUTORVE (u_new1,u_new2)
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF
CALL UPDATE
CALL FROMRVETOU (u_new1,u_new2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Final step in SSPRK(5,4) !
CALL BACKUPCONS (u_new1, u_temp1, u_new2, u_temp2)
CALL SPATIAL (u_new1, u_new2)

IF(RUNDM_flag == 1) THEN
	DO i = imin1, imax1
		DO j = 1, length_step_part_1
			u_new1 (i, j) = 0.517231671970585D0 * u2_dm (i, j) + 0.096059710526147D0 * u3_dm (i, j) + 0.386708617503269D0 * u_new1 (i, j) &
				+ 0.063692468666290D0 * dt * l3_dm(i, j) + 0.226007483236906D0 * dt * l1 (i, j)
		END DO
	END DO
END IF
DO i = imin2, imax2
	DO j = 1, length_step_part_2
		u_new2 (i, j) = 0.517231671970585D0 * u2_nm (i, j) + 0.096059710526147D0 * u3_nm (i, j) + 0.386708617503269D0 * u_new2 (i, j) &
			+ 0.063692468666290D0 * dt * l3_nm(i, j) + 0.226007483236906D0 * dt * l2 (i, j)
	END DO
END DO

! Update dx !
if(movinggriddm_flag == 1) then
	dx1 = 0.517231671970585D0 * dx1_2 + 0.096059710526147D0 * dx1_3 + 0.386708617503269D0 * dx1 &
		+ 0.063692468666290D0 * dt * a1_3 + 0.226007483236906D0 * dt * vel1_max * dx1 / radius1
	call getgrid_DM
endif
if(movinggridnm_flag == 1) then
	dx2 = 0.517231671970585D0 * dx2_2 + 0.096059710526147D0 * dx2_3 + 0.386708617503269D0 * dx2 &
		+ 0.063692468666290D0 * dt * a2_3 + 0.226007483236906D0 * dt * vel2_max * dx2 / radius2
	call getgrid_NM
endif

CALL BOUNDARY2D_DM (u_new1)
CALL BOUNDARY2D_NM (u_new2)
CALL FROMUTORVE (u_new1, u_new2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! adjust the atmospheric motion !
if(sponge_flag == 1) THEN
	call getsponge
END IF

! We update abar and zbar !
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF

! Update temperature
if(helmeos_flag == 1) THEN
	CALL FINDHELMTEMP
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the update in time that is seperated from the hydro equation !
! This is the simpe operator splitting technique !
IF(flame_flag == 1) THEN

	! We need to propagate the deflagration !
	CALL UPDATE_FLAME_RADIUS
	IF (fusion_flag == 1) THEN
	
		! This does the Carbon burning
		IF(carburn_flag == 1) THEN
			CALL BURN_PHASE1
		END IF

		! This do the O- and Si- burning
		IF(advburn_flag == 1) THEN
			CALL BURN_PHASE2
		END IF

		! Update the AZbar and temperature accordingly
		CALL FIND_AZBAR()
		CALL FINDHELMTEMP

		! For completely burnt zone, check if NSE applies
		IF(convert_nse_flag == 1) THEN
			CALL NSE
		END IF

		! Check if the change of isotope perserve the sum
		CALL CHECKXISOTOPE

		! Update Abar and Zbar and temperature again
		CALL Find_AZBAR()
		CALL FINDHELMTEMP

   	ENDIF

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Section for adjusting atmospheric density !
IF(fixrhonm_flag == 1) THEN
	rhoaold = rho2_a
	rho2_a = min(min(maxval(rho2), rho2_c)*rhofac_2, rho2_a, minval(rho2))

	! Assign atmoshperic epsilon !
	IF(nm_epsilon == 1 .AND. fermieosnm_flag == 0) THEN
		epsilon2_a = k2 * rho2_a ** (gamma2 - 1.0E0_DP) / (gamma2 - 1.0E0_DP)
	ELSEIF(helmeos_flag == 1) THEN
		CALL HELMEOS_RtoE (rho2_a, temp2_a, abar2_ini, zbar2_ini, ye2_ini, epsilon2_a, dummy)
	ELSE
		CALL EOSEPSILON (rho2_a, dummy, epsilon2_a, 2)
	END IF

	! Adjust density !
	DO j = -4, length_step_2 + 5
		IF(rho2(j) == rhoaold) THEN
			rho2(j) = rho2_a
		END IF
	END DO
END IF

! For DM !
IF(fixrhodm_flag == 1) THEN
	rhoaold = rho1_a
	rho1_a = min(min(maxval(rho1), rho1_c)*rhofac_1, rho1_a)

	! Assign atmoshperic epsilon !
	CALL EOSEPSILON (rho1_a, dummy, epsilon1_a, 1)

	! Adjust density !
	DO j = -4, length_step_1 + 5
		IF(rho1(j) == rhoaold) THEN
			rho1(j) = rho1_a
		END IF
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! At last Check if the density reach atmospheric density !
If (test_model == 0) THEN
	CALL CHECKRHO
END IF

! We update abar and zbar !
IF(xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF

! Do the update on the remaining variables correspondingly !
! The integer indicate that whether we should calculate the gravitational potential !
CALL UPDATE

! Once again update the numerical value of u since we did some changes in density !
CALL FROMRVETOU (u_new1, u_new2)
CALL BACKUPCONS (u_new1, u_old1, u_new2, u_old2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE