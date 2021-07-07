!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine solve for the initial hydrostatic equilibrium star !
! assuming a two fluid formalism. We assume the newtonian gravity    !
! and is solving for the initial density profile using RK-5 method   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETRHO_2F
USE DEFINITION
USE NUCLEAR_MODULE
IMPLICIT NONE

! The number of hydrostatic equation and the improved accuracy !
INTEGER, PARAMETER :: no_of_eq_ini = 6, more = 10 ** ini_acc

! The extra array arising from the extra accuracy !
INTEGER, PARAMETER :: length_morestep = max(length_step_1, length_step_2) * more

! The smaller length step arising from the extra accuracy !
REAL (DP), PARAMETER :: dxmore = min(dx1_ini, dx2_ini) / DBLE (more)

! Dummy variables !
REAL (DP) :: dummy, rho1_min, rho2_min

! Integer parameters !
INTEGER :: i, j, k, r_grid_more1, r_grid_more2

! Dummy arrays for NM and DM density to be interpolated !	
REAL (DP), DIMENSION (-4 : length_morestep + 5) :: den1, den2

! Temporal distance arrays !
REAL (DP), DIMENSION (-4 : length_morestep + 5) :: rtemp

! This is necessary for any finite temperature EOS since the deriative of pressure plays a role !
REAL (DP), DIMENSION (-4 : length_morestep + 5) :: dpdrho2_temp

! The pressure, density, at center, atmosphere and the dummy varaiables used to stored them !
REAL (DP) :: p1_c, p2_c, p1_a, p2_a, x, ini_rho1, ini_rho2, ini_p1, ini_p2

! Variables essential in the RK-5 method !
REAL (DP), DIMENSION (1 : no_of_eq_ini) :: y_zero, y_one, y_two, y_three, y_four, y_five
REAL (DP), DIMENSION (1 : no_of_eq_ini) :: der_one, der_two, der_three, der_four, der_five, der_new

! Assigning the variables essential in the RK-5 method !
REAL (DP), DIMENSION (1 : no_of_eq_ini, -4 : length_morestep + 5) :: y	

! Quantity related to chemical composition if !
! you are using variable compositions of star !
REAL (DP), DIMENSION (total_ion) :: xiso_ini

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign minimum density !
rho1_min = 1.0D-4*rho1_a
rho2_min = 1.0D-4*rho2_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We read the EOS table for DM !
IF (fermieosdm_flag == 1) THEN
	OPEN (UNIT = 99, FILE = 'EOS_Table1.eos', STATUS = 'OLD')
	DO i = 1, eoslineno
		READ (99, *) eostable1 (i, 1), eostable1 (i, 2)
	END DO
	eosline1 = 1
	CLOSE (99)
END IF

! We read the EOS table for NM !
IF (fermieosnm_flag == 1) THEN
	OPEN (UNIT = 100, FILE = 'EOS_Table2.eos', STATUS = 'OLD')
	DO i = 1, eoslineno
		READ (100, *) eostable2 (i, 1), eostable2 (i, 2)
	END DO
	CLOSE (100)
	eosline2 = 1
END IF

! We convert density at center and atmosphere to the corresponding pressure !
CALL GETRHO_EOSRTOP (p1_c, rho1_c, gs1, mb1, me1, ye1, 1)
CALL GETRHO_EOSRTOP (p1_a, rho1_min, gs1, mb1, me1, ye1, 1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We choose which table is needed to read according to the chosen EOS !
IF(helmeos_flag == 1) THEN

	! We assign initial chemical composition accordingly !
	! If you are using helmeos EOS in construction !
	! It is because the realistic EOS not only depends on !
	! Ye but also Abar and Zbar !	
	xiso_ini = 0.0D0
	xiso_ini(che4) = xhe4_ini
	xiso_ini(cc12) = xc12_ini
	xiso_ini(co16) = xo16_ini
	xiso_ini(cne20) = xne20_ini
	xiso_a = xiso_ini

	! Now convert the composition into mean atomic and mass number !
	CALL private_helmeos_azbar(xiso_ini, abar2_ini, zbar2_ini, ye2_ini)

	! assign atmospheric electron fraction !
	ye2_a = ye2_ini
	
	! Read the helmeos eos table for interpolation !
	CALL read_helm_table()

	! We convert density at center and atmosphere to the corresponding pressure !
	! Note that at center the deriative of density should be zero !
	CALL HELMEOS_RtoP(rho2_c, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, p2_c, dpdrho2_temp(0), dummy, dummy, dummy, dummy)
	CALL HELMEOS_RtoP(rho2_min, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, p2_a,  dummy, dummy, dummy, dummy, dummy)

ELSEIF (ccsneos_flag == 1) THEN

	! CCSN !
        CALL GETRHO_CCSNRTOP (p2_c, rho2_c)
        CALL GETRHO_CCSNRTOP (p2_a, rho2_min)

ELSE

	! We convert density at center and atmosphere to the corresponding pressure !
	CALL GETRHO_EOSRTOP (p2_c, rho2_c, gs2, mb2, me2, ye2_old, 2)
	CALL GETRHO_EOSRTOP (p2_a, rho2_min, gs2, mb2, me2, ye2_old, 2)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! assign distance !
DO j = 1, length_morestep 
	rtemp(j) = dxmore*DBLE(j - 1/2)
END DO

! We assign the value of y at center !
y (1, 0) = (4.0D0/3.0D0)*pi_old*rtemp(1)**3*rho1_c
y (2, 0) = p1_c - 0.5D0*y(1,0)*rho1_c/rtemp(1)
y (3, 0) = rho1_c
y (4, 0) = (4.0D0/3.0D0)*pi_old*rtemp(1)**3*rho2_c
y (5, 0) = p2_c - 0.5D0*y(4,0)*rho2_c/rtemp(1)
y (6, 0) = rho2_c

! initialize the integer parameter !
r_grid_more1 = length_morestep
r_grid_more2 = length_morestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this is the main structure of RK-5 method !
! this should not be unfamiliar so I would skip most the explaination !
DO j = 0, length_morestep - 1
	
	! Update the value of x and y !
	DO i = 1, no_of_eq_ini
		y_zero (i) = y (i, j)
	END DO

	x = (DBLE (j) + 0.5D0) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the first step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_2F (der_one, x, y_zero, no_of_eq_ini)
	
	! We artificially make the value to zero to avoid singularity at center !
	IF (j == 0) THEN
		der_one (2) = 0.0E0_DP
		der_one (5) = 0.0E0_DP
	END IF

	! If the density reach atmospheric values, no changes in all the quantity !
	IF (y_zero (3) == rho1_min) THEN
		der_one (1) = 0.0E0_DP
		der_one (2) = 0.0E0_DP
	END IF
	IF (y_zero (6) == rho2_min) THEN
		der_one (4) = 0.0E0_DP
		der_one (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_one (i) = y_zero (i) + (1.0E0_DP / 4.0E0_DP) * dxmore * der_one (i)
	END DO

	! We determine whether the pressure reached atmospheric pressure !
	IF (j < r_grid_more1) THEN
		ini_p1 = y_one (2)
		IF (ini_p1 <= p1_a) THEN
			r_grid_more1 = j
			ini_p1 = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
	END IF
	IF (j < r_grid_more2) THEN
		ini_p2 = y_one (5)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eostable1, eosline1, 1)

	if(helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, ini_rho2, ini_rho2)
        elseif (ccsneos_flag == 1) then
                call GETRHO_CCSNPTOR (ini_p2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	! Assign density according to pressure !
	If (ini_p1 == p1_a) THEN
		y_one (3) = rho1_min
	ELSE
		y_one (3) = ini_rho1
	END IF
	If (ini_p2 == p2_a) THEN
		y_one (6) = rho2_min
	ELSE
		y_one (6) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (1.0E0_DP / 4.0E0_DP) * dxmore
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the second step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_2F (der_two, x, y_one, no_of_eq_ini)

	IF (y_one (3) == rho1_min) THEN
		der_two (1) = 0.0E0_DP
		der_two (2) = 0.0E0_DP
	END IF
	IF (y_one (6) == rho2_min) THEN
		der_two (4) = 0.0E0_DP
		der_two (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_two (i) = y_zero (i) + (1.0E0_DP / 8.0E0_DP) * dxmore * der_one (i) &
				+ (1.0E0_DP / 8.0E0_DP) * dxmore * der_two (i)
	END DO

	IF (j < r_grid_more1) THEN
		ini_p1 = y_two (2)
		IF (ini_p1 <= p1_a) THEN
			r_grid_more1 = j
			ini_p1 = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
	END IF
	IF (j < r_grid_more2) THEN
		ini_p2 = y_two (5)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eostable1, eosline1, 1)

	if(helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, ini_rho2, ini_rho2)
        elseif (ccsneos_flag == 1) then
                call GETRHO_CCSNPTOR (ini_p2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p1 == p1_a) THEN
		y_two (3) = rho1_min
	ELSE
		y_two (3) = ini_rho1
	END IF
	If (ini_p2 == p2_a) THEN
		y_two (6) = rho2_min
	ELSE
		y_two (6) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (1.0E0_DP / 4.0E0_DP) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the third step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_2F (der_three, x, y_two, no_of_eq_ini)
	
	IF (y_two (3) == rho1_min) THEN
		der_three (1) = 0.0E0_DP
		der_three (2) = 0.0E0_DP
	END IF
	IF (y_two (6) == rho2_min) THEN
		der_three (4) = 0.0E0_DP
		der_three (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_three (i) = y_zero (i) - (1.0E0_DP / 2.0E0_DP) * dxmore * der_two (i) &
			+ dxmore * der_three (i)
	END DO

	IF (j < r_grid_more1) THEN
		ini_p1 = y_three (2)
		IF (ini_p1 <= p1_a) THEN
			r_grid_more1 = j
			ini_p1 = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
	END IF
	IF (j < r_grid_more2) THEN
		ini_p2 = y_three (5)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eostable1, eosline1, 1)

	if(helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, ini_rho2, ini_rho2)
        elseif (ccsneos_flag == 1) then
                call GETRHO_CCSNPTOR (ini_p2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p1 == p1_a) THEN
		y_three (3) = rho1_min
	ELSE
		y_three (3) = ini_rho1
	END IF
	If (ini_p2 == p2_a) THEN
		y_three (6) = rho2_min
	ELSE
		y_three (6) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (1.0E0_DP / 2.0E0_DP) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the forth step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_2F (der_four, x, y_three, no_of_eq_ini)
	
	IF (y_three (3) == rho1_min) THEN
		der_four (1) = 0.0E0_DP
		der_four (2) = 0.0E0_DP
	END IF
	IF (y_three (6) == rho2_min) THEN
		der_four (4) = 0.0E0_DP
		der_four (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_four (i) = y_zero (i) + (3.0E0_DP / 1.6E1_DP) * dxmore * der_one (i) &
			+ (9.0E0_DP / 1.6E1_DP) * dxmore * der_four (i) 
	END DO

	IF (j < r_grid_more1) THEN
		ini_p1 = y_three (2)
		IF (ini_p1 <= p1_a) THEN
			r_grid_more1 = j
			ini_p1 = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
	END IF
	IF (j < r_grid_more2) THEN
		ini_p2 = y_three (5)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eostable1, eosline1, 1)

	if(helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, ini_rho2, ini_rho2)
        elseif (ccsneos_flag == 1) then
                call GETRHO_CCSNPTOR (ini_p2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p1 == p1_a) THEN
		y_four (3) = rho1_min
	ELSE
		y_four (3) = ini_rho1
	END IF
	If (ini_p2 == p2_a) THEN
		y_four(6) = rho2_min
	ELSE
		y_four (6) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (3.0E0_DP / 4.0E0_DP) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the fifth step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_2F (der_five, x, y_four, no_of_eq_ini)
	
	IF (y_four (3) == rho1_min) THEN
		der_five (1) = 0.0E0_DP
		der_five (2) = 0.0E0_DP
	END IF
	IF (y_four (6) == rho2_min) THEN
		der_five (4) = 0.0E0_DP
		der_five (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_five (i) = y_zero (i) - (3.0E0_DP / 7.0E0_DP) * dxmore * der_one (i) &
			+ (2.0E0_DP / 7.0E0_DP) * dxmore * der_two (i) &
			+ (1.2E1_DP / 7.0E0_DP) * dxmore * der_three (i) &
			- (1.2E1_DP / 7.0E0_DP) * dxmore * der_four (i) &
			+ (8.0E0_DP / 7.0E0_DP) * dxmore * der_five (i)
	END DO

	IF (j < r_grid_more1) THEN
		ini_p1 = y_three (2)
		IF (ini_p1 <= p1_a) THEN
			r_grid_more1 = j
			ini_p1 = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
	END IF
	IF (j < r_grid_more2) THEN
		ini_p2 = y_three (5)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eostable1, eosline1, 1)

	if(helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, ini_rho2, ini_rho2)
        elseif (ccsneos_flag == 1) then
                call GETRHO_CCSNPTOR (ini_p2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p1 == p1_a) THEN
		y_five (3) = rho1_min
	ELSE
		y_five (3) = ini_rho1
	END IF
	If (ini_p2 == p2_a) THEN
		y_five (6) = rho2_min
	ELSE
		y_five (6) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the final step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_2F (der_new, x, y_five, no_of_eq_ini)

	IF (y_five (3) == rho1_min) THEN
		der_new (1) = 0.0E0_DP
		der_new (2) = 0.0E0_DP
	END IF
	IF (y_five (6) == rho2_min) THEN
		der_new (4) = 0.0E0_DP
		der_new (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y (i, j + 1) = y_zero (i) + (7.0E0_DP * der_one (i) + 3.2E1_DP * der_three (i) & 
				+ 1.2E1_DP * der_four (i) + 3.2E1_DP * der_five (i) & 
				+ 7.0E0_DP * der_new (i)) * dxmore / 9.0E1_DP
	END DO

	IF (j < r_grid_more1) THEN
		ini_p1 = y (2, j + 1)
		IF (ini_p1 <= p1_a) THEN
			r_grid_more1 = j
			ini_p1 = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
	END IF
	IF (j < r_grid_more2) THEN
		ini_p2 = y (5, j + 1)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eostable1, eosline1, 1)

	if(helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, ini_rho2, ini_rho2)
        elseif (ccsneos_flag == 1) then
                call GETRHO_CCSNPTOR (ini_p2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p1 == p1_a) THEN
		y (3, j + 1) = rho1_min
	ELSE
		y (3, j + 1) = ini_rho1
	END IF
	If (ini_p2 == p2_a) THEN
		y (6, j + 1) = rho2_min
	ELSE
		y (6, j + 1) = ini_rho2
	END IF

	if(helmeos_flag == 1) then
		call HELMEOS_RtoP(y(6, j+1), temp2_ini, abar2_ini, zbar2_ini, ye2_ini, dummy, dpdrho2_temp(j+1), dummy, dummy, dummy, dummy)
	END IF
	
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We assign atmospheric density !
DO j = r_grid_more1 + 1, length_morestep
	y (3, j) = rho1_min
END DO
DO j = r_grid_more2 + 1, length_morestep
	y (6, j) = rho2_min
END DO

! Assign density !
den1 = y (3, :)
den2 = y (6, :)

! Now we assign the value of density and mass to the original array !
DO j = 1, length_step_2
	IF(r2(j) > rtemp(length_morestep)) THEN
		rho2(j) = rho2_min
	ELSE	
		DO  k = 1, length_morestep
			IF(rtemp(k) == r2(j)) THEN
				rho2(j) = den2(k)
				EXIT
			ELSEIF(rtemp(k) > r2(j)) THEN
				CALL AKIMA(rtemp(k-3), rtemp(k-2), rtemp(k-1), rtemp(k), rtemp(k+1), rtemp(k+2), & 
				den2(k-3), den2(k-2), den2(k-1), den2(k), den2(k+1), den2(k+2), r2(j), rho2(j))
				EXIT
			END IF
		END DO
	END IF
END DO
DO j = 1, length_step_1
	IF(r1(j) > rtemp(length_morestep)) THEN
		rho1(j) = rho1_min
	ELSE
		DO  k = 1, length_morestep
			IF(rtemp(k) > r1(j)) THEN
				rho1(j) = den1(k)
				EXIT
			ELSEIF(rtemp(k) > r1(j)) THEN
				CALL AKIMA(rtemp(k-3), rtemp(k-2), rtemp(k-1), rtemp(k), rtemp(k+1), rtemp(k+2), & 
				den1(k-3), den1(k-2), den1(k-1), den1(k), den1(k+1), den1(k+2), r1(j), rho1(j))
				EXIT
			END IF
		END DO
	END IF
END DO

! Check the density !
DO j = 1, length_step_2
	IF(rho2(j) < rho2_a) THEN
		rho2(j) = rho2_a
	END IF
END DO
DO j = 1, length_step_1
	IF(rho1(j) < rho1_a) THEN
		rho1(j) = rho1_a
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
m_cell = 0.0D0

! Assign mass faced coordinate !
DO j = 1, length_step_2
	DO k = 1, j
		If (rho2(k) > rho2_a) THEN
			m_cell (j) = m_cell(j) + vol2(k) * rho2(k)
		ELSE
			CYCLE
		END IF
	END DO
ENDDO

! Assign mass centered coordinate by interpolation !
DO j = 1, length_step_2
	m_r(j) = 0.5D0*(m_cell (j) + m_cell (j-1))
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We copy the density to ghost shell !
CALL BOUNDARY1D_DM (rho1, even)
CALL BOUNDARY1D_NM (rho2, even)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We choose which table is needed to read according to the chosen EOS !
IF(helmeos_flag == 1) THEN

	! Now of course, we need boundary condition for composition !
	! We first assign the appropriate value to the arrays !
	Do i = 1, length_step_2
		if (rho2(i) > 0.0D0) then	
			abar2(i) = abar2_ini
			zbar2(i) = zbar2_ini
			xiso(:,i) = xiso_ini(:)
			temp2(i) = temp2_ini
			ye2(i) = ye2_ini
		END IF
	END DO

	! Assign initial electron fraction !
	ye_ini = ye2_ini

	! Copy the arrays to the ghost shell !
	CALL BOUNDARY2D_X (xiso)
	CALL BOUNDARY1D_NM (abar2, even)       
	CALL BOUNDARY1D_NM (zbar2, even)
	CALL BOUNDARY1D_NM (temp2, even)
	CALL BOUNDARY1D_NM (ye2, even) 

END IF 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine solve for the initial hydrostatic equilibrium star !
! assuming a one fluid formalism. We assume the newtonian gravity    !
! and is solving for the initial density profile using RK-Butcher    !
! fifth order method   						     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETRHO_1F
USE DEFINITION
USE nuclear_module
IMPLICIT NONE

! The number of hydrostatic equation and the improved accuracy !
INTEGER, PARAMETER :: no_of_eq_ini = 3, more = 10 ** ini_acc

! The extra array arising from the extra accuracy !
INTEGER, PARAMETER :: length_morestep = length_step_2 * more

! The smaller length step arising from the extra accuracy !
REAL (DP), PARAMETER :: dxmore = dx2_ini / DBLE (more) 

! Dummy variables !
REAL (DP) :: dummy, rho2_min

! Integer parameters !
INTEGER :: i, j, k, r_grid_more2

! The pressure, density, at center, atmosphere and the dummy varaiables used to stored them !
REAL (DP) ::  p2_c, p2_a, x, ini_rho2, ini_p2

! Dummy arrays for NM and DM density to be interpolated !	
REAL (DP), DIMENSION (-4 : length_morestep + 5) :: den2

! Temporal distance arrays !
REAL (DP), DIMENSION (-4 : length_morestep + 5) :: rtemp

! Variables essential in the RK-5 method !
REAL (DP), DIMENSION (1 : no_of_eq_ini) :: y_zero, y_one, y_two, y_three, y_four, y_five
REAL (DP), DIMENSION (1 : no_of_eq_ini) :: der_one, der_two, der_three, der_four, der_five, der_new

! Assigning the variables essential in the RK-5 method !
REAL (DP), DIMENSION (1 : no_of_eq_ini, -4 : length_morestep + 5) :: y		

! This is necessary for any EOS since the deriative of pressure plays a role !
REAL (DP), DIMENSION (-4 : length_morestep + 5) :: dpdrho2_temp

! Quantity related to chemical composition if !
! you are using variable compositions of star !
REAL (DP), DIMENSION (total_ion) :: xiso_ini

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign minmium density !
rho2_min = 1.0D-4*rho2_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We read the EOS table for NM !
IF (fermieosnm_flag == 1) THEN
	OPEN (UNIT = 100, FILE = 'EOS_Table2.eos', STATUS = 'OLD')
	DO i = 1, eoslineno
		READ (100, *) eostable2 (i, 1), eostable2 (i, 2)
	END DO
	CLOSE (100)
	eosline2 = 1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We choose which table is needed to read according to the chosen EOS !
IF(helmeos_flag == 1) THEN

	! We assign initial chemical composition accordingly !
	! If you are using helmeos EOS in construction !
	! It is because the realistic EOS not only depends on !
	! Ye but also Abar and Zbar !	
	xiso_ini = 0.0D0
	xiso_ini(che4) = xhe4_ini
	xiso_ini(cc12) = xc12_ini
	xiso_ini(co16) = xo16_ini
	xiso_ini(cne20) = xne20_ini
	xiso_a = xiso_ini

	! Now convert the composition into mean atomic and mass number !
	CALL private_helmeos_azbar(xiso_ini, abar2_ini, zbar2_ini, ye2_ini)

	! assign atmospheric electron fraction !
	ye2_a = ye2_ini
	
	! Read the helmeos eos table for interpolation !
	CALL read_helm_table()	

	! We convert density at center and atmosphere to the corresponding pressure !
	! Note that at center the deriative of density should be zero !
	CALL HELMEOS_RtoP(rho2_c, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, p2_c, dpdrho2_temp(0), dummy, dummy, dummy, dummy)
	CALL HELMEOS_RtoP(rho2_min, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, p2_a, dummy, dummy, dummy, dummy, dummy)

ELSEIF(ccsneos_flag == 1) THEN
	
	! CCSN EOS !
        CALL GETRHO_CCSNRTOP (p2_c, rho2_c)
        CALL GETRHO_CCSNRTOP (p2_a, rho2_a)

ELSE

	! We convert density at center and atmosphere to the corresponding pressure !
	CALL GETRHO_EOSRTOP (p2_c, rho2_c, gs2, mb2, me2, ye2_old, 2)
	CALL GETRHO_EOSRTOP (p2_a, rho2_min, gs2, mb2, me2, ye2_old, 2)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! assign distance !
DO j = 1, length_morestep 
	rtemp(j) = dxmore*DBLE(j - 1/2)
END DO

! We assign the value of y at center !
y (1, 0) = (4.0D0/3.0D0)*pi_old*rtemp(1)**3*rho2_c
y (2, 0) = p2_c - 0.5D0*y (1, 0)*rho2_c/rtemp(1)
y (3, 0) = rho2_c

! initialize the integer parameter !
r_grid_more2 = length_morestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this is the main structure of RK-5 method !
! this should not be unfamiliar so I would skip most the explaination !
DO j = 0, length_morestep - 1	
	
	! Update the value of x and y !
	DO i = 1, no_of_eq_ini
		y_zero (i) = y (i, j)
	END DO

	x = (DBLE (j) + 0.5D0) * dxmore
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the first step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_1F (der_one, x, y_zero, no_of_eq_ini)
	
	! We artificially make the value to zero to avoid singularity at center !
	IF (j == 0) THEN
		der_one (2) = 0.0E0_DP
	END IF

	! If the density reach atmospheric values, no changes in all the quantity !
	IF (y_zero (3) == rho2_min) THEN
		der_one (1) = 0.0E0_DP
		der_one (2) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_one (i) = y_zero (i) + (1.0E0_DP / 4.0E0_DP) * dxmore * der_one (i)
	END DO

	! We determine whether the pressure reached atmospheric pressure !
	IF (j < r_grid_more2) THEN
		ini_p2 = y_one (2)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	if(helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, ini_rho2, ini_rho2)
        elseif (ccsneos_flag == 1) then
                call GETRHO_CCSNPTOR (ini_p2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF
	
	! Assign density according to pressure !
	If (ini_p2 == p2_a) THEN
		y_one (3) = rho2_min
	ELSE
		y_one (3) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (1.0E0_DP / 4.0E0_DP) * dxmore
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the second step in RK-5 !	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_1F (der_two, x, y_one, no_of_eq_ini)

	IF (y_one (3) == rho2_min) THEN
		der_two (1) = 0.0E0_DP
		der_two (2) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_two (i) = y_zero (i) + (1.0E0_DP / 8.0E0_DP) * dxmore * der_one (i) &
				+ (1.0E0_DP / 8.0E0_DP) * dxmore * der_two (i)
	END DO

	IF (j < r_grid_more2) THEN
		ini_p2 = y_two (2)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	if(helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, ini_rho2, ini_rho2)
        elseif (ccsneos_flag == 1) then
                call GETRHO_CCSNPTOR (ini_p2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p2 == p2_a) THEN
		y_two (3) = rho2_min
	ELSE
		y_two (3) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (1.0E0_DP / 4.0E0_DP) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the third step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_1F (der_three, x, y_two, no_of_eq_ini)

	IF (y_two (3) == rho2_min) THEN
		der_three (1) = 0.0E0_DP
		der_three (2) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_three (i) = y_zero (i) - (1.0E0_DP / 2.0E0_DP) * dxmore * der_two (i) &
			+ dxmore * der_three (i)
	END DO

	IF (j < r_grid_more2) THEN
		ini_p2 = y_three (2)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	if(helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, ini_rho2, ini_rho2)
        elseif (ccsneos_flag == 1) then
                call GETRHO_CCSNPTOR (ini_p2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p2 == p2_a) THEN
		y_three (3) = rho2_min
	ELSE
		y_three (3) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (1.0E0_DP / 2.0E0_DP) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the forth step in RK-5 !	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_1F (der_four, x, y_three, no_of_eq_ini)

	IF (y_three (3) == rho2_min) THEN
		der_four (1) = 0.0E0_DP
		der_four (2) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_four (i) = y_zero (i) + (3.0E0_DP / 1.6E1_DP) * dxmore * der_one (i) &
			+ (9.0E0_DP / 1.6E1_DP) * dxmore * der_four (i) 
	END DO

	IF (j < r_grid_more2) THEN
		ini_p2 = y_two (2)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	if(helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, ini_rho2, ini_rho2)
        elseif (ccsneos_flag == 1) then
                call GETRHO_CCSNPTOR (ini_p2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p2 == p2_a) THEN
		y_four (3) = rho2_min
	ELSE
		y_four (3) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + (3.0E0_DP / 4.0E0_DP) * dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the fifth step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_1F (der_five, x, y_four, no_of_eq_ini)

	IF (y_four (3) == rho2_min) THEN
		der_five (1) = 0.0E0_DP
		der_five (2) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_five (i) = y_zero (i) - (3.0E0_DP / 7.0E0_DP) * dxmore * der_one (i) &
			+ (2.0E0_DP / 7.0E0_DP) * dxmore * der_two (i) &
			+ (1.2E1_DP / 7.0E0_DP) * dxmore * der_three (i) &
			- (1.2E1_DP / 7.0E0_DP) * dxmore * der_four (i) &
			+ (8.0E0_DP / 7.0E0_DP) * dxmore * der_five (i)
	END DO

	IF (j < r_grid_more2) THEN
		ini_p2 = y_three (2)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	if(helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, ini_rho2, ini_rho2)
        elseif (ccsneos_flag == 1) then
                call GETRHO_CCSNPTOR (ini_p2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p2 == p2_a) THEN
		y_five (3) = rho2_min
	ELSE
		y_five (3) = ini_rho2
	END IF

	x = (DBLE (j) + 0.5D0) * dxmore + dxmore

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the final step in RK-5 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL INI_DER_1F (der_new, x, y_five, no_of_eq_ini)

	IF (y_five (3) == rho2_min) THEN
		der_new (1) = 0.0E0_DP
		der_new (2) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y (i, j + 1) = y_zero (i) + (7.0E0_DP * der_one (i) + 3.2E1_DP * der_three (i) & 
				+ 1.2E1_DP * der_four (i) + 3.2E1_DP * der_five (i) & 
				+ 7.0E0_DP * der_new (i)) * dxmore / 9.0E1_DP
	END DO

	IF (j < r_grid_more2) THEN
		ini_p2 = y (2, j + 1)
		IF (ini_p2 <= p2_a) THEN
			r_grid_more2 = j
			ini_p2 = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
	END IF

	if(helmeos_flag == 1) then
		call HELMEOS_PtOR(ini_p2, temp2_ini, abar2_ini, zbar2_ini, ye2_ini, ini_rho2, ini_rho2)
        elseif (ccsneos_flag == 1) then
                call GETRHO_CCSNPTOR (ini_p2, ini_rho2)
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eostable2, eosline2, 2)
	ENDIF

	If (ini_p2 == p2_a) THEN
		y (3, j + 1) = rho2_min
	ELSE
		y (3, j + 1) = ini_rho2
	END IF

	if(helmeos_flag == 1) then
		call HELMEOS_RtoP(y(3, j+1), temp2_ini, abar2_ini, zbar2_ini, ye2_ini, dummy, dpdrho2_temp(j+1), dummy, dummy, dummy, dummy)
	END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We assign atmospheric density !
DO j = r_grid_more2 + 1, length_morestep
	y (3, j) = rho2_min
END DO

! Assign density !
den2 = y (3, :)

! Now we assign the value of density and mass to the original array !
DO j = 1, length_step_2
	IF(r2(j) > rtemp(length_morestep)) THEN
		rho2(j) = rho2_min
	ELSE
		DO  k = 1, length_morestep
			IF(rtemp(k) == r2(j)) THEN
				rho2(j) = den2(k)
				EXIT
			ELSEIF(rtemp(k) > r2(j)) THEN
				CALL AKIMA(rtemp(k-3), rtemp(k-2), rtemp(k-1), rtemp(k), rtemp(k+1), rtemp(k+2), & 
				den2(k-3), den2(k-2), den2(k-1), den2(k), den2(k+1), den2(k+2), r2(j), rho2(j))
				EXIT
			END IF
		END DO
	END IF
END DO

! Check the density !
DO j = 1, length_step_2
	IF(rho2(j) < rho2_a) THEN
		rho2(j) = rho2_a
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
m_cell = 0.0D0

! Assign mass faced coordinate !
DO j = 1, length_step_2
	DO k = 1, j
		If (rho2(k) > rho2_a) THEN
			m_cell (j) = m_cell(j) + vol2(k) * rho2(k)
		ELSE
			CYCLE
		END IF
	END DO
ENDDO

! Assign mass centered coordinate by interpolation !
DO j = 1, length_step_2
	m_r(j) = 0.5D0*(m_cell (j) + m_cell (j-1))
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We copy the density to ghost shell !
CALL BOUNDARY1D_NM (rho2, even)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We choose which table is needed to read according to the chosen EOS !
IF(helmeos_flag == 1) THEN

	! Now of course, we need boundary condition for composition !
	! We first assign the appropriate value to the arrays !
	! Only if user wants a variable composition of the star !
	! Now of course, we need boundary condition for composition !
	! We first assign the appropriate value to the arrays !
	Do i = 1, length_step_2
		if (rho2(i) > 0.0D0) then	
			abar2(i) = abar2_ini
			zbar2(i) = zbar2_ini
			xiso(:,i) = xiso_ini(:)
			temp2(i) = temp2_ini
			ye2(i) = ye2_ini
		END IF
	END DO

	! Assign initial electron fraction !
	ye_ini = ye2_ini

	! Copy the arrays to the ghost shell !
	CALL BOUNDARY2D_X (xiso)
	CALL BOUNDARY1D_NM (abar2, even)       
	CALL BOUNDARY1D_NM (zbar2, even)
	CALL BOUNDARY1D_NM (temp2, even) 
	CALL BOUNDARY1D_NM (ye2, even)  

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE