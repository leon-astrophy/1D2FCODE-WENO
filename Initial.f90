!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine solve the initial star model to be simulate !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INITIAL
USE DEFINITION
USE FLAMETABLE_MODULE
USE NUSPEC_MODULE
USE WENO_MODULE
USE NUCLEAR_MODULE
USE FLAME_MODULE
USE ECAPTABLE_MODULE
USE NU_MODULE
IMPLICIT NONE

! Integer parameter !
INTEGER :: j

! Backup parameter !
n_backup = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize dx !
dx1 = dx1_ini
dx2 = dx2_ini

! Initialize dt !
dt = min(cfl * dx1, cfl * dx2)

! Initialize all global quantity !
mass1 = 0.0D0
mass2 = 0.0D0
energy1 = 0.0D0
energy2 = 0.0D0

! Initialize !
energy_input = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) 'We shall now build the initial model for simulation'
WRITE (*,*) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We need to get the polytropic index !
WRITE(*,*) 'Get Polytropic Index'
CALL GETPOLY

! Use fermigas eos according to users will !
IF (fermieosnm_flag == 1 .OR. fermieosdm_flag == 1) THEN
	WRITE (*,*) 
	WRITE (*,*) 'We use exact fermi gas eos'
	WRITE (*,*) 
	WRITE (*,*) 'Build EOS table'
	WRITE (*,*) 
	CALL EOSTABLE
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we initialize, or allocate, other variables !
! We build arrays for hydro variables !
WRITE(*,*) 'Build hydro variables'
WRITE (*,*) 
CALL BUILDHYDRO

! We build the variables neccessary for system of hyperbolic equation !
WRITE(*,*) 'Build WENO variables'
WRITE (*,*) 
call BUILDWENO

! We build the distance and volume !
WRITE(*,*) 'Get distance volume variables'
WRITE (*,*) 
call GetGrid_NM
IF(DM_flag == 1) THEN
	call GetGrid_DM
END IF

! We build variables neccessary for the chemical compositon !
If(helmeos_flag == 1) THEN
	WRITE(*,*) 'Build nuclear composition variables'
	WRITE (*,*) 
	call buildnuclear
END IF

! We build variables for neutrino emisstivity !
If(nuspec_flag == 1) THEN
	WRITE(*,*) 'Build neutrino emisstivity variables'
	WRITE (*,*) 
	call BUILDNU
END IF

! We build variables neccessary for the deflagration !
IF(flame_flag == 1) THEN
	WRITE(*,*) 'Build flame variables'
	CALL BUILDFLAME
	WRITE(*,*)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize all the variables, arrays that related to isotopes composition !
! This is needed before solving for the hydrostatic star !
IF(helmeos_flag == 1) THEN
	WRITE(*,*) 'Initialize isotopes composition variables'
	CALL InitializeISO
	WRITE (*,*)
END IF

! Nse table !
IF(convert_nse_flag == 1) THEN
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!WRITE(*,*) 'Read NSE Table'
	!CALL read_nse_table
	!WRITE (*,*)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	WRITE(*,*) 'Read Electron Capture Table'
	CALL readEcapRate
	WRITE (*,*)
END IF

! Flame table !
If(flame_flag == 1) THEN
	WRITE(*,*) 'Read Flame And Deton Table'
	CALL read_flame_table
	CALL read_deton_table
	WRITE (*,*)
END IF

! Neutrino emisstivity !
IF(nuspec_flag == 1) THEN
   	WRITE(*,*) 'Read neutirno coefficient'
	WRITE (*,*)
	CALL READCOFF
   	WRITE(*,*) 'Done reading coefficient'
   	WRITE(*,*)
	WRITE(*,*) 'Initialize neutrino spectra variables'
	CALL GETNUSPEC
	WRITE(*,*)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

If(test_model == 0) THEN

	! Solve for the initial density, velocity and epsilon !
	! We determine whether the user wants a 1F or 2F star !
	WRITE(*,*) 'Build a real star'
	WRITE (*,*) 
	If (DM_flag == 1) THEN
		CALL GETRHO_2F
	ELSE
		CALL GETRHO_1F
	END IF
	
	! solve for the initial velocity and epsilon !
	CALL GETVEL

	! Override polytropic index !
	IF(newpolyeosnm_flag == 0) THEN
		gamma2 = 1.3334D0
	END IF
	
	! For core collapse supernovae !
	IF (corecollapse_flag == 1) THEN
		gammac_1 = 1.32E0_DP
		Ec_1 = 1.46925E15_DP * pressure / density ** gammac_1
		kc_1 = (gammac_1 - 1.0D0) * Ec_1
		Ec_2 = (gammac_1 - 1.0D0) / (gammac_2 - 1.0D0) * Ec_1 * rhoc_b ** (gammac_1 - gammac_2)
		Ec_3 = (gammac_2 - gammac_1) / (gammac_2 - 1.0D0) * Ec_1 * rhoc_b ** (gammac_1 - 1.0D0)
		kc_2 = Ec_2 * (gammac_2 - 1.0D0)
	END IF

	! Get internal energy !
	CALL GETEPSILON

	! Override polytropic index !
	IF(newpolyeosnm_flag == 0) THEN
		gamma2 = (4.0D0/3.0D0)
	END IF

	! We initialize the deflagration !
	IF (flame_flag == 1) THEN
		CALL GETFLAME
	
		! For central ignited model no detonation !
		founddeton_flag = 0
	END IF

	! Check whether the density reached atmospheric density !
	CALL CHECKRHO

	! We assign the abar, zbar because there are some changes !
	! Due to the ignition or initialization of deflagration !
	IF (xisotran_flag == 1) THEN
		CALL find_AZbar()
	END IF

ELSE

	! We initialize a hydrotest !
	CALL GETRHO_SIMPLE

	! We do not need a atmosphere !
	rho1_a = 0.0E0_DP
	rho2_a = 0.0E0_DP

END IF	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the mass !
CALL FINDMASS

! Find the NM radius !
Do j = 1, length_step_2
	If (rho2 (j) == rho2_a) THEN
		r_grid2 = j - 1
		EXIT
	END IF
END DO

! Find the DM radius !
If (DM_flag == 1) THEN
	Do j = 1, length_step_1
		If (rho1 (j) == rho1_a) THEN
			r_grid1 = j - 1
			EXIT
		END IF
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Report some stellar parameters !
WRITE(*,*) 'Report initial model parameters'
WRITE(*,*)
WRITE(*,*) '--------------------'
WRITE(*,*) 'Stellar parameters :'
WRITE(*,*) '--------------------'
WRITE(*,*) 'NM Central Density ', rho2_c, ' Code Unit'
WRITE(*,*) 'NM Mass ', mass2, ' Solar Mass'
WRITE(*,*) 'NM Radius ', r2(r_grid2), ' Code Unit'
WRITE(*,*) 'NM Radius ', (r2(r_grid2)/length/1.0D3), ' km'
If (DM_flag == 1) Then
	WRITE(*,*) 'There is a DM component!'
	WRITE(*,*) 'DM Central Density ', rho1_c, ' Code Unit'
	WRITE(*,*) 'DM Mass ', mass1, ' Solar Mass'
	WRITE(*,*) 'DM Radius ', r1(r_grid1), ' Code Unit'
	WRITE(*,*) 'DM Radius ', (r1(r_grid1)/length/1.0D3), ' km'
END IF
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) 'Finish construct initial model'
WRITE (*,*) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign the initial energy input from deflagration/detonation
IF (fusion_flag == 1) THEN
		
	! Now we initialize the energy input due to deflagration !
	CALL FLAME_INI

endif

! We check the density again because there may be some change in density !
CALL CHECKRHO

! We find again the azbar, xiso after checkrho because there may be some change !
IF (xisotran_flag == 1) THEN
	CALL FIND_AZBAR()
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Miscellaneous setting for nse flag !
If (convert_nse_flag == 1) THEN
	nse_flag = 0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We see whether the user wants a sedov test !
! Choose the test that you want !
IF (sedov_flag == 1) THEN
	WRITE(*,*) 'Welcome to Sedov Explosions Problem'
	WRITE (*,*) 
	CALL sedov
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update the variables !
WRITE(*,*) 'Do Update'
CALL UPDATE
WRITE (*,*)
WRITE(*,*) 'Done initial update'
WRITE (*,*)  

! Convert to conservative variables !
WRITE(*,*) 'Build conservative variables'
WRITE (*,*) 
CALL FROMRVETOU (u_new1, u_new2)
CALL BACKUPCONS (u_new1, u_old1, u_new2, u_old2)
WRITE(*,*) 'Done building initial conservative variables'
WRITE (*,*)

! Write finish message !
WRITE(*,*) 'Finish initial...'
WRITE (*,*) 

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialize the polytropic index and exponent for !
! Polytropic equation of state 					   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETPOLY
USE DEFINITION 
IMPLICIT NONE

If (newpolyeosnm_flag == 1) THEN

	! assign newtonian polytropic index and exponent for NM !
	k2 = (3.0E0_DP ** (2.0E0_DP/3.0E0_DP) * pi_old ** (4.0E0_DP/3.0E0_DP) * &
		h_bar ** (2.0E0_DP) * ye2_old ** (5.0E0_DP/3.0E0_DP) ) &
		/ (5.0E0_DP * me2 * mb2 ** (5.0E0_DP/3.0E0_DP) )

	gamma2 = (5.0E0_DP/3.0E0_DP)


ELSE If (ccsneos_flag == 1) THEN

	kc_1 = 4.64025719E-1_DP
	gammac_1 = (4.0E0_DP/3.0E0_DP)
	gammac_2 = 2.5E0_DP
	gammac_th = 1.5E0_DP
	rhoc_b = 3.238E-4_DP !5.8284E-4_DP

	Ec_1 = (kc_1 / (gammac_1 - 1.0D0))
	Ec_2 = (gammac_1 - 1.0D0) / (gammac_2 - 1.0D0) * Ec_1 * rhoc_b ** (gammac_1 - gammac_2)
	Ec_3 = (gammac_2 - gammac_1) / (gammac_2 - 1.0D0) * Ec_1 * rhoc_b ** (gammac_1 - 1.0D0)
	kc_2 = Ec_2 * (gammac_2 - 1.0D0)

ELSE

	! assign relativistic polytropic index and exponent for NM !
	k2 = (3.0E0_DP ** (1.0E0_DP/3.0E0_DP) * pi_old ** (2.0E0_DP/3.0E0_DP) * &
		h_bar * ye2_old ** (4.0E0_DP/3.0E0_DP) ) &
		/ (4.0E0_DP * mb2 ** (4.0E0_DP/3.0E0_DP) )

	gamma2 = (4.0E0_DP/3.0E0_DP)

END IF

If (ccsneos_flag == 1) THEN
	WRITE (*,*) '-------------------------------'
	WRITE (*,*) 'report polytropic index for NM:'
	WRITE (*,*) '-------------------------------'
	WRITE (*,*) 'kc_1', kc_1
	WRITE (*,*) 'gammac_1', gammac_1
	WRITE (*,*) 
ELSE
	WRITE (*,*) '-------------------------------'
	WRITE (*,*) 'report polytropic index for NM:'
	WRITE (*,*) '-------------------------------'
	WRITE (*,*) 'k2', k2
	WRITE (*,*) 'gamma2', gamma2
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! assign the polytopic index and exponent for DM if it exist !
If (DM_flag == 1) THEN
	k1 = (3.0E0_DP ** (2.0E0_DP/3.0E0_DP) * pi_old ** (4.0E0_DP/3.0E0_DP) * & 
		h_bar ** (2.0E0_DP) * ye1 ** (5.0E0_DP/3.0E0_DP) ) &
		/ (5.0E0_DP * me1 * mb1 ** (5.0E0_DP/3.0E0_DP) )

	gamma1 = (5.0E0_DP/3.0E0_DP)

	WRITE (*,*) '-------------------------------'
	WRITE (*,*) 'report polytropic index for DM:'
	WRITE (*,*) '-------------------------------'
	WRITE (*,*) 'k1', k1
	WRITE (*,*) 'gamma1', gamma1
END IF

END SUBROUTINE