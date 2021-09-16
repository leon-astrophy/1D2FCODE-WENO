!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the parameter files containing constants necessary for the hydro simulation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define double precision !
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND (15, 307)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Conversion between units !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Physical constants to be as one !
REAL (DP), PARAMETER :: gconst = 6.67430D-8
REAL (DP), PARAMETER :: clight = 2.99792458D10
REAL (DP), PARAMETER :: solar = 1.98847D33

! Conversion between units !
REAL (DP), PARAMETER :: length = (clight**2)/(solar*gconst)
REAL (DP), PARAMETER :: mass = (1.0D0/solar)
REAL (DP), PARAMETER :: time = (clight**3)/(solar*gconst)

! Derived conversion !
REAL (DP), PARAMETER :: density = (mass/length**3)
REAL (DP), PARAMETER :: epsilon = (1.0D0/clight**2)
REAL (DP), PARAMETER :: h_bar = (1.054571817D-27)*(length**2*mass/time)
REAL (DP), PARAMETER :: pressure = density*epsilon
REAL (DP), PARAMETER :: qdot = pressure/time

! We use GK as default temperature unit !
REAL (DP), PARAMETER :: temperature = 1.0D-9

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section is related to the basic parameters governing the simulation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Boundary condition for the simulation box !
! The first one is inner boundary and second is outer !
! 0 = periodic, 1 = reflecting, 2 = outgoing !
INTEGER, PARAMETER :: boundary_flag(2) = (/1,2/)

! Physical length (dimensionless) of the DM simulation box !
REAL (DP), PARAMETER :: total_length_1 = 6.0E3_DP

! Physical length (dimensionless) of the NM simulation box !
REAL (DP), PARAMETER :: total_length_2 = 6.0E3_DP

! Value of DM spatial grid size dx !
REAL (DP), PARAMETER :: dx1_ini = 6.0E0_DP

! Value of NM spatial grid size dx !
REAL (DP), PARAMETER :: dx2_ini = 6.0E0_DP

! The total number of array stored by each DM variables !
INTEGER, PARAMETER :: length_step_1 = INT (total_length_1 / dx1_ini) 

! The total number of array stored by each NM variables !
INTEGER, PARAMETER :: length_step_2 = INT (total_length_2 / dx2_ini) 

! The minimum number of array to simulate the whole star !
! This would be change accordingly when you switch on Checkstep flag !
INTEGER	:: length_step_part_1 = length_step_1
INTEGER	:: length_step_part_2 = length_step_2

! Value of CFL number that govern the stability condition (dt / dx) !				
REAL (DP) :: cfl = 2.0E-1_DP

! Physical time (dimensionless) that stop the program !
REAL (DP) :: total_time = 4.0D5

! Physical time (dimensionless) interval for each output !	
REAL (DP), PARAMETER :: output_time = 1.0E2_DP

! Physical time (dimensionless) interval for each profile output !
REAL (DP) :: output_time_profile1 = 1.0E1_DP

! Physical time (dimensionless) interval for each profile output !
REAL (DP) :: output_time_profile2 = 3.0E1_DP

! The maximum number of iteration for Rungekutta method !
INTEGER, PARAMETER :: total_time_step = 1000000

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section governs the initial condition of the hydrostatic star ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			

! Initial central density of dark matter component, do not set to zero !
REAL (DP), PARAMETER :: rho1_c = 3.0E8_DP*density

! Initial central density of normal matter component, do not set to zero !			
REAL (DP), PARAMETER :: rho2_c = 3.0E9_DP*density

! The muplicity factor to the dm central density which obtains the atmoshperic density !
REAL (DP), PARAMETER :: rhofac_1 = 1.0E-5_DP

! The muplicity factor to the nm central density which obtains the atmoshperic density !
REAL (DP), PARAMETER :: rhofac_2 = 1.0E-5_DP

! Density of the dark matter at the atmosphere, do not set to zero !
REAL (DP) :: rho1_a = rho1_c * rhofac_1

! Density of the normal matter at the atmosphere, do not set to zero !			
REAL (DP) :: rho2_a = rho2_c * rhofac_2

! Pressure of the dark matter at the atmosphere !
REAL (DP), PARAMETER :: p1a = 0.0E0_DP

! Pressure of the normal matter at the atmosphere !			
REAL (DP), PARAMETER :: p2a = 0.0E0_DP

! Initial velocity of dark matter, it can be set to zero !			
REAL (DP), PARAMETER :: ini_vel1 = 0.0E0_DP

! Initial velocity of normal matter, it can be set to zero !
REAL (DP), PARAMETER :: ini_vel2 = 0.0E0_DP

! Velocity of the dark matter at the atmosphere, it can be set to zero !			
REAL (DP) :: vel1_a = 0.0E0_DP

! Velocity of the normal matter at the atmosphere, it can be set to zero !
REAL (DP) :: vel2_a = 0.0E0_DP             

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This section governs the physics of the EOS of NM or DM    !
! CAUTION : We assumed an ideal completely degenerate fermi  !
! gas EOS. To change the EOS, you need to input the required !
! parameters by yourself. For example, temperature           !     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			

! Baryonic mass for dark matter !
REAL (DP), PARAMETER :: mb1 = 1.78266192D-25*mass

! Baryonic mass for normal matter !			
REAL (DP), PARAMETER :: mb2 = 1.66053906660D-24*mass

! Fermionic mass (dark matter particles) for dark matter !
REAL (DP), PARAMETER :: me1 = 1.78266192D-25*mass

! Fermionic mass (electrons) for normal matter !			
REAL (DP), PARAMETER :: me2 = 9.1093837015D-28*mass

! Dark matter fraction for dark matter !
REAL (DP), PARAMETER :: ye1 = 1.0E0_DP

! Electron fraction for normal matter !
! Note that in principle it can altered !
! Depending on the elements that you choose !
REAL (DP) :: ye2_old = 5.0E-1_DP

! Multiplicity factor for dark matter !	
REAL (DP), PARAMETER :: gs1 = 2.0E0_DP

! Multiplicity factor for normal matter				
REAL (DP), PARAMETER :: gs2 = 2.0E0_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section related to the variables in the levelset surface !
! Tracking module. Select the position of initial levelset 	!
! Surface. Remark : 1 code length = 1.4766839 km		!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Surface where you wanted to enclosed DM mass is !
REAL (DP), PARAMETER :: surf1 = 150.D0*length

! Surface where you wanted to enclosed NM mass is !
REAL (DP), PARAMETER :: surf2 = 150.D0*length

! Surface where you wanted to initialize deflagration !
REAL (DP), PARAMETER :: defla2 = 150.D0*length

! Surface where you wanted to enclosed DM mass is !
REAL (DP), PARAMETER :: masssurf1 = 0.1E0_DP

! Surface where you wanted to enclosed NM mass is !
REAL (DP), PARAMETER :: masssurf2 = 0.1E0_DP

! Surface where you wanted to initialize deflagration !
REAL (DP), PARAMETER :: massdefla2 = 0.01E0_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section governs the chemical composition the hydrostatic star !
! The composition Xi is defined by rhoi/rho, the fraction of density !
! occupied by that elements for each desity at a radial distance     !
! Caution : The composition should sum up to 1, also, once you       !
! assume variables composition, you can only use helmeos	     ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

! This includes !
! 1. Initial He-4 mass fraction
! 2. Initial C-12 mass fraction
! 3. Initial O-16 mass fraction
! 4. Initial metallicity (Assume fixed)
! 5. Initial electron fraction
! 6. Atmospheric electron fraction
REAL (DP), PARAMETER :: xhe4_ini  = 0.00D0                        
REAL (DP), PARAMETER :: xc12_ini  = 0.50D0                       
REAL (DP), PARAMETER :: xo16_ini = 0.50D0                       
REAL (DP), PARAMETER :: xne20_ini  = 0.00D0
REAL (DP), PARAMETER :: metalz  = 0.00D0
REAL (DP) :: ye_ini = 0.5D0 
REAL (DP) :: ye2_a = 0.5D0                        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section is the extra physics related to the !
! Helmeos finite temperature EOS. If you hate it   !
! Please feel free to delete them.                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Inital temperature (assume isothermal) of the star !
REAL (DP), PARAMETER :: temp2_ini = 1.0E-1_DP

! atmospheric temperature of the star !
REAL (DP), PARAMETER :: temp2_a = 1.0E-1_DP

! maximum temperature allowed in the star !
REAL (DP), PARAMETER :: temp2_max = 7.0D1

! minimum tempearture allowed in the star !
REAL (DP), PARAMETER :: temp2_min = 1.0D-4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Limiter of electron fraction			  !
! Note: Setting them according to your EOS choice !
! temp_max = Maximum Ye allowed			  !
! temp_min = Minimum Ye allowed			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL (DP), PARAMETER :: ye_max = 0.5D0
REAL (DP), PARAMETER :: ye_min = 0.2D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section related to parameters that are essential in the spatial discritization !
! Poisson gravity solver and initial star model construction. We do not recommend you !
! to change any of it without solid understanding to the knowledge behind them        !                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Spacial structure: 0 = 1-D Cartesian; 1 = cylindrical symmetric; 2 = spherical symmetric !	
INTEGER, PARAMETER :: sp_dim_i = 2

! Presence of gravity. 1 = present of gravity if sp_dim_i = 2. 0 = absent of gravity if sp_sim_i = 2 !
! CAUTION : w_gravity_i MUST be set to 0 if sp_dim_i = 0 or 1 !
INTEGER, PARAMETER :: w_gravity_i = 1

! Accuracy of the initial date if it is generated !
! (Do NOT change this unless you know what you are doing) !
INTEGER, PARAMETER :: ini_acc = 2

! Parameter govern which internal energy is used !
REAL (DP), PARAMETER :: eta_1 = 1.0D-3
REAL (DP), PARAMETER :: eta_2 = 1.0D-1

! Maximum run time in relaxation of the potential !
! (Do NOT change this unless you know what you are doing) !
INTEGER, PARAMETER :: relax_max = 100000

! Tolerance in relaxation of the potential !
! (Do NOT change this unless you know what you are doing) !			
REAL (DP), PARAMETER :: tolerance = 3.0E-8_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section is some logical flag to determine whether some extra feature should be turned on or not !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine whether to solve the  energy equation     !
! CAUTION: Turn off energy equation only for cold EOS !	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 0 = Choose not to solve NM energy equation !
INTEGER, PARAMETER :: nm_epsilon = 1

! 0 = Choose not to solve NM  equation !
INTEGER, PARAMETER :: dual_energy = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control the reconstruction method and its variants !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This choose to use the reconstruction schemes !
INTEGER, PARAMETER :: TVD_flag = 0
INTEGER, PARAMETER :: WENO_flag = 1
INTEGER, PARAMETER :: MP5_flag = 0
INTEGER, PARAMETER :: PPM_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Choose the desire NM Riemann solvers !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This choose to use the reconstruction schemes !
INTEGER, PARAMETER :: HLLC_flag = 0
INTEGER, PARAMETER :: HLL_flag = 0
INTEGER, PARAMETER :: LF_flag = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control the simulation box size and grid size !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Flag for checking the simulation box !
! 1 = Only simulate the box area with matter !
! other meshes without matter is ignored !
INTEGER :: checkstepdm_flag = 0
INTEGER :: checkstepnm_flag = 1

! Flag for turning on moving-grid !
INTEGER :: movinggriddm_flag = 0
INTEGER :: movinggridnm_flag = 0

! Flag for finding the moving-grid !
INTEGER, PARAMETER :: found_movinggriddm_flag = 0
INTEGER, PARAMETER :: found_movinggridnm_flag = 1

! Flag for fixing the atmoshpere !
INTEGER, PARAMETER :: fixrhodm_flag = 0
INTEGER, PARAMETER :: fixrhonm_flag = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control the atmospheric motion or parameter !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This flag creates a sponge near the atmosphere to ensure a stable star !
INTEGER, PARAMETER :: sponge_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control the initial profile or model !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1 = explosion test !
! Note that the test can not be done for ideal Fermi Gas EOS !
INTEGER, PARAMETER :: sedov_flag = 0

! Hydro Test Model. Choose from 1, 2, 3 ,4, 5, 0 Means not running hydro test !
INTEGER, PARAMETER :: test_model = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control the DM existence and its motion !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The presence of dark matter component !
! 1 = presence, 0 = absent !
INTEGER, PARAMETER ::  DM_flag = 0

! Whether dark matter component are movable !
! 1 = movable, 0 = stationary !
! CAUTION - If you set the absent of DM !
! MAKE SURE to set this flag to 0 !
INTEGER, PARAMETER ::  RUNDM_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control the levelset function for DM and NM !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1 = use level set method to trace the NM !
INTEGER, PARAMETER :: levelset_flag_NM = 0

! 1 = use level set method to trace the DM !
! CAUTION - If you set the absent of DM !
! MAKE SURE to set this flag to 0 !
INTEGER, PARAMETER :: levelset_flag_DM = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control the extra outputed quantity other than hydro variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1 = output quantity related to potential !
INTEGER, PARAMETER :: outputpotential = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control the nuclear reaction and thermonuclear explosion !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1 = allow a propagation of deflagration and nuclear reaction !
INTEGER, PARAMETER :: fusion_flag = 1

! 1 = allow a transition to detonation !
INTEGER, PARAMETER :: deton_flag = 1

! 1 = allow a the three step burning progress !
INTEGER, PARAMETER :: flame_flag = 1

! Flag for allowing 1st step burning input for level set 1 & 2
! 1 = Allow energy input by carbon burning
INTEGER, PARAMETER :: carburn_flag = 1

! Flag for allowing 2nd step burning input for level set 1 & 2
! 1 = Allow energy input by advanced burning
INTEGER, PARAMETER :: advburn_flag = 1  
 
! Flag for allowing final burning input for level set 1 & 2
! 1 = Allow energy input by NSE evolution
INTEGER, PARAMETER :: convert_nse_flag = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for tracer particles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1 = With tracer particle to record the thermodynamics history
INTEGER, PARAMETER :: tracer_flag = 0

! Physical time interval for each tracer profiles
REAL (DP) :: output_PPTtime = 2.5D3				
REAL (DP) :: output_PPTtime_last = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for more frequent output of profile after DDT !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Physical time interval for each ddt profiles
REAL (DP) :: output_ddt = 1.0D4
REAL (DP) :: output_ddt_last = 0.0D0
			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control the neutrino, electron and isotope microphysics !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1 = Switch on the neutrino spectra and emisstivity calculator
INTEGER, PARAMETER :: nuspec_flag = 0

! 1 = allow variable and advection of chemical composition !
INTEGER, PARAMETER :: xisotran_flag = 1

! Flag for allowing the electron fraction to be transported !
! 1 = Allow Ye to be advected
INTEGER, PARAMETER :: etran_flag = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control the equation of state for DM and NM !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1 = use Fermi Gas EOS For DM !
INTEGER, PARAMETER :: fermieosdm_flag = 1

! 1 = use Fermi Gas EOS For NM !
INTEGER, PARAMETER :: fermieosnm_flag = 0

! 1 = use newtonian polytropic EOS For NM !
INTEGER, PARAMETER :: newpolyeosnm_flag = 0

! 1 = use helmholtz EOS For NM !
INTEGER, PARAMETER :: helmeos_flag = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for core collapse supernovae !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1 = use GR potential !
INTEGER, PARAMETER :: gr_flag = 0

! 1 = use EOS for core collapse !
INTEGER, PARAMETER :: ccsneos_flag = 0

! 1 = initiate core collapse !
INTEGER, PARAMETER :: corecollapse_flag = 0
