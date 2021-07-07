!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the definition module contains the key variables used in the hydro code !				  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE DEFINITION
IMPLICIT NONE
SAVE
INCLUDE "Parameter.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Mathematical constants and physical constants !
REAL (DP), PARAMETER :: pi_old = 3.1415926535897932384626433832795E0_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Number of conservative hydro - equation solved by the code !
INTEGER :: no_of_eq

! The min/max number of equation for DM hydro problem !
INTEGER :: imin1, imax1

! The min/max number of equation for NM hydro problem !
INTEGER :: imin2, imax2

! The min/max number of scalar equation for DM hydro problem !
INTEGER :: iminsca1, imaxsca1

! The min/max number of scalar equation for NM hydro problem !
INTEGER :: iminsca2, imaxsca2

! The indexs for conservative equation of NM and DM !
INTEGER :: irho1, ivel1, itau1, ieps1
INTEGER :: irho2, ivel2, itau2, ieps2, iye2

! Boundary condition flag for the primitive variables !
INTEGER, PARAMETER :: even = 0, odd = 1
INTEGER :: foundnan
INTEGER :: bfac(50)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The number of iteration in the RK method !
INTEGER, PARAMETER :: time_step = 100000000 

! The time step dt and step size dx. Note that it is variable !
REAL (DP) :: dt, dx1, dx2

! The global time of simulation !
REAL (DP) :: global_time  

! The lastest output time for discrete data and profile !
REAL (DP) :: last_outputtime, last_outputtime_profile1, last_outputtime_profile2
REAL (DP) :: proindex1, proindex2

! Related to the geometry of the simulation !
REAL (DP), PARAMETER :: sp_dim = DBLE (sp_dim_i)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Integer related to initial NM radius, backup procedure and NM eos table !
! Storing NM and DM EOS table information !
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: eostable1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: eostable2

! Integer related to initial NM radius and DM eos table !
INTEGER :: r_grid1, eosline1

! Line number of EOS table !
INTEGER :: eoslineno = 28800
INTEGER :: r_grid2, n_backup, eosline2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Variables for finding gravitational potential !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: phinm, phidm, phinm_1, phinm_2, phidm_1, phidm_2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: phip_nm, phip_dm, phipface_1, phipface_2

! Variables for effective GR potential !
REAL (DP), ALLOCATABLE, DIMENSION (:) ::  phigr

! DM hydrodynamic variables. I guess they are easy to understand !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: rho1, vel1, epsilon1, p1, dpdrho1, dpdepsilon1, dpdr1, cs1

! NM primitive variables variables. I guess they are easy to understand !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: rho2, vel2, epsilon2, p2, dpdrho2, dpdepsilon2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: epsp2, rhoe2, et2, bige2, dpdr2, temp2, ye2, chempo2, cs2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Variables related to rungekutta time evolution !
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: u_temp1, u_old1, u_new1, u2_dm, u3_dm, l3_dm, l1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: u_temp2, u_old2, u_new2, u2_nm, u3_nm, l3_nm, l2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Constant in front of complete degenerate fermi gas EOS !
REAL (DP) :: a_max1, a_max2
REAL (DP) :: b_max1, b_max2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Other necessary variables for NM !
REAL (DP) :: epsilon2_a, k2, gamma2, dlfmmo2
REAL (DP) :: mass2, energy2, gravenergy2, kinenergy2, intenergy2

! Other necessary variables for DM !
REAL (DP) :: epsilon1_a, k1, gamma1, dlfmmo1
REAL (DP) :: mass1, energy1, gravenergy1, kinenergy1, intenergy1

! Total energy !
REAL (DP) :: total_energy, energy_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! distance and volume array !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: r1, r2, vol1, vol2

! mass array for NM !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: m_r, m_cell

! left and rigth face 1/4 from cell enter !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: r1L, r1R, r1F, r2L, r2R, r2F

! face value and differential volume !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: dvol1, dvol2

! Maximum velocity and radius !
REAL (DP) :: vel1_max, radius1, boundary1
REAL (DP) :: vel2_max, radius2, boundary2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! CCSN !
REAL (DP) :: kc_1, kc_2, Ec_1, Ec_2, Ec_3
REAL (DP) :: gammac_1, gammac_2, gammac_th, rhoc_b
REAL (DP), ALLOCATABLE, DIMENSION (:) :: p_p, p_th, eps_p, eps_th
REAL (DP), ALLOCATABLE, DIMENSION (:) :: dpdrho_p, dpdrho_th, dpdepsilon_p, dpdepsilon_th

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Output control !
LOGICAL :: output_file = .false.

END MODULE DEFINITION