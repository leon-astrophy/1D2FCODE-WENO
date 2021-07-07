!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains all the tools for reading the !
! deflagration and detonation energy release and     !
! chemical composition of the resultant		     !
! Written by Leung Shing Chi in 2016		     !
!						     !	
! To add your table, do the following:		     !
! 1. Set the density/temp limit of your table	     !
! 2. Give the size of your table		     !
! 3. Give a subroutine to read your table	     !
! 4. Give a subroutine to do the interpolation	     !
! 5. Patch it to Helmeos_module for the energy input !
! 6. Set up in Initial for the initialization	     !
! 7. Set up flag for your table in WENO2D_module     !
!  						     !
! This module contains the following subroutines:    !
! 1. subroutine read_flame_table		     !
! 2. subroutine read_deton_table		     !
! 3. subroutine readtable_flamestate		     !
! 4. subroutine readtable_flameenergy		     !
! 5. subroutine readtable_detonenergy		     !
! 6. subroutine readtable_detonvel		     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE FLAMETABLE_MODULE
USE NUCLEAR_MODULE
USE DEFINITION
IMPLICIT NONE

! New set of table
! Size of deflagration table
INTEGER, PARAMETER :: temp_rowno1 = 0
INTEGER, PARAMETER :: den_rowno1 = 30

! Deflagration quantities
REAL (DP), DIMENSION(0:den_rowno1, 0:temp_rowno1):: ashtable_eps
REAL (DP), DIMENSION(0:den_rowno1, 0:temp_rowno1, 9):: ashtable_state
REAL (DP), DIMENSION(0:den_rowno1, 0:temp_rowno1, 1:total_ion):: ashtable_xiso

! Size of detonation table
INTEGER, PARAMETER :: temp_rowno2 = 0
INTEGER, PARAMETER :: den_rowno2 = 30

! Detonation quantities
REAL (DP), DIMENSION(0:den_rowno2, 0:temp_rowno2):: dettable_eps, dettable_vel
REAL (DP), DIMENSION(0:den_rowno2, 0:temp_rowno2, 9):: dettable_state
REAL (DP), DIMENSION(0:den_rowno2, 0:temp_rowno2, 1:total_ion):: dettable_xiso

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine read in the input deflagration ash table   !	
! Written by Leung Shing Chi in 2016			     !
! The ashtable contains the following:			     !
! The preheated state (rho, temp and epsilon)		     !
! The postheated state (rho, temp and epsilon)		     !
! The post-burn state (rho, temp, epsilon)		     !
! The energy release by comparing different stages	     !
! The final chemical composition			     !
! Notice that 50% C and 50% O is assumed as fuel composition !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_FLAME_TABLE
USE NUCLEAR_MODULE
USE DEFINITION
IMPLICIT NONE
  
! Integer variables !
integer :: i, j, k 
  
! Dummy variable !
REAL (DP) :: dummy

! open the deflagration table !
open(unit=500, file='./Library/def_table_7iso_2.0e9.dat',action='read')

do i = 0, den_rowno1
	do j = 0, temp_rowno1

		! read the deflagration table !
		read(500,*) (ashtable_state(i,j,k), k=1,9), ashtable_eps(i,j), dummy, dummy, dummy, &
		(ashtable_xiso(i,j,k), k = 1, total_ion)

		! Convert it to code unit !
		ashtable_eps(i,j) = ashtable_eps(i,j) / 9.0D20
		if(ashtable_state(i,j,1) < 5.0D7) then
			ashtable_eps(i,j) = 3.0871D-4
			ashtable_xiso(i,j,:) = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0/)
		endif
		ashtable_state(i,j,1) = ashtable_state(i,j,1)*density
		ashtable_state(i,j,2) = ashtable_state(i,j,2)*temperature
		ashtable_state(i,j,3) = ashtable_state(i,j,3)*epsilon
		ashtable_state(i,j,4) = ashtable_state(i,j,4)*density 
		ashtable_state(i,j,5) = ashtable_state(i,j,5)*temperature  
		ashtable_state(i,j,6) = ashtable_state(i,j,6)*epsilon
		ashtable_state(i,j,7) = ashtable_state(i,j,7)*density
		ashtable_state(i,j,8) = ashtable_state(i,j,8)*temperature
		ashtable_state(i,j,9) = ashtable_state(i,j,9)*epsilon

	enddo
enddo

close(500)

end subroutine read_flame_table

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine does the interpolation to get the ash density !
! by deflagration						!
! Written by Leung Shing Chi in 2016				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readtable_flamestate(temp_in, rho_in, rho_out, temp_out, eps_out)
USE NUCLEAR_MODULE
USE DEFINITION
IMPLICIT NONE

! The position of the target grid !
integer :: i_table, j_table

! Input data !  
real (DP), INTENT(IN) :: temp_in, rho_in

! Intermediate data !      
real (DP) :: di, dj

! Output data !
real (DP), INTENT(OUT) :: rho_out, temp_out, eps_out

! To a good approximation, there is no temperature dependence
j_table = 0
dj = 0.0D0

! But there is density dependence
i_table = INT((LOG10(rho_in/density) - 7.0D0) / 0.1D0)
di = ((LOG10(rho_in/density) - 7.0D0) - DBLE(i_table) * 0.1D0) / 0.1D0

! Do the bilinear interpolation for the ash density
rho_out = ashtable_state(i_table, j_table,7) + &
	di * (ashtable_state(i_table + 1, j_table,7) - ashtable_state(i_table, j_table,7))

temp_out = ashtable_state(i_table, j_table,8) + &
	      di * (ashtable_state(i_table + 1, j_table,8) - ashtable_state(i_table, j_table,8))

eps_out = ashtable_state(i_table, j_table,9) + &
	     di * (ashtable_state(i_table + 1, j_table,9) - ashtable_state(i_table, j_table,9))

end subroutine readtable_flamestate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine does the interpolation to get the energy release !
! by deflagration Written by Leung Shing Chi in 2016		   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readtable_flameenergy(temp_in, rho_in, fusionrate_out, x_isotope_out)
USE NUCLEAR_MODULE
USE DEFINITION
IMPLICIT NONE

! The position of the target grid !
integer :: i_table, j_table

! Input data !
real (DP), INTENT(IN) :: temp_in, rho_in 

! Intermediate data !
real (DP) :: di, dj
 
! Output data !
real (DP), INTENT(OUT) :: fusionrate_out
real (DP), INTENT(OUT), dimension(total_ion) :: x_isotope_out

! To a good approximation, there is no temperature dependence !
j_table = 0
dj = 0.0D0

! But there is temperature dependence !
i_table = INT((LOG10(rho_in/density) - 7.0D0) / 0.1D0) 
di = ((LOG10(rho_in/density) - 7.0D0) - DBLE(i_table) * 0.1D0) / 0.1D0

! Do the bilinear interpolation for the reaction energy
fusionrate_out = ashtable_eps(i_table, j_table) + &
	di * (ashtable_eps(i_table + 1, j_table) - ashtable_eps(i_table, j_table)) 

! Do the bilinear interpolation for the reaction composition
x_isotope_out = ashtable_xiso(i_table, j_table, :) + &
	di * (ashtable_xiso(i_table + 1, j_table, :) - ashtable_xiso(i_table, j_table,:))

end subroutine readtable_flameenergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine read in the input detonation ash table     !    
! Written by Leung Shing Chi in 2016 			     !
! The dettable contains the following:			     !
! The preheated state (rho, temp and epsilon)		     !
! The postheated state (rho, temp and epsilon)		     !
! The post-burn state (rho, temp, epsilon and xiso)	     !
! The energy release at different stages		     !
! The velocity of the detonation			     !
! The final chemical composition			     !
! Notice that 50% C and 50% O is assumed as fuel composition !   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_DETON_TABLE
USE NUCLEAR_MODULE
USE DEFINITION
IMPLICIT NONE

! Dummy variables  
integer :: i, j, k

! Dummy variables
REAL (DP) :: dummy

! For deton data
open(unit=500, file='./Library/det_table_7iso_2.5e9.dat',action='read')

! Now do the loop !
do i = 0, den_rowno2
	do j = 0, temp_rowno2

	! Read the table
	read(500,*) (dettable_state(i,j,k), k=1,9), dettable_vel(i,j), dettable_eps(i,j), &
		      dummy, dummy, dummy, (dettable_xiso(i,j,k), k = 1, total_ion)

	! Convert to code unit
	dettable_vel(i,j) = dettable_vel(i,j)/clight
	dettable_eps(i,j) = dettable_eps(i,j)*epsilon
	if(dettable_state(i,j,1) < 5.0D7) then
		dettable_eps(i,j) = 3.0871D-4        
		dettable_xiso(i,j,:) = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0/)
	endif
	dettable_state(i,j,1) = dettable_state(i,j,1)*density
	dettable_state(i,j,2) = dettable_state(i,j,2)*temperature
	dettable_state(i,j,3) = dettable_state(i,j,3)*epsilon
	dettable_state(i,j,4) = dettable_state(i,j,4)*density
	dettable_state(i,j,5) = dettable_state(i,j,5)*temperature
	dettable_state(i,j,6) = dettable_state(i,j,6)*epsilon
	dettable_state(i,j,7) = dettable_state(i,j,7)*density
	dettable_state(i,j,8) = dettable_state(i,j,8)*temperature
	dettable_state(i,j,9) = dettable_state(i,j,9)*epsilon

	enddo
enddo

close(500)

end subroutine read_deton_table

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine interpolates the table and finds the energy by detonation !						     
! Written by Leung Shing Chi in 2016					    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readtable_detonenergy(temp_in, rho_in, fusionrate_out, x_isotope_out)
USE NUCLEAR_MODULE
USE DEFINITION
IMPLICIT NONE

! Position of the target grid
integer :: i_table, j_table

! Input data
real (DP) :: temp_in, rho_in

! Intermediate data
real (DP) :: di, dj

! Output data    
real (DP) ::fusionrate_out
real (DP), dimension(total_ion) :: x_isotope_out
     
! Do the check
if(rho_in > 1.0D6*density .and. rho_in < 1.0D9*density) then

	! To a good approximation, there is no temperature dependence 
	j_table = 0
	dj = 0.0D0
   
	! But there is density dependence
	i_table = INT((LOG10(rho_in/density) - 6.0D0) / 0.1D0)
	di = ((LOG10(rho_in/density) - 6.0D0) - DBLE(i_table) * 0.1D0) / 0.1D0

	! Interpolate to find the reaction energy
	fusionrate_out = dettable_eps(i_table, j_table) + &
			di * (1.0D0 - dj) * (dettable_eps(i_table + 1, j_table) - dettable_eps(i_table, j_table)) + &
			dj * (1.0D0 - di) * (dettable_eps(i_table, j_table + 1) - dettable_eps(i_table, j_table)) + &
			di * dj * (dettable_eps(i_table + 1, j_table + 1) - dettable_eps(i_table, j_table))

	! Interpolate to find the ash composition
	x_isotope_out = dettable_xiso(i_table, j_table, :) + &
			di * (1.0D0 - dj)*(dettable_xiso(i_table + 1, j_table, :) - dettable_xiso(i_table, j_table,:)) + &
			dj * (1.0D0 - di)*(dettable_xiso(i_table, j_table + 1, :) - dettable_xiso(i_table, j_table,:)) + &
			di * dj * (dettable_xiso(i_table + 1, j_table + 1, :) - dettable_xiso(i_table, j_table,:))

endif

end subroutine readtable_detonenergy 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine does the interpolation to get the front velocity  of detonation !
! Written by Leung Shing Chi in 2016						  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readtable_detonvel(rho_in, vel_out)
USE NUCLEAR_MODULE
USE DEFINITION
IMPLICIT NONE

! Position of the target grid
integer :: i_table, j_table

! Input data
real (DP) :: temp_in, rho_in

! Intermediate data
real (DP) :: di, dj

! Output data 
real (DP) :: vel_out

! To a good approximation, there is no temperature dependence
j_table = 0
dj = 0.0D0
                 
! But there is density dependence
i_table = INT((LOG10(rho_in/density) - 6.0D0) / 0.1D0)
di = ((LOG10(rho_in/density) - 6.0D0) - DBLE(i_table) * 0.1D0) / 0.1D0

! DO the density check
! Note: This may be cancelled because this check
! is done somewhere before entering this subroutine
if(rho_in > 1.0D6*density .and. rho_in < 1.0D9*density) then

	! Do the bulinear Interpolation 
	vel_out = dettable_vel(i_table, j_table) + &
		di * (1.0D0 - dj) * (dettable_vel(i_table + 1, j_table) - dettable_vel(i_table, j_table)) + &
		dj * (1.0D0 - di) * (dettable_vel(i_table, j_table + 1) - dettable_vel(i_table, j_table)) + &
		di * dj * (dettable_vel(i_table + 1, j_table + 1) - dettable_vel(i_table, j_table))

else

	! No detonation propagation when out of the density
	vel_out = 0.0D0

endif

end subroutine readtable_detonvel

end module flametable_module