!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialize the system of hyperbolic conservation !
! laws by setting up index for conservative variables              !
! If you add more physics, you  need to be aware whether more 	   !
! equation is needed to be solved 				   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDWENO
USE DEFINITION
USE WENO_MODULE
USE FLAME_MODULE
USE NUCLEAR_MODULE
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dummy integer !
INTEGER :: dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We initialize the no_of_eq to 0 !
no_of_eq = 0

! Initialize imax imin !
imin1 = 0
imax1 = 0
imin2 = 0
imax2 = 0

! Initialize max/min scalar equation !
iminsca1 = 0
iminsca2 = 0
imaxsca1 = 0
imaxsca2 = 0

! We now find out how many conservative equation is needed to be solved !
! We first do the DM part !
! We do the DM part only if the DM is precense and movable !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (RUNDM_flag == 1) THEN

	! DM density conservative equation is needed to be solved !
	! increase no of eq !
	no_of_eq = no_of_eq + 1

	! increase the maximum number of equation needed to be solved !
	imin1 = no_of_eq

	! assign the order of conservative equation accordingly !
	irho1 = no_of_eq

	! boundary factor of density related to boundary condition !
	bfac(no_of_eq) = 1

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Now the procdeure are more or less the same, I will no quote everyone of them !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! DM velocity conservative equation is needed to be solved !
	no_of_eq = no_of_eq + 1
	imax1 = no_of_eq
	ivel1 = no_of_eq 
	bfac(no_of_eq) = -1

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! DM epsilon conservative equation is needed to be solved !
	! Mutted for cold equation of state !
	!no_of_eq = no_of_eq + 1
	!imax1 = no_of_eq
	!itau1 = no_of_eq
	!bfac(no_of_eq) = 1
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Dummy integer !
	dummy = no_of_eq

	! Minimum scalar equation !
	IF(imaxsca1 > 0) THEN
		iminsca1 = dummy + 1
	END IF
	
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we turn on to the NM part !
! So now the NM conservative equation start at order no_of_eq + 1 !
imin2 = no_of_eq + 1

! NM density conservative equation is needed to be solved !
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
irho2 = no_of_eq
bfac(no_of_eq) = 1

! NM velocity conservative equation is needed to be solved !
no_of_eq = no_of_eq + 1
imax2 = no_of_eq
ivel2 = no_of_eq 
bfac(no_of_eq) = -1

! NM epsilon conservative equation is needed to be solved !
IF(nm_epsilon == 1) THEN
	no_of_eq = no_of_eq + 1
	imax2 = no_of_eq
	itau2 = no_of_eq
	bfac(no_of_eq) = 1
END IF

! NM internal energy equation is needed to be solved !
IF(dual_energy == 1) THEN
	no_of_eq = no_of_eq + 1
	imax2 = no_of_eq
	ieps2 = no_of_eq
	bfac(no_of_eq) = 1
END IF

! Dummy integer !
dummy = no_of_eq

! Set up the conservative equation for electron fraction !
IF(Etran_flag == 1) then 
	no_of_eq = no_of_eq + 1 
	imaxsca2 = no_of_eq                 
	imax2 = no_of_eq
	iye2 = no_of_eq
	bfac(no_of_eq) = 1
ENDIF

! We check whether there is a deflagration !
IF(flame_flag == 1) THEN
	no_of_eq = no_of_eq + 1
	imaxsca2 = no_of_eq     
	imax2 = no_of_eq
	iscaG = no_of_eq
	bfac(no_of_eq) = 1
END IF

! Detonation levelset !
If (deton_flag == 1) THEN
	no_of_eq = no_of_eq + 1
	imaxsca2 = no_of_eq     
	imax2 = no_of_eq
	iscaG2 = no_of_eq
	bfac(no_of_eq) = 1
END IF

! Set up the conservative equation for chemical composition !
! CAUTION : If you change the number of elements !
! Please also change this section to include the elements you want !
IF(xisotran_flag == 1) THEN

	! Helium-4
	no_of_eq = no_of_eq + 1
	imaxsca2 = no_of_eq     
	imax2 = no_of_eq
	ihe4 = no_of_eq
	bfac(no_of_eq) = 1

	! Carbon-12
	no_of_eq = no_of_eq + 1
	imaxsca2 = no_of_eq     
	imax2 = no_of_eq
	ic12 = no_of_eq   
	bfac(no_of_eq) = 1
   
	! Oxygen-16
	no_of_eq = no_of_eq + 1
	imaxsca2 = no_of_eq     
	imax2 = no_of_eq
	io16 = no_of_eq   
	bfac(no_of_eq) = 1

	! Neon-20
	no_of_eq = no_of_eq + 1
	imaxsca2 = no_of_eq     
	imax2 = no_of_eq
	ine20 = no_of_eq   
	bfac(no_of_eq) = 1

	! Magnesium-24
	no_of_eq = no_of_eq + 1
	imaxsca2 = no_of_eq     
	imax2 = no_of_eq
	img24 = no_of_eq   
	bfac(no_of_eq) = 1

	! Silicon-28
	no_of_eq = no_of_eq + 1
	imaxsca2 = no_of_eq     
	imax2 = no_of_eq
	isi28 = no_of_eq   
	bfac(no_of_eq) = 1

	! Nickel-56
	no_of_eq = no_of_eq + 1
	imaxsca2 = no_of_eq     
	imax2 = no_of_eq
	ini56 = no_of_eq   
	bfac(no_of_eq) = 1

ENDIF

! Minimum scalar equation !
IF(imaxsca2 > 0) THEN
	iminsca2 = dummy + 1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Write out the quantity !
WRITE (*,*) 'There are ', no_of_eq, 'of equation to be solved'
WRITE (*,*) 'Minimum DM equation index is ', imin1
WRITE (*,*) 'Minimum DM scalar equation index is ', iminsca1
WRITE (*,*) 'Maximum DM equation index is ', imax1
WRITE (*,*) 'Maximum DM scalar equation index is ', imaxsca1
WRITE (*,*) 'Minimum NM equation index is ', imin2
WRITE (*,*) 'Minimum NM scalar equation index is ', iminsca2
WRITE (*,*) 'Maximum NM equation index is ', imax2
WRITE (*,*) 'Maximum NM scalar equation index is ', imaxsca2
WRITE (*,*) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we allocate the u, and l accordingly, they will be used in rungekutta evolution !
ALLOCATE(u_temp2(imin2 : imax2, -4 : length_step_2 + 5))
ALLOCATE(u_old2(imin2 : imax2, -4 : length_step_2 + 5))
ALLOCATE(u_new2(imin2 : imax2, -4 : length_step_2 + 5))
ALLOCATE(u2_nm(imin2 : imax2, -4 : length_step_2 + 5))
ALLOCATE(u3_nm(imin2 : imax2, -4 : length_step_2 + 5))
ALLOCATE(l3_nm(imin2 : imax2, -4 : length_step_2 + 5))
ALLOCATE(l2(imin2 : imax2, -4 : length_step_2 + 5))

! For existence of DM !
ALLOCATE(u_temp1(imin1 : imax1, -4 : length_step_1 + 5))
ALLOCATE(u_old1(imin1 : imax1, -4 : length_step_1 + 5))
ALLOCATE(u_new1(imin1 : imax1, -4 : length_step_1 + 5))
ALLOCATE(u2_dm(imin1 : imax1, -4 : length_step_1 + 5))
ALLOCATE(u3_dm(imin1 : imax1, -4 : length_step_1 + 5))
ALLOCATE(l3_dm(imin1 : imax1, -4 : length_step_1 + 5))
ALLOCATE(l1(imin1 : imax1, -4 : length_step_1 + 5))

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEALLOCATE SSPRK variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DESTROYWENO
USE DEFINITION
IMPLICIT NONE

! Now we deallocate !
DEALLOCATE(u_temp2)
DEALLOCATE(u_old2)
DEALLOCATE(u_new2)
DEALLOCATE(u2_nm)
DEALLOCATE(u3_nm)
DEALLOCATE(l3_nm)
DEALLOCATE(l2)

! For existence of DM !
DEALLOCATE(u_temp1)
DEALLOCATE(u_old1)
DEALLOCATE(u_new1)
DEALLOCATE(u2_dm)
DEALLOCATE(u3_dm)
DEALLOCATE(l3_dm)
DEALLOCATE(l1)

END SUBROUTINE