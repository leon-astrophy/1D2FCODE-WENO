!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the NUCLEAR composition module that used to caluculate         ! 
! everything, such as atomic mass, binding energy, mean atomic number... !
! if you assume a variable chemical composition. Naturally, it should be !
! merged with the helmeos EOS module, but for the code development, I    !
! seperated them out to make things easier                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE NUCLEAR_MODULE
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for isotope composition !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Number of isotope in the star !
INTEGER, PARAMETER :: total_ion = 7

! Threshold mass fractions !
REAL, PARAMETER :: xthres = 1.0D-150

! Initial chemical composition of the star !
! We focus on mean atomic number and mean mass number !
REAL (DP) :: abar2_ini, zbar2_ini, ye2_ini

! Atmospheric composition of the star !
REAL (DP), dimension(total_ion) :: xiso_a

! Labeling the chemical composition !
! We assume only 7 isotope which will change manually !
! They are helium, carbon, oxygen, neon, magnesium, sillicon, nickel56 !
integer :: ihe4, ic12, io16, ine20, img24, isi28, ini56
integer :: che4, cc12, co16, cne20, cmg24, csi28, cni56

! name of the chemical elements !
character(len=5) :: ionam(1:total_ion)

! These are the nuclear properties of the chemical elements !
real (selected_real_kind(15,307)), dimension(1:total_ion) :: aion, zion, bion, nion, mion, wion

! We turn to some important arrays !
! Mean atomic mass !
REAL (DP), allocatable, dimension (:) :: abar2

! Mean atomic number !
real (DP), allocatable, dimension (:) :: zbar2

! Sum of isotope masses !
real (DP), allocatable, dimension (:) :: xmass_sum

! Chemical composition !
REAL (DP), allocatable, DIMENSION (:,:) :: xiso

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section of NSE table 1 !
!!!!!!!!!!!!!!!!!!!!!!!!!!
! Table size !
INTEGER, PARAMETER :: temp_rowno_nse = 70
INTEGER, PARAMETER :: den_rowno_nse = 30

! Binding energy and composition of NSE
REAL (DP), DIMENSION(0:den_rowno_nse, 0:temp_rowno_nse):: nsetable_binde
REAL (DP), DIMENSION(0:den_rowno_nse, 0:temp_rowno_nse, 1:total_ion):: nsetable_xiso

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section of NSE table 2 !
!!!!!!!!!!!!!!!!!!!!!!!!!!
! Table size         
INTEGER, PARAMETER :: ent_rowno_nse2 = 50
INTEGER, PARAMETER :: ye_rowno_nse2 = 26
INTEGER, PARAMETER :: den_rowno_nse2 = 50

! New binding energy and composition of NSE with a larger network
REAL (DP), DIMENSION(0:den_rowno_nse2, 0:ye_rowno_nse2+1, 3):: nsetable2_head
REAL (DP), DIMENSION(0:den_rowno_nse2, 0:ye_rowno_nse2+1, 0:ent_rowno_nse2+1, 6):: nsetable2_binde

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section of NSE table 3 !      
!!!!!!!!!!!!!!!!!!!!!!!!!!                 
! Table size !
INTEGER, PARAMETER :: temp_rowno_nse3 = 48     
INTEGER, PARAMETER :: ye_rowno_nse3 = 122      
INTEGER, PARAMETER :: den_rowno_nse3 = 23

! New binding energy and composition of NSE with a larger network
REAL (DP), DIMENSION(0:den_rowno_nse3+1, 0:ye_rowno_nse3+1, 0:temp_rowno_nse3+1, 6):: nsetable3_binde

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flag for being in NSE state !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nse_flag = 0 means in C-burning 
! nse_flag = 1 means in O- and Mg- burning
! nse_flag = 2 means in NSE burning
INTEGER  , allocatable, DIMENSION (:) :: nse_flag

!!!!!!!!!!!!!!!!!!!!!!!!
! Section for neutrino !
!!!!!!!!!!!!!!!!!!!!!!!!
! Some global quantities
! 1. Total non-thermal neutrino loss
! 2. Luminsoity
! 3. Luminosity by burning, def. and det.
REAL (DP) :: total_ecap_nu_qdot
REAL (DP) :: lumino
REAL (DP) :: lumino_burn, lumino_flame, lumino_deton

! MAss burned by deflgration/detonation
REAL (DP) :: burn_mass

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine allocates all the neccessary arrays !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDNUCLEAR
USE DEFINITION
IMPLICIT NONE

! allocate mean atomic number and mean atomic mass !
allocate(abar2(-4:length_step_2 + 5))
allocate(zbar2(-4:length_step_2 + 5))

! allocate nse flag !
allocate(nse_flag(-4:length_step_2 + 5))

! sum of isotope mass !
allocate(xmass_sum(total_ion))

! allocate mean atomic composition !
allocate(xiso(total_ion,-4:length_step_2 + 5))

END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates all the neccessary arrays !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DESTROYNUCLEAR
USE DEFINITION
IMPLICIT NONE

! Deallocate mean atomic number and mean atomic mass !
Deallocate(abar2) 
Deallocate(zbar2)

! Deallocate nse flag !
Deallocate(nse_flag) 

! sum of isotope mass !
deallocate(xmass_sum)

! Deallocate mean atomic composition !
Deallocate(xiso)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads the NSE table which assumed a NSE	!
! Of Si28 <=> Ni56 + alpha If sufficient temperature is reached !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_NSE_TABLE
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
INTEGER :: i, j, k
REAL (DP) :: dummy

! Open the nse table !
OPEN(unit=500, file='./Library/nse_table_7iso.dat',action='read')

! Now read the nse table !
DO i = 0, den_rowno_nse, 1
	DO j = 0, temp_rowno_nse, 1
		READ(500,*) dummy, dummy, nsetable_binde(i,j), (nsetable_xiso(i,j,k), k = 1, total_ion)
		nsetable_binde(i,j) = nsetable_binde(i,j) / 9.0D20
	ENDDO
ENDDO

CLOSE(500)

END SUBROUTINE read_nse_table

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine computes the energy release and the composition of !
! isotopes when nuclear statistical equilibrium is reached 	     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETNSESTATE(rho_in, temp_in, xiso_nse_out)
USE DEFINITION
IMPLICIT NONE

! Input parameter !
REAL (DP), INTENT(IN) :: rho_in, temp_in

! Output parameter !
REAL (DP), INTENT(OUT) :: xiso_nse_out(1:total_ion)

! Dummy varaibles !
REAL (DP) :: log10rho

! Integer parameter !
INTEGER :: rho_grid, temp_grid, ye_grid

! Real parameter !
REAL (DP) :: rho_dgrid, temp_dgrid, ye_dgrid

! assign the logarithm of density !
log10rho = LOG10(rho_in/density)                

IF(log10rho > 7.0D0) THEN            ! Minimum = 5.3
	if(temp_in > 5.0D0) THEN          ! Minimum = 3.0 and maximum 11.0
		rho_grid = INT((log10rho - 7.0D0) / 0.1D0)
		temp_grid = INT((temp_in - 4.0D0) / 0.1D0)

		rho_dgrid = (log10rho - (7.0D0 + (DBLE(rho_grid) * 0.1D0))) / 0.1D0
		temp_dgrid = (temp_in - (4.0D0 + (DBLE(temp_grid) * 0.1D0))) / 0.1D0

		IF(rho_grid >= 30) THEN
			rho_grid = 30         
			rho_dgrid = 0.0D0
		ENDIF

		IF(temp_grid >= 70) THEN
			temp_grid = 70
			temp_dgrid = 0.0D0                    
		ENDIF

		xiso_nse_out(:) =  nsetable_xiso(rho_grid, temp_grid, :) + &
			rho_dgrid * temp_dgrid * &
			(nsetable_xiso(rho_grid+1, temp_grid+1, :) - nsetable_xiso(rho_grid, temp_grid, :)) + & 
			rho_dgrid * (1.0D0 - temp_dgrid) * &
			(nsetable_xiso(rho_grid+1, temp_grid, :) - nsetable_xiso(rho_grid, temp_grid, :)) + &
			(1.0D0 - rho_dgrid) * temp_dgrid * &
			(nsetable_xiso(rho_grid, temp_grid+1, :) - nsetable_xiso(rho_grid, temp_grid, :)) 
	ENDIF
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculate all the nuclear variables for the elements ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INITIALIZEISO
USE DEFINITION, ONLY: DP
IMPLICIT NONE

! Integer parameter !
INTEGER :: i
 
! this is the conversion between different units of energy !
! Include MeV to erg and MeV to gram !
REAL (DP) :: 		mev2erg,mev2gr
REAL (DP), PARAMETER :: 	ev2erg = 1.60217648740d-12

! speed of light !
REAL (DP), PARAMETER :: 	clight  = 2.99792458d10

! avogardo number !
REAL (DP), PARAMETER :: 	avo     = 6.0221367d23

! Of course, the conversion is fixed by physical constant !
PARAMETER        (mev2erg = ev2erg*1.0d6, mev2gr  = mev2erg/clight**2)

burn_mass = 0.0D0

! set the id numbers of the elements !
! If you add more elements, this needed to be changed !
che4  = 1
cc12  = 2
co16  = 3
cne20 = 4
cmg24 = 5
csi28 = 6
cni56 = 7

! set the names of the elements
ionam(che4)  = 'he4 '
ionam(cc12)  = 'c12 '
ionam(co16)  = 'o16 '
ionam(cne20) = 'ne20'
ionam(cmg24) = 'mg24'
ionam(csi28) = 'si28'
ionam(cni56) = 'ni56'

! set the number of nucleons in the element !
aion(che4)  = 4.0d0
aion(cc12)  = 12.0d0
aion(co16)  = 16.0d0
aion(cne20) = 20.0d0
aion(cmg24) = 24.0d0
aion(csi28) = 28.0d0
aion(cni56) = 56.0d0

! set the number of protons in the element !
zion(che4)  = 2.0d0
zion(cc12)  = 6.0d0
zion(co16)  = 8.0d0
zion(cne20) = 10.0d0
zion(cmg24) = 12.0d0
zion(csi28) = 14.0d0
zion(cni56) = 28.0d0

! set the binding energy of (MeV) the element !
bion(che4)  =  28.29603d0
bion(cc12)  =  92.16294d0
bion(co16)  = 127.62093d0
bion(cne20) = 160.64788d0
bion(cmg24) = 198.25790d0
bion(csi28) = 236.53790d0
bion(cni56) = 484.00300d0

! set the number of neutrons and mass !
DO i = 1, total_ion
	nion(i) = aion(i) - zion(i)
ENDDO

! mass in gram of each isotope !
! the calculation is total neutron mass + total proton mass !
! minus the binding energy !
DO i = 1, total_ion
	mion(i) = nion(i)*1.67492721184d-24 + zion(i)*1.67262163783d-24 - bion(i)*mev2gr
ENDDO

! molar mass !      
DO i = 1, total_ion
	wion(i) = avo * mion(i)             
ENDDO

! a common approximation !   
DO i = 1, total_ion
	wion(i) = aion(i)       
ENDDO

! Also, initial the network for NSE !
CALL init_network

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine compute the binding energy for     !
! a given composition. The subroutine takes          !
! chemical composition as input and return the       !
! binding energy. Written by Leung Shing Chi in 2016 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_binde(x_in, binde_out)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Input chemical composition !
! noted that the input must be an array !
real (DP) :: x_in(1:total_ion)

! Output binding energy 
real (DP) :: binde_out

! Dummy variable
integer :: i

! Binding energy
real (DP) :: binde_sum

! Initialiation
binde_sum = 0.0D0
  
! Sum all the isotopes
do i = 1, total_ion
	binde_sum = binde_sum + x_in(i) / mion(i) * bion(i)
enddo

! Change to code unit
binde_out = 1.782661907D-27 * binde_sum

end subroutine compute_binde

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
! This subroutine checks if the local isotope mass fraction sum !
! to one. If not, normalize it. 				!
! Written by Leung Shing Chi in 2016.				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkxisotope
USE DEFINITION
IMPLICIT NONE

! dummy variables !
integer :: i, j
 
! local variables !
real (DP):: xiso_sum

! First check if any grid has any unusual mass fraction !
do i = 1, total_ion
	do j = 1, length_step_part_2
            	if(xiso(i,j) < xthres) THEN
			xiso(i,j) = xthres
		END IF
 	enddo
enddo

! Then do the sum and normalization !
do j = 1, length_step_part_2
       xiso_sum = 0.0D0
       		do i = 1, total_ion
         		xiso_sum = xiso_sum + xiso(i,j)
        	enddo
        xiso(:,j) = xiso(:,j) / xiso_sum
enddo

! Copy to ghost shell !
CALL BOUNDARY2D_X()

end subroutine checkxisotope

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
! This subroutine normalize the mass fractions inputed from the !
! reconstruction steps for Riemann solvers. Then it compute the !
! AZBAR used to compute the edge state speed of sound 		!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EDGEAZBAR(xiso_in, abar_out, zbar_out)
USE DEFINITION
IMPLICIT NONE

! dummy variables !
integer :: i, j

! input !
REAL (DP), INTENT(IN), DIMENSION (1 : total_ion) :: xiso_in

! output !
REAL (DP), INTENT(OUT) :: abar_out, zbar_out

! local variables !
REAL (DP), DIMENSION (1 : total_ion) :: xiso_local
real (DP) :: xiso_sum, dummy

! Assign !
xiso_local = xiso_in

! First check if any grid has any unusual mass fraction !
do i = 1, total_ion
	if(xiso_local(i) < xthres) THEN
		xiso_local(i) = xthres
	END IF
enddo

! Then do the sum and normalization !
xiso_sum = 0.0D0
do i = 1, total_ion
	xiso_sum = xiso_sum + xiso_local(i)
enddo
xiso_local(:) = xiso_local(:) / xiso_sum

! Then find the azbar !
CALL PRIVATE_HELMEOS_AZBAR(xiso_local(:), abar_out, zbar_out, dummy)

end subroutine EDGEAZBAR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                        
! This subroutine finds the mean mass number and mean atomic  !
! number of all grids. Written by Leung Shing Chi	      !	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine find_AZbar()
USE DEFINITION
IMPLICIT NONE

! Dummy variables !
integer :: i

! Dummy variable
real (DP) :: dummy

! Simply call the subroutine to help you
! find the abar and zbar directly
do i = 1, length_step_2
 	CALL PRIVATE_HELMEOS_AZBAR(xiso(:,i), abar2(i), zbar2(i), dummy)
enddo

! Copy the result to ghost cells
CALL BOUNDARY1D_NM (abar2, even)
CALL BOUNDARY1D_NM (zbar2, even)

end subroutine find_AZbar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                     
! This subroutine bridges the NUCLEAR_NET subroutine which ! 
! aims at finding the abar and zbar.			   !
! Written by Leung Shin Chi in 2016.			   !
! For more details about nuclear net, refer Timmes (2000b) !
! or his website.                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine private_helmeos_azbar(xiso_in, abar_out, zbar_out, electron_frac_out)
IMPLICIT NONE

! Input quantities
Real (selected_real_kind(15,307)), dimension(total_ion):: xiso_in

! Output quantitites
Real (selected_real_kind(15,307)) :: abar_out, zbar_out, electron_frac_out

! Dummy variables
Real (selected_real_kind(15,307)) :: dummy
Real (selected_real_kind(15,307)), dimension(total_ion) :: dummy_x

! The traditional way
!ymass(1:ionmax) = xmass(1:ionmax)/wion(1:ionmax)      
!wbar  = 1.0d0/sum(ymass(1:ionmax))
!sum1  = sum(aion(1:ionmax)*ymass(1:ionmax))
!abar2  = wbar * sum1
!electron_frac2  = sum(zion(1:ionmax)*ymass(1:ionmax))
!zbar2  = wbar * ye

! Call the preset subroutine
call private_azbar(xiso_in,aion,zion,wion,total_ion,dummy_x,abar_out,zbar_out,dummy,electron_frac_out,dummy)

end subroutine private_helmeos_azbar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                             
! This subroutine is obtained from Timmes' nuclear network !
! program to calculate the abar and zbar for a given       !
! chemical composition in a flexible structure             !
! Merged by Leung Shing Chi in 2016                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine private_azbar(xmass,aion,zion,wion,ionmax,ymass,abar,zbar,wbar,ye,nxcess)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! this routine calculates composition variables

! input:
! mass fractions               = xmass(1:ionmax)  dimensionless
! number of nucleons           = aion(1:ionmax)   dimensionless
! charge of nucleus            = zion(1:ionmax)   dimensionless
! atomic weight or molar mass  = wion(1:ionmax)    g/mole
! number of isotopes           = ionmax
!
! output:
! molar abundances        = ymass(1:ionmax)   mole/g
! mean number of nucleons = abar              dimensionless
! mean nucleon charge     = zbar              dimensionless
! mean weight             = wbar              g/mole
! electron fraction       = ye                mole/g
! neutron excess          = xcess

! declare the pass
integer          ionmax
real (selected_real_kind(15,307)), dimension(1:total_ion) :: xmass,aion,zion,wion,ymass
real (selected_real_kind(15,307)) :: abar,zbar,wbar,ye,nxcess

! local variables !
real (selected_real_kind(15,307)) asum,sum1

! molar abundances !
ymass(1:total_ion) = xmass(1:total_ion)/wion(1:total_ion)

! mean numner of nucleons
      abar  = 1.0d0/sum(ymass(1:ionmax))
      wbar  = abar

! mean charge !
ye  = sum(zion(1:total_ion)*ymass(1:total_ion))
zbar  = abar * ye

! neutron excess
nxcess = sum1 - 2.0d0 * ye

return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the total mass of !
! individual isopte			  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDISOTOPEMASS
USE DEFINITION
IMPLICIT NONE

! dummy variables !
INTEGER :: i, j

! Initialization !
xmass_sum = 0.0D0
  
! Do the integration
Do i = 1, total_ion
	DO j = 1, length_step_part_2

		! Check whether the mass fraction of isotopes reached zero !
		IF (rho2 (j) > rho2_a .AND. xiso(i, j) > xthres) THEN
			xmass_sum (i) = xmass_sum (i) + vol2 (j) * rho2 (j) * xiso(i, j)
		END IF	

	ENDDO
ENDDO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine open the files for outputing the isotope    !
! composition where we expect it should be changing with time !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENXISOFILE
USE DEFINITION
IMPLICIT NONE
 
! Integer variables !
INTEGER ::  fileno_len
CHARACTER (len = 256) :: fileno
 
WRITE (fileno, *) n_backup / time_step + 1
fileno = ADJUSTL (fileno)
fileno_len = LEN_TRIM (fileno)

! Isotope mass fractions profile !
OPEN (UNIT = 501, FILE = './Outfile/Isotope/Star_WENO_XHe4_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 502, FILE = './Outfile/Isotope/Star_WENO_XC12_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 503, FILE = './Outfile/Isotope/Star_WENO_XO16_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 504, FILE = './Outfile/Isotope/Star_WENO_XNe20_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 505, FILE = './Outfile/Isotope/Star_WENO_XMg24_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 506, FILE = './Outfile/Isotope/Star_WENO_XSi28_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 507, FILE = './Outfile/Isotope/Star_WENO_XNi56_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')

! Isotope mass fractions !
OPEN (UNIT = 521, FILE = './Outfile/Isotope/Star_WENO_MassHe4_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 522, FILE = './Outfile/Isotope/Star_WENO_MassC12_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 523, FILE = './Outfile/Isotope/Star_WENO_MassO16_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 524, FILE = './Outfile/Isotope/Star_WENO_MassNe20_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 525, FILE = './Outfile/Isotope/Star_WENO_MassMg24_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 526, FILE = './Outfile/Isotope/Star_WENO_MassSi28_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 527, FILE = './Outfile/Isotope/Star_WENO_MassNi56_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')

! Central Properties And Others
OPEN (UNIT = 531, FILE = './Outfile/Isotope/Star_WENO_CentralABar_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 532, FILE = './Outfile/Isotope/Star_WENO_CentralZBar_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')

! Energy loss
OPEN (UNIT = 541, FILE = './Outfile/Isotope/Star_WENO_EcapNeutrinoEnergyLoss_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine close the files for outputing the isotope   !
! composition where we expect it should be changing with time !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CLOSEXISOFILE
USE DEFINITION
IMPLICIT NONE

INTEGER :: j

Do j = 501, 507
	CLOSE(j)
END DO

Do j = 521, 527
	CLOSE(j)
END DO

Do j = 531, 532
	CLOSE(j)
END DO

Do j = 541, 541
	CLOSE(j)
END DO
 
END SUBROUTINE
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do the output of  the isotope               !
! composition where we expect it should be changing with time !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OUTPUTXISOPROFILE(time)
USE DEFINITION
IMPLICIT NONE

! The input time !
REAL(DP) :: time

! integer parameter !
INTEGER :: j 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Output xiso profile !
WRITE (501, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (501, 702) r2(j), xiso(1, j)
END DO
WRITE (501, *)

WRITE (502, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (502, 702) r2(j), xiso(2, j)
END DO
WRITE (502, *)

WRITE (503, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (503, 702) r2(j), xiso(3, j)
END DO
WRITE (503, *)

WRITE (504, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (504, 702) r2(j), xiso(4, j)
END DO
WRITE (504, *)

WRITE (505, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (505, 702) r2(j), xiso(5, j)
END DO
WRITE (505, *)

WRITE (506, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (506, 702) r2(j), xiso(6, j)
END DO
WRITE (506, *)

WRITE (507, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (507, 702) r2(j), xiso(7, j)
END DO
WRITE (507, *)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

701 FORMAT (F33.15, ES33.15)
702 FORMAT (F33.15, E26.15E3)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do the output of  the chemical              !
! composition where we expect it should be changing with time !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OUTPUTXISO(time)
USE DEFINITION
IMPLICIT NONE

! The input time !
REAL(DP) :: time

! Integer !
INTEGER :: n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Also find the neutrino loss !
CALL FINDISOTOPEMASS

! Isotope mass fractions !
WRITE (521, 701) time, xmass_sum(1)
WRITE (522, 701) time, xmass_sum(2)
WRITE (523, 701) time, xmass_sum(3)
WRITE (524, 701) time, xmass_sum(4)
WRITE (525, 701) time, xmass_sum(5)
WRITE (526, 701) time, xmass_sum(6)
WRITE (527, 701) time, xmass_sum(7)

! Central Properties !
WRITE (531, 701) time, abar2(1)
WRITE (532, 701) time, zbar2(1)     

! Energies loss !
WRITE (541, 701) time, total_ecap_nu_qdot / dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

701 FORMAT (F33.15, ES33.15)

END SUBROUTINE

END MODULE