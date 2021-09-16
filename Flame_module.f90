!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the deflagration module that is used to simulate deflagration !
! or detonation in 1D type Ia supernova explosion. In this module  	!
! we use levelset method the capturing the propagation of deflagration  !
! front (which we call flame in short).					!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE FLAME_MODULE
USE NUCLEAR_MODULE
USE DEFINITION
IMPLICIT NONE

! flag for finding detonation !
INTEGER :: founddeton_flag

! time when finding the detonation transition !
REAL (DP) :: found_deton_time

! flame and detonation density !
REAL (DP) :: flame_rho, deton_rho

! flame front quantity !
REAL (DP) :: flame_temp, flame_abar, flame_zbar, flame_ye

! deton front quantity !
REAL (DP) :: deton_temp, deton_abar, deton_zbar, deton_ye

! level set function for NM and DM!
REAL (DP), allocatable, dimension(:) :: scaG, scaG2

! The fraction occupied by the deflagration                                    
real (DP), allocatable, dimension(:) :: flame_ratio, flame_ratio_old 
real (DP), allocatable, dimension(:) :: deton_ratio, deton_ratio_old
real (DP), allocatable, dimension(:) :: burn_ratio

! The type of grid depending on the geometry !
integer, allocatable, dimension(:) :: flamegrid_flag, detongrid_flag

! Energy input by deflagration, deontation and burning
real (DP), allocatable, dimension (:) :: flame_qdot, deton_qdot, burn_qdot

! the integer required as determining extra conservative equation for levelset !
INTEGER :: iscaG, iscaG2

! Nuclear reaction density range
real (DP), parameter :: rho2_burn_max = 1.0D7*density
real (DP), parameter :: rho2_burn_min = 1.0D6*density

! Deflagration density range
real (DP), parameter :: rho2_flame_max = 1.0D10*density 
real (DP), parameter :: rho2_flame_min = 1.0D7*density

! Detonation density range
real (DP), parameter :: rho2_deton_max = 1.0D9*density
REAL (DP), parameter :: rho2_deton_min = 1.0D6*density 

! Mass ash !
REAL (DP) :: mass_ash

! We go on to the subroutine for levelset method for capturing deflagration !

contains
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine allocate the variables used in the levelset method !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDFLAME
USE DEFINITION
IMPLICIT NONE

! Allocate variables for deflagrations !	
allocate(scaG(-4:length_step_2 + 5))
allocate(flame_ratio(-4:length_step_2 + 5))
allocate(flame_ratio_old(-4:length_step_2 + 5))
allocate(flamegrid_flag(-4:length_step_2 + 5))
allocate(flame_qdot(-4:length_step_2+5))

! For detonations !
allocate(scaG2(-4:length_step_2 + 5))
allocate(deton_ratio(-4:length_step_2 + 5))
allocate(deton_ratio_old(-4:length_step_2 + 5))
allocate(burn_ratio(-4:length_step_2 + 5))
allocate(detongrid_flag(-4:length_step_2 + 5))
allocate(deton_qdot(-4:length_step_2+5))
allocate(burn_qdot(-4:length_step_2+5))

end subroutine buildFlame

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocate the variables used in the levelset method !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DESTROYFLAME
USE DEFINITION
IMPLICIT NONE

! Deallocate NM variables !
deallocate(scaG)
deallocate(scaG2)
deallocate(flame_ratio)
deallocate(flame_ratio_old)
deallocate(deton_ratio)
deallocate(deton_ratio_old)
deallocate(burn_ratio)
deallocate(flamegrid_flag)
deallocate(detongrid_flag)
deallocate(flame_qdot)
deallocate(deton_qdot)
deallocate(burn_qdot)

end subroutine destroyFlame

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine open suitable files related to nuclear burning !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENFILEFLAME
USE DEFINITION
IMPLICIT NONE
 
! Integer variables !
INTEGER ::  fileno_len
CHARACTER (len = 256) :: fileno
 
WRITE (fileno, *) n_backup / time_step + 1
fileno = ADJUSTL (fileno)
fileno_len = LEN_TRIM (fileno)

! This part is quite easy so no detail explaination !
OPEN (UNIT = 400, FILE = './Outfile/Flame/Star_WENO_BurnRatio_'//fileno(1:fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 401, FILE = './Outfile/Flame/Star_WENO_FlameRatio_'//fileno(1:fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 402, FILE = './Outfile/Flame/Star_WENO_FlameRatioOld_'//fileno(1:fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 403, FILE = './Outfile/Flame/Star_WENO_FlamegridFlag_'//fileno(1:fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 404, FILE = './Outfile/Flame/Star_WENO_DetonRatio_'//fileno(1:fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 405, FILE = './Outfile/Flame/Star_WENO_DetonRatioOld_'//fileno(1:fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 406, FILE = './Outfile/Flame/Star_WENO_DetongridFlag_'//fileno(1:fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 407, FILE = './Outfile/Flame/Star_WENO_NSEFlag_'//fileno(1:fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 408, FILE = './Outfile/Flame/Star_WENO_ScaG_'//fileno(1:fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 409, FILE = './Outfile/Flame/Star_WENO_ScaG2_'//fileno(1:fileno_len)//'.dat', STATUS = 'REPLACE')

! Global !
OPEN (UNIT = 411, FILE = './Outfile/Flame/Star_WENO_Luminosity_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 412, FILE = './Outfile/Flame/Star_WENO_Ash_'//fileno(1:fileno_len)//'.dat', STATUS = 'REPLACE')
	
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine close suitable files related to nuclear burning !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CLOSEFILEFLAME
USE DEFINITION
IMPLICIT NONE

! Integer !
Integer :: i 

! This part is quite easy so no detail explaination !
DO i = 400, 409
	CLOSE (i)
END DO
DO i = 411, 412
	CLOSE (i)
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine output the quantity that is related to the level set !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OUTPUTFLAME(time)
USE DEFINITION
USE NUCLEAR_MODULE
IMPLICIT NONE

! reak parameter !
REAL (DP) :: time

! Call subroutine !
CALL FINDASH
CALL findluminosity

! Write the output !
WRITE (411, 702) time, lumino, lumino_flame, lumino_deton, lumino_burn, burn_mass
WRITE (412, 701) global_time, mass_ash

! Format !
701 FORMAT (F33.15, ES33.15)
702 FORMAT (F33.15, 20ES33.15)
 
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine output the quantity that is related to the level set !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OUTPUTPROFILEFLAME(time)
USE DEFINITION
IMPLICIT NONE
 
! real parameter !
REAL (DP) :: time
 
! integer parameter !
INTEGER :: j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Related to flame !
WRITE (400, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (400, 701) r2(j), burn_ratio(j)
END DO
WRITE (400, *)

WRITE (401, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (401, 701) r2(j), flame_ratio(j)
END DO
WRITE (401, *)

WRITE (402, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (402, 701) r2(j), flame_ratio_old(j)
END DO
WRITE (402, *)

WRITE (403, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (403, 701) r2(j), DBLE(flamegrid_flag(j))
END DO
WRITE (403, *)

WRITE (404, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (404, 701) r2(j), deton_ratio(j)
END DO
WRITE (404, *)

WRITE (405, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (405, 701) r2(j), deton_ratio_old(j)
END DO
WRITE (405, *)

WRITE (406, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (406, 701) r2(j), DBLE(detongrid_flag(j))
END DO
WRITE (406, *)

WRITE (407, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (407, 701) r2(j), DBLE(nse_flag(j))
END DO
WRITE (407, *)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Level Set !
WRITE (408, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (408, 701) r2(j), ScaG(j)
END DO
WRITE (408, *)

WRITE (409, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (409, 701) r2(j), ScaG2(j)
END DO
WRITE (409, *)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

701 FORMAT (F33.15, ES33.15)
 
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the total luminosity !
! by summing up all energy-generating source !
! Merged by Leung Shing Chi in 2016	     !
! Updated by Leung Shing CHi in 2017	     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDLUMINOSITY
USE DEFINITION
IMPLICIT NONE

! dummy variables
INTEGER :: j

! Initilization
lumino_flame = 0.0D0
lumino_deton = 0.0D0
lumino_burn = 0.0D0
   
! Do the integration
DO j = 1, length_step_part_2, 1
	lumino_flame = lumino_flame + vol2(j) * rho2(j) * flame_qdot(j)
ENDDO

DO j = 1, length_step_part_2, 1
	lumino_deton = lumino_deton + vol2(j) * rho2(j) * deton_qdot(j)
ENDDO

DO j = 1, length_step_part_2, 1
	lumino_burn = lumino_burn + vol2(j) * rho2(j) * burn_qdot(j)
ENDDO

! Divide dt to get the time-rate
lumino_flame = lumino_flame / dt
lumino_deton = lumino_deton / dt
lumino_burn = lumino_burn / dt

! Sum them up
lumino = (lumino_flame + lumino_burn + lumino_deton) 

END SUBROUTINE findluminosity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialize the level set function using signed distance function ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFLAME
USE DEFINITION
IMPLICIT NONE

! integer variables !
INTEGER :: j

! Real variables !
REAL(DP) :: m_defla, x2

! Temporal !
REAL (DP) :: dist1, dist2

! Initilization of flame ratio !
flame_ratio = 0.0D0
flame_ratio_old = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The central ignited model !

! Locate the initial deflagration point by interpolation !
m_defla = 0.02D0
DO j = 1, length_step_2
	If(m_cell(j) > m_defla) THEN
		CALL LINEAR(m_cell(j-1), m_cell(j), r2F(j-1), r2F(j), m_defla, x2) 
		EXIT
	END IF
END DO

! Locate the initial deflagration near star center !					   
DO j = 1, length_step_2
	scaG(j) =  r2(j) - x2
END DO

! We need the updated flame ratio !
DO j = 1, length_step_2
	IF(scaG(j) > 0.0D0) THEN
		flame_ratio (j) = 0.0D0
	ELSEIF(scaG(j) < 0.0D0) THEN
		flame_ratio (j) = 1.0D0
	END IF
END DO
DO j = 1, length_step_2
	IF(x2 < r2F(j)) THEN
		flame_ratio (j) = (x2 - r2F(j - 1))/(dx2)
		EXIT
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Copy the data to ghost cells !
CALL boundary1D_NM(scaG,even)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reinitialize the level set so that it remains signed distance !
! As long as the initial level set is signed - distance, no reinitialization is needed !
!call ReiniFlame
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Backup initial result
flame_ratio_old = flame_ratio

END SUBROUTINE GetFlame

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialize the level set function using signed distance function ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETDETON
USE DEFINITION
IMPLICIT NONE

! integer variables !
INTEGER :: j

! Real variables !
REAL(DP) :: m_deton, x2

! Initilization of flame ratio !
deton_ratio_old = 0.0D0
deton_ratio = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The central ignited model !

! Locate the initial deflagration point by interpolation !
m_deton = 0.02D0
DO j = 1, length_step_2
	If(m_cell(j) > m_deton) THEN
		CALL LINEAR(m_cell(j-1), m_cell(j), r2F(j-1), m_deton, x2) 
		EXIT
	END IF
END DO

! Locate the initial deflagration near star center !					   
DO j = 1, length_step_2
	scaG2(j) =  r2(j) - x2
END DO

! We need the updated flame ratio !
DO j = 1, length_step_2
	IF(scaG2(j) > 0.0D0) THEN
		deton_ratio (j) = 0.0D0
	ELSEIF(scaG2(j) < 0.0D0) THEN
		deton_ratio (j) = 1.0D0
	END IF
END DO
DO j = 1, length_step_2
	IF(x2 < r2F(j)) THEN
		deton_ratio (j) = (x2 - r2F(j - 1))/(dx2)
		EXIT
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Copy the data to ghost cells !
CALL boundary1D_NM(scaG2,even)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reinitialize the level set so that it remains signed distance !
! As long as the initial level set is signed - distance, no reinitialization is needed !
!call Reinideton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Backup initial result
deton_ratio_old = deton_ratio

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine reinitialize the level set so to maintain its property of being signed distance function !
! this reinitialization process follows from reinecke et. al (1999)					   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE REINIFLAME
USE DEFINITION
IMPLICIT NONE

! Output array of intersection point
REAL (DP), dimension(-4 : length_step_2 + 5) :: flame_position_out
REAL (DP), dimension(-4 : length_step_2 + 5) :: flame_position_in

! Flame distance
REAL (DP), dimension(-4 : length_step_2 + 5) :: flame_distance_out
REAL (DP), dimension(-4 : length_step_2 + 5) :: flame_distance_in

! Flame distance
REAL (DP), dimension(-4 : length_step_2 + 5) :: dist_out
REAL (DP), dimension(-4 : length_step_2 + 5) :: dist_in

! Temporal variables !
REAL (DP) :: last_distance, distance, flame_distance

! Flame front !
integer :: flame_out, flame_in

! integer variables !
integer :: k_in, k_out

! integer variables !
integer :: j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
flame_in = 0 
flame_out = 0
flame_position_in = 0.0D0
flame_position_out = 0.0D0
flame_distance_in = 0.0D0
flame_distance_out = 0.0D0

! Set large flame distance !
dist_in = 1.0D10
dist_out = 1.0D10

! Get the flame position !
DO j = 1, length_step_2 - 1
	IF (scaG(j) < 0.0D0 .AND. scaG(j+1) > 0.0D0) then
		flame_out = flame_out + 1
		CALL LINEAR (scaG(j), scaG(j+1), r2(j), r2(j+1), 0.0D0, flame_position_out(flame_out))
	ELSEIF (scaG(j) > 0.0D0 .AND. scaG(j+1) < 0.0D0) then
		flame_in = flame_in + 1
		CALL LINEAR (scaG(j), scaG(j+1), r2(j), r2(j+1), 0.0D0, flame_position_in(flame_in))
	END IF
END DO

! For the outward flame front !
IF(flame_out > 0) THEN
	DO j = 1, length_step_2
		flame_distance = 0.0D0
		last_distance = 10.0D0*r2(length_step_2)
		DO k = 1, flame_out
			distance = ABS(flame_position_out (k) - r2(j))
  	    		IF ((flame_out > 1 .and. distance <= last_distance) .or. flame_out == 1) then
	       			flame_distance = distance 
   	       			last_distance = distance
				k_out = k
			END IF
		END DO
		dist_out (j) = flame_distance
		flame_distance_out (j) = flame_position_out (k_out)
	END DO
END IF
	
! For the inward flame front !
IF(flame_in > 0) THEN
	DO j = 1, length_step_2
		flame_distance = 0.0D0
		last_distance = 10.0D0*r2(length_step_2)
		DO k = 1, flame_in
			distance = ABS(flame_position_in (k) - r2(j))
  	    		IF ((flame_in > 1 .and. distance <= last_distance) .or. flame_in == 1) then
	       			flame_distance = distance 
   	       			last_distance = distance
				k_in = k
			END IF
		END DO
		dist_in (j) = flame_distance
		flame_distance_in (j) = flame_position_in (k_in)
	END DO
END IF

! Reinitialize !
DO j = 1, length_step_2
	
	! Find the minimum between the two !
	IF(dist_in (j) > dist_out (j)) THEN
		IF(r2(j) > flame_distance_out (j)) THEN
			scaG(j) = dist_out (j)
		ELSEIF(r2(j) < flame_distance_out (j)) THEN
			scaG(j) = - dist_out (j)
		END IF
	ELSEIF(dist_in (j) < dist_out (j)) THEN
		IF(r2(j) > flame_distance_in (j)) THEN
			scaG(j) = - dist_in (j)
		ELSEIF(r2(j) < flame_distance_in (j)) THEN
			scaG(j) = dist_in (j)
		END IF
	END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Copy the results to ghost cells !
CALL BOUNDARY1D_NM(scaG, even)
 
END SUBROUTINE ReiniFlame

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine reinitialize the level set so to maintain its property of being signed distance function !
! this reinitialization process follows from reinecke et. al (1999)					   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE REINIDETON
USE DEFINITION
IMPLICIT NONE

! Output array of intersection point
REAL (DP), dimension(-4 : length_step_2 + 5) :: flame_position_out
REAL (DP), dimension(-4 : length_step_2 + 5) :: flame_position_in

! Flame distance
REAL (DP), dimension(-4 : length_step_2 + 5) :: flame_distance_out
REAL (DP), dimension(-4 : length_step_2 + 5) :: flame_distance_in

! Flame distance
REAL (DP), dimension(-4 : length_step_2 + 5) :: dist_out
REAL (DP), dimension(-4 : length_step_2 + 5) :: dist_in

! Temporal variables !
REAL (DP) :: last_distance, distance, flame_distance

! Flame front !
integer :: flame_out, flame_in

! integer variables !
integer :: k_in, k_out

! integer variables !
integer :: j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
flame_in = 0 
flame_out = 0
flame_position_in = 0.0D0
flame_position_out = 0.0D0
flame_distance_in = 0.0D0
flame_distance_out = 0.0D0

! Set large flame distance !
dist_in = 1.0D10
dist_out = 1.0D10

! Get the flame position !
DO j = 1, length_step_2 - 1
	IF (scaG2(j) < 0.0D0 .AND. scaG2(j+1) > 0.0D0) then
		flame_out = flame_out + 1
		CALL LINEAR (scaG2(j), scaG2(j+1), r2(j), r2(j+1), 0.0D0, flame_position_out(flame_out))
	ELSEIF (scaG2(j) > 0.0D0 .AND. scaG2(j+1) < 0.0D0) then
		flame_in = flame_in + 1
		CALL LINEAR (scaG2(j), scaG2(j+1), r2(j), r2(j+1), 0.0D0, flame_position_in(flame_in))
	END IF
END DO

! For the outward flame front !
IF(flame_out > 0) THEN
	DO j = 1, length_step_2
		flame_distance = 0.0D0
		last_distance = 10.0D0*r2(length_step_2)
		DO k = 1, flame_out
			distance = ABS(flame_position_out (k) - r2(j))
  	    		IF ((flame_out > 1 .and. distance <= last_distance) .or. flame_out == 1) then
	       			flame_distance = distance 
   	       			last_distance = distance
				k_out = k
			END IF
		END DO
		dist_out (j) = flame_distance
		flame_distance_out (j) = flame_position_out (k_out)
	END DO
END IF
	
! For the inward flame front !
IF(flame_in > 0) THEN
	DO j = 1, length_step_2
		flame_distance = 0.0D0
		last_distance = 10.0D0*r2(length_step_2)
		DO k = 1, flame_in
			distance = ABS(flame_position_in (k) - r2(j))
  	    		IF ((flame_in > 1 .and. distance <= last_distance) .or. flame_in == 1) then
	       			flame_distance = distance 
   	       			last_distance = distance
				k_in = k
			END IF
		END DO
		dist_in (j) = flame_distance
		flame_distance_in (j) = flame_position_in (k_in)
	END DO
END IF

! Reinitialize !
DO j = 1, length_step_2
	
	! Find the minimum between the two !
	IF(dist_in (j) > dist_out (j)) THEN
		IF(r2(j) > flame_distance_out (j)) THEN
			scaG2(j) = dist_out (j)
		ELSEIF(r2(j) < flame_distance_out (j)) THEN
			scaG2(j) = - dist_out (j)
		END IF
	ELSEIF(dist_in (j) < dist_out (j)) THEN
		IF(r2(j) > flame_distance_in (j)) THEN
			scaG2(j) = - dist_in (j)
		ELSEIF(r2(j) < flame_distance_in (j)) THEN
			scaG2(j) = dist_in (j)
		END IF
	END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Copy the results to ghost cells !
CALL BOUNDARY1D_NM(scaG2, even)
 
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the fraction of grid !
! occupied by one the deflagration 		  !
! Written by Leung Shing Chi in 2016		  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE COMPUTE_FLAMERATIO()
USE DEFINITION
IMPLICIT NONE

! Dummy variables !
integer :: j, k

! The location of zero level set!
Real(DP) :: x_in, x_out

! Temporal variables !
REAL (DP) :: dummy_in, dummy_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialization the flame ratio to zero !
flame_ratio = 0.0D0

! We need the updated flame ratio !
DO j = 1, length_step_2
	IF(scaG(j) > 0.0D0) THEN
		flame_ratio (j) = 0.0D0
	ELSEIF(scaG(j) < 0.0D0) THEN
		flame_ratio (j) = 1.0D0
	END IF
END DO

! locate the deflagration crossing point by interpolation ! 
DO j = 1, length_step_2 - 1
	If (scaG(j) < 0.0D0 .AND. scaG(j+1) > 0.0E0_DP) then
		CALL LINEAR (scaG(j), scaG(j+1), r2(j), r2(j+1), 0.0D0, x_out)
		DO k = 1, length_step_2
			IF(r2F(k) > x_out) THEN
				flame_ratio(k) = (x_out - r2F(k-1))/dx2
				EXIT
			END IF
		END DO
	END IF
	If (scaG(j) > 0.0D0 .AND. scaG(j+1) < 0.0E0_DP) then
		CALL LINEAR (scaG(j), scaG(j+1), r2(j), r2(j+1), 0.0D0, x_in)
		DO k = 1, length_step_2
			IF(r2F(k) > x_in) THEN
				flame_ratio(k) = (r2F(k) - x_in)/dx2
				EXIT
			END IF
		END DO
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! According to the summed fraction, find the !
! grid type according to the geometry !
do j = 1, length_step_2 

	! Get the effective flame ratio !
	! We dont allow flame_ratio to be reduced !
	! I.e. Turning ash into fuel !
	flame_ratio(j) = MIN(MAX(flame_ratio(j), flame_ratio_old(j)), 1.0D0 - deton_ratio(j))

 	if(flame_ratio(j) == 0.0D0) then

 		! Completely empty
 		flamegrid_flag(j) = 0

 	elseif(flame_ratio(j) == 1.0D0) then

 		! Completely filled
 		flamegrid_flag(j) = 1

 	else
	
 		! Partially filled	
 		flamegrid_flag(j) = 2

	endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Copy the results to ghost cells
CALL BOUNDARY1D_NM(flame_ratio,even)

end subroutine compute_flameratio

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the fraction of grid !
! occupied by one the deflagration 		  !
! Written by Leung Shing Chi in 2016		  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE COMPUTE_DETONRATIO()
USE DEFINITION
IMPLICIT NONE

! Dummy variables !
integer :: j, k

! The location of zero level set!
Real(DP) :: x_in, x_out

! Temporal variables !
REAL (DP) :: dummy_in, dummy_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialization the flame ratio to zero !
deton_ratio = 0.0D0

! We need the updated flame ratio !
DO j = 1, length_step_2
	IF(scaG2(j) > 0.0D0) THEN
		deton_ratio (j) = 0.0D0
	ELSEIF(scaG2(j) < 0.0D0) THEN
		deton_ratio (j) = 1.0D0
	END IF
END DO

! locate the deflagration crossing point by interpolation ! 
DO j = 1, length_step_2 - 1
	If (scaG2(j) < 0.0D0 .AND. scaG2(j+1) > 0.0E0_DP) then
		CALL LINEAR (scaG2(j), scaG2(j+1), r2(j), r2(j+1), 0.0D0, x_out)
		DO k = 1, length_step_2
			IF(r2F(k) > x_out) THEN
				deton_ratio(k) = (x_out - r2F(k-1))/dx2
				EXIT
			END IF
		END DO
	END IF
	If (scaG2(j) > 0.0D0 .AND. scaG2(j+1) < 0.0E0_DP) then
		CALL LINEAR (scaG2(j), scaG2(j+1), r2(j), r2(j+1), 0.0D0, x_in)
		DO k = 1, length_step_2
			IF(r2F(k) > x_in) THEN
				deton_ratio(k) = (r2F(k) - x_in)/dx2
				EXIT
			END IF
		END DO
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! According to the summed fraction, find the !
! grid type according to the geometry !
do j = 1, length_step_2

	! Get the effective flame ratio !
	! We dont allow flame_ratio to be reduced !
	! I.e. Turning ash into fuel !
	deton_ratio(j) = MIN(MAX(deton_ratio(j), deton_ratio_old(j)), 1.0D0 - flame_ratio (j))

 	if(deton_ratio(j) == 0.0D0) then

 		! Completely empty
 		detongrid_flag(j) = 0

 	elseif(deton_ratio(j) == 1.0D0) then

 		! Completely filled
 		detongrid_flag(j) = 1

 	else
	
 		! Partially filled	
 		detongrid_flag(j) = 2

	endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Copy the results to ghost cells
CALL BOUNDARY1D_NM(deton_ratio,even)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the flame propagation velocity by finding ! 
! the local sound speed. The deflagration should propagate at 	  !
! speed as suggested by timmers et. al 1992			  !
! Written by Leung Shing Chi in 2016				  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_flame_velocity(j_in, cs_out)
USE DEFINITION
USE FLAMETABLE_MODULE
USE NUCLEAR_MODULE
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: j_in

! Integer !
INTEGER :: j

! Transition Density !
REAL (DP), PARAMETER :: rhot = 2.2D0 * 1.6190D-11

! Distance !
REAL (DP) :: dist1, dist2

! the local sound speed and atwood number !
REAL (DP) :: cs_out, flame_vel

! If deton_flag is on, then check for DDT !
if(deton_flag == 1 .AND. founddeton_flag == 0) then

	! Typing DDT required density, see khokhlov et al. 1997 !
	if(rho2(j_in) <= rhot .and. rho2(j_in) > rho2_deton_min .AND. rho2(j_in) <= rho2_deton_max) then

		! Write the relevant output !
		write(*,*)
		write(*,*) 'DDT Occured at = ', global_time
		write(*,*) 'Found deton at position ', r2(j_in)
		write(*,*)

		! Set deton ratio !
		deton_ratio_old (:) = 0.0D0

		! Plant the level set
                do j = 1, length_step_2, 1   
			dist1 = r2(j) - (r2(j_in) + 10.0D0)
			dist2 = (r2(j_in) - 10.0D0) - r2(j)
                        scaG2(j) = max(dist1, dist2)
		enddo

		! Compute deton ratio !
		CALL compute_detonratio()
        
		! Boundary !
                CALL boundary1D_NM (scaG2,even)

		! switch the found deton flag
 		founddeton_flag = 1

		! Set the found deton flag		
		found_deton_time = global_time
		output_ddt_last = found_deton_time

		! Set output file !
		output_file = .true.

	END IF

END IF

! Find flame speed !
IF(rho2(j_in) > rho2_flame_min .AND. rho2(j_in) < rho2_flame_max) then

	! Use deflagration speed only flame within the deflagration density range can propagate !
	CALL HELMEOS_CS(rho2(j_in), temp2(j_in), abar2(j_in), zbar2(j_in), ye2(j_in), flame_vel)
	cs_out = 0.06D0*flame_vel

ELSE
	cs_out = 0.0D0
END IF

! Format !
701 FORMAT (F33.15, 4ES33.15)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do the propagation of level-set due !
! to its normal-flow (propagation of flame in the     !
! supernovae context)				      !
! Written by Leung SHing CHi in 2016		      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE UPDATE_SCAG
USE DEFINITION
IMPLICIT NONE
      
! dummy variables !
integer :: j, k

! local sound speed !
real (DP) :: cs_grid
      
! initialization
cs_grid = 0.0D0

! Update the level-set !
do j = 1, length_step_2
	
	! First find the grid where it is partially burnt !
	IF(flamegrid_flag(j) > 1) THEN
		call update_flame_velocity(j, cs_grid) 
	ELSE
		cs_grid = 0.0D0
	END IF

	! Update level-set !
 	scaG(j) = scaG(j) - dt * cs_grid

enddo

! Fill the level set ghost grid
CALL BOUNDARY1D_NM(scaG, even)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do the propagation of level-set due !
! to its normal-flow (propagation of flame in the     !
! supernovae context)				      !
! Written by Leung SHing CHi in 2016		      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE UPDATE_SCAG2
USE DEFINITION
USE FLAMETABLE_MODULE
USE NUCLEAR_MODULE
IMPLICIT NONE
      
! dummy variables !
integer :: j

! local sound speed !
real (DP) :: cs_grid
      
! initialization
cs_grid = 0.0D0

! Update the level-set !
do j = 1, length_step_2

	! First find the grid where it is partially burnt !
	IF(detongrid_flag(j) > 1) THEN
		
		! Only flame within the deflagration density range can propagate !
		If(rho2(j) > rho2_flame_min .AND. rho2(j) < rho2_flame_max) then
			CALL HELMEOS_CS(rho2(j), temp2(j), abar2(j), zbar2(j), ye2(j), cs_grid)
		else
			cs_grid = 0.0D0
		endif 

	ELSE
		cs_grid = 0.0D0
	END IF
	
	! Update level-set !
 	scaG2(j) = scaG2(j) - dt * cs_grid

enddo

! Fill the level set ghost grid
CALL BOUNDARY1D_NM(scaG2, even)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine collects all the subroutines !
! for updating the level-set		       !
! Written by Leung Shing Chi in 2016           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE UPDATE_FLAME_RADIUS
USE DEFINITION
IMPLICIT NONE

! We first calculate the updated flame ratio !
! After the flame has been advected by the fluid !
CALL compute_flameratio()  

! For detonation !
IF(founddeton_flag == 1) THEN
	CALL compute_detonratio()  
END IF

! Backup the data
flame_ratio_old = flame_ratio 

! For detonation !
IF(founddeton_flag == 1) THEN
	deton_ratio_old = deton_ratio 
END IF

! Update the level-set !
CALL update_scaG

! For detonation !
IF(founddeton_flag == 1) THEN
	CALL update_scaG2
END IF

! We need to reinitialize !
CALL REINIFLAME

! For detonation !
IF(founddeton_flag == 1) THEN
	CALL REINIDETON
END IF

! comput the geometry for the deflagration !
CALL compute_flameratio

! For detonation !
IF(founddeton_flag == 1) THEN
	CALL compute_detonratio
END IF

! Conmpute the total local fraction
burn_ratio (:) = flame_ratio (:) + deton_ratio (:)

end subroutine update_flame_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This initialize all the flame, deton and burn energy release !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE QDOT_INI
USE DEFINITION 
IMPLICIT NONE

flame_qdot (:) = 0.0D0
deton_qdot (:) = 0.0D0
burn_qdot (:) = 0.0D0

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the initial input  !
! of energy due to deflagration/detonation !
! Written by Leung Shing Chi in 2016	   !
! Updated by Leung Shing Chi in 2017	   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FLAME_INI()
USE DEFINITION
USE NUCLEAR_MODULE
USE FLAMETABLE_MODULE
IMPLICIT NONE

! integer variables !
INTEGER :: j

! Local variables !
REAL (DP):: rho_mid, temp_mid, rho_ash, eps_ash, eps_old

! Local dummy !
REAL (DP) :: dummy
REAL (DP) :: flame_mid

! Ash composition !
REAL (DP), dimension(total_ion) :: x_burn

! Initilzation
flame_qdot(:) = 0.0D0

DO j = 1, length_step_part_2

	! Store to local variables
	rho_mid = rho2 (j) 
	temp_mid = temp2 (j) 
   	eps_old = epsilon2(j)

	! Energy injection and change in chemical composition are done !
	! to grids which are first assumed to be burnt !
	IF(rho_mid > rho2_flame_min .and. flame_ratio(j) > 0.0D0) THEN
		CALL readtable_flameenergy(temp_mid, rho_mid, flame_mid, x_burn(:))
		CALL readtable_flamestate(temp_mid, rho_mid, rho_ash, dummy, eps_ash)
		epsilon2 (j) = (1.0D0 - flame_ratio(j)) * epsilon2 (j) + flame_ratio (j) * eps_ash
		rho2(j) = (1.0D0 - flame_ratio(j)) * rho_mid + flame_ratio(j) * rho_ash
		xiso (:, j) = (1.0D0 - flame_ratio (j)) * xiso(:, j) + flame_ratio (j) * x_burn(:) 
		burn_mass = burn_mass + rho_mid * vol2(j)
	endif

	! assign energy input !
 	energy_input = energy_input + (rho2(j)*epsilon2(j) - rho_mid*eps_old)*vol2(j)

ENDDO

! Copy to new results to ghost grids !
CALL BOUNDARY1D_NM(rho2, even)
CALL BOUNDARY2D_X(xiso)
CALL BOUNDARY1D_NM(epsilon2, even)

END SUBROUTINE flame_ini

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the energy release !
! By assuming 12C --> 24Mg			!
! Written by Leung Shing Chi in 2016		!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BURN_PHASE1
USE DEFINITION
USE NUCLEAR_MODULE
IMPLICIT NONE  

! Integer variables !
INTEGER :: j

! Local variables !
REAL (DP) :: rho_mid, temp_mid, vol_mid, eps_old
REAL (DP) :: rho_burn_min, rho_burn_max

! Change of fraction for level sets !
REAL (DP) :: x1, x2

! local flame energy, defined by Dqbar_C12 * xiso_ini(c12) + Dqbar_Ne22 * metallicity
! i.e. 6.2762E-4 * (xiso_ini(c12) - Z/2) + 1.9527E-4 * Z
! Example: Z = 0.00, flame_ene = 3.13810E-4
! Example: Z = 0.02, flame_ene = 3.11439E-4
! Example: Z = 0.04, flame_ene = 3.09068E-4
! Example: Z = 0.06, flame_ene = 3.06698E-4
! Example: Z = 0.08, flame_ene = 3.04327E-4
! Example: Z = 0.10, flame_ene = 3.01956E-4
real (DP), parameter :: flame_ene = 3.11439D-4

! fuel and Ash composition !
real (DP), dimension(total_ion) :: x_fuel1 = (/0.0D0, xc12_ini, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
real (DP), dimension(total_ion) :: x_ash1 = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, xc12_ini, 0.0D0, 0.0D0/)

! local deton energy, defined by Dqbar_C12 * xiso_ini(c12) + Dqbar_Ne22 * metallicity
! i.e. 6.2762E-4 * (xiso_ini(c12) - Z/2) + 1.9527E-4 * Z
real (DP), parameter :: deton_ene = 3.11439D-4  ! He to Ne

! First, initialize the grid, to erase all past data
flame_qdot(:) = 0.0D0
deton_qdot(:) = 0.0D0

! Then, find the energy release due to nuclear fusion
do j = 1, length_step_part_2
	rho_mid = rho2 (j)
	temp_mid = temp2 (j)
	vol_mid = vol2 (j)
	eps_old = epsilon2(j)

	! Just compute the change of fraction occupies !
	! by deflagration or detonation !
	x1 = MIN(flame_ratio(j) - flame_ratio_old(j), xiso(cc12,j) / 0.51D0)
	x2 = MIN(deton_ratio(j) - deton_ratio_old(j), xiso(cc12,j) / 0.51D0)
	
	! Remember to switch flame_ini as well !
	! When there is a change in area fraction !
	! Change the chemical composition accordingly !
	! And inject temperature !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! For no deton level-set approach, burning density is changing !
	!if(x1 > 0.0D0 .and. rho_mid > rho2_flame_min .and. rho_mid < rho2_flame_max) then
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(x1 > 0.0D0 .and. rho_mid > rho2_flame_min .and. rho_mid < rho2_flame_max) then
		flame_qdot(j) = flame_ene * x1
		epsilon2 (j) = epsilon2(j) + flame_qdot(j)
		xiso (:,j) = xiso (:,j) - x1 * x_fuel1(:) + x1 * x_ash1(:) 
		burn_mass = burn_mass + x1 * rho_mid * vol_mid  
	endif 

	! Repeat the same procedure for detonation
	if(x2 > 0.0D0 .and. rho_mid > rho2_deton_min .and. rho_mid < rho2_deton_max) then
		deton_qdot(j) = deton_ene * x2 
		epsilon2(j) = epsilon2(j) + deton_qdot(j) 
		xiso (:,j) = xiso (:,j) - x2 *  x_fuel1(:) + x2 * x_ash1(:)
		burn_mass = burn_mass + x2 * rho_mid * vol_mid
	endif
	energy_input = energy_input + rho_mid*(epsilon2(j) - eps_old)*vol2(j)
enddo

! Copy the results to ghost cells
call boundary1D_NM(epsilon2, even)
call boundary2D_X

end subroutine burn_phase1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the energy release !
! By assuming 16O + 24Mg --> 28Si		!
! Written by Leung Shing Chi in 2016		!
! Updated by LEung Shing Chi in 2017		!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BURN_PHASE2
USE DEFINITION
USE NUCLEAR_MODULE
IMPLICIT NONE  

! integer variables !
integer :: j          

! Local variables !
real (DP) :: rho_mid, temp_mid, eps_old

! Change of fraction for both level sets !
real (DP) :: x1, x2

! local flame energy, defined by !
! Energy release form 16O  --> 28Si: Dqbar_O16 * (xiso_ini(o16) - Z/2)) !
! Energy release from 24Mg --> 28Si: Dqbar_Mg24 * (xiso_ini(c12) + Z/2) !
! i.e. 5.1004E-4 * (xiso_ini(o16) - Z/2) !
! i.e. 2.2037E-4 * (xiso_ini(c12) + Z/2) !
! Example: !
real (DP) :: burn_ene1 = 5.1004D-4   ! Energy release form 16O  --> 28Si
real (DP) :: burn_ene2 = 2.0237D-4   ! Energy release from 24Mg --> 28Si  

! NQSE timescale !
real (DP) :: nqse_burntime, nqse_burntime2, temp_nse, fnse

! fuel and Ash composition !
real (DP), dimension(total_ion) :: x_mid
real (DP), dimension(total_ion) :: x_fuel1 = (/0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
real (DP), dimension(total_ion) :: x_fuel2 = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0/)
real (DP), dimension(total_ion) :: x_ash = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0/)

! local deton energy            
real (DP) :: deton_ene 

! First, initialize the grid, to erase all past data
burn_qdot(:) = 0.0D0   

! Then, find the energy release due to nuclear fusion
do j = 1, length_step_part_2

	! Pass the data to local variables
	rho_mid = rho2 (j)
	temp_mid = temp2 (j)
	x_mid (:) = xiso (:, j)
	eps_old = epsilon2(j)
	if(nse_flag(j) == 2) cycle

	! Only the flame_ratio or deton_ratio = 1 can burn
	! Meaning that completely C-burnt cells can start
	! Notice the MAX nature of flame_ratio or deton_ratio
	if(burn_ratio(j) == 1.0D0 .and. rho_mid > 2.43D-11) then
		
		! When there is a change in area fraction
		! Change the chemical composition accordingly
		! And inject temperature
		! Calculate the mimimum energy timescale !
		temp_nse = -8.6803D0 + 1.8787D0 * LOG10(rho_mid/density)
		nqse_burntime = EXP(182.06D0 / temp_nse - 46.054D0) / 4.9282D-6

		! Just compute the change of fraction by considering
		! the maximum energy release equals to the amount of fuel remained
		x1 = x_mid(co16)
		x2 = x_mid(cmg24)
		if(dt > nqse_burntime) then
	   
			! Complete burning occurs
			burn_qdot(j) = burn_ene1 * x1 + burn_ene2 * x2
			epsilon2(j) = epsilon2(j) + burn_qdot(j)
			xiso(:,j) = x_mid(:) - x1 * x_fuel1(:) - x2 * x_fuel2(:) + (x1 + x2) * x_ash(:)

			! If it is completely burnt, then go to the next burning phase
			nse_flag(j) = 2

		else

			fnse =  -1.0604D0 * (dt / nqse_burntime)**2 + 2.0604D0 * (dt / nqse_burntime)
			! InComplete burning occurs
			burn_qdot(j) = (burn_ene1 * x1 + burn_ene2 * x2) * fnse
			epsilon2(j) = epsilon2(j) + burn_qdot(j)
			xiso(:,j) = x_mid(:) + (-x1 * x_fuel1(:) - x2 * x_fuel2(:) + (x1 + x2) * x_ash(:)) * fnse
			nse_flag(j) = 1

		endif
	endif
	energy_input = energy_input + rho_mid*(epsilon2(j) - eps_old)*vol2(j)

enddo  

! Copy the results to ghost cells
call boundary1D_NM(epsilon2, even)
call boundary2D_X      
       
end subroutine burn_phase2  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the 3rd phase initial burning to Nuclear statistical equilibrium !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BURN_PHASE3_INI
USE DEFINITION                               
USE NUCLEAR_MODULE
USE ECAPTABLE_MODULE
IMPLICIT NONE                                   

! Dummy variables !
integer :: i, j

! parameters !
integer, PARAMETER :: imax = 5000

! tolerance !
REAL (DP), PARAMETER :: tor = 1.0D-12

! Dummy variables !
real (DP) :: dummy, tempold, fnse

! Local variables !
real (DP) :: abar_mid, zbar_mid, ye_mid, eps_old

! Input variables !
real (DP) :: temp_beg     
real (DP) :: eps_beg

! Trial variables !
real (DP) :: temp_mid                        
real (DP) :: eps_mid         
real (DP) :: rho_mid

! Change in temperature !
real (DP) :: dtemp_mid 

! Number of successful trial !                 
integer  :: count_digit

! Ecap rate and neutrino energy loss !
real (DP) :: ecaprate_mid
real (DP) :: eneurate_mid

! Initial and Expected chemical composition !
real (DP), dimension(total_ion) :: x_mid
real (DP), dimension(total_ion) :: x_burn

! Binding energy !
real (DP) :: binde_bf, binde_af, deps_nuc

! Energy balance equation !
real (DP) :: check_e, check_e_last

! Timescale for burning to NSE !
real (DP) :: nse_burntime, temp_nse

! Initialize ecap neutrino energy lost !
total_ecap_nu_qdot = 0.0D0
        
do j = 1, length_step_part_2                 

	! Do the checking first !              
	if(temp2 (j) < 5.0D0) cycle
	if(ye2 (j) < ye_min) cycle
	if(ye2 (j) > ye_max) cycle
	if(rho2 (j) < rho2_flame_min) cycle
	if(burn_ratio (j) /= 1.0D0) cycle

	! Special test if multi-stage burning is used !
	if(nse_flag (j) == 1) cycle

	! If they pass then start the iteration !
	nse_flag (j) = 2                	

	! The density is not changed, so !
	! the trial value is also the !
	! final value !
	rho_mid = rho2 (j)

	! Give a trial temperature !
	temp_beg = temp2 (j)                  
	eps_beg = epsilon2 (j)
	eps_old = epsilon2 (j)
	
	! Also give the trial composition !
	abar_mid = abar2 (j)
	zbar_mid = zbar2 (j)
	ye_mid = ye2 (j)                
	x_mid (:) = xiso (:, j)
	
	If(temp2 (j) < 15.00) then
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Initially no electron capture !
		! Scheme for electron capture
		! Get the Ecap rate and energy loss
		!call getecaprate(rho_mid, temp_beg, ye2 (j), ecaprate_mid, eneurate_mid)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		ecaprate_mid = 0.0D0
		eneurate_mid = 0.0d0

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Initially no electron capture !
		! Update the Ye rate
		! Note: No iteration is done here because
		! the temperature sensitivity of Ecap
		! rate is much smaller than NSE composition
		!ye_mid = ye2 (j) + ecaprate_mid * dt 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

		! Compute the binding energy of the 
		! initial composition
		call compute_binde(x_mid, binde_bf)                

		! Prepare for the search of temperature
		count_digit = 0
		temp_mid = temp_beg
		tempold = temp_beg
		dtemp_mid = 0.01D0 * temp_beg          
		
		do i = 0, imax   

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! No electron capture initially !
			!call getecaprate(rho_mid, temp_mid, ye2 (j), ecaprate_mid, eneurate_mid)
			!ye_mid = ye2 (j) + ecaprate_mid * dt
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Now do the search
			! Get the trial NSE state by the trial temp
			call getnsestate(rho_mid, temp_mid, x_burn)
			!call getnse(rho_mid, temp_mid, x_burn)				

			! Compute the trial binding energy
			call compute_binde(x_burn, binde_af)

			! Calculate the trial abar and zbar
			call private_helmeos_azbar(x_burn, abar_mid, zbar_mid, dummy)

			! Get the trial epsilon
			call HELMEOS_RtoE(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid, dummy)

			! Calculate the binding energy change
			deps_nuc = binde_af - binde_bf

			! Check if the energy is balanced
			check_e = eps_mid - eps_beg - deps_nuc - 8.32696D-4 * ecaprate_mid * dt + eneurate_mid * dt

			! Make sure you go to the right direction of dtemp
			if(i == 0) then
				if(check_e > 0.0D0) dtemp_mid = -dtemp_mid
			endif

			! Use shooting method    
			if(check_e_last * check_e < 0.0D0 .and. i /= 0) then
				temp_mid = temp_mid - dtemp_mid
				dtemp_mid = dtemp_mid * 0.1D0
				temp_mid = temp_mid + dtemp_mid  
				count_digit = count_digit + 1
			else
				temp_mid = temp_mid + dtemp_mid
				check_e_last = check_e
			endif

			! My version, I need the temperature to converge !
			dummy = abs((temp_mid-tempold)/tempold)
			If (dummy < tor) THEN
				EXIT
			ELSE
				tempold = temp_mid
			END IF
			If(i == imax) THEN
				WRITE (*,*) 'imax reached', i, j
			END IF

			! Check if the nse solver fail to converge or reached low temp !
			if(i == imax .AND. temp_mid < 5.0D0) then
				temp_mid = temp_beg        
				eps_mid = eps_beg               
				call getnsestate(rho_mid, temp_mid, x_burn)
				!call getnse(rho_mid, temp_mid, x_burn)	
				exit
			endif            

		enddo

		! Find the minimum nse burn time scale !
		temp_nse = -8.6803D0 + 1.8787D0 * LOG10(rho_mid/density)
            	nse_burntime = EXP(196.02D0 / temp_nse - 41.646D0) / 4.9282D-6
		if(dt > nse_burntime) then

			! When things are done, move the refined
			! trial results to the outpuit
			temp2 (j) = temp_mid       
			epsilon2 (j) = eps_mid
			ye2 (j) = ye_mid
			burn_qdot(j) = burn_qdot(j) + binde_af - binde_bf                
			xiso (:, j) = x_burn (:) 

		else

			! When things are partially done, use
			! linear interpolation
			fnse =  -1.0604D0 * (dt / nse_burntime)**2 + 2.0604D0 * (dt / nse_burntime)
			temp_mid = temp_beg + (temp_mid - temp_beg) * fnse
			x_burn (:) = x_mid (:) + (x_burn (:) - x_mid (:)) * fnse
			call compute_binde(x_burn, binde_af)
			call private_helmeos_azbar(x_burn, abar_mid, zbar_mid, dummy)	
			call HELMEOS_RtoE(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid, dummy)	   

			! Now patch the result
			temp2 (j) = temp_mid
 			epsilon2 (j) = eps_mid
			ye2 (j) = ye_mid    
			xiso (:, j) = x_burn (:)
			burn_qdot(j) = burn_qdot(j) + binde_af - binde_bf

		endif

	else   
     
		! If temperature is too high, ignore NSE
		! composition change and focus on E-capture
		! For detonation, no NSE energy is accounted
		temp_mid = temp2 (j)
		call getecaprate(rho_mid, temp_beg, ye_mid, ecaprate_mid, eneurate_mid)
		ye2 (j) = ye_mid + ecaprate_mid * dt
		epsilon2 (j) = eps_beg + 8.32696D-4 * ecaprate_mid * dt - eneurate_mid * dt
		call invert_helm_ed(epsilon2 (j), rho2 (j), abar2 (j), zbar2 (j), ye2 (j), temp_mid, temp2 (j), epsilon2(j))
		burn_qdot(j) = 0.0D0

	endif

	! Find the total neutrino loss by electron capture !
	total_ecap_nu_qdot = total_ecap_nu_qdot + vol2 (j) * rho_mid * (eneurate_mid * dt)
      	energy_input = energy_input + rho_mid*(epsilon2(j) - eps_old)*vol2(j)

enddo

! Copy the results to ghost cells
call boundary1D_NM(temp2, even)                 
call boundary1D_NM(epsilon2, even)
call boundary1D_NM(ye2, even)
call boundary2D_X()                

end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the isotope composition  !
! due to nuclear statistical equilibrium (NSE)   !
! Written by Leung Shing Chi in 2016, but 	 !
! assuming the binding energy change due to 	 !
! the composition difference causes gain/loss in !
! local energy density				 !
! Written by Leung Shing Chi in 2016		 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE NSE
USE DEFINITION                               
USE NUCLEAR_MODULE
USE ECAPTABLE_MODULE
IMPLICIT NONE                                   

! Dummy variables !
integer :: i, j

! parameters !
integer, PARAMETER :: imax = 5000

! tolerance !
REAL (DP), PARAMETER :: tor = 1.0D-12

! Dummy variables !
real (DP) :: dummy, tempold

! Local variables !
real (DP) :: abar_mid, zbar_mid, ye_mid, eps_old

! Input variables !
real (DP) :: temp_beg     
real (DP) :: eps_beg

! Trial variables !
real (DP) :: temp_mid                        
real (DP) :: eps_mid         
real (DP) :: rho_mid

! Change in temperature !
real (DP) :: dtemp_mid 

! Number of successful trial !                 
integer  :: count_digit

! Ecap rate and neutrino energy loss !
real (DP) :: ecaprate_mid
real (DP) :: eneurate_mid

! Initial and Expected chemical composition !
real (DP), dimension(total_ion) :: x_mid
real (DP), dimension(total_ion) :: x_burn

! Binding energy !
real (DP) :: binde_bf, binde_af, deps_nuc

! Energy balance equation !
real (DP) :: check_e, check_e_last

! Timescale for burning to NSE !
real (DP) :: nse_burntime, temp_nse, fnse

! Initialize ecap neutrino energy lost !
total_ecap_nu_qdot = 0.0D0
        
do j = 1, length_step_part_2                 

	! Do the checking first !              
	if(temp2 (j) < 5.0D0) cycle
	if(ye2 (j) < ye_min) cycle
	if(ye2 (j) > ye_max) cycle
	if(rho2 (j) < rho2_flame_min) cycle
	if(burn_ratio (j) /= 1.0D0) cycle
	
	! Special test if multi-stage burning is used !
	if(nse_flag (j) == 1) cycle
	
	! If they pass then start the iteration !
	nse_flag (j) = 2                	
	
	! The density is not changed, so !
	! the trial value is also the !
	! final value !
	rho_mid = rho2 (j)
	
	! Give a trial temperature !
	temp_beg = temp2 (j)                  
	eps_beg = epsilon2 (j)
	eps_old = epsilon2 (j)
	
	! Also give the trial composition !
	abar_mid = abar2 (j)
	zbar_mid = zbar2 (j)
	ye_mid = ye2 (j)                
	x_mid (:) = xiso (:, j)
	
	If(temp2 (j) < 15.00) then

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Scheme for electron capture
		! Get the Ecap rate and energy loss
		!call getecaprate(rho_mid, temp_beg, ye2 (j), ecaprate_mid, eneurate_mid)
		! Update the Ye rate
		! Note: No iteration is done here because
		! the temperature sensitivity of Ecap
		! rate is much smaller than NSE composition
		!ye_mid = ye2 (j) + ecaprate_mid * dt   
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

		! Compute the binding energy of the 
		! initial composition
		call compute_binde(x_mid, binde_bf)                

		! Prepare for the search of temperature
		count_digit = 0
		temp_mid = temp_beg
		tempold = temp_beg
		dtemp_mid = 0.01D0 * temp_beg          

		do i = 0, imax   

			! My patch !
			call getecaprate(rho_mid, temp_mid, ye2 (j), ecaprate_mid, eneurate_mid)
			ye_mid = ye2 (j) + ecaprate_mid * dt

			! Now do the search
			! Get the trial NSE state by the trial temp
			call getnsestate(rho_mid, temp_mid, x_burn)
			!call getnse(rho_mid, temp_mid, x_burn)
			
			! Compute the trial binding energy
			call compute_binde(x_burn, binde_af)

			! Calculate the trial abar and zbar
			call private_helmeos_azbar(x_burn, abar_mid, zbar_mid, dummy)

			! Get the trial epsilon
			call HELMEOS_RtoE(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid, dummy)

			! Calculate the binding energy change
			deps_nuc = binde_af - binde_bf

			! Check if the energy is balanced
			check_e = eps_mid - eps_beg - deps_nuc - 8.32696D-4 * ecaprate_mid * dt + eneurate_mid * dt

			! Make sure you go to the right direction of dtemp
			if(i == 0) then
				if(check_e > 0.0D0) dtemp_mid = -dtemp_mid
			endif

			! Use shooting method    
			if(check_e_last * check_e < 0.0D0 .and. i /= 0) then
				temp_mid = temp_mid - dtemp_mid
				dtemp_mid = dtemp_mid * 0.1D0
				temp_mid = temp_mid + dtemp_mid  
				count_digit = count_digit + 1
			else
				temp_mid = temp_mid + dtemp_mid
				check_e_last = check_e
			endif

			! My version, I need the temperature to converge !
			dummy = abs((temp_mid-tempold)/tempold)
			If (dummy < tor) THEN
				EXIT
			ELSE
				tempold = temp_mid
			END IF
			If(i == imax) THEN
				WRITE (*,*) 'imax reached', i, j
			END IF

			! Check if the nse solver fail to converge or reached low temp !
			if(i == imax .AND. temp_mid < 5.0D0) then
				WRITE (*,*) 'temp lower than 5GK', i, imax
				temp_mid = temp_beg        
				eps_mid = eps_beg                 
				call getnsestate(rho_mid, temp_mid, x_burn)
				!call getnse(rho_mid, temp_mid, x_burn)	
				exit
			endif            

		enddo

		! Find the minimum nse burn time scale !
		temp_nse = -8.6803D0 + 1.8787D0 * LOG10(rho_mid/density)
            	nse_burntime = EXP(196.02D0 / temp_nse - 41.646D0) / 4.9282D-6
		if(dt > nse_burntime) then

			! When things are done, move the refined
			! trial results to the outpuit
			temp2 (j) = temp_mid       
			epsilon2 (j) = eps_mid
			ye2 (j) = ye_mid
			burn_qdot(j) = burn_qdot(j) + binde_af - binde_bf                
			xiso (:, j) = x_burn (:) 

		else

			! When things are partially done, use
			! linear interpolation
			fnse =  -1.0604D0 * (dt / nse_burntime)**2 + 2.0604D0 * (dt / nse_burntime)
			temp_mid = temp_beg + (temp_mid - temp_beg) * fnse
			x_burn (:) = x_mid (:) + (x_burn (:) - x_mid (:)) * fnse
			call compute_binde(x_burn, binde_af)
			call private_helmeos_azbar(x_burn, abar_mid, zbar_mid, dummy)	
			call HELMEOS_RtoE(rho_mid, temp_mid, abar_mid, zbar_mid, ye_mid, eps_mid, dummy)	   

			! Now patch the result
			temp2 (j) = temp_mid
 			epsilon2 (j) = eps_mid
			ye2 (j) = ye_mid    
			xiso (:, j) = x_burn (:)
			burn_qdot(j) = burn_qdot(j) + binde_af - binde_bf

		endif

	else   
     
		! If temperature is too high, ignore NSE
		! composition change and focus on E-capture
		! For detonation, no NSE energy is accounted
		temp_mid = temp2 (j)
		call getecaprate(rho_mid, temp_mid, ye_mid, ecaprate_mid, eneurate_mid)
		xiso(:,j) = x_burn(:)
		burn_qdot(j) = 0.0D0

	endif

	! Find the total neutrino loss by electron capture !
	total_ecap_nu_qdot = total_ecap_nu_qdot + vol2 (j) * rho_mid * (eneurate_mid * dt)
      	energy_input = energy_input + rho_mid*(epsilon2(j) - eps_old)*vol2(j)

enddo

! Copy the results to ghost cells
call boundary1D_NM(temp2, even)                 
call boundary1D_NM(epsilon2, even)
call boundary1D_NM(ye2, even)
call boundary2D_X()                

end subroutine NSE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE