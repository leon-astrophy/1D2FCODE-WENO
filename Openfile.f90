!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine open files that contains essential !
! Hydro variables of the simulation, such as density !
! velocity, specific energy density, etc......       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENFILEHYDRO
USE DEFINITION
IMPLICIT NONE

! Integer variables !
INTEGER ::  fileno_len
CHARACTER (len = 256) :: fileno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE (fileno, *) n_backup / time_step + 1
fileno = ADJUSTL (fileno)
fileno_len = LEN_TRIM (fileno)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We open specific files only if DM is precense !
IF (DM_flag == 1) THEN
	
	! DM Profile !
	OPEN (UNIT = 100, FILE = './Outfile/Hydro/Star_WENO_MassCoordinate_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 101, FILE = './Outfile/Hydro/Star_WENO_Density_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 102, FILE = './Outfile/Hydro/Star_WENO_Epsilon_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 103, FILE = './Outfile/Hydro/Star_WENO_Pressure_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')

	! DM Velocity is needed only if DM is movable !
	IF (RUNDM_flag == 1) THEN
		OPEN (UNIT = 104, FILE = './Outfile/Hydro/Star_WENO_Velocity_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
		OPEN (UNIT = 105, FILE = './Outfile/Hydro/Star_WENO_TotalEnergyDensity_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
		OPEN (UNIT = 106, FILE = './Outfile/Hydro/Star_WENO_AsymptoticVelocity_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	END IF
			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! DM Log !
	OPEN (UNIT = 111, FILE = './Outfile/Hydro/Star_WENO_CentralDensity_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 112, FILE = './Outfile/Hydro/Star_WENO_Mass_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 113, FILE = './Outfile/Hydro/Star_WENO_Energy_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 114, FILE = './Outfile/Hydro/Star_WENO_MechEnergy_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 115, FILE = './Outfile/Hydro/Star_WENO_IntEnergy_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 116, FILE = './Outfile/Hydro/Star_WENO_TotalEnergy_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! NM Profile !
OPEN (UNIT = 200, FILE = './Outfile/Hydro/Star_WENO_MassCoordinate_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 201, FILE = './Outfile/Hydro/Star_WENO_Density_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 202, FILE = './Outfile/Hydro/Star_WENO_Epsilon_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 203, FILE = './Outfile/Hydro/Star_WENO_Pressure_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 204, FILE = './Outfile/Hydro/Star_WENO_Velocity_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 205, FILE = './Outfile/Hydro/Star_WENO_TotalEnergyDensity_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 206, FILE = './Outfile/Hydro/Star_WENO_AsymptoticVelocity_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')

! We output temperature if using finite temp EOS !
If (helmeos_flag == 1) THEN
	OPEN (UNIT = 207, FILE = './Outfile/Hydro/Star_WENO_Temperature_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
END IF

! We output electron fraction !
If (etran_flag == 1) THEN
	OPEN (UNIT = 208, FILE = './Outfile/Hydro/Star_WENO_Ye_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
END IF

! We output the files related to gravity only if you allowed so !
IF (outputpotential == 1) THEN
	IF(RUNDM_flag == 1) THEN
		OPEN (UNIT = 209, FILE = './Outfile/Hydro/Star_WENO_Potential_DM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
	END IF
	OPEN (UNIT = 210, FILE = './Outfile/Hydro/Star_WENO_Potential_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! NM Log !
OPEN (UNIT = 211, FILE = './Outfile/Hydro/Star_WENO_CentralDensity_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 212, FILE = './Outfile/Hydro/Star_WENO_Mass_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 213, FILE = './Outfile/Hydro/Star_WENO_Energy_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 214, FILE = './Outfile/Hydro/Star_WENO_MechEnergy_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 215, FILE = './Outfile/Hydro/Star_WENO_IntEnergy_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')

! We output temperature if using finite temp EOS !
If (helmeos_flag == 1) THEN
	OPEN (UNIT = 216, FILE = './Outfile/Hydro/Star_WENO_CentralTemperature_NM_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
END IF

! We output electron fraction !
If (etran_flag == 1) THEN
	OPEN (UNIT = 217, FILE = './Outfile/Hydro/Star_WENO_CentralYe_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine close files that contains essential !
! Hydro variables of the simulation, such as density  !
! velocity, specific energy density, etc......        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CLOSEFILEHYDRO
USE DEFINITION
IMPLICIT NONE

! Integer variables !
INTEGER :: m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Close DM output files !
IF (DM_flag == 1) THEN
	
	! Profiles And Log !
	DO m = 100, 103
		CLOSE (m)
	END DO
	DO m = 111, 106
		CLOSE (m)
	END DO

	! For RUNDM !
	IF (RUNDM_flag == 1) THEN
		DO m = 104, 106
			CLOSE (m)
		END DO
	END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Close NM output files !
DO m = 200, 206
	CLOSE (m)
END DO
DO m = 211, 215
	CLOSE (m)
END DO

! For others flag !
If (helmeos_flag == 1) THEN
	CLOSE (207)
	CLOSE (216)
END IF

If (etran_flag == 1) THEN
	CLOSE (208)
	CLOSE (217)
END IF

IF (outputpotential == 1) THEN
	IF(RUNDM_flag == 1) THEN
		CLOSE (209)
	END IF
	CLOSE (210)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE