!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine output the profiles, such as density !
! pressure, mass, epsilon ...etc                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OUTPUTPROFILEHYDRO (time)
USE DEFINITION
USE WENO_MODULE
USE NUCLEAR_MODULE
USE FLAME_MODULE
IMPLICIT NONE

! The input time !
REAL(DP), INTENT(IN) :: time

! Integer parameter !
INTEGER :: j

! Dummy !
REAL(DP) :: dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We output NM profiles !
! NM density !
WRITE (201, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (201, 701) r2(j), rho2 (j)
END DO
WRITE (201, *)

! NM epsilon !
WRITE (202, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (202, 701) r2(j), epsilon2 (j)
END DO
WRITE (202, *)

! NM pressure !
WRITE (203, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (203, 701) r2(j), p2 (j)
END DO
WRITE (203, *)

! NM velocity !
WRITE (204, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	WRITE (204, 701) r2(j), vel2 (j)
END DO
WRITE (204, *)

! NM Total Energy Density !
WRITE (205, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	dummy = epsilon2(j) + 0.5*(vel2(j)**2 + phinm_2(j)) + phinm_1(j)
	If (rho2(j) > rho2_a) THEN
		WRITE (205, 701) r2(j), dummy
	ELSE	
		WRITE (205, 701) r2(j), 0.0D0
	END IF
END DO
WRITE (205, *)

! NM Asymptotic Velocity !
WRITE (206, *) '"Time = ', time
DO j = -4, length_step_2 + 5
	dummy = vel2(j)**2 + phinm_2(j) + 2.0D0*phinm_1(j)
	If (rho2(j) > rho2_a) THEN
		WRITE (206, 701) r2(j), dummy
	ELSE	
		WRITE (206, 701) r2(j), 0.0D0
	END IF
END DO
WRITE (206, *)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Temperature For finite temperature EOS !
If (helmeos_flag == 1) THEN
	WRITE (207, *) '"Time = ', time
	DO j = -4, length_step_2 + 5
		WRITE (207, 701) r2(j), temp2 (j)
	END DO
	WRITE (207, *)
END IF

! NM electron fraction !
If (etran_flag == 1) THEN
	WRITE (208, *) '"Time = ', time
	DO j = -4, length_step_2 + 5
		WRITE (208, 701) r2(j), ye2 (j)
	END DO
	WRITE (208, *)
END IF

! We output the profiles related to gravity !
IF (outputpotential == 1) THEN

	! gravitational potential deriative !
	WRITE (210, *) '"Time = ', time
	DO j = -4, length_step_2 + 5
		WRITE (210, 701) r2(j), phinm (j)
	END DO
	WRITE (210, *)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We output the DM profile only if the DM is precense and movable !
IF (DM_flag == 1) THEN

	! DM density !
	WRITE (101, *) '"Time = ', time
	DO j = -4, length_step_1 + 5
			WRITE (101, 701) r1(j), rho1 (j)
	END DO
	WRITE (101, *)

	! DM epsilon !
	WRITE (102, *) '"Time = ', time
	DO j = -4, length_step_1 + 5
		WRITE (102, 701) r1(j), epsilon1 (j)
	END DO
	WRITE (102, *)

	! DM pressure !
	WRITE (103, *) '"Time = ', time
	DO j = -4, length_step_1 + 5
		WRITE (103, 701) r1(j), p1 (j)
	END DO
	WRITE (103, *)

	! output velocity only if DM is movable !
	IF(RUNDM_flag == 1) THEN
		WRITE (104, *) '"Time = ', time
		DO j = -4, length_step_1 + 5
			WRITE (104, 701) r1(j), vel1 (j)
		END DO
		WRITE (104, *)

		! DM Total Energy Density !
		WRITE (105, *) '"Time = ', time
		DO j = -4, length_step_1 + 5
			dummy = epsilon1(j) + 0.5D0 *(vel1(j)**2 + phidm_1(j)) + phidm_2(j)
			IF (rho1(j) > rho1_a) THEN
				WRITE (105, 701) r1(j), dummy
			ELSE	
				WRITE (105, 701) r1(j), 0.0D0
			END IF
		END DO
		WRITE (105, *)

		! DM Asymptotic Velocity !
		WRITE (106, *) '"Time = ', time
		DO j = -4, length_step_1 + 5
			dummy = vel1(j)**2 + phidm_1(j) + 2.0D0*phidm_2(j)
			IF (rho1(j) > rho1_a) THEN
				WRITE (106, 701) r1(j), dummy
			ELSE	
				WRITE (106, 701) r1(j), 0.0D0
			END IF
		END DO
		WRITE (106, *)

		IF (outputpotential == 1) THEN
			! gravitational potential !
			WRITE (209, *) '"Time = ', time
			DO j = -4, length_step_1 + 5
				WRITE (209, 702) r1(j), phidm (j)
			END DO
			WRITE (209, *)
		END IF
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Special Section For Mass Coordinate !
CALL MASSCOORDINATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

701 FORMAT (F33.15, ES33.15)
702 FORMAT (F33.15, 3ES33.15)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine output the discrete data at each time step !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OUTPUTHYDRO (time)
USE DEFINITION
USE NU_MODULE
USE WENO_MODULE
USE FLAME_MODULE
USE NUSPEC_MODULE
IMPLICIT NONE

! Real parameter !
REAL (DP) :: time, x1, x2

! Integer parameter !
INTEGER :: j

! Temporaily variables !
REAL (DP) :: centralrho1, centralrho2, centraltemp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We masses and energies !
CALL FINDEPS
CALL FINDMASS
CALL FINDENERGY

! We find the central density !
CALL FINDCENTRALDENSITY(centralrho1, centralrho2)

! We find the central temperature 
If (helmeos_flag == 1) THEN
	CALL FINDCENTRALTEMP(centraltemp2)
END IF

! Section for neutrino !
If(nuspec_flag == 1) THEN
	CALL findneutrinoloss
	CALL FindTotalNeutrinoLoss
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We output the DM data only if the DM is precense !
IF (DM_flag == 1) THEN
	WRITE (111, 702) time, centralrho1, maxval(rho1)
	WRITE (112, 701) time, mass1
	WRITE (113, 701) time, energy1
	WRITE (114, 702) time, kinenergy1, gravenergy1
	WRITE (115, 701) time, intenergy1
	WRITE (116, 701) time, total_energy
END IF

! We output NM data !
WRITE (211, 702) time, centralrho2, maxval(rho2)
WRITE (212, 701) time, mass2
WRITE (213, 702) time, energy2, energy_input
WRITE (214, 702) time, kinenergy2, gravenergy2
WRITE (215, 701) time, intenergy2
If (helmeos_flag == 1) THEN
	WRITE (216, 702) time, centraltemp2, maxval(temp2)
END IF
If (etran_flag == 1) THEN
	WRITE (217, 702) time, ye2(1), minval(ye2)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

701 FORMAT (F33.15, ES33.15)
702 FORMAT (F33.15, 2ES33.15)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write out the parameter for the ease of plotting graphs !				
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PARAMETER(ppt_in)
USE DEFINITION 
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: ppt_in

! Open file !
OPEN (UNIT = 721, FILE = './Outfile/Parameter.dat', STATUS = 'REPLACE')

! WRITE !
WRITE (721,*) total_length_1, dx1_ini
WRITE (721,*) total_length_2, dx2_ini
WRITE (721,*) ppt_in, ppt_in

! Close file!
CLOSE (721)

END SUBROUTINE