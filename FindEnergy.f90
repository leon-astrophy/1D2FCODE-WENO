!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine find the total energy (internal + Mechanical) !
! of the star assuming a spherical symmetric geometry		!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDENERGY
USE DEFINITION
IMPLICIT NONE

! Integer Parameter !
INTEGER :: j

! Intermediate variables !
REAL (DP) :: energy_1, energy_2

! Initialize !
energy1 = 0.0E0_DP
energy2 = 0.0E0_DP
energy_1 = 0.0E0_DP
energy_2 = 0.0E0_DP
intenergy1 = 0.0E0_DP
intenergy2 = 0.0E0_DP
kinenergy1 = 0.0E0_DP
kinenergy2 = 0.0E0_DP
gravenergy1 = 0.0E0_DP
gravenergy2 = 0.0E0_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The sum of energy is the sum of internal + mechanical !
DO j = 1, length_step_part_2
	IF (rho2(j) > rho2_a) THEN
		intenergy2 = intenergy2 + vol2(j) * rho2(j) * epsilon2(j)
		kinenergy2 = kinenergy2 + vol2(j) * rho2(j) * 5.0E-1_DP*vel2(j)**(2.0E0_DP)
		gravenergy2 = gravenergy2 + vol2(j) * rho2(j) * (5.0E-1_DP*phinm_2(j) + phinm_1(j))
		energy_2 = energy_2 + vol2(j) * rho2(j) * (5.0E-1_DP*(vel2(j)**(2.0E0_DP) + phinm(j)) + epsilon2(j))		
	ENDIF
END DO
	
! We sum up the energy contribution !
energy2 = kinenergy2 + gravenergy2 + intenergy2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the DM case if the DM is precense (not necessarily movable) !
IF (DM_flag == 1) THEN
	IF (RUNDM_flag == 1) THEN
		DO j = 1, length_step_part_1
			IF (rho1(j) > rho1_a) THEN
				intenergy1 = intenergy1 + vol1(j) * rho1(j) * epsilon1(j)
				kinenergy1 = kinenergy1 + vol1(j) * rho1(j) * 5.0E-1_DP*vel1(j)**(2.0E0_DP)
				gravenergy1 = gravenergy1 + vol1(j) * rho1(j) * (5.0E-1_DP*phidm_1(j) + phidm_2(j))
				energy_1 = energy_1 + vol1(j) * rho1(j) * (5.0E-1_DP*(vel1(j)**(2.0E0_DP) + phidm(j)) + epsilon1(j))			
			ENDIF
		END DO

		! We sum up the energy contribution !
		energy1 = kinenergy1 + gravenergy1 + intenergy1
	ELSE
		DO j = 1, length_step_part_1
			IF (rho1(j) > rho1_a) THEN
				intenergy1 = intenergy1 + vol1(j) * rho1(j) * epsilon1(j)
				gravenergy1 = gravenergy1 + vol1(j) * rho1(j) * (5.0E-1_DP*phidm_1(j) + phidm_2(j))
				energy_1 = energy_1 + vol1(j) * rho1(j) * (5.0E-1_DP * phidm(j) + epsilon1(j))	
			ENDIF
		END DO

		! We sum up the energy contribution !
		energy1 = gravenergy1 + intenergy1
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sum the total energy !
total_energy = energy_1 + energy_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE