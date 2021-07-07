!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine find the total mass of fluid assuming a spherical symmetry !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDMASS
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
INTEGER :: j, grid_j

! Dummy variables !
REAL (DP) :: integrand 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize normal matter mass !
mass2 = 0.0E0_DP

DO j = 1, length_step_part_2
	If (rho2(j) > rho2_a) THEN
		mass2 = mass2 + vol2(j) * rho2 (j)
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now do the DM case only if the users wants dm component !
if (DM_flag == 1) Then
	mass1 = 0.0E0_DP
	DO j = 1, length_step_part_1
		If (rho1(j) > rho1_a) THEN
			mass1 = mass1 + vol1(j) * rho1 (j)
		END IF
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine find the total ejected mass and their	      !
! corresponding energy at the end of the supernova simulation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EJECTEDMASS
USE NUCLEAR_MODULE
USE NU_MODULE
USE DEFINITION
IMPLICIT NONE

! Integer Parameter !
INTEGER :: j

! Real Parameter !
REAL (DP) :: ejected_1
REAL (DP) :: ejected_2
REAL (DP) :: iron_peak
REAL (DP) :: dummy
REAL (DP) :: c2erg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
ejected_1 = 0.0D0
ejected_2 = 0.0D0
iron_peak = 0.0D0
c2erg = 5.59429D-55

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sum all NM ejected mass and ejected iron peak elements !
DO j = 1, length_step_part_2
	IF (rho2(j) > rho2_a) THEN
		! If the total energy of NM is positive, then count only the DM potential !
		If(energy2 >= 0.0D0) THEN
			dummy = 5.0E-1_DP*vel2(j)**(2.0E0_DP) + phinm_1 (j) + epsilon2(j)
		ELSE
			dummy = 5.0E-1_DP*(vel2(j)**(2.0E0_DP) + phinm_2 (j)) + phinm_1(j) + epsilon2(j)
		END IF
		IF (dummy >= 0.0E0_DP) THEN
			ejected_2 = ejected_2 + vol2(j) * rho2 (j)
			iron_peak = iron_peak + vol2(j) * rho2 (j) * xiso(cni56, j)
		END IF
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sum all DM ejected if they are movable !
If (RUNDM_flag == 1) THEN
	DO j = 1, length_step_part_1
		IF (rho1(j) > rho1_a) THEN
			! If the total energy of DM is positive, then count only the NM potential !
			If(energy1 >= 0.0D0) THEN
				dummy = 5.0E-1_DP*vel1(j)**(2.0E0_DP) + phidm_2 (j) + epsilon1(j)
			ELSE
				dummy = 5.0E-1_DP*(vel1(j)**(2.0E0_DP) + phidm_1 (j)) + phidm_2(j) + epsilon1(j)
			END IF
			IF (dummy >= 0.0E0_DP) THEN
				ejected_1 = ejected_1 + vol1(j) * rho1 (j)
			END IF
		END IF
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Output the ejected mass !	
OPEN (UNIT = 999, FILE = './Outfile/Star_WENO_Ejecta.dat', STATUS = 'REPLACE')
WRITE(999,*) '-------------------'
WRITE(999,*) 'EXPLOSION RESULTS :'
WRITE(999,*) '-------------------'
WRITE(999,*) 'NM Bound Mass ', mass2 - ejected_2, ' Solar Mass'
WRITE(999,*) 'NM Ejected Mass ', ejected_2, ' Solar Mass'
WRITE(999,*) 'NM Ejected Ni56 ', iron_peak, ' Solar Mass'
IF(NUSPEC_flag == 1) THEN
	WRITE(999,*) 'Time Integrated Neutrino Loss ', SUM(total_loss)/c2erg, ' Erg'
	WRITE(999,*) 'Time Integrated Neutrino Production ', ' Pair Neutrino ', ' Plasma Neutrino '
	DO j = 1, 5 
		WRITE(999,*) 'Number ', total_pair(j), total_plas(j), j, 'MeV'
	END DO
END IF
IF(RUNDM_flag ==1) THEN
	WRITE(999,*) 'DM Bound Mass ', mass1 - ejected_1, ' Solar Mass'
	WRITE(999,*) 'DM Ejected Mass ', ejected_1, ' Solar Mass'
END IF
CLOSE (999)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

701 FORMAT (A, 5ES33.15, A)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine find the mass coordinate with respect to     !
! the radial coordinate, as an input to the light curve solver !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MASSCOORDINATE
USE DEFINITION
IMPLICIT NONE

! Integer Parameter !
INTEGER :: j, k

! Temporial array to store mass !
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: m_face1
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: m_center1
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: m_face2
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: m_center2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
m_face2 = 0.0D0
m_center2 = 0.0D0

! DO NM first Assign mass faced coordinate !
DO j = 1, length_step_2
	DO k = 1, j
		If (rho2(k) > rho2_a) THEN
			m_face2 (j) = m_face2(j) + vol2(k) * rho2(k)
		ELSE
			CYCLE
		END IF
	END DO 
ENDDO

! Assign mass centered coordinate by interpolation !
DO j = 1, length_step_2
	m_center2(j) = 0.5D0*(m_face2 (j) + m_face2 (j-1))
ENDDO

! Output !	
WRITE (200, *) '"Time = ', global_time
DO j = -4, length_step_2 + 5
	WRITE (200, 701) r2(j), m_center2(j)
END DO
WRITE (200, *)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do DM !
If(DM_flag == 1) THEN

	! Initialize !
	m_face1 = 0.0D0
	m_center1 = 0.0D0

	! DO DM mass coordinate !
	DO j = 1, length_step_1
		DO k = 1, j
			If (rho1(k) > rho1_a) THEN
				m_face1 (j) = m_face1(j) + vol1(k) * rho1(k)
			ELSE
				CYCLE
			END IF
		END DO 
	ENDDO

	! Assign mass centered coordinate by interpolation !
	DO j = 1, length_step_1
		m_center1(j) = 0.5D0*(m_face1 (j) + m_face1 (j-1))
	ENDDO

	! Output !	
	WRITE (100, *) '"Time = ', global_time
	DO j = -4, length_step_1 + 5
		WRITE (100, 701) r1(j), m_center1(j)
	END DO
	WRITE (100, *)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

701 FORMAT (F33.15, ES33.15)

END SUBROUTINE