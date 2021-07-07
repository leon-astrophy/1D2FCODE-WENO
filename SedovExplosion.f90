!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine find the position where the enclosed mass !
! is large enough for the initial perturbation in epsilon   !
! or velocity so that the star could be "exploded" without  !
! being code crashed. We usually choose M ~ 0.01 to 0.1     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDRADIUS(n2)
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
INTEGER :: j

! Real parameter !
REAL (DP) :: temp

! Output parameter !
Integer, INTENT(OUT) :: n2

! We find the radius of the star !
DO j = 1, length_step_part_2
	If (rho2(j) <= rho2_a) THEN
		temp = DBLE(j)
		EXIT
	end if
END DO

! assign the output !
n2 = INT(temp/5.0E0_DP)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do a similar blast test as suggested	     !	
! by Sedov. We include a sudden amount of energy around the  !
! star center, the energy should be large enough to overcome !
! gravitational energy binding. Note that for sedov test, we !
! assume a helmeos EOS so make sure you choose EOS correctly !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SEDOV
USE DEFINITION
IMPLICIT NONE

! integer parameter !
INTEGER :: j, temp, n2

! Real Parameter !
REAL (DP) :: energy_in

! We need to first find the gravitational potential 1
If (sp_dim_i == 2) THEN
	CALL FINDPOTENTIAL
END IF

! We find the mass !
CALL FINDMASS

! We find the gravitational binding energy !
CALL FINDENERGY

! We need to find the radius where we apply the perturbation !
CALL FINDRADIUS(n2)
	
! We need to increase the applied energy to explode the star !
! we only include for NM !
If( DM_flag == 1) THEN
	energy_in = abs(5.0D0 * (energy1 + energy2))
	WRITE (*,*) energy1 + energy2, 'total energy'
ELSE
	energy_in = abs(5.0D0 * energy2)
	WRITE (*,*) energy2, 'total energy'
END IF

WRITE (*,*) energy_in/4.0E0_DP, 'input energy'
WRITE (*,*) n2, r2(n2), 'radius of sedov'

! We need to update the internal energy of the first n2 grid !
Do j = 1, n2
	epsilon2(j) = epsilon2(j) + (3.0E0_DP * energy_in * (r2(n2) - r2(j)))&
		 / (rho2(j) * 4.0E0_DP * pi_old * (r2(n2))**(4.0E0_DP)) 	
END DO

! We copy the results to ghost shell !
CALL BOUNDARY1D_NM (epsilon2,even)

END SUBROUTINE