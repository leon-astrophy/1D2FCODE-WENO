!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine update the pressure source, temperature !
! and find the gravitaional potential for next RK step	  !                          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE UPDATE
USE DEFINITION
USE NUCLEAR_MODULE
IMPLICIT NONE

! Integer parameter !
INTEGER :: i, j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
! Since we assume spherical symmetry, the potential !
! Can only be found if you set to spherical symmetry !
! If would be nice if you add the other coordinate version !
CALL FINDPOTENTIAL

! For finite temperautere EOS !
IF(helmeos_flag == 1) THEN 
	CALL FINDHELMTEMP
END IF

! Find Pressure !
CALL FINDPRESSURE

! Find internal energy rhoe !
IF(dual_energy == 1) THEN
	CALL FINDRHOE
END IF

! Find sound speed ! 
CALL FINDSOUNDSPEED

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine backup the conservative variables for the !
! next temporal evolution				    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BACKUPCONS (v1, u1, v2, u2)
USE DEFINITION
USE WENO_MODULE
IMPLICIT NONE

! The numerical value of v in the WENO reconstruction !
REAL (DP), INTENT (IN), DIMENSION (imin1 : imax1, -4 : length_step_1 + 5) :: v1
REAL (DP), INTENT (IN), DIMENSION (imin2 : imax2, -4 : length_step_2 + 5) :: v2

! The numerical value of u in the hyperbolic convervation law !
REAL (DP), INTENT (OUT), DIMENSION (imin1 : imax1, -4 : length_step_1 + 5) :: u1
REAL (DP), INTENT (OUT), DIMENSION (imin2 : imax2, -4 : length_step_2 + 5) :: u2

! Integer parameter !
INTEGER :: i, j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! These are the neccesarry variables that need to be updated each time !
! Assign the value of v !
IF(RUNDM_flag == 1) THEN
	DO i = imin1, imax1
		DO j = -4, length_step_1 + 5
			u1 (i, j) = v1 (i, j)
		END DO
	END DO
END IF

DO i = imin2, imax2
	DO j = -4, length_step_2 + 5
		u2 (i, j) = v2 (i, j)
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE