!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determine the alpha, maximum       !
! eigenvalue of the jacobian matrix or system of     !
! conservation laws, alpha will be supplied to the   !				         
! riemann solver of lax-freidrich		     !
! The eigenvalue, according to Toro's book, are      !
! u, u-c u+c					     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ALPHASPLIT (alpha)
USE DEFINITION
USE WENO_MODULE
IMPLICIT NONE

! integer parameter !
INTEGER :: i, j

! the alpha that we want !
REAL (DP), INTENT (OUT) :: alpha(no_of_eq)

! effective local speed for DM and NM !
REAL (DP) :: lambda(2)

! The maximum local effective speed for DM and NM!
REAL (DP) :: lambda_max(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We initialize the maximum local speed !
lambda_max = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we do the NM part !   
DO j = 1, length_step_part_2    
	lambda(2) = ABS(vel2(j)) + cs2(j) 
	lambda_max(2) = MAX(lambda(2), lambda_max(2))
END DO

! We assign the value of alpha to NM equaton !
Do j = imin2, imax2
	alpha(j) = lambda_max(2)
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we do the DM part only if the DM is precense and movable !
IF (RUNDM_flag == 1) THEN

	! We find the maximum local speed for DM first !
	DO j = 1, length_step_part_1
		lambda(1) = ABS(vel1(j)) + cs1(j) 
		lambda_max(1) = MAX(lambda(1), lambda_max(1))
	END DO

	! We assign the value of alpha to DM equation !
	Do j = imin1, imax1
		alpha(j) = lambda_max(1)
	END DO

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE