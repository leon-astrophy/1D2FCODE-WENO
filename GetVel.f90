!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine assign the initial velocity to both nm and dm !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETVEL
USE DEFINITION
IMPLICIT NONE

! Integer Parameter !
INTEGER :: j

! Get velocity !
DO j = 1, length_step_2
	vel2(j) = 0.0D0
END DO

! Copy the value of velocity to ghost shell !
CALL BOUNDARY1D_NM (vel2, odd)

! For DM !
If (rundm_flag == 1) THEN
	DO j = 1, length_step_1
		vel1(j) = 0.0E0_DP
	END DO
	CALL BOUNDARY1D_DM (vel1, odd)
END IF

END SUBROUTINE