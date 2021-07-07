!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine serves as a check on !
! Whether computed values are Nan      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CHECKNAN
USE DEFINITION
USE IEEE_ARITHMETIC
IMPLICIT NONE

! integer parameter !
Integer :: i, j

! Initialize logic flag !
foundnan = 0

! We find out whether there are Nan in NM density profile !
DO j = 1, length_step_part_2
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!If (isnan(rho2(j)) .eqv. .TRUE.) THEN
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO i = imin2, imax2
		if(ieee_is_nan(u_new2(i,j))) then
			foundnan = foundnan + 1
			EXIT
		END IF
	END DO
END DO

! Do the same for DM density profile !
If (RUNDM_flag == 1) THEN
	DO j = 1, length_step_part_1
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!If (isnan(rho1(j)) .eqv. .TRUE.) THEN
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DO i = imin1, imax1
			if(ieee_is_nan(u_new1(i,j))) then
				foundnan = foundnan + 1
				EXIT
			END IF
		END DO
	END DO
END IF

END SUBROUTINE