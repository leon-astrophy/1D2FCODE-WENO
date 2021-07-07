!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine set the boundary conditions for the nuclear isotopes !			  
! for the meaning of boundary flag, please refer to parameter.h	       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY2D_X()
USE DEFINITION
USE NUCLEAR_MODULE
IMPLICIT NONE

! integer variables !
integer :: i, j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We first do the left boundary !
IF(boundary_flag(1) == 0) THEN
	DO j = 1, 5
		xiso(:,1-j) = xiso(:,length_step_2 + 1 - j)
	ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
	DO j = 1, 5
		xiso(:,1-j) = xiso(:,j)
	ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
	DO j = 1, 5
		xiso(:,1-j) = xiso(:,1)
	ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We then do the right boundary !
IF(boundary_flag(2) == 0) THEN
	DO j = 1, 5
		xiso(:,length_step_2 + j) = xiso(:,j)
	ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
	DO j = 1, 5
		xiso(:,length_step_2 + j) = xiso(:,length_step_2 + 1 - j)
	ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
	DO j = 1, 5
		xiso(:,length_step_2 + j) = xiso(:,length_step_2)                  
	ENDDO
ENDIF

END SUBROUTINE boundary2D_X