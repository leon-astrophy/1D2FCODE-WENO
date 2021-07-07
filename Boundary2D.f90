!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the boundary condition for the conservative variables !
! For the meaning of each boundary flag (0, 1, 2), Please refer !
! to the parameter.h file. This is for the DM boundary 2D       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY2D_DM (u)
USE DEFINITION
IMPLICIT NONE

! The temporaily array used in this subroutine for boundary condition !
REAL (DP), INTENT (INOUT), DIMENSION (imin1 : imax1, -4 : length_step_1 + 5) :: u

! Integer variable !
INTEGER :: i, j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! do the inner boundary !
IF(boundary_flag(1) == 0) THEN
	DO i = imin1, imax1
		DO j = 1, 5
			u(i, 1 - j) = u(i, length_step_1 + 1 - j)
		ENDDO
	ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
	DO i = imin1, imax1
		DO j = 1, 5
			u(i, 1 - j) = bfac(i) * u(i, j)
		ENDDO
	ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
	DO i = imin1, imax1
		DO j = 1, 5
			u(i, 1 - j) = u(i , 1)
		ENDDO
	ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the outer boundary !
IF(boundary_flag(2) == 0) THEN
	DO i = imin1, imax1
		DO j = 1, 5
			u(i,length_step_1 + j) = u(i,j)
		ENDDO
	ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
	DO i = imin1, imax1
		DO j = 1, 5
			u(i, length_step_1 + j) = bfac(i) * u(i, length_step_1 + 1 - j)
		ENDDO
	ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
	DO i = imin1, imax1
		DO j = 1, 5
			u(i,length_step_1 + j) = u(i,length_step_1)                  
		ENDDO
	ENDDO
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This boundary subroutine is for the NM conservative variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY2D_NM (u)
USE DEFINITION
IMPLICIT NONE

! The temporaily array used in this subroutine for boundary condition !
REAL (DP), INTENT (INOUT), DIMENSION (imin2 : imax2, -4 : length_step_2 + 5) :: u

! Integer variable !
INTEGER :: i, j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! do the inner boundary !
IF(boundary_flag(1) == 0) THEN
	DO i = imin2, imax2
		DO j = 1, 5
			u(i, 1 - j) = u(i, length_step_2 + 1 - j)
		ENDDO
	ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
	DO i = imin2, imax2
		DO j = 1, 5
			u(i, 1 - j) = bfac(i) * u(i, j)
		ENDDO
	ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
	DO i = imin2, imax2
		DO j = 1, 5
			u(i, 1 - j) = u(i , 1)
		ENDDO
	ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the outer boundary !
IF(boundary_flag(2) == 0) THEN
	DO i = imin2, imax2
		DO j = 1, 5
			u(i,length_step_2 + j) = u(i,j)
		ENDDO
	ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
	DO i = imin2, imax2
		DO j = 1, 5
			u(i, length_step_2 + j) = bfac(i) * u(i, length_step_2 + 1 - j)
		ENDDO
	ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
	DO i = imin2, imax2
		DO j = 1, 5
			u(i,length_step_2 + j) = u(i,length_step_2)                  
		ENDDO
	ENDDO
ENDIF

END SUBROUTINE