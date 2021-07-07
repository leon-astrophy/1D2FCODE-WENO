!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the boundary condition for perimitive variables !
! such as density, velocity, potential, pressure...       !
! In this version, I choose to use the full boundary      !
! box for boundary 1D, this is for DM primitive variable  !					  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY1D_DM (array, sign)
USE DEFINITION
IMPLICIT NONE

! the temporaily array for boundary condition !
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_1 + 5) :: array

! the input sign !
INTEGER, INTENT (IN) :: sign

! integer variables !
INTEGER :: j, fac

! Set up the parity factor according to the input sign
IF(sign == 0) THEN
	fac = 1
ELSEIF(sign == 1) THEN
	fac = -1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We assign the boundary values to inner boundary !
IF(boundary_flag(1) == 0) THEN
	DO j = 1, 5
		array(1 - j) = array(length_step_1 + 1 - j)
	ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
	DO j = 1, 5
		array(1 - j) = fac * array(j)
	ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
	DO j = 1, 5
		array(1 - j) = array(1)
	ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We assign the boundary conditions to outer boundary !
IF(boundary_flag(2) == 0) THEN
	DO j = 1, 5
		array(length_step_1 + j) = array(j)
	ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
	DO j = 1, 5
		array(length_step_1 + j) = fac * array(length_step_1 + 1 - j)
	ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
	DO j = 1, 5
		array(length_step_1 + j) = array(length_step_1)                  
	ENDDO
ENDIF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is boundary subroutines for DM primitive variable !					  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY1D_NM (array, sign)
USE DEFINITION
IMPLICIT NONE

! the temporaily array for boundary condition !
REAL (DP), INTENT (INOUT), DIMENSION (-4 : length_step_2 + 5) :: array

! the input sign !
INTEGER, INTENT (IN) :: sign

! integer variables !
INTEGER :: j, fac

! Set up the parity factor according to the input sign
IF(sign == 0) THEN
	fac = 1
ELSEIF(sign == 1) THEN
	fac = -1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We assign the boundary values to inner boundary !
IF(boundary_flag(1) == 0) THEN
	DO j = 1, 5
		array(1 - j) = array(length_step_2 + 1 - j)
	ENDDO
ELSEIF(boundary_flag(1) == 1) THEN
	DO j = 1, 5
		array(1 - j) = fac * array(j)
	ENDDO
ELSEIF(boundary_flag(1) == 2) THEN
	DO j = 1, 5
		array(1 - j) = array(1)
	ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We assign the boundary conditions to outer boundary !
IF(boundary_flag(2) == 0) THEN
	DO j = 1, 5
		array(length_step_2 + j) = array(j)
	ENDDO
ELSEIF(boundary_flag(2) == 1) THEN
	DO j = 1, 5
		array(length_step_2 + j) = fac * array(length_step_2 + 1 - j)
	ENDDO
ELSEIF(boundary_flag(2) == 2) THEN
	DO j = 1, 5
		array(length_step_2 + j) = array(length_step_2)                  
	ENDDO
ENDIF

END SUBROUTINE