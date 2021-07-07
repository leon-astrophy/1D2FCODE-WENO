!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine sets up the position of each grid			 !
! Notice if dx is changed, this subroutine needs			 !
! to be called every time						 !
! This subroutine also finds the local volume to save computing time	 !
! Note: Here the Cart. grid style is used, meaning that the 		 !
! r- and z- position is a constant along a constant z and r respectively !
! Written by Leung Shing Chi in 2016					 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETGRID_DM
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER  :: j

DO j = -4, length_step_1 + 5
	r1L (j) = (DBLE(j) - 0.75D0)*dx1
	r1R (j) = (DBLE(j) - 0.25D0)*dx1
ENDDO
DO j = -4, length_step_1 + 5
	r1(j) = (DBLE(j) - 0.5D0) * dx1
ENDDO
DO j = -4, length_step_1 + 5
	r1F(j) = (DBLE(j)) * dx1
ENDDO
DO j = -4, length_step_1 + 5
	IF(sp_dim /= 0.0D0) THEN
		dvol1(j) = dx1*r1(j)**(sp_dim_i)
	ELSE 
		dvol1(j) = dx1
	END IF
END DO
IF(sp_dim_i == 0) THEN
	DO j = -4, length_step_1 + 5
		vol1(j) = dx1
	ENDDO
ELSEIF(sp_dim_i == 1) THEN
	DO j = -4, length_step_1 + 5
		vol1(j) = 2.0D0 * pi_old * (r1F(j)**2 - r1F(j-1)**2)
	ENDDO
ELSEIF(sp_dim_i == 2) THEN
	DO j = -4, length_step_1 + 5   
		vol1(j) = (4.0D0/3.0D0) * pi_old * (r1F(j)**3 - r1F(j-1)**3)
	ENDDO
ENDIF

END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GETGRID_NM
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER  :: j

! Get the face position
DO j = -4, length_step_2 + 5
	r2L (j) = (DBLE(j) - 0.75D0)*dx2
	r2R (j) = (DBLE(j) - 0.25D0)*dx2
ENDDO

! Get the r-position
DO j = -4, length_step_2 + 5
	r2(j) = (DBLE(j) - 0.5D0) * dx2
ENDDO

! Get the cell boundary position 
DO j = -4, length_step_2 + 5
	r2F(j) = (DBLE(j)) * dx2
ENDDO

! Get the volume elements !
IF(sp_dim_i == 0) THEN
	DO j = -4, length_step_2 + 5
		vol2(j) = dx2
	ENDDO
ELSEIF(sp_dim_i == 1) THEN
	DO j = -4, length_step_2 + 5
		vol2(j) = 2.0D0 * pi_old * (r2F(j)**2 - r2F(j-1)**2)
	ENDDO
ELSEIF(sp_dim_i == 2) THEN
	DO j = -4, length_step_2 + 5   
		vol2(j) = (4.0D0/3.0D0) * pi_old * (r2F(j)**3 - r2F(j-1)**3)
	ENDDO
ENDIF

! Differential volume !
DO j = -4, length_step_2 + 5
	IF(sp_dim /= 0.0D0) THEN
		dvol2(j) = dx2*r2(j)**(sp_dim_i)
	ELSE 
		dvol2(j) = dx2
	END IF
END DO

END SUBROUTINE 