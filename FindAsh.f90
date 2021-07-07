!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the mass of matter !
! which is burnt by deflagration 		!
! rates and assumes NSE composition 		!
! Written by Leung Shing Chi in 2016 		!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDASH
USE DEFINITION
USE FLAME_MODULE
IMPLICIT NONE

! Dummy variables
INTEGER :: j

! Initialization
mass_ash = 0.0D0

DO j = 1, length_step_part_2, 1

	! Just sum the mass up when the grid is burnt
	! no matter partially or completely
	mass_ash = mass_ash + vol2(j) * rho2 (j) * flame_ratio(j) 

ENDDO

END SUBROUTINE findash