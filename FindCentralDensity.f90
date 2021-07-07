!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine find the central density of NM and DM by interpolation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDCENTRALDENSITY(centralrho1, centralrho2)
USE DEFINITION
IMPLICIT NONE

! Real output !
REAL (DP), INTENT(OUT) :: centralrho1, centralrho2

! Interpolate DM central density !
! We do this only if the users wants DM component !
If (DM_flag == 1) THEN
	centralrho1 = rho1(1) 
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Old interpolation scheme muted !
!((1.5E0_DP * dx1) ** 2 * rho1 (1) - (5.0E-1_DP * dx1) ** 2 * rho1 (2)) &
!/ ((1.5E0_DP * dx1) ** 2 - (5.0E-1_DP * dx1) ** 2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Interpolate NM central density !
centralrho2 = rho2(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Old interpolation scheme muted !
!((1.5E0_DP * dx2) ** 2 * rho2 (1) - (5.0E-1_DP * dx2) ** 2 * rho2 (2)) &
!/ ((1.5E0_DP * dx2) ** 2 - (5.0E-1_DP * dx2) ** 2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine find the central temperature of NM !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDCENTRALTEMP(centraltemp2)
USE DEFINITION
IMPLICIT NONE

! Real output !
REAL (DP), INTENT(OUT) :: centraltemp2

! Interpolate NM central density !
centraltemp2 = temp2(1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Old interpolation scheme muted !
!((1.5E0_DP * dx2) ** 2 * temp2 (1) - (5.0E-1_DP * dx2) ** 2 * temp2 (2)) &
!/ ((1.5E0_DP * dx2) ** 2 - (5.0E-1_DP * dx2) ** 2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

END SUBROUTINE