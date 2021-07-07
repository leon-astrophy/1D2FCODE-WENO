!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine find the maximum time step dt that !
! statisfy the CFL condition of stability	     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDDT
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
integer :: i

! Real parameter !
! Local effective speed for DM and NM !
REAL(DP) :: lambda_locdm, lambda_locnm 

! Local maximum effective speed !
REAL(DP) :: lambda_loc 

! Local minimum time step for DM and NM !
REAL(DP) :: dt_dm, dt_nm 

! local time step !
REAL(DP) :: dt_1, dt_2, dt_3, dt_4, dt_5, dt_temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize by setting an arbitrarily large number !
dt_dm = 1.0E5_DP
dt_nm = 1.0E5_DP
dt_1 = 1.0E5_DP
dt_2 = 1.0E5_DP
dt_3 = 1.0E5_DP
dt_4 = 1.0E5_DP
dt_5 = 1.0E5_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the find dt for DM case !
! We do it only if DM is precense !
! for movable DM, we also include the fluid velocity !
IF (DM_flag == 1) THEN
	IF (RUNDM_flag == 1) THEN
		DO i=1, length_step_1
			IF(rho1(i) > rho1_a) THEN
				lambda_locdm = abs(vel1(i)) + cs1(i) 
				dt_dm = MIN(cfl * dx1 / lambda_locdm, dt_dm)
			END IF
		END DO	

	! or else, we include only the local sound speed !
	ELSE
		DO i=1, length_step_1		
			IF(rho1(i) > rho1_a) THEN		
				lambda_locdm = cs1(i) 
				dt_dm = MIN(cfl * dx1 / lambda_locdm, dt_dm)
			END IF
		END DO
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the find dt for NM case !
! Only grid with density above atmospheric density is counted !
DO i=1, length_step_2
	IF(rho2(i) > rho2_a) THEN
		lambda_locnm = abs(vel2(i)) + cs2(i) 
		dt_nm = MIN(cfl * dx2 / lambda_locnm, dt_nm)
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We choose the minimum dt accordingly !
dt_temp = MIN(dt_nm, dt_dm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We adjust dt to match the data output time !
!IF(global_time + dt_temp - last_outputtime >= output_time) THEN
!	dt_1 = MIN(output_time - (global_time - last_outputtime), dt_temp)
!ENDIF
! We adjust dt to match the data output time !
!IF(global_time + dt_temp - last_outputtime_profile1 >= output_time_profile1 * 1.0E1_DP ** (proindex1)) THEN
!	dt_2 = MIN(output_time_profile1 * 1.0E1_DP ** (proindex1) - (global_time - last_outputtime_profile1), dt_temp)
!ENDIF
! We adjust dt to match the data output time !
!IF(global_time + dt_temp - last_outputtime_profile2 >= output_time_profile2 * 1.0E1_DP ** (proindex2)) THEN
!	dt_3 = MIN(output_time_profile2 * 1.0E1_DP ** (proindex2) - (global_time - last_outputtime_profile2), dt_temp)
!ENDIF
! Check if we reach the maximum time !
!IF(global_time + dt_temp >= total_time) THEN
!	dt_5 = MIN(total_time - global_time, dt_temp)
!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set the output !
dt = MIN(dt_1, dt_2, dt_3, dt_4, dt_5, dt_temp)

END SUBROUTINE