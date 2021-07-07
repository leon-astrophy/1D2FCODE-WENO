!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine update the pressure profile and their deriative    !
! once the density profile is being updated through rungekutta time  !
! evolution. It is being used in every time step, do not confused it !
! with subroutine GETRHOEOSRTOP                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDPRESSURE
USE DEFINITION
USE NUCLEAR_MODULE
IMPLICIT NONE

! Integer parameter !
INTEGER :: j

! Dummy variables !
REAL (DP) :: dummy, dxdrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the DM case first !
! This would be done only if users wants DM component !
IF (DM_flag == 1) THEN
	IF (fermieosdm_flag == 1) Then
		DO j = 1, length_step_1
			CALL FERMIMO (dlfmmo1, rho1 (j), 1)
			CALL FINDDXDRHO (dxdrho, rho1 (j), 1)
			IF (dlfmmo1 <= 1.0E-2_DP) THEN
				p1 (j) = a_max1*small_pressure(dlfmmo1)
			ELSE
				p1 (j) = a_max1*large_pressure(dlfmmo1)
			END IF

			! Derivative !
			dpdrho1 (j) = a_max1*dxdrho*dpdx(dlfmmo1)
			dpdepsilon1 (j) = 0.0E0_DP
		END DO
	ELSE
		DO j = 1, length_step_1
			p1 (j) = k1 * rho1(j) ** gamma1
			dpdrho1 (j) = k1 * gamma1 * rho1(j) ** (gamma1 - 1.0D0)
			dpdepsilon1 (j) = 0.0D0
		END DO	
	END IF
	CALL BOUNDARY1D_DM (p1, even)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The following steps are more or less similar , so no repeat !
IF (helmeos_flag /= 1) THEN
	IF (fermieosnm_flag == 1) THEN
		DO j = 1, length_step_2
			CALL FERMIMO (dlfmmo2, rho2 (j), 2)
			CALL FINDDXDRHO (dxdrho, rho2 (j), 2)
			IF (dlfmmo2 <= 1.0E-2_DP) THEN
				p2 (j) = a_max2*small_pressure(dlfmmo2)
			ELSE
				p2 (j) = a_max2*large_pressure(dlfmmo2)
			END IF

			! Derivative !
			dpdrho2 (j) = a_max2*dxdrho*dpdx(dlfmmo2)
		END DO
	ELSEIF (ccsneos_flag == 1) THEN
                DO j = 1, length_step_2

			! polytrope epsilon !
			IF (rho2(j) < rhoc_b) THEN
				eps_p (j) = Ec_1 * rho2(j) ** (gammac_1 - 1.0E0_DP)
			ELSE
				eps_p (j) = Ec_2 * rho2(j) ** (gammac_2 - 1.0E0_DP) + Ec_3
			END IF

			! thermal epsilon !
			eps_th (j) = epsilon2 (j) - eps_p (j)

			! polytrope pressure !
			IF (rho2(j) < rhoc_b) THEN
				p_p (j) = kc_1 * rho2(j) ** (gammac_1)
			ELSE
				p_p (j) = kc_2 * rho2(j) ** (gammac_2)
			END IF

			! thermal pressure !
			p_th (j) = rho2(j) * eps_th (j) * (gammac_th - 1.0E0_DP)

			! pressure derivatives !
			p2 (j) = p_p (j) + p_th (j)
			!IF(p2(j) < 0.0D0) THEN
			!	p2(j) = p_p (j)
			!END IF
			IF (rho2(j) < rhoc_b) THEN
				dpdrho2 (j) = kc_1 * gammac_1 * rho2(j) ** (gammac_1 - 1.0E0_DP) + eps_th(j) * (gammac_th - 1.0E0_DP)
			ELSE
				dpdrho2 (j) = kc_2 * gammac_2 * rho2(j) ** (gammac_2 - 1.0E0_DP) + eps_th(j) * (gammac_th - 1.0E0_DP)
			END IF

			! epsilon derivatives !
			dpdepsilon2 (j) = rho2(j) * (gammac_th - 1.0E0_DP)

		END DO
	ELSE	
		IF(nm_epsilon == 0) THEN
			DO j = 1, length_step_2
				p2 (j) = k2 * rho2(j) ** gamma2
				dpdrho2 (j) = k2 * gamma2 * rho2(j) ** (gamma2 - 1.0D0)
				dpdepsilon2 (j) = 0.0D0
			END DO	
		ELSE
			DO j = 1, length_step_2
				p2 (j) = rho2(j) * epsilon2(j) * (gamma2 - 1.0E0_DP) 
				dpdrho2 (j) = epsilon2(j) * (gamma2 - 1.0E0_DP)
				dpdepsilon2 (j) = rho2(j) * (gamma2 - 1.0E0_DP)
			END DO
		END IF
	END IF
ELSE
	DO j = 1, length_step_2
		CALL HELMEOS_RtoP(rho2 (j), temp2(j), abar2(j), zbar2(j), ye2(j), p2 (j), dpdrho2 (j), dpdepsilon2 (j), chempo2(j))
	END DO	
END IF

! Copy to boundary !
CALL BOUNDARY1D_NM (p2, even)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
	real(DP) function large_pressure(x)
	implicit none
	real(DP) :: x
	large_pressure = x*SQRT(x**2 + 1.0D0)*(2.0D0*x**2 - 3.0D0) + 3.0D0*log(x + SQRT(x**2 + 1.0D0))
	end function

	real(DP) function dpdx(x)
	implicit none
	real(DP) :: x
	dpdx = 8.0D0*x**4/SQRT(x**2 + 1.0D0)
	end function

	real(DP) function small_pressure(x)
	implicit none
	real(DP) :: x
	small_pressure = 1.6D0*x**5 - (4.0D0/7.0D0)*x**7 + (1.0D0/3.0D0)*x**9 - (5.0D0/2.2D1)*x**11 & 
			+ (3.5D1/2.08D2)*x**13 - (2.1D1/1.6D2)*x**15 + (2.31D2/2.176D3)*x**17
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the soundspeed of DM or NM !                              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDSOUNDSPEED
USE DEFINITION
USE NUCLEAR_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: j

! For DM !
IF (DM_flag == 1) THEN
	DO j=1, length_step_1
		cs1(j) = sqrt(dpdrho1(j)+dpdepsilon1(j)*p1(j)/rho1(j)**(2.0E0_DP))
	END DO
	CALL BOUNDARY1D_DM (cs1, even)
ENDIF

! For NM !
If(helmeos_flag == 1) THEN
	DO j=1, length_step_2
		CALL HELMEOS_CS(rho2(j), temp2(j), abar2(j), zbar2(j), ye2(j), cs2(j))
	END DO
ELSE
	DO j=1, length_step_2
		cs2(j) = sqrt(dpdrho2(j)+dpdepsilon2(j)*p2(j)/rho2(j)**(2.0E0_DP))
	END DO
END IF
CALL BOUNDARY1D_NM (cs2, even)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the pressure gradient !                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDDPDR
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: j

! Loop !
DO j = 0, length_step_2
	dpdr2 (j) = (p2L(j) - p2R(j-1))/dx2
END DO

END SUBROUTINE