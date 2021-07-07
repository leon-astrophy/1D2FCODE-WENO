!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the specific energy density epsilon      !
! which is necessary for hydro dynamics evolution. We assume     !
! completely degenerate fermi gas EOS for DM and either          !
! finite temperature EOS or completely degenerate EOS for NM     !
! If you want to use your own EOS, you need to take care of this !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETEPSILON
USE DEFINITION
USE NUCLEAR_MODULE
IMPLICIT NONE

! Integer variable !
INTEGER :: j

! Dummy variable !
REAL (DP) :: dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the DM case first !
! This part is done only if the users wants DM component
IF (DM_flag == 1) THEN
	IF (fermieosdm_flag == 1) THEN
		DO j = -4, length_step_1 + 5

			! We need to get the dimensionless fermi momentum for further calculation !
			CALL FERMIMO (dlfmmo1, rho1 (j), 1)        
			IF (dlfmmo1<=1.0E-2_DP) THEN
				epsilon1 (j) = a_max1*small_energy(dlfmmo1)/rho1(j)
			ELSE
				epsilon1 (j) = a_max1*large_energy(dlfmmo1)/rho1(j)
			END IF
		END DO
	ELSE
		DO j = -4, length_step_1 + 5
			epsilon1 (j) = k1 * rho1(j) ** (gamma1 - 1.0E0_DP) / (gamma1 - 1.0E0_DP)
		END DO
	END IF

	! We assign the atmospheric specific energy density !
	epsilon1_a = epsilon1(length_step_1)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the NM case !
! We see whether finite temperature EOS is used !
IF(helmeos_flag == 1) THEN
	DO j = -4, length_step_2 + 5
		CALL HELMEOS_RtoE(rho2 (j), temp2(j), abar2(j), zbar2(j), ye2(j), epsilon2 (j), dummy)
	END DO
ELSEIF (fermieosnm_flag == 1) THEN
	DO j = -4, length_step_2 + 5

		! We need to get the dimensionless fermi momentum for further calculation !
		CALL FERMIMO (dlfmmo2, rho2 (j), 2)        
		IF (dlfmmo2<=1.0E-2_DP) THEN
			epsilon2 (j) = a_max2*small_energy(dlfmmo2)/rho2(j)
		ELSE
			epsilon2 (j) = a_max2*large_energy(dlfmmo2)/rho2(j)
		END IF

	END DO
ELSEIF (ccsneos_flag == 1) THEN
	DO j = -4, length_step_2 + 5
		epsilon2 (j) = kc_1 * rho2(j) ** (gammac_1 - 1.0E0_DP) / (gammac_1 - 1.0E0_DP)
	END DO
ELSE
	! The following steps are more or less the same !
	DO j = -4, length_step_2 + 5
		epsilon2 (j) = k2 * rho2(j) ** (gamma2 - 1.0E0_DP) / (gamma2 - 1.0E0_DP)
	END DO

END IF

! We assign the atmospheric specific energy density !
epsilon2_a = epsilon2(length_step_2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
	real(DP) function large_energy(x)
	implicit none
	real(DP) :: x
	large_energy = 3.0D0*x*SQRT(x**2 + 1.0D0)*(1.0D0 + 2.0D0*x**2) - 3.0D0*log(x + SQRT(x**2 + 1.0D0)) - 8.0D0*x**3
	end function

	real(DP) function small_energy(x)
	implicit none
	real(DP) :: x
	small_energy = 8.0D0*x**3 + (1.2D1/5.0D0)*x**5 - (3.0D0/7.0D0)*x**7 + (1.0D0/6.0D0)*x**9 - (1.5D1/1.76D2)*x**11 & 
			+ (2.1D1/4.16D2)*x**13 - (2.1D1/6.40D2)*x**15 + (9.9D1/4.352D3)*x**17 - 8.0D0*x**3
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the internal energy for cold EOS !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDEPS
USE DEFINITION
IMPLICIT NONE

! Integer variable !
INTEGER :: j

! Dummy variable !
REAL (DP) :: dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the DM case first !
IF (DM_flag == 1) THEN
	IF (fermieosdm_flag == 1) THEN
		DO j = -4, length_step_1 + 5

			! We need to get the dimensionless fermi momentum for further calculation !
			CALL FERMIMO (dlfmmo1, rho1 (j), 1)        
			IF (dlfmmo1<=1.0E-2_DP) THEN
				epsilon1 (j) = a_max1*small_energy(dlfmmo1)/rho1(j)
			ELSE
				epsilon1 (j) = a_max1*large_energy(dlfmmo1)/rho1(j)
			END IF

		END DO
	ELSE
		DO j = -4, length_step_1 + 5
			epsilon1 (j) = k1 * rho1(j) ** (gamma1 - 1.0E0_DP) / (gamma1 - 1.0E0_DP)
		END DO
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the NM case !
IF (NM_epsilon == 0) THEN
	IF (fermieosnm_flag == 1) THEN
		DO j = -4, length_step_2 + 5

			! We need to get the dimensionless fermi momentum for further calculation !
			CALL FERMIMO (dlfmmo2, rho2 (j), 2)        
			IF (dlfmmo2<=1.0E-2_DP) THEN
				epsilon2 (j) = a_max2*small_energy(dlfmmo2)/rho2(j)
			ELSE
				epsilon2 (j) = a_max2*large_energy(dlfmmo2)/rho2(j)
			END IF

		END DO
	ELSE
		! The following steps are more or less the same !
		DO j = -4, length_step_2 + 5
			epsilon2 (j) = k2 * rho2(j) ** (gamma2 - 1.0E0_DP) / (gamma2 - 1.0E0_DP)
		END DO
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
	real(DP) function large_energy(x)
	implicit none
	real(DP) :: x
	large_energy = 3.0D0*x*SQRT(x**2 + 1.0D0)*(1.0D0 + 2.0D0*x**2) - 3.0D0*log(x + SQRT(x**2 + 1.0D0)) - 8.0D0*x**3
	end function

	real(DP) function small_energy(x)
	implicit none
	real(DP) :: x
	small_energy = 8.0D0*x**3 + (1.2D1/5.0D0)*x**5 - (3.0D0/7.0D0)*x**7 + (1.0D0/6.0D0)*x**9 - (1.5D1/1.76D2)*x**11 & 
			+ (2.1D1/4.16D2)*x**13 - (2.1D1/6.40D2)*x**15 + (9.9D1/4.352D3)*x**17 - 8.0D0*x**3
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update the rhoe for dual energy formalism !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDRHOE
USE DEFINITION 
IMPLICIT NONE

! Integer !
Integer :: j

DO j = -4, length_step_2 + 5
	rhoe2(j) = rho2(j)*epsilon2(j)
END DO

CALL BOUNDARY1D_NM (rhoe2,even)

END SUBROUTINE