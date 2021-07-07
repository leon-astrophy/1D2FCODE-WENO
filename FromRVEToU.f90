!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine convert the primitive variables (like pressure, density...) !    
! to conservative variables in the system of hyperbolic equation              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FROMRVETOU (u1, u2)
USE DEFINITION
USE FLAME_MODULE
USE NUCLEAR_MODULE
IMPLICIT NONE

! The real parameter u in the hypperbolic equation !
REAL (DP), INTENT (OUT), DIMENSION (imin1 : imax1, -4 : length_step_1 + 5) :: u1
REAL (DP), INTENT (OUT), DIMENSION (imin2 : imax2, -4 : length_step_2 + 5) :: u2

! Integer parameter !
INTEGER :: j

! Converting natural variables to conservative variables !
! We include the hydro equation for DM only if DM precense or DM is movable !
! Mutted energy equation for cold DM-EOS !
If (RUNDM_flag == 1) then
	DO j = -4, length_step_1 + 5
		u1 (irho1, j) = rho1 (j)
		u1 (ivel1, j) = rho1 (j) * vel1 (j)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! No need to solve energy equation for cold DM equation of states !
		!u (itau1, j) = rho1 (j) * epsilon1 (j) + 5.0E-1_DP * rho1 (j) * vel1 (j) ** 2
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	END DO
END IF

! we do the same thing for NM !
DO j = -4, length_step_2 + 5
	u2 (irho2, j) = rho2 (j) 
	u2 (ivel2, j) = rho2 (j) * vel2(j)
END DO

! For Energy equation !
IF(nm_epsilon == 1) THEN 
	DO j = -4, length_step_2 + 5
		u2 (itau2, j) = rho2 (j) * epsilon2 (j) + 5.0E-1_DP * rho2 (j) * vel2 (j) ** 2
	END DO
END IF
IF(dual_energy == 1) THEN 
	DO j = -4, length_step_2 + 5
		u2 (ieps2, j) = rho2 (j) * epsilon2 (j)
	END DO
END IF

! There is an extra equation if we allows deflagration !
IF(flame_flag == 1) THEN
	DO j = -4, length_step_2 + 5
		u2 (iscaG,j) = scaG(j) * rho2(j)
	END DO
END IF

! detonation level set !
IF(deton_flag == 1) THEN
	DO j = -4, length_step_2 + 5
		u2 (iscaG2,j) = scaG2(j) * rho2(j)
	END DO
END IF

! There are extra equations if we allow electron fractions to be advected !
IF(etran_flag == 1) THEN
	DO j = -4, length_step_2 + 5 
		u2 (iye2,j) = ye2 (j) * rho2 (j)
	END DO
END IF

! There are extra equations if we allow isotopes to be advected !
IF (xisotran_flag == 1) THEN
	DO j = -4, length_step_2 + 5 
		u2 (ihe4,j) = xiso(che4,j) * rho2(j)
		u2 (ic12,j) = xiso(cc12,j) * rho2(j)
		u2 (io16,j) = xiso(co16,j) * rho2(j)
		u2 (ine20,j) = xiso(cne20,j) * rho2(j)
		u2 (img24,j) = xiso(cmg24,j) * rho2(j)
		u2 (isi28,j) = xiso(csi28,j) * rho2(j)
		u2 (ini56,j) = xiso(cni56,j) * rho2(j)
	END DO
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine convert the conservative variables to primitive variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FROMUTORVE (u1, u2)
USE DEFINITION
USE flame_module
use nuclear_module
IMPLICIT NONE

! The numerical u in the hyperbolic PDE !
REAL (DP), INTENT (IN), DIMENSION (imin1 : imax1, -4 : length_step_1 + 5) :: u1
REAL (DP), INTENT (IN), DIMENSION (imin2 : imax2, -4 : length_step_2 + 5) :: u2

! Dummy variables !
REAL (DP) :: dummy_max

! Integer !
INTEGER :: j

! We assign the values to natural variables !
! We include the hydro equation for DM only if DM precense and DM is movable !
! Mutted energy equation for cold DM EOS !
If (RUNDM_flag == 1) then
	DO j = -4, length_step_1 + 5
		rho1 (j) = u1 (irho1, j)
		vel1 (j) = u1 (ivel1, j) / rho1 (j)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! No energy equation for cold DM-EOS !
		!epsilon1 (j) = u1 (itau1, j) / rho1 (j) - 5.0E-1_DP * vel1 (j) ** 2
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	END DO
END IF
	
! We do the same for NM !
DO j = -4, length_step_2 + 5
	rho2 (j) = u2 (irho2, j) 
	vel2 (j) = u2 (ivel2, j) / rho2 (j) 
END DO

! For NM energy equation !
IF(dual_energy == 1) THEN
	DO j = -4, length_step_2 + 5
		rhoe2 (j) = u2 (ieps2, j)
		bige2 (j) = u2 (itau2, j) / rho2 (j)
		et2 (j) = u2 (itau2, j) / rho2 (j)  - 5.0E-1_DP * vel2 (j) ** 2
	END DO

	! Determine the epsilon for pressure equation !
	DO j = -4, length_step_2 + 5
		If(et2 (j) > eta_1*bige2 (j)) THEN
			epsp2(j) = et2(j)
		ELSE
			epsp2(j) = u2 (ieps2, j) / rho2 (j)
		END IF
	END DO

	! Determine the epsilon for epsilon equation !
	DO j = -3, length_step_2 + 4
		!dummy_max = max(rho2(j-1)*bige2 (j-1),rho2(j)*bige2 (j),rho2(j+1)*bige2 (j+1))
		If(et2 (j) > 1.0D-4*bige2 (j)) THEN
			epsilon2(j) = et2 (j)
		ELSE
			epsilon2(j) = u2 (ieps2, j) / rho2 (j)
		END IF
	END DO

	! Special setting for boundary values !
	CALL BOUNDARY1D_NM (epsilon2, even)
ELSE
	IF(nm_epsilon == 1) THEN
		DO j = -4, length_step_2 + 5
			epsilon2 (j) = u2 (itau2, j) / rho2 (j) - 5.0E-1_DP * vel2 (j) ** 2
			epsp2 (j) = epsilon2 (j)
		END DO
	END IF
END IF

! There is an extra equation if we allows deflagration !
IF(flame_flag == 1) THEN
	DO j = -4, length_step_2 + 5
		scaG(j) = u2 (iscaG,j) / rho2(j)
	END DO
END IF

! Detonation level-set !
!IF(deton_flag == 1) THEN
!	DO j = -4, length_step_2 + 5
!		scaG2(j) = u2 (iscaG2,j) / rho2(j)
!	END DO
!END IF

! Convert the electron fraction !
IF(etran_flag == 1) THEN
	DO j = -4, length_step_2 + 5
		ye2 (j) = u2 (iye2,j) / rho2 (j) 

		! Ensure the electron fraction lies in suitable range !
		ye2 (j) = MIN(MAX(ye2 (j), ye_min), ye_max)
	END DO
ENDIF

! Convert the isotope composition !
IF (xisotran_flag == 1) THEN
	DO j = -4, length_step_2 + 5
		xiso (che4, j) = u2 (ihe4, j) / rho2(j)
		xiso (cc12, j) = u2 (ic12, j) / rho2(j)
		xiso (co16, j) = u2 (io16, j) / rho2(j)
		xiso (cne20, j) = u2 (ine20, j) / rho2(j)
		xiso (cmg24, j) = u2 (img24, j) / rho2(j)
		xiso (csi28, j) = u2 (isi28, j) / rho2(j)
		xiso (cni56, j) = u2 (ini56, j) / rho2(j)
	END DO

	! since we do a conversion, we need to check whether abundance sum to 1 !
	CALL checkXisotope
END IF
	
END SUBROUTINE