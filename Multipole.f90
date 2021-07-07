!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find Gravitational Potential Using Multipole Method !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MULTIPOLE_DM
USE DEFINITION
IMPLICIT NONE

! Potential at left and right faces !
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: phi1L, phi1R

! Integer ! 
integer :: j, k

! temporaily arrays !
REAL (DP) :: tmp1, tmp2, sum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the potential at left 1/4 from cell center !
DO j = 1, length_step_1
tmp1 = 0.0D0
tmp2 = 0.0D0
	DO k = 1, length_step_1
		! Contribution to the intter multipole moments !
		If(r1(k) < r1L(j)) THEN
			If(rho1(k) > rho1_a) THEN
				tmp1 = tmp1 + vol1(k) * rho1(k)
			ELSE
				tmp1 = tmp1
			END IF
		
		! For outer multipole moments !
		ELSEIF(r1(k) > r1L(j)) THEN
			If(rho1(k) > rho1_a) THEN
				tmp2 = tmp2 + vol1(k) * rho1(k) / r1(k)
			ELSE
				tmp2 = tmp2
			END IF
		END IF
	
	END DO

	! Assign potential !
	phi1L(j) =  - (tmp1/r1L(j) + tmp2)

END DO

! Find the potential at right 1/4 from cell center !
DO j = 1, length_step_1
tmp1 = 0.0D0
tmp2 = 0.0D0
	DO k = 1, length_step_1
	
		! Contribution to the intter multipole moments !
		If(r1(k) < r1R(j)) THEN
			If(rho1(k) > rho1_a) THEN
				tmp1 = tmp1 + vol1(k) * rho1(k)
			ELSE
				tmp1 = tmp1
			END IF
		
		! For outer multipole moments !
		ELSEIF(r1(k) > r1R(j)) THEN
			If(rho1(k) > rho1_a) THEN
				tmp2 = tmp2 + vol1(k) * rho1(k) / r1(k)
			ELSE
				tmp2 = tmp2
			END IF
		END IF
	
	END DO

	! Assign potential !
	phi1R(j) =  - (tmp1/r1R(j) + tmp2)

END DO

! Do the averaging !
DO j = 1, length_step_1
	phidm_1 (j) = 0.5D0*(phi1R(j) + phi1L(j))
END DO

! Copy to boundary !
CALL BOUNDARY1D_DM (phidm_1,even)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find Gravitational Potential Using Multipole Method !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MULTIPOLE_NM
USE DEFINITION
IMPLICIT NONE

! Potential at left and right faces !
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: phi2L, phi2R

! Integer ! 
integer :: j, k

! temporaily arrays !
REAL (DP) :: tmp1, tmp2, sum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the potential at left 1/4 from cell center !
DO j = 1, length_step_2
tmp1 = 0.0D0
tmp2 = 0.0D0
	DO k = 1, length_step_2
		! Contribution to the intter multipole moments !
		If(r2(k) < r2L(j)) THEN
			If(rho2(k) > rho2_a) THEN
				tmp1 = tmp1 + vol2(k) * rho2(k)
			ELSE
				tmp1 = tmp1
			END IF
		
		! For outer multipole moments !
		ELSEIF(r2(k) > r2L(j)) THEN
			If(rho2(k) > rho2_a) THEN
				tmp2 = tmp2 + vol2(k) * rho2(k) / r2(k)
			ELSE
				tmp2 = tmp2
			END IF
		END IF
	
	END DO

	! Assign potential !
	phi2L(j) =  - (tmp1/r2L(j) + tmp2)

END DO

! Find the potential at right 1/4 from cell center !
DO j = 1, length_step_2
tmp1 = 0.0D0
tmp2 = 0.0D0
	DO k = 1, length_step_2
	
		! Contribution to the intter multipole moments !
		If(r2(k) < r2R(j)) THEN
			If(rho2(k) > rho2_a) THEN
				tmp1 = tmp1 + vol2(k) * rho2(k)
			ELSE
				tmp1 = tmp1
			END IF
		
		! For outer multipole moments !
		ELSEIF(r2(k) > r2R(j)) THEN
			If(rho2(k) > rho2_a) THEN
				tmp2 = tmp2 + vol2(k) * rho2(k) / r2(k)
			ELSE
				tmp2 = tmp2
			END IF
		END IF
	
	END DO

	! Assign potential !
	phi2R(j) =  - (tmp1/r2R(j) + tmp2)

END DO

! Do the averaging !
DO j = 1, length_step_2
	phinm_2 (j) = 0.5D0*(phi2R(j) + phi2L(j))
END DO

! Copy to boundary !
CALL BOUNDARY1D_NM (phinm_2,even)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Map the potential between grids !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INTERPOLATION
USE DEFINITION 
IMPLICIT NONE

! Integer !
Integer :: j, k

! Find the mass !
CALL FINDMASS

! Map only when DM is precense !
IF (DM_flag == 1) THEN

	! Map from DM grid to NM grid !
	DO j = 1, length_step_2
		IF(r2(j) > r1(length_step_1)) THEN
			phinm_1(j) = - mass1/r2(j)
		ELSE
			DO k = 1, length_step_1
				IF(r1(k) == r2(j)) THEN
					phinm_1(j) = phidm_1(k)
					EXIT
				ELSEIF(r1(k) > r2(j)) THEN
					CALL AKIMA(r1(k-3), r1(k-2), r1(k-1), r1(k), r1(k+1), r1(k+2), & 
					phidm_1(k-3), phidm_1(k-2), phidm_1(k-1), phidm_1(k), phidm_1(k+1), phidm_1(k+2), r2(j), phinm_1(j))
					EXIT
				END IF
			END DO
		END IF
	END DO

	! Map from NM grid to DM grid !
	DO j = 1, length_step_1
		IF(r1(j) > r2(length_step_2)) THEN
			phidm_2(j) = - mass2/r1(j)
		ELSE
			DO k = 1, length_step_2
				IF(r2(k) == r1(j)) THEN
					phidm_2(j) = phinm_2(k)
					EXIT
				ELSEIF(r2(k) > r1(j)) THEN
					CALL AKIMA(r2(k-3), r2(k-2), r2(k-1), r2(k), r2(k+1), r2(k+2), & 
					phinm_2(k-3), phinm_2(k-2), phinm_2(k-1), phinm_2(k), phinm_2(k+1), phinm_2(k+2), r1(j), phidm_2(j))
					EXIT
				END IF
			END DO
		END IF
	END DO

	! Sum the total potential !
	DO j = 1, length_step_1
		phidm (j) = phidm_1 (j) + phidm_2 (j)
	END DO
	DO j = 1, length_step_2
		phinm (j) = phinm_1 (j) + phinm_2 (j)	
	END DO

ELSE

	! No interpolation !
	phinm = phinm_2
	
END IF

! Map to boundary !
CALL BOUNDARY1D_NM (phinm, even)

! Compute the gravitational force using 4th order center difference !
DO j = 1, length_step_2
	phip_nm (j) = (- phinm (j + 2) + 8.0E0_DP * phinm (j + 1) - 8.0E0_DP * phinm (j - 1) + phinm (j - 2)) / (1.2E1_DP * dx2)
END DO
CALL BOUNDARY1D_NM (phip_nm, odd)

! NM Gravitational force at cell interface !
DO j = 0, length_step_2
	phipface_2 (j) = (phinm (j + 1) - phinm (j)) / dx2
END DO

! For DM !
IF(DM_flag == 1) THEN
	CALL BOUNDARY1D_DM (phidm, even)
	! Compute the gravitational force using 4th order center difference !
	DO j = 1, length_step_1
		phip_dm (j) = (- phidm (j + 2) + 8.0E0_DP * phidm (j + 1) - 8.0E0_DP * phidm (j - 1) + phidm (j - 2)) / (1.2E1_DP * dx1)
	END DO
	CALL BOUNDARY1D_DM (phip_dm, odd)

	! DM Gravitational force at cell interface !
	DO j = 0, length_step_1
		phipface_1 (j) = (phidm (j + 1) - phidm (j)) / dx1
	END DO
END IF

END SUBROUTINE
 
!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate GR potential !
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GRPOTENTIAL
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: j, k

! Integrated tov mass !
REAL (DP), DIMENSION (-4:length_step_2 + 5) :: mtov

! Cell-interace tov mass !
REAL (DP), DIMENSION (-4:length_step_2 + 5) :: tovface

! Lorentz factor !
REAL (DP), DIMENSION (-4:length_step_2 + 5) :: lorentz

! Integrand to be integrated !
REAL (DP), DIMENSION (-4:length_step_2 + 5) :: integrand

! Initialize !
integrand = 0.0D0
lorentz = 0.0D0
tovface = 0.0D0
phigr = 0.0D0
mtov = 0.0D0

! Find TOV mass lying at cell-interace !
DO j = 1, length_step_2
	DO k = 1, j
		IF(rho2(k) > rho2_a) THEN
			tovface(j) = tovface(j) + vol2(k)*rho2(k)*(1.0D0 + epsilon2(k))
		ELSE
			CONTINUE
		END IF
	END DO
END DO

DO j = 1, length_step_2

	! Interpolate them to cell center !
	mtov(j) = 0.5D0*(tovface(j-1) + tovface(j))

	! Get lorentz factor !
	IF(rho2(j) > rho2_a) THEN
		lorentz(j) = SQRT(1.0D0 + vel2(j)**2 - 2.0D0*mtov(j)/r2(j))
	ELSE
		lorentz(j) = SQRT(1.0D0 - 2.0D0*mtov(j)/r2(j))
	END IF

	! Form the integrand !
	IF(rho2(j) > rho2_a) THEN
		integrand(j) = (mtov(j)/4.0D0/pi_old + r2(j)**3*p2(j))*(rho2(j) + rho2(j)*epsilon2(j) + p2(j))/rho2(j)/lorentz(j)**2/r2(j)**2
	ELSE
		integrand(j) = mtov(j)/4.0D0/pi_old/lorentz(j)**2/r2(j)**2
	END IF

END DO

! Get the GR potential !
DO j = 1, length_step_2
	phigr(j) = 0.5D0*(integrand(j) + integrand(length_step_2))
 	DO k = j + 1, length_step_2 - 1
		phigr(j) = phigr(j) + integrand(k)
	END DO
	phigr(j) = phigr(j)*dx2*(-4.0D0*pi_old) + 0.5D0*log(1.0D0 - 2.0D0*mtov(length_step_2)/r2(length_step_2))
END DO

! Copy to boundary !
CALL BOUNDARY1D_NM (phigr,even)

END SUBROUTINE