!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This files contain all the riemann solvers available for !
! simulating hydrodynamics				   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Riemann_Module
USE DEFINITION
IMPLICIT NONE

! Left and right hydro-states for DM !
REAL (DP), ALLOCATABLE, DIMENSION(:) :: rho1R, rho1L, p1R, p1L, vel1R, vel1L, eps1R, eps1L, vf1R, vf1L

! Left and right hydro-states for NM !
REAL (DP), ALLOCATABLE, DIMENSION(:) :: rho2R, rho2L, p2R, p2L, vel2R, vel2L, eps2R, eps2L, vf2R, vf2L

! For finite temperature EOS !
REAL (DP), ALLOCATABLE, DIMENSION(:) :: temp2L, temp2R, rhoe2L, rhoe2R, abar2L, abar2R, zbar2L, zbar2R

! Speed of sound !
REAL (DP), ALLOCATABLE, DIMENSION(:) :: cs1L, cs1R
REAL (DP), ALLOCATABLE, DIMENSION(:) :: cs2L, cs2R

! Left and right scalars for DM !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: sca1, sca1L, sca1R

! Left and right scalars for NM !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: sca2, sca2L, sca2R

! Left and right fluxes, conserved quantity !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: fluxL1, fluxR1, uL1, uR1
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: fluxL2, fluxR2, uL2, uR2

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign left and right states and fluxes for riemann problem !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDRIEMANN
USE DEFINITION 
IMPLICIT NONE

! Left and right fluxes, conserved quantity !
IF(RUNDM_flag == 1) THEN
	ALLOCATE(rho1L(-4 : length_step_1 + 5))
	ALLOCATE(rho1R(-4 : length_step_1 + 5))
	ALLOCATE(p1L(-4 : length_step_1 + 5))
	ALLOCATE(p1R(-4 : length_step_1 + 5))
	ALLOCATE(vel1L(-4 : length_step_1 + 5))
	ALLOCATE(vel1R(-4 : length_step_1 + 5))
	ALLOCATE(eps1L(-4 : length_step_1 + 5))
	ALLOCATE(eps1R(-4 : length_step_1 + 5))
	ALLOCATE(cs1L(-4 : length_step_1 + 5))
	ALLOCATE(cs1R(-4 : length_step_1 + 5))
	ALLOCATE(vf1L(-4 : length_step_1 + 5))
	ALLOCATE(vf1R(-4 : length_step_1 + 5))
	ALLOCATE(fluxL1(imin1:imax1, -4 : length_step_1 + 5))
	ALLOCATE(fluxR1(imin1:imax1, -4 : length_step_1 + 5))
	ALLOCATE(uL1(imin1:imax1, -4 : length_step_1 + 5))
	ALLOCATE(uR1(imin1:imax1, -4 : length_step_1 + 5))
	ALLOCATE(sca1(iminsca1 : imaxsca1, -4 : length_step_1 + 5))
	ALLOCATE(sca1L(iminsca1 : imaxsca1, -4 : length_step_1 + 5))
	ALLOCATE(sca1R(iminsca1 : imaxsca1, -4 : length_step_1 + 5))
END IF

! NM !
ALLOCATE(rho2L(-4 : length_step_2 + 5))
ALLOCATE(rho2R(-4 : length_step_2 + 5))
ALLOCATE(p2L(-4 : length_step_2 + 5))
ALLOCATE(p2R(-4 : length_step_2 + 5))
ALLOCATE(vel2L(-4 : length_step_2 + 5))
ALLOCATE(vel2R(-4 : length_step_2 + 5))
ALLOCATE(eps2L(-4 : length_step_2 + 5))
ALLOCATE(eps2R(-4 : length_step_2 + 5))
ALLOCATE(vf2L(-4 : length_step_2 + 5))
ALLOCATE(vf2R(-4 : length_step_2 + 5))
ALLOCATE(cs2L(-4 : length_step_2 + 5))
ALLOCATE(cs2R(-4 : length_step_2 + 5))
ALLOCATE(fluxL2(imin2:imax2, -4 : length_step_2 + 5))
ALLOCATE(fluxR2(imin2:imax2, -4 : length_step_2 + 5))
ALLOCATE(uL2(imin2:imax2, -4 : length_step_2 + 5))
ALLOCATE(uR2(imin2:imax2, -4 : length_step_2 + 5))
ALLOCATE(sca2(iminsca2 : imaxsca2, -4 : length_step_2 + 5))
ALLOCATE(sca2L(iminsca2 : imaxsca2, -4 : length_step_2 + 5))
ALLOCATE(sca2R(iminsca2 : imaxsca2, -4 : length_step_2 + 5))

! Finite temperature EOS !
IF(helmeos_flag == 1) THEN
	ALLOCATE(temp2L(-4 : length_step_2 + 5))
	ALLOCATE(temp2R(-4 : length_step_2 + 5))
	ALLOCATE(rhoe2L(-4 : length_step_2 + 5))
	ALLOCATE(rhoe2R(-4 : length_step_2 + 5))
	ALLOCATE(abar2L(-4 : length_step_2 + 5))
	ALLOCATE(abar2R(-4 : length_step_2 + 5))
	ALLOCATE(zbar2L(-4 : length_step_2 + 5))
	ALLOCATE(zbar2R(-4 : length_step_2 + 5))
END IF

END 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign left and right states and fluxes for riemann problem !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DESTROYRIEMANN
USE DEFINITION 
IMPLICIT NONE

! Left and right fluxes, conserved quantity !
IF(RUNDM_flag == 1) THEN
	DEALLOCATE(rho1L)
	DEALLOCATE(rho1R)
	DEALLOCATE(p1L)
	DEALLOCATE(p1R)
	DEALLOCATE(vel1L)
	DEALLOCATE(vel1R)
	DEALLOCATE(eps1L)
	DEALLOCATE(eps1R)
	DEALLOCATE(cs1L)
	DEALLOCATE(cs1R)
	DEALLOCATE(fluxL1)
	DEALLOCATE(fluxR1)
	DEALLOCATE(uL1)
	DEALLOCATE(uR1)
	DEALLOCATE(sca1)
	DEALLOCATE(sca1L)
	DEALLOCATE(sca1R)
	DEALLOCATE(vf1L)
	DEALLOCATE(vf1R)
END IF

! NM !
DEALLOCATE(rho2L)
DEALLOCATE(rho2R)
DEALLOCATE(p2L)
DEALLOCATE(p2R)
DEALLOCATE(vel2L)
DEALLOCATE(vel2R)
DEALLOCATE(eps2L)
DEALLOCATE(eps2R)
DEALLOCATE(vf2L)
DEALLOCATE(vf2R)
DEALLOCATE(cs2L)
DEALLOCATE(cs2R)
DEALLOCATE(fluxL2)
DEALLOCATE(fluxR2)
DEALLOCATE(uL2)
DEALLOCATE(uR2)
DEALLOCATE(sca2)
DEALLOCATE(sca2L)
DEALLOCATE(sca2R)

! Finite temperature EOS !
IF(helmeos_flag == 1) THEN
	DEALLOCATE(temp2L)
	DEALLOCATE(temp2R)
	DEALLOCATE(rhoe2L)
	DEALLOCATE(rhoe2R)
	DEALLOCATE(abar2L)
	DEALLOCATE(abar2R)
	DEALLOCATE(zbar2L)
	DEALLOCATE(zbar2R)
END IF

END 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LF_DM (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j

! output !
REAL (DP), INTENT (OUT), DIMENSION (imin1 : imax1, -4 : length_step_1 + 5) :: flux_out

! the alpha in the LF flux !
REAL (DP) :: alpha(1 : no_of_eq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find maximum eigenvalues !
CALL ALPHASPLIT (alpha)

! For DM !
DO i = imin1, imax1
	DO j = -1, length_step_part_1 + 1
		flux_out(i, j) = 0.5D0 * (fluxL1 (i, j) + fluxR1 (i, j) - alpha(i) * (uR1 (i, j) - uL1 (i,j)))
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The LF flux, see for example. shu et al. for the original WENO paper !					    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LF_NM (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j

! output !
REAL (DP), INTENT (OUT), DIMENSION (imin2 : imax2, -4 : length_step_2 + 5) :: flux_out

! the alpha in the LF flux !
REAL (DP) :: alpha(1 : no_of_eq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find maximum eigenvalues !
CALL ALPHASPLIT (alpha)

! For NM !
DO i = imin2, imax2
	DO j = -1, length_step_part_2 + 1
		flux_out(i, j) = 0.5D0 * (fluxL2 (i, j) + fluxR2 (i, j) - alpha(i) * (uR2 (i, j) - uL2 (i,j)))
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The HLL riemann solver, see Toro. 2009 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HLL_DM (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! output !
REAL (DP), INTENT (OUT), DIMENSION (imin1 : imax1, -4 : length_step_1 + 5) :: flux_out

! some local array !
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: sL, sR

! For more general signal speed !
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: ubar

! For more general signal speed !
REAL (DP), DIMENSION (-4 : length_step_1 + 5) :: dbar, dbar2, neta2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
DO j = -1, length_step_part_1 + 1
	! Get roe averaged quantity !
	ubar(j) = compute_roe(vel1L(j),vel1R(j), rho1L(j), rho1R(j))

	! Get more general signal speed parameter !
	neta2(j) = 0.5D0*sqrt(rho1R(j)*rho1L(j))/(sqrt(rho1L(j)) + sqrt(rho1R(j)))**2
	dbar2(j) = compute_roe(cs1L(j)**2,cs1R(j)**2,rho1L(j),rho1R(j)) + neta2(j)*(vel1R(j) - vel1L(j))**2
	dbar(j) = sqrt(dbar2(j))
END DO

! Get signal speed !
DO j = -1, length_step_part_1 + 1
	sL(j) = min(vel1L(j) - cs1L(j), ubar(j) - dbar(j))
	sR(j) = max(vel1R(j) + cs1R(j), ubar(j) + dbar(j))
END DO

! Find the flux !
DO j = -1, length_step_part_1 + 1
	IF(sL(j) >= 0.0D0) THEN
		DO i = imin1, imax1
			flux_out(i,j) = fluxL1(i,j)
		END DO
	ELSEIF(sL(j) < 0.0D0 .AND. sR(j) > 0.0D0) THEN
		DO i = imin1, imax1
			flux_out(i,j) = compute_fluxhll(fluxL1(i,j),fluxR1(i,j),uL1(i,j),uR1(i,j),sL(j),sR(j))
		END DO
	ELSEIF(sR(j) <= 0.0D0) THEN
		DO i = imin1, imax1
			flux_out(i,j) = fluxR1(i,j)
		END DO
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	real(DP) function compute_fluxhll(yl,yr,xl,xr,sl,sr)
	implicit none
	real(DP) :: yl, yr, xl, xr, sl, sr
	compute_fluxhll = (sr*yl-sl*yr+sl*sr*(xr-xl))/(sr-sl)
	end function

	real(DP) function compute_roe(xl,xr,rhol,rhor)
	implicit none
	real(DP) :: xl,xr,rhol,rhor
	compute_roe = (sqrt(rhol)*xl + sqrt(rhor)*xr)/(sqrt(rhol) + sqrt(rhor))
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The HLL riemann solver, see Toro. 2009 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HLL_NM (flux_out)
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! output !
REAL (DP), INTENT (OUT), DIMENSION (imin2 : imax2, -4 : length_step_2 + 5) :: flux_out

! some local array !
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: sL, sR

! For more general signal speed !
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: ubar

! For more general signal speed !
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: dbar, dbar2, neta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For NM !
DO j = -1, length_step_part_2 + 1
	! Get roe averaged quantity !
	ubar(j) = compute_roe(vel2L(j),vel2R(j), rho2L(j), rho2R(j))

	! Get more general signal speed parameter !
	neta(j) = 0.5D0*sqrt(rho2R(j)*rho2L(j))/(sqrt(rho2L(j)) + sqrt(rho2R(j)))**2
	dbar2(j) = compute_roe(cs2L(j)**2,cs2R(j)**2,rho2L(j),rho2R(j)) + neta(j)*(vel2R(j) - vel2L(j))**2
	dbar(j) = sqrt(dbar2(j))
END DO
DO j = -1, length_step_part_2 + 1
	sL(j) = min(vel2L(j) - cs2L(j), ubar(j) - dbar(j))
	sR(j) = max(vel2R(j) + cs2R(j), ubar(j) + dbar(j))
END DO

! Find the flux !
DO j = -1, length_step_part_2 + 1
	IF(sL(j) >= 0.0D0) THEN
		DO i = imin2, imax2
			flux_out(i,j) = fluxL2(i,j)
		END DO
	ELSEIF(sL(j) < 0.0D0 .AND. sR(j) > 0.0D0) THEN
		DO i = imin2, imax2
			flux_out(i,j) = compute_fluxhll(fluxL2(i,j),fluxR2(i,j),uL2(i,j),uR2(i,j),sL(j),sR(j))
		END DO
	ELSEIF(sR(j) <= 0.0D0) THEN
		DO i = imin2, imax2
			flux_out(i,j) = fluxR2(i,j)
		END DO
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	real(DP) function compute_fluxhll(yl,yr,xl,xr,sl,sr)
	implicit none
	real(DP) :: yl, yr, xl, xr, sl, sr
	compute_fluxhll = (sr*yl-sl*yr+sl*sr*(xr-xl))/(sr-sl)
	end function

	real(DP) function compute_roe(xl,xr,rhol,rhor)
	implicit none
	real(DP) :: xl,xr,rhol,rhor
	compute_roe = (sqrt(rhol)*xl + sqrt(rhor)*xr)/(sqrt(rhol) + sqrt(rhor))
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE