!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine take the conservative variables as input   !
! and calculate the flux function for each conservative      !
! equation. The flux function is being further split using   !
! global LF flux splitting scheme. They are input into 	     !
! WENO reconstruction scheme to calculate the numerical flux !
! at both left and right hand side of at each boundary cell  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SPATIAL (u1, u2)
USE DEFINITION
USE RIEMANN_MODULE
USE MP5_MODULE
USE PPM_MODULE
USE WENO_MODULE
USE FLAME_MODULE
USE NUCLEAR_MODULE
IMPLICIT NONE

! integer parameter !
INTEGER :: i, j

! real parameter !
REAL (DP) :: rho_min_dm, rho_min_nm

! the flux, geometric and other source terms !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: sa1, sb1, sa2, sb2

! The flux and flux differences !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: flux1, dfdx1, flux2, dfdx2

! the conservative variables ! 
REAL (DP), INTENT (IN), DIMENSION (imin1 : imax1, -4 : length_step_1 + 5) :: u1
REAL (DP), INTENT (IN), DIMENSION (imin2 : imax2, -4 : length_step_2 + 5) :: u2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We need to allocate the neccesarry variables !
ALLOCATE(sa1(imin1 : imax1, -4 : length_step_1 + 5))
ALLOCATE(sb1(imin1 : imax1, -4 : length_step_1 + 5))
ALLOCATE(flux1(imin1 : imax1, -4 : length_step_1 + 5))
ALLOCATE(dfdx1(imin1 : imax1, -4 : length_step_1 + 5))

ALLOCATE(sa2(imin2 : imax2, -4 : length_step_2 + 5))
ALLOCATE(sb2(imin2 : imax2, -4 : length_step_2 + 5))
ALLOCATE(flux2(imin2 : imax2, -4 : length_step_2 + 5))
ALLOCATE(dfdx2(imin2 : imax2, -4 : length_step_2 + 5))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the frame velocity !
IF(movinggriddm_flag == 1) THEN
	CALL FINDFRAMEVEL_DM
END IF
IF(movinggridnm_flag == 1) THEN
	CALL FINDFRAMEVEL_NM
END IF

! For riemann solver !
CALL BUILDRIEMANN

! Initialize source term !
IF(RUNDM_flag == 1) THEN
	flux1(:,:) = 0.0D0
	dfdx1(:,:) = 0.0D0
	sa1(:,:) = 0.0D0
	sb1(:,:) = 0.0D0
END IF

! Initialize fluxes !
flux2(:,:) = 0.0D0
dfdx2(:,:) = 0.0D0
sa2(:,:) = 0.0D0
sb2(:,:) = 0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign scalars !
IF(RUNDM_flag == 1) THEN
	IF(iminsca1 > 0) THEN
		DO i = iminsca1, imaxsca1
			DO j = -4, length_step_part_1 + 5
				sca1(i,j) = u1(i,j)/u1(irho1,j)
			END DO
		END DO
	END IF
END IF
IF(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO j = -4, length_step_part_2 + 5
			sca2(i,j) = u2(i,j)/u2(irho2,j)
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Then reconstruct primitive variables !
IF(MP5_flag == 1) THEN
	CALL RECONSTRUCT_MP5
ELSEIF(PPM_flag == 1) THEN
	CALL RECONSTRUCT_PPM
ELSEIF(WENO_flag == 1) THEN
	CALL RECONSTRUCT_WENO
ELSEIF(TVD_flag == 1) THEN
	CALL RECONSTRUCT_TVD
END IF

! Then compute the conservative variables and flux function !
CALL BUILDSTATES

! Choose an appropriate riemann solver !
IF(RUNDM_flag == 1) THEN
	IF(HLL_flag == 1) THEN
		CALL HLL_DM(flux1(imin1:imax1,:))
	ELSEIF(LF_flag == 1) THEN
		CALL LF_DM(flux1(imin1:imax1,:))
	END IF
END IF

! For NM !
IF(HLL_flag == 1) THEN
	CALL HLL_NM(flux2(imin2:imax2,:))
ELSEIF(LF_flag == 1) THEN
	CALL LF_NM(flux2(imin2:imax2,:))
END IF

! Find pressure gradient !
CALL FINDDPDR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign DM source term !
If (RUNDM_flag == 1) then

	! Threshold density !
	rho_min_dm = 1.1E0_DP * rho1_a	

	! Geometric source for momentum equation !
	DO j = 1, length_step_part_1
		sa1 (ivel1, j) = - p1(j)
	END DO

	! Gravitational source terms !
	DO j = 1, length_step_part_1
	   	IF(rho1(j) > rho_min_dm) THEN
               		sb1 (ivel1, j) = rho1(j)*phip_dm (j)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! No DM internal energy !
               		!sb1 (itau1, j) = rho1(j)*vel1(j)*phip_dm (j)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            	ENDIF
	   	IF(rho1(j) <= rho_min_dm) THEN
               		sb1 (ivel1, j) = 0.0D0
               		sb1 (itau1, j) = 0.0D0
            	ENDIF	  	    
	ENDDO

	! Moving grid !
	IF(movinggriddm_flag == 1) THEN
		DO j = 1, length_step_part_1
			DO i = imin1, imax1
				sb1 (i,j) = sb1 (i,j) + u1 (i,j) * (3.0D0 * vel1_max / radius1)
			ENDDO
		ENDDO
	END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we do the same thing for NM !
! Set NM density threshold !
rho_min_nm = 1.1E0_DP * rho2_a	

! Geometric source for momentum equation !
DO j = 1, length_step_part_2
	sa2 (ivel2, j) = - p2(j)
END DO

! Geometric and gravitational source terms !
DO j = 1, length_step_part_2
	IF(rho2(j) > rho_min_nm) THEN
		sb2 (ivel2, j) = rho2(j)*phip_nm (j)
		If(nm_epsilon == 1) THEN
			sb2 (itau2, j) = rho2(j)*vel2(j)*phip_nm (j)
		END IF
	ENDIF	
	IF(rho2(j) <= rho_min_nm) THEN
		sb2 (ivel2, j) = 0.0D0
		sb2 (itau2, j) = 0.0D0
	END IF
ENDDO

! Dual energy !
IF(dual_energy == 1) THEN
	DO j = 1, length_step_part_2
		sb2 (ieps2, j) = - vel2(j)*dpdr2(j)
	END DO
END IF 

! Moving grid !
IF(movinggridnm_flag == 1) THEN
	DO j = 1, length_step_part_2
		DO i = imin2, imax2
			sb2 (i,j) = sb2 (i,j) + u2 (i,j) * (3.0D0 * vel2_max / radius2)
		ENDDO
	ENDDO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! now we calculate the flux difference at each ith grid dfdx !
IF (RUNDM_flag == 1) THEN

	! Initialize !
	l1(:,:) = 0.0D0

	DO i = imin1, imax1
		DO j = 1, length_step_part_1
	
			! we find the flux difference at the jth cell center !
			If(sp_dim_i /= 0) THEN
				dfdx1 (i, j) = (r1F(j)**DBLE(sp_dim)*flux1 (i, j) - r1F(j-1)**DBLE(sp_dim)*flux1 (i, j - 1)) / dvol1(j)
			ELSE
				dfdx1 (i, j) = (flux1 (i, j) - flux1 (i, j - 1)) / dx1
			END IF		

		END DO

		! we sum this up with the source term to get the L !	
		DO j = 1, length_step_part_1
			l1 (i, j) = - dfdx1 (i, j) - sp_dim*sa1 (i, j)/r1(j) - sb1 (i, j) 
		END DO

	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
l2(:,:) = 0.0D0

! DO for NM !
DO i = imin2, imax2
	DO j = 1, length_step_part_2
	
		! we find the flux difference at the jth cell center !
		If(sp_dim_i /= 0) THEN
			dfdx2 (i, j) = (r2F(j)**DBLE(sp_dim)*flux2 (i, j) - r2F(j-1)**DBLE(sp_dim)*flux2 (i, j - 1)) / dvol2(j)
		ELSE
			dfdx2 (i, j) = (flux2 (i, j) - flux2 (i, j - 1)) / dx2
		END IF		

	END DO

	! we sum this up with the source term to get the L !	
	DO j = 1, length_step_part_2
		l2 (i, j) = - dfdx2 (i, j) - sp_dim*sa2 (i, j)/r2(j) - sb2 (i, j) 
	END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For Riemann solver !
CALL DESTROYRIEMANN

! We need to deallocate the neccesarry variables !
DEALLOCATE(sa1)
DEALLOCATE(sb1)
DEALLOCATE(flux1)
DEALLOCATE(dfdx1)

DEALLOCATE(sa2)
DEALLOCATE(sb2)
DEALLOCATE(flux2)
DEALLOCATE(dfdx2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine construct the conservative variables and flux !
! function and left and right at each cell boundaries		!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDSTATES
USE DEFINITION
USE RIEMANN_MODULE
IMPLICIT NONE

! Integer !
INTEGER :: i, j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build conserved variables !
If (RUNDM_flag == 1) then
	DO j = 0, length_step_part_1
		uL1 (irho1, j) = rho1L (j)
		uL1 (ivel1, j) = rho1L (j) * vel1L (j)
		uR1 (irho1, j) = rho1R (j)
		uR1 (ivel1, j) = rho1R (j) * vel1R (j)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! No DM internal energy !
		!uL1 (itau1, j) = rho1L (j) * eps1L (j) + 5.0E-1_DP * rho1L (j) * vel1L (j) ** 2
		!uR1 (itau1, j) = rho1R (j) * eps1R (j) + 5.0E-1_DP * rho1R (j) * vel1R (j) ** 2
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	END DO

	! DM scalar equation !
	If(iminsca1 > 0) THEN
		DO i = iminsca1, imaxsca1
			DO j = 0, length_step_part_1
				uL1 (i,j) = rho1L (j) * sca1L(i,j)
				uR1 (i,j) = rho1R (j) * sca1R(i,j)
			END DO
		END DO
	END IF
END IF

! Do the same for NM !
DO j = 0, length_step_part_2
	uL2 (irho2, j) = rho2L (j)
	uL2 (ivel2, j) = rho2L (j) * vel2L (j)
	uR2 (irho2, j) = rho2R (j)
	uR2 (ivel2, j) = rho2R (j) * vel2R (j)
END DO

! NM energy equation !
IF(nm_epsilon == 1) THEN
	DO j = 0, length_step_part_2
		uL2 (itau2, j) = rho2L (j) * eps2L (j) + 5.0E-1_DP * rho2L (j) * vel2L (j) ** 2
		uR2 (itau2, j) = rho2R (j) * eps2R (j) + 5.0E-1_DP * rho2R (j) * vel2R (j) ** 2
	END DO
END IF
IF(dual_energy == 1) THEN 
	DO j = 0, length_step_part_2
		uL2 (ieps2, j) = rhoe2L(j)
		uR2 (ieps2, j) = rhoe2R(j)
	END DO
END IF

! NM scalar equation !
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO j = 0, length_step_part_2
			uL2 (i,j) = rho2L (j) * sca2L(i,j)
			uR2 (i,j) = rho2R (j) * sca2R(i,j)
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build fluxes For DM !
If (RUNDM_flag == 1) then
   	DO i = imin1, imax1, 1
		DO j = 0, length_step_part_1
			fluxL1 (i, j) = uL1 (i, j) * vel1L (j)
			fluxR1 (i, j) = uR1 (i, j) * vel1R (j)
		ENDDO
	ENDDO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
   	! Add the presusre work done term to the energy equation             
       	!DO j = 0, length_step_part_1
	!	fluxL1 (itau1, j) = fluxL1 (itau1, j) + p1L(j) * vel1L(j)
	!	fluxR1 (itau1, j) = fluxR1 (itau1, j) + p1R(j) * vel1R(j)
      	!ENDDO   
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

   	! Add the pressure term to momentum equation !
      	DO j = 0, length_step_part_1
         	fluxL1 (ivel1, j) = fluxL1 (ivel1, j) + p1L(j)
		fluxR1 (ivel1, j) = fluxR1 (ivel1, j) + p1R(j)
      	ENDDO

	! Extra flux term for moving grid !
	If(movinggriddm_flag == 1) THEN
		DO i = imin1, imax1, 1
			DO j = 0, length_step_part_1
         			fluxL1 (i, j) = fluxL1 (i, j) - uL1(i,j)*vf1L(j)
				fluxR1 (i, j) = fluxR1 (i, j) - uR1(i,j)*vf1R(j)
      			ENDDO
		END DO
	END IF
END IF

! For NM !
DO i = imin2, imax2, 1
	DO j = 0, length_step_part_2
		fluxL2 (i, j) = uL2 (i, j) * vel2L (j)
		fluxR2 (i, j) = uR2 (i, j) * vel2R (j)
	ENDDO
ENDDO

! Add the presusre work done term to the energy equation             
IF(nm_epsilon == 1) THEN   
	DO j = 0, length_step_part_2
		fluxL2 (itau2, j) = fluxL2 (itau2, j) + p2L(j) * vel2L(j)
		fluxR2 (itau2, j) = fluxR2 (itau2, j) + p2R(j) * vel2R(j)
	ENDDO 
END IF  
IF(dual_energy == 1) THEN    
	DO j = 0, length_step_part_2
		fluxL2 (ieps2, j) = fluxL2 (ieps2, j) + p2L(j) * vel2L(j)
		fluxR2 (ieps2, j) = fluxR2 (ieps2, j) + p2R(j) * vel2R(j)
	ENDDO 
END IF  	

! Add the pressure term to momentum equation !
DO j = 0, length_step_part_2
	fluxL2 (ivel2, j) = fluxL2 (ivel2, j) + p2L(j)
	fluxR2 (ivel2, j) = fluxR2 (ivel2, j) + p2R(j)
ENDDO

! Extra flux term for moving grid !
If(movinggridnm_flag == 1) THEN
	DO i = imin2, imax2, 1
		DO j = 0, length_step_part_2
         		fluxL2 (i, j) = fluxL2 (i, j) - uL2(i,j)*vf2L(j)
			fluxR2 (i, j) = fluxR2 (i, j) - uR2(i,j)*vf2R(j)
      		ENDDO
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the specific internal energy !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOSEPSILON (den, pre, eps, type)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT (IN) :: type

! Input density !
REAL (DP), INTENT (IN) :: den, pre

! Output value !
REAL (DP), INTENT (OUT) :: eps

! Extra patches by Ivan !
REAL (DP) :: p_poly, eps_poly, eps_thermal

! Fermi-momentum !
REAL (DP) :: fermi

! For DM Output !
IF(type == 1) THEN
	IF(fermieosdm_flag == 1) THEN
		CALL FERMIMO (fermi, den, 1)
		IF (fermi<=1.0E-2_DP) THEN
			eps = a_max1*small_energy(fermi)/den
		ELSE
			eps = a_max1*large_energy(fermi)/den
		END IF
	ELSE
		eps = k1 * den ** (gamma1 - 1.0E0_DP) / (gamma1 - 1.0E0_DP)
	END IF

! For NM !
ELSEIF(type == 2) THEN
	IF(fermieosnm_flag == 1) THEN
		CALL FERMIMO (fermi, den, 2)
		IF (fermi<=1.0E-2_DP) THEN
			eps = a_max2*small_energy(fermi)/den
		ELSE
			eps = a_max2*large_energy(fermi)/den
		END IF
	ELSEIF (ccsneos_flag == 1) THEN
		IF (den < rhoc_b) THEN
			p_poly = kc_1 * den ** (gammac_1)
			eps_poly = Ec_1 * den ** (gammac_1 - 1.0E0_DP)
		ELSE
			p_poly = kc_2 * den ** (gammac_2)
			eps_poly = Ec_2 * den ** (gammac_2 - 1.0E0_DP) + Ec_3
		END IF
		eps_thermal = (pre - p_poly)/(den*(gammac_th - 1.0E0_DP))
		eps = eps_poly + eps_thermal
	ELSE
		IF(nm_epsilon == 0) THEN
			eps = k2 * den ** (gamma2 - 1.0E0_DP) / (gamma2 - 1.0E0_DP)
		ELSE
			eps = pre/den/(gamma2 - 1.0D0)
		END IF
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the speed of sound !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOSSOUNDSPEED(pre, den, eps, cs_out, type)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT(IN) :: type

! Input density !
REAL (DP), INTENT (IN) :: pre, den, eps

! Output value !
REAL (DP), INTENT (OUT) :: cs_out

! Local real variables !
REAL (DP) :: fermi, dpdden, dpdeps, dxdrho

! Extra patches by Ivan !
REAL (DP) :: p_poly, eps_poly, eps_thermal

! We do the DM case first !
IF (type == 1) Then
	IF(fermieosdm_flag == 1) THEN
		CALL FERMIMO (fermi, den, 1)
		CALL FINDDXDRHO (dxdrho, den, 1)
		dpdden = a_max1*dxdrho*dpdx(fermi)
		dpdeps = 0.0E0_DP
		cs_out = sqrt(dpdden+dpdeps*pre/den**(2.0E0_DP))
	ELSE
		dpdden = k1 * gamma1 * den ** (gamma1 - 1.0D0)
		dpdeps = 0.0D0
		cs_out = sqrt(dpdden+dpdeps*pre/den**(2.0E0_DP))
	END IF

! For NM !
ELSEIF(type == 2) THEN
	IF (fermieosnm_flag == 1) THEN
		CALL FERMIMO (fermi, den, 2)
		CALL FINDDXDRHO (dxdrho, den, 2)
		dpdden = a_max2*dxdrho*dpdx(fermi)
		dpdeps = 0.0E0_DP
		cs_out = sqrt(dpdden+dpdeps*pre/den**(2.0E0_DP))
	ELSEIF (ccsneos_flag == 1) THEN
		IF (den < rhoc_b) THEN
			p_poly = kc_1 * den ** (gammac_1)
			eps_poly = Ec_1 * den ** (gammac_1 - 1.0E0_DP)
		ELSE
			p_poly = kc_2 * den ** (gammac_2)
			eps_poly = Ec_2 * den ** (gammac_2 - 1.0E0_DP) + Ec_3
		END IF
		eps_thermal = (pre - p_poly)/(den*(gammac_th - 1.0E0_DP))
		IF (den < rhoc_b) THEN
			dpdden = kc_1 * gammac_1 * den ** (gammac_1 - 1.0E0_DP) + eps_thermal * (gammac_th - 1.0E0_DP)
		ELSE
			dpdden = kc_2 * gammac_2 * den ** (gammac_2 - 1.0E0_DP) + eps_thermal * (gammac_th - 1.0E0_DP)
		END IF
		dpdeps = 0.0D0
		cs_out = sqrt(dpdden+dpdeps*pre/den**(2.0E0_DP))
	ELSE
		IF(nm_epsilon == 0) THEN
			dpdden = k2 * gamma2 * den ** (gamma2 - 1.0D0)
			dpdeps = 0.0D0
			cs_out = sqrt(dpdden+dpdeps*pre/den**(2.0E0_DP))
		ELSE
			dpdden = eps * (gamma2 - 1.0E0_DP)
			dpdeps = den * (gamma2 - 1.0E0_DP)
			cs_out = sqrt(dpdden+dpdeps*pre/den**(2.0E0_DP))
		END IF
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
	real(DP) function dpdx(x)
	implicit none
	real(DP) :: x
	dpdx = 8.0D0*x**4/SQRT(x**2 + 1.0D0)
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the first adiabatic index for gas given certain !
! input of pressure, density and speed of sound 			!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOSGAMMA(pre, den, cs, gam_out)
USE DEFINITION 
IMPLICIT NONE

! Input density !
REAL (DP), INTENT (IN) :: pre, den, cs

! Output value !
REAL (DP), INTENT (OUT) :: gam_out

! Assign !
gam_out = den*cs**2/pre

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the suitable frame velocity v_f appears in !
! v = v_f * r/R so that we can keep track on the whole star with   !
! co-expanding grid 						   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDFRAMEVEL_NM
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
INTEGER :: j

! My patch !
vel2_max = maxval(vel2)
radius2 = r2(length_step_2)

! Adjust grid velocity !
IF(boundary2 < 0.90D0*r2(length_step_2)) THEN
	vel2_max = vel2_max*0.5D0
ELSE
        vel2_max = vel2_max
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the suitable frame velocity v_f appears in !
! v = v_f * r/R so that we can keep track on the whole star with   !
! co-expanding grid 						   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDFRAMEVEL_DM
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
INTEGER :: j

! My patch !
vel1_max = maxval(vel1)
radius1 = r1(length_step_1)

! Adjust grid velocity !
IF(boundary1 < 0.90D0*r1(length_step_1)) THEN
	vel1_max = vel1_max*0.5D0
ELSE
        vel1_max = vel1_max
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find the outermost radius of the star !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDBOUNDARY_DM
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Temporal integer !
INTEGER :: jr

! Initialize !
jr = 1

! For NM !
DO j = 1, length_step_1 - 1, 1
      if(rho1(j) > rho1_a .and. rho1(j+1) <= rho1_a) then
	jr = max(j, jr)
      endif                 
ENDDO
boundary1 = r1(jr)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find the outermost radius of the star !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDBOUNDARY_NM
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Temporal integer !
INTEGER :: jr

! Initialize !
jr = 1

! For NM !
DO j = 1, length_step_2 - 1, 1
      if(rho2(j) > rho2_a .and. rho2(j+1) <= rho2_a) then
	jr = max(j, jr)
      endif                 
ENDDO
boundary2 = r2(jr)

END SUBROUTINE