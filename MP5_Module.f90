!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the MP5 subroutine that contains all necessary tool for   !
! reconstructing interface values of primitive variables that       !
! will be used as input into the riemann solvers		    !
! for details theoretical aspect of MP5, please refer to the        !
! original MP5 paper Suresh, A., and H. T. Huynh. 		    !
!"Accurate monotonicity-preserving schemes with Runge–Kutta time   !
!stepping." Journal of Computational Physics 136.1 (1997): 83-99.   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MP5_MODULE
USE DEFINITION
IMPLICIT NONE

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculate the required variables on the left hand !
! side of each i + 1/2 call boundary using MP5 method. The same     !
! goes for the right hand side flux, with a symmetric reflection    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RECONSTRUCT_MP5
USE RIEMANN_MODULE
USE FLAME_MODULE
USE NUCLEAR_MODULE
USE DEFINITION  
IMPLICIT NONE

! integer parameter !
INTEGER :: i, j, in

! Dummy !
REAL (DP) :: left, right, dummy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do reconstruction steps for DM !
IF(RUNDM_flag == 1) THEN
	
	! Reconstruct hydrodynamics variables at cell boundaries for DM !
	DO j = -1, length_step_part_1 + 1
		CALL MP5 (j, p1(j-2), p1(j-1), p1(j), p1(j+1), p1(j+2), p1R(j-1), p1L(j))
		CALL MP5 (j, rho1(j-2), rho1(j-1), rho1(j), rho1(j+1), rho1(j+2), rho1R(j-1), rho1L(j))
		CALL MP5 (j, vel1(j-2), vel1(j-1), vel1(j), vel1(j+1), vel1(j+2), vel1R(j-1), vel1L(j))
	END DO
	
	! For DM scalars !
	If(iminsca1 > 0) THEN
		DO i = iminsca1, imaxsca1
			DO j = -1, length_step_part_1 + 1
				CALL MP5 (j, sca1(i,j-2), sca1(i,j-1), sca1(i,j), sca1(i,j+1), sca1(i,j+2), sca1R(i,j-1), sca1L(i,j))
			END DO
		END DO
	END IF

	! Get epsilon at boundary !
	DO j = -1, length_step_part_1 + 1
		CALL EOSEPSILON(rho1R(j), p1R(j), eps1R(j), 1)
		CALL EOSEPSILON(rho1L(j), p1L(j), eps1L(j), 1)
	END DO

	! Get speed of sound at boundary !
	DO j = -1, length_step_part_1 + 1
		CALL EOSSOUNDSPEED(p1R(j), rho1R(j), eps1R(j), cs1R(j), 1)
		CALL EOSSOUNDSPEED(p1L(j), rho1L(j), eps1L(j), cs1L(j), 1)
	END DO

	If(movinggriddm_flag == 1) THEN
		DO j = -1, length_step_part_1 + 1
			vf1R(j) = vel1_max*r1F(j)/radius1
			vf1L(j) = vel1_max*r1F(j)/radius1
		END DO
	END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reconstruct hydrodynamics variables at cell boundaries for NM !
DO j = -1, length_step_part_2 + 1
	CALL MP5 (j, p2(j-2), p2(j-1), p2(j), p2(j+1), p2(j+2), p2R(j-1), p2L(j))
	CALL MP5 (j, rho2(j-2), rho2(j-1), rho2(j), rho2(j+1), rho2(j+2), rho2R(j-1), rho2L(j))
	CALL MP5 (j, vel2(j-2), vel2(j-1), vel2(j), vel2(j+1), vel2(j+2), vel2R(j-1), vel2L(j))
END DO

! Do extra reconstuctions for dual energy !
IF (dual_energy == 1) THEN
	DO j = -1, length_step_part_2 + 1
		CALL MP5 (j, rhoe2(j-2), rhoe2(j-1), rhoe2(j), rhoe2(j+1), rhoe2(j+2), rhoe2R(j-1), rhoe2L(j))
	END DO
END IF

! For NM scalars !
If(iminsca2 > 0) THEN
	DO i = iminsca2, imaxsca2
		DO j = -1, length_step_part_2 + 1
			CALL MP5 (j, sca2(i,j-2), sca2(i,j-1), sca2(i,j), sca2(i,j+1), sca2(i,j+2), sca2R(i,j-1), sca2L(i,j))
		END DO
	END DO
END IF

! Get epsilon at boundary !
IF (dual_energy == 1) THEN
	DO j = -1, length_step_part_2 + 1
		eps2R(j) = rhoe2R(j)/rho2R(j)
		eps2L(j) = rhoe2L(j)/rho2L(j)
	END DO
ELSE
	DO j = -1, length_step_part_2 + 1
		CALL EOSEPSILON(rho2R(j), p2R(j), eps2R(j), 2)
		CALL EOSEPSILON(rho2L(j), p2L(j), eps2L(j), 2)
	END DO
END IF

! Get speed of sound at boundary !
IF(helmeos_flag == 1) THEN
	DO j = -1, length_step_part_2 + 1
		CALL MP5 (j, cs2(j-2), cs2(j-1), cs2(j), cs2(j+1), cs2(j+2), cs2R(j-1), cs2L(j))
	END DO
ELSE
	DO j = -1, length_step_part_2 + 1
		CALL EOSSOUNDSPEED(p2R(j), rho2R(j), eps2R(j), cs2R(j), 2)
		CALL EOSSOUNDSPEED(p2L(j), rho2L(j), eps2L(j), cs2L(j), 2)
	END DO
END IF

! No reconstruction for frame velocity since we assume analytic continous form !
If(movinggridnm_flag == 1) THEN
	DO j = -1, length_step_part_2 + 1
		vf2R(j) = vel2_max*r2F(j)/radius2
		vf2L(j) = vel2_max*r2F(j)/radius2
	END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculate the required variables on the left hand !
! side of each i + 1/2 call boundary using MP5 method. The same     !
! goes for the right hand side flux, with a symmetric reflection    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MP5 (i, vm2, vm1, v, vp1, vp2, vm_out, vp_out)
USE DEFINITION 
IMPLICIT NONE

! Input integer !
INTEGER, INTENT(IN) :: i

! The input into the subroutine, including conservative variable and input flux function !
REAL (DP), INTENT (IN) :: vm2, vm1, v, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL (DP), INTENT (OUT) :: vm_out, vp_out

! Interpolation coefficients ! 
REAL (DP) :: wrm2, wrm1, wrc, wrp1, wrp2
REAL (DP) :: wlm2, wlm1, wlc, wlp1, wlp2

! Intermediate left and right state !
REAL (DP) :: tmpl, tmpr, vl, vr, ppml, ppmr

! Max and min reconstructed variables !
REAL (DP) :: vminr, vmaxr, vminl, vmaxl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get coefficients !
wrm2 = 1.0D0/3.0D1
wrm1 = - 1.3D1/6.0D1
wrc = 4.7D1/6.0D1
wrp1 = 9.0D0/2.0D1
wrp2 = - 1.0D0/2.0D1
wlm2 = - 1.0D0/2.0D1
wlm1 = 9.0D0/2.0D1
wlc = 4.7D1/6.0D1
wlp1 = - 1.3D1/6.0D1
wlp2 = 1.0D0/3.0D1

! First interpolate using 5th order interpolation
tmpr = (wrm2 * vm2 + wrm1 * vm1 + wrc * v + wrp1 * vp1 + wrp2 * vp2)
tmpl = (wlm2 * vm2 + wlm1 * vm1 + wlc * v + wlp1 * vp1 + wlp2 * vp2)

! Find the MP-limited left and right state !
CALL MPLIMITER (tmpr, vm2, vm1, v, vp1, vp2, vr, vminr, vmaxr)
CALL MPLIMITER (tmpl, vp2, vp1, v, vm1, vm2, vl, vminl, vmaxl)

! Assign output !
vm_out = vl
vp_out = vr 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Further reconsturct usnig HMP5 scheme !
!CALL HMP5(vm2, vm1, v, vp1, vp2, ppmr, ppml)

! Further limit the face value using PPM face value !
!vm_out = 0.5D0*(vr + ppmr) - 0.5D0*sign(1.0D0,(vr - vminr)*(vr - vmaxr))*(vr - ppmr)
!vp_out = 0.5D0*(vl + ppml) - 0.5D0*sign(1.0D0,(vl - vminl)*(vl - vmaxl))*(vl - ppml)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine did a classic PPM reconstruction as an input to the 			       !
! hybrid PPM-MP5 scheme stated in the paper :						       !
! "Hybrid monotonicity-preserving piecewise parabolic method for compressible Euler equations" !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HMP5 (vm2, vm1, v, vp1, vp2, vhalfp_out, vhalfm_out)
USE DEFINITION
IMPLICIT NONE

! The input into the subroutine, including conservative variable and input flux function !
REAL (DP), INTENT (IN) :: vm2, vm1, v, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL (DP), INTENT (OUT) :: vhalfp_out, vhalfm_out

! Temporal variables !
REAL (DP) :: vr, vl, dv, dvp, dvm, delv, v6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Apply monotonicty contraints !
If((vp1 - v)*(vm1 - v) < 0.0D0) THEN
	dv = min(0.5D0*abs(vp1 - vm1),2.0D0*abs(v - vm1),2.0D0*abs(vp1 - v))*sign(1.0D0,0.5D0*(vp1 - vm1))
ELSE
	dv = 0.0D0
END IF
If((vp2 - vp1)*(v - vp1) < 0.0D0) THEN
	dvp = min(0.5D0*abs(vp2 - v),2.0D0*abs(vp1 - v),2.0D0*abs(vp2 - vp1))*sign(1.0D0,0.5D0*(vp2 - v))
ELSE
	dvp = 0.0D0
END IF
If((v - vm1)*(vm2 - vm1) < 0.0D0) THEN
	dvm = min(0.5D0*abs(v - vm2),2.0D0*abs(vm1 - vm2),2.0D0*abs(v - vm1))*sign(1.0D0,0.5D0*(v - vm2))
ELSE
	dvm = 0.0D0
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Interpolate !
vr = v + 0.5D0*(vp1 - v) + (dv - dvp)/6.0D0
vl = vm1 + 0.5D0*(v - vm1) + (dvm - dv)/6.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign some parameters !
delv = vr - vl
v6 = 6.0D0*(v - 0.5D0*(vr + vl))

! Apply further constraints !
If(abs(delv) < abs(v6)) THEN

	! Case by case !
	IF((vr - v)*(vl - v) >= 0.0D0) THEN
		vr = v
		vl = v
	ELSEIF ((vr - vl)*(v - 0.5D0*(vr + vl)) > ((vr - vl)**2)/6.0D0) THEN
		vl = 3.0D0*v - 2.0D0*vr
	ELSEIF ((vr - vl)*(v - 0.5D0*(vr + vl)) < ((vr - vl)**2)/6.0D0) THEN
		vr = 3.0D0*v - 2.0D0*vl
	END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign output !
vhalfp_out = vr
vhalfm_out = vl

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is a non-linear function of limiter that takes input !
! as f(vdag, vm2, vm1, v0, vp1, vp2) so a MP limited output is given   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MPLIMITER (vdag, vm2, vm1, v0, vp1, vp2, vlim, vmin_out, vmax_out)
USE DEFINITION
IMPLICIT NONE

! Some important parameter !
REAL (DP) :: alpha = 4.0D0 !4.0E0_DP
REAL (DP) :: tor = 0.0D0 !1.0E-10_DP

! Temporal variables !
REAL (DP) :: mmd2, mmd4, med

! The input into the subroutine, including conservative variable and input flux function !
REAL (DP), INTENT (IN) :: vdag, vm2, vm1, v0, vp1, vp2

! The output of the subroutine, the flux at cell boundary !
REAL (DP), INTENT (OUT) :: vlim, vmin_out, vmax_out

! Some intermediate variables !
REAL (DP) :: vmp, dm1, d, dp1, vul, vav, dm4p12, dm4m12
REAL (DP) :: vmd, vlc, vmin, vmax, rec, deltap, deltam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate vmp !
CALL minmod2(vp1 - v0, alpha*(v0 - vm1), mmd2)
vmp = v0 + mmd2

! Determine whether we use the interpolated value !
If ((vdag - v0)*(vdag - vmp) <= tor) THEN
	vlim = vdag
ELSE
	dm1 = vm2 + v0 - 2.0E0_DP * vm1
	d = vm1 + vp1 - 2.0E0_DP * v0
	dp1 = v0 + vp2 - 2.0E0_DP * vp1

	CALL minmod4(4.0E0_DP * d - dp1, 4.0E0_DP * dp1 - d, d, dp1, mmd4)
	dm4p12 = mmd4
	
	CALL minmod4(4.0E0_DP * dm1 - d, 4.0E0_DP * d - dm1, dm1, d, mmd4)
	dm4m12 = mmd4	

	vul = v0 + alpha * (v0 - vm1)
	vav = 0.5E0_DP * (v0 + vp1)
	vmd = vav - 0.5E0_DP * dm4p12
	vlc = v0 + 0.5E0_DP * (v0 - vm1) + (4.0E0_DP)/(3.0E0_DP) * dm4m12
	vmin = max(min(v0, vp1, vmd), min(v0, vul, vlc))
	vmax = min(max(v0, vp1, vmd), max(v0, vul, vlc))
	vmin_out = vmin
	vmax_out = vmax
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Finally modify the edges value !
	!CALL MEDIAN(vdag, vmin, vmax, med)
	!vlim = med
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Unknown origins of modifying the face values !
	CALL MINMOD2(vmax - vdag, vmin - vdag, mmd2)
	vlim = vdag + mmd2
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the minmod function (also limiter) introduced in the paper !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MINMOD2(x, y, out)
USE DEFINITION
IMPLICIT NONE

! input variables !
REAL (DP), INTENT (IN) :: x, y

! output variables !
REAL (DP), INTENT (OUT) :: out

! assign the output !
out = 0.5E0_DP*(sign(1.0E0_DP,x) + sign(1.0E0_DP,y))*min(abs(x), abs(y))

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the minmod function (also limiter) introduced in the paper !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MINMOD4(w, x, y, z, out)
USE DEFINITION
IMPLICIT NONE

! input variables !
REAL (DP), INTENT (IN) :: w, x, y, z

! output variables !
REAL (DP), INTENT (OUT) :: out

! assign the output !
out = 0.125E0_DP * ( sign(1.0E0_DP,w) + sign(1.0E0_DP,x) ) &
		* abs( (sign(1.0E0_DP,w) + sign(1.0E0_DP,y)) &
		* (sign(1.0E0_DP,w) + sign(1.0E0_DP,z)) ) &
		* min( abs(w), abs(x), abs(y), abs(z) )

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the median function introduced in the paper !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MEDIAN(x ,y, z, out)
USE DEFINITION
IMPLICIT NONE

! input variables !
REAL (DP), INTENT (IN) :: x, y, z

! output variables !
REAL (DP), INTENT (OUT) :: out

! intermediate variables !
REAL (DP) :: mmd2

CALL minmod2( y - x, z - x, mmd2)

out = x + mmd2

END SUBROUTINE

END MODULE