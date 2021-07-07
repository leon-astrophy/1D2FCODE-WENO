!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine supply a sponge on the stellar !
! surface to avoid spurious oscillation due to	 !
! the imperfect initial model			 !
! Written by Leung Shing Chi in 2016		 !
! based on the Zingale 2002			 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETSPONGE
USE DEFINITION
IMPLICIT NONE

! Dummy variables
integer :: i, j

! Position factor defined in the paper
real (DP) :: r_tp, r_md, r_sp, rad_dist
real (DP) :: f_damp

! Threshold density and damping factor
real (DP) :: rho_sponge = 1.0D-13
real (DP) :: kap_sponge = 4.9282D-3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We do the NM component first !
! Give r_tp and r_sp an arbitrarily large number 
rho_sponge = rho2(1) * 0.01D0
r_md = 1.0D10
r_sp = 1.0D10
i = 0

! First find r_sp according to the prescription in paper
do j = 1, length_step_2, 1
	rad_dist = r1(j)

	if(rho2(j) <= rho_sponge) then
		r_sp = MIN(rad_dist, r_sp)
	endif
	if(rho2(j) <= rho2_a) then
		r_md = MIN(rad_dist, r_md)
	endif
enddo
r_tp = 2.0D0 * r_md - r_sp

! Then assign a relaxed veloicty according
! to the damping strength
do j = 1, length_step_2, 1
	rad_dist = r2(j)

	if(rad_dist > r_tp) then
		f_damp = 1.0D0
	elseif(rad_dist <= r_tp .and. rad_dist > r_sp) then
		f_damp = 0.5D0 * (1.0D0 - DCOS(pi_old * (rad_dist - r_sp)/(r_tp - r_sp)))
	else
		f_damp = 0.0D0
	endif

	vel2(j) = vel2(j) / (1.0D0 + kap_sponge * f_damp * dt)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we do the DM if the DM is presence and movable !
If (RUNDM_flag == 1) THEN

	! Give r_tp and r_sp an arbitrarily large number 
	rho_sponge = rho1(1) * 0.01D0
	r_md = 1.0D10
	r_sp = 1.0D10
	i = 0

	! First find r_sp according to the prescription in paper
	do j = 1, length_step_1, 1
		rad_dist = r1(j)

   		if(rho1(j) <= rho_sponge) then
   	  		r_sp = MIN(rad_dist, r_sp)
   		endif

   		if(rho1(j) <= rho1_a) then
		r_md = MIN(rad_dist, r_md)
		endif
	enddo
	r_tp = 2.0D0 * r_md - r_sp

	! Then assign a relaxed veloicty according
	! to the damping strength
	do j = 1, length_step_1, 1
		rad_dist = r1(j)

		if(rad_dist > r_tp) then
			f_damp = 1.0D0
		elseif(rad_dist <= r_tp .and. rad_dist > r_sp) then
			f_damp = 0.5D0 * (1.0D0 - DCOS(pi_old * (rad_dist - r_sp)/(r_tp - r_sp)))
		else
			f_damp = 0.0D0
		endif

		vel1(j) = vel1(j) / (1.0D0 + kap_sponge * f_damp * dt)
	end do

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine