!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This are the Rienmann problem test for typical 1D hydro code !
! Please refer to Ka Wing's paper for details 		       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETRHO_SIMPLE
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
integer :: j

! Set the polytropic index !
gamma2 = 1.4E0_DP

if(test_model == 1) then
	do j=1,length_step_2
		if(j <= length_step_2*3/10) then
			rho2(j) = 1.0E0_DP
			vel2(j) = 0.75E0_DP
			p2(j) = 1.0E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		else
			rho2(j) = 0.125E0_DP
			vel2(j) = 0.0E0_DP
			p2(j) = 0.1E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
      		endif
	enddo

	call boundary1d_NM (rho2,even)
	call boundary1d_NM (vel2,odd)
	call boundary1d_NM (p2,even)
	call boundary1d_NM (epsilon2,even)
elseif(test_model == 2) then
	do j=1,length_step_2
		if(j <= length_step_2*5/10) then
			rho2(j) = 1.0E0_DP
			vel2(j) = -2.0E0_DP
			p2(j) = 0.4E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		else
			rho2(j) = 1.0E0_DP
			vel2(j) = 2.0E0_DP
			p2(j) = 0.4E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		endif
	enddo
	call boundary1d_NM (rho2,even)
	call boundary1d_NM (vel2,odd)
	call boundary1d_NM (p2,even)
	call boundary1d_NM (epsilon2,even)
elseif(test_model == 4) then
	do j=1,length_step_2
		if(j <= length_step_2*4/10) then
			rho2(j) = 5.99924E0_DP
			vel2(j) = 19.5975E0_DP
			p2(j) = 460.894E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		else
			rho2(j) = 5.99242E0_DP
			vel2(j) = -6.19633E0_DP
			p2(j) = 46.0950E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		endif
	enddo

	call boundary1d_NM (rho2,even)
	call boundary1d_NM (vel2,odd)
	call boundary1d_NM (p2,even)
	call boundary1d_NM (epsilon2,even)
elseif(test_model == 5) then
	do j=1,length_step_2
		if(j <= length_step_2*8/10) then
			rho2(j) = 1.0E0_DP
			vel2(j) = -19.59745E0_DP
			p2(j) = 1000.0E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		else
			rho2(j) = 1.0E0_DP
			vel2(j) = -19.59745E0_DP
			p2(j) = 0.01E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		endif
	enddo

	call boundary1d_NM (rho2,even)
	call boundary1d_NM (vel2,odd)
	call boundary1d_NM (p2,even)
	call boundary1d_NM (epsilon2,even)
elseif(test_model == 3) then
	do j=1,length_step_2
		if(j <= length_step_2*5/10) then
			rho2(j) = 1.0E0_DP
			vel2(j) = 0.0E0_DP
			p2(j) = 1000.0E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		else
			rho2(j) = 1.0E0_DP
			vel2(j) = 0.0E0_DP
			p2(j) = 0.01E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		endif
	enddo

	call boundary1d_NM (rho2,even)
	call boundary1d_NM (vel2,odd)
	call boundary1d_NM (p2,even)
	call boundary1d_NM (epsilon2,even)
elseif(test_model == 6) then
	do j=1,length_step_2
		if(j <= length_step_2*5/10) then
			rho2(j) = 1.4E0_DP
			vel2(j) = 0.0E0_DP
			p2(j) = 1.0E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		else
			rho2(j) = 1.0E0_DP
			vel2(j) = 0.0E0_DP
			p2(j) = 1.0E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		endif
	enddo

	call boundary1d_NM (rho2,even)
	call boundary1d_NM (vel2,odd)
	call boundary1d_NM (p2,even)
	call boundary1d_NM (epsilon2,even)
elseif(test_model == 7) then
	do j=1,length_step_2
		if(j <= length_step_2*5/10) then
			rho2(j) = 1.4E0_DP
			vel2(j) = 0.1E0_DP
			p2(j) = 1.0E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		else
			rho2(j) = 1.0E0_DP
			vel2(j) = 0.1E0_DP
			p2(j) = 1.0E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		endif
	enddo

	call boundary1d_NM (rho2,even)
	call boundary1d_NM (vel2,odd)
	call boundary1d_NM (p2,even)
	call boundary1d_NM (epsilon2,even)
elseif(test_model == 8) then
	do j=1,length_step_2
		if(j <= length_step_2*1/10) then
			rho2(j) = 1.0E0_DP
			vel2(j) = 0.0E0_DP
			p2(j) = 1000.0E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		elseif(j >= length_step_2*9/10) then
			rho2(j) = 1.0E0_DP
			vel2(j) = 0.0E0_DP
			p2(j) = 100.0E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		else
			rho2(j) = 1.0E0_DP
			vel2(j) = 0.0E0_DP
			p2(j) = 0.01E0_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		endif
	enddo

	call boundary1d_NM (rho2,even)
	call boundary1d_NM (vel2,odd)
	call boundary1d_NM (p2,even)
	call boundary1d_NM (epsilon2,even)
elseif(test_model == 9) then
	do j=1,length_step_2
		if(j <= 4) then
			rho2(j) = 1.0E0_DP
			vel2(j) = 0.0E0_DP
			p2(j) = (3.0D0 * (gamma2 - 1.0D0) * 1.0D0)/ &
				(4.0D0 * pi_old * r2(4) ** 3)
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		else
			rho2(j) = 1.0E0_DP
			vel2(j) = 0.0E0_DP
			p2(j) = 1.0E-5_DP
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		endif
	enddo

	call boundary1d_NM (rho2,even)
	call boundary1d_NM (vel2,odd)
	call boundary1d_NM (p2,even)
	call boundary1d_NM (epsilon2,even)
elseif(test_model == 10) then
	do j=1,length_step_2
		if(j < length_step_2/8) then
			rho2(j) = 3.857143D0
			vel2(j) = 2.629369D0
			p2(j) = 31.0D0/3.0D0
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		else
			rho2(j) = 1.0E0_DP + 0.2D0*sin(8.0D0*(DBLE(j) - 0.5D0)*dx2)
			vel2(j) = 0.0E0_DP
			p2(j) = 1.0D0
			epsilon2(j) = p2(j) / rho2(j) / (gamma2-1.0E0_DP)
		endif
	enddo

	call boundary1d_NM (rho2,even)
	call boundary1d_NM (vel2,odd)
	call boundary1d_NM (p2,even)
	call boundary1d_NM (epsilon2,even)
else
	stop "no such test model"
END IF

endsubroutine