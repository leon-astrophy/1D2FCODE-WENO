!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the tracer particle module that records the		!
! thermodynamics history of each particule trajectory		!
! The information will be used to input as data post processing !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE PPT_MODULE
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Number of partcles !
integer, parameter :: nop = 1000

! Thermodynamics history of particles !
real (DP), allocatable, dimension(:) :: posold_p
real (DP), allocatable, dimension(:) :: pos_p
real (DP), allocatable, dimension(:) :: pos2_p
real (DP), allocatable, dimension(:) :: pos3_p
real (DP), allocatable, dimension(:) :: rho_p
real (DP), allocatable, dimension(:) :: vel_p
real (DP), allocatable, dimension(:) :: vel3_p
real (DP), allocatable, dimension(:) :: eps_p
real (DP), allocatable, dimension(:) :: temp_p
real (DP), allocatable, dimension(:) :: phi_p
real (DP), allocatable, dimension(:) :: mass_p
real (DP), allocatable, dimension(:) :: mend_p

! Grid variables !
INTEGER :: grid_p(1:nop)
REAL (DP), allocatable, dimension(:) :: dgrid_p

! Determine whether partilce leave the computational domain !
logical, allocatable, dimension(:) :: leavebox_p

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutione allocate arrays for tracer particles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDPPT
IMPLICIT NONE

! Allocate arrays !
allocate(pos_p(1:nop))
allocate(pos2_p(1:nop))
allocate(pos3_p(1:nop))
allocate(posold_p(1:nop))
allocate(vel_p(1:nop))
allocate(vel3_p(1:nop))    
allocate(rho_p(1:nop))
allocate(eps_p(1:nop))        
allocate(temp_p(1:nop))  
allocate(phi_p(1:nop))
allocate(mass_p(1:nop))
allocate(mend_p(1:nop))
allocate(leavebox_p(1:nop))
allocate(dgrid_p(1:nop))

end subroutine buildPPT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocate arrays for tracer particles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DESTROYPPT
IMPLICIT NONE

! Deallocate arrays !
deallocate(pos_p)
deallocate(posold_p)
deallocate(vel_p)
deallocate(rho_p)
deallocate(eps_p)
deallocate(temp_p)
deallocate(phi_p)
deallocate(mass_p)
deallocate(mend_p)
deallocate(leavebox_p)
deallocate(dgrid_p)

end subroutine destroyPPT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutione output files related to tracer particles !
! trajectory's histroy					    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine outputPPT(n)
USE DEFINITION, ONLY : GLOBAL_TIME
IMPLICIT NONE
   
! integer parameter !
integer :: i, n
integer :: sqrt_nop

! Related to output 
character(len=256) :: fileno
integer :: fileno_len

! Assign related to output file !
write(fileno,*) n
fileno = ADJUSTL(fileno)
fileno_len = LEN_TRIM(fileno)
sqrt_nop = DSQRT(DBLE(nop))

! Open the file 
open(unit=601,file='./Outfile/Tracer/Star_WENO_PPT_'//fileno(1:fileno_len)//'.dat', STATUS = 'REPLACE')

! Write the output file !
write(601,*) n, global_time
write(601,*) nop
do i = 1, nop
	write(601,100) i, rho_p(i), temp_p(i), pos_p(i), pos_p(i), vel_p(i), vel_p(i), phi_p(i), eps_p(i), mass_p(i), mend_p(i), pos_p(i), leavebox_p(i)
enddo 
write(601,*)

! Close the file !
close(601)

100 format (I6, 11ES16.8, L3)

end subroutine outputPPT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This subroutine initialize tracer particles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETPPT
IMPLICIT NONE

! Write out message 
write(*,*) 'Now initialize tracer particles'

! initialize tracer particles !
call getpos_p
call getgrid_p
call getnse_p
call getdata_p

end subroutine getppt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine evolve the tracer particles trajectory !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EVOLVE_P_1ST
USE DEFINITION, ONLY : dt
IMPLICIT NONE

! Integer parameter !
INTEGER :: i

! Back up old position !
posold_p = pos_p      
 
! Get velocity !
call getgrid_p
call getvel_p

! Update the position !
DO i = 1, nop
	pos_p(i) = posold_p(i) + vel_p(i) * 0.391752226571890D0 * dt
ENDDO

! Assign appropriate boundary condition for particles !
DO i = 1, nop
	IF(pos_p(i) < 0.0D0) pos_p(i) = -pos_p(i)
ENDDO

end subroutine evolve_p_1st

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine evolve the tracer particles trajectory !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EVOLVE_P_2ND
USE DEFINITION, ONLY : dt
IMPLICIT NONE

! Integer parameter !
INTEGER :: i

! Get velocity !
call getgrid_p
call getvel_p

! Update the position !
DO i = 1, nop
	pos_p(i) = 0.444370493651235D0*posold_p(i)+0.555629506348765D0*pos_p(i)+0.368410593050371D0*dt*vel_p(i)
	pos2_p(i) = pos_p(i)
ENDDO

! Assign appropriate boundary condition for particles !
DO i = 1, nop
	IF(pos_p(i) < 0.0D0) pos_p(i) = -pos_p(i)
ENDDO

end subroutine evolve_p_2nd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine evolve the tracer particles trajectory !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EVOLVE_P_3RD
USE DEFINITION, ONLY : dt
IMPLICIT NONE

! Integer parameter !
INTEGER :: i

! Get velocity !
call getgrid_p
call getvel_p

! Update the position !
DO i = 1, nop
	pos_p(i) = 0.620101851488403D0*posold_p(i)+0.379898148511597D0*pos_p(i)+0.251891774271694D0*dt*vel_p(i)
      pos3_p(i) = pos_p(i)
      vel3_p(i) = vel_p(i)
ENDDO

! Assign appropriate boundary condition for particles !
DO i = 1, nop
	IF(pos_p(i) < 0.0D0) pos_p(i) = -pos_p(i)
ENDDO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine evolve the tracer particles trajectory !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EVOLVE_P_4TH
USE DEFINITION, ONLY : dt
IMPLICIT NONE

! Integer parameter !
INTEGER :: i

! Get velocity !
call getgrid_p
call getvel_p

! Update the position !
DO i = 1, nop
	pos_p(i) = 0.178079954393132D0*posold_p(i)+0.821920045606868D0*pos_p(i)+0.544974750228521D0*dt*vel_p(i)
ENDDO

! Assign appropriate boundary condition for particles !
DO i = 1, nop
	IF(pos_p(i) < 0.0D0) pos_p(i) = -pos_p(i)
ENDDO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine evolve the tracer particles trajectory !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EVOLVE_P_5TH
USE DEFINITION, ONLY : dt
IMPLICIT NONE

! Integer parameter !
INTEGER :: i

! Get velocity !
call getgrid_p
call getvel_p

! Update the position !
DO i = 1, nop
	pos_p(i) = 0.517231671970585D0*pos2_p(i)+0.096059710526147D0*pos3_p(i)+0.386708617503269D0*pos_p(i) &
		+ 0.063692468666290D0*dt*vel3_p(i) + 0.226007483236906D0*dt*vel_p(i)
ENDDO

! Assign appropriate boundary condition for particles !
DO i = 1, nop
	IF(pos_p(i) < 0.0D0) pos_p(i) = -pos_p(i)
ENDDO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine gets the final thermodynamic data of each particles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EVOLVE_P_FINAL
IMPLICIT NONE

! Integer parameter !
INTEGER :: i

! Get data !
call getgrid_p
call getdata_p

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutione get the tracer particles initial position !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETPOS_P
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
INTEGER :: i, j

! Real parameter !
REAL (DP) :: pos_r
REAL (DP) :: target_mass
REAL (DP) :: mass_shell

! Assign variables !
mass_shell = mass2*0.99999D0 / DBLE(nop)

! Do the loop !
do i = 1, nop

	! We check whether the enclosed mass exceeds the mass shell limit !
	target_mass = DBLE(i) * mass_shell
		do j = 1, length_step_2
			if(m_cell(j) > target_mass) then

				! Use interpolation subroutine !
				mass_p(i) = target_mass
				CALL LINEAR (m_cell(j-1), m_cell(j), r2F(j-1), r2F(j), target_mass, pos_r)
				exit

			END IF
		enddo

	! Now, assign the position of the tracers !
	! We want constant amount of tracers per unit volume !
	! For each mass shell (length in 1D, area in 2D) !
	pos_p(i) = pos_r

enddo

END SUBROUTINE getpos_p


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the grid of tracer particles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETGRID_P
USE DEFINITION
IMPLICIT NONE

! Integer !
integer :: i

! Real variables !
REAL (DP) :: half_dx
REAL (DP) :: ratio, grid_p_raw

! For estimating homologous expansion
REAL (DP) :: rad

! Assign  !
half_dx = dx2 * 0.5D0

! Loop !
DO i = 1, nop, 1
 
	! Assume the hydro quantity is defined on the grid center
	grid_p(i) = INT((pos_p(i) + half_dx) / dx2) + 1
	dgrid_p(i) = (pos_p(i) - (DBLE(grid_p(i)) - 1.5D0) * dx2) / dx2	    

ENDDO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the trajectory history of each tracer particles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETDATA_P
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
INTEGER :: i, j

! Related to grid variables
INTEGER :: grid_x
REAL (DP) :: dgrid_x

! Do the loop !
DO i = 1, nop

      ! Assign !
      grid_x = grid_p(i)
      dgrid_x = dgrid_p(i)

	! Condition !
	IF(leavebox_p(i) .eqv. .false.) THEN
         	if(grid_x > 0 .and. grid_x < length_step_2) then

	    		! Interpolation !
      	    		temp_p(i) = temp2(grid_x-1) * (1.0D0 - dgrid_x) + temp2(grid_x) * dgrid_x
      	    		rho_p(i) = rho2(grid_x-1) * (1.0D0 - dgrid_x) + rho2(grid_x) * dgrid_x
      	    		phi_p(i) = phinm(grid_x-1) * (1.0D0 - dgrid_x) + phinm(grid_x) * dgrid_x

		ENDIF
	ELSE
	    leavebox_p(i) = .true.

      	    if(grid_x > length_step_2) then
      	 	grid_x = length_step_2
      	 	dgrid_x = 0.0D0
      	    endif

            temp_p(i) = temp2(grid_x)
            rho_p(i) = rho2(grid_x)
	    phi_p(i) = phinm(grid_x)
	END IF

ENDDO

END SUBROUTINE getdata_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine get the velocity of each particles 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getvel_p
use definition
implicit none

! Grid variables !
integer :: i
integer :: grid_x
real (DP) :: dgrid_x

! Related to homologus expansion
REAL (DP) :: rad

! Assign 
do i = 1, nop, 1

      grid_x = grid_p(i)
      dgrid_x = dgrid_p(i)

      IF(leavebox_p(i) == .false.) THEN
  
         IF(grid_x > 0 .and. grid_x < length_step_2) THEN

            ! Interpolation !
      	    vel_p(i) = vel2(grid_x-1) * (1.0D0 - dgrid_x) + vel2(grid_x) * dgrid_x

         ELSE

            leavebox_p(i) = .true.
            rad = pos_p(i)
            vel_p(i) = vel_p(i) * pos_p(i) / rad

         ENDIF

      ENDIF

ENDDO

END SUBROUTINE getvel_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the logical flag for whether particles leave the box !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETNSE_P
IMPLICIT NONE

! Integer parameter !
integer :: i

! Initialize !
do i = 1, nop
	leavebox_p(i) = .false.
enddo

end subroutine getnse_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the final mass of tracer particles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETFINAL_P
USE DEFINITION
IMPLICIT NONE

! Integer parameter !
integer :: i, j, k

! Temporial array to store mass !
REAL (DP), DIMENSION (-4 : length_step_2 + 5) :: m_temp

! Initialize temporal mass !
m_temp = 0.0E0_DP

! Integrate to get the enclosed mass !
DO j = 1, length_step_2

	! First, integrate !
	DO k = 1, j
		If (rho2(k) > rho2_a) THEN
			m_temp(j) = m_temp(j) + vol2(k) * rho2(k)
		ELSE
			CYCLE
		END IF
	END DO

END DO

! Get the cell centered velocity values !
do i = 1, nop

	! We find the the mass by integraiton !
	IF(leavebox_p(i) .eqv. .false.) THEN
		IF(pos_p(i) > 0 .and. pos_p(i) < r2F(length_step_2)) THEN
			DO j = 1, length_step_2
				If(r2F(j) > pos_p(i)) THEN
					! linearly interpolated 
					CALL Linear(r2F(j-1), r2F(j), m_temp(j-1), m_temp(j), pos_p(i), mend_p(i))

					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					! Old liner interpolation is not used !
					!mend_p(i) = mass_temp(j-1) + (mass_temp(j) - mass_temp(j-1))* & 
					! 	(pos_p(i) - r2(j-1))/(r2(j) - r2(j-1))
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					EXIT 
				END IF
			END DO
		ENDIF
	END IF
ENDDO

END SUBROUTINE

end module PPT_module
