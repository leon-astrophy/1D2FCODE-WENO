!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D Two Fluid Hydrodynamics simulation Code	   !
! Can also be used as simulating type 1a supernova !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM STAR_WENO
USE WENO_MODULE
USE DEFINITION
USE NUCLEAR_MODULE
USE NU_MODULE
USE FLAME_MODULE
USE PPT_MODULE
IMPLICIT NONE

! PPT !	
integer :: filecount_PPT

! Dummy integer !
INTEGER :: n, j, exno

! Real variables !
! Initialize global simulation time !
global_time = 0.0E0_DP

! Initialize latest output time !
last_outputtime = 0.0E0_DP
last_outputtime_profile1 = 0.0E0_DP
last_outputtime_profile2 = 0.0E0_DP
proindex1 = 0.0E0_DP
proindex2 = 0.0E0_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!'
WRITE(*,*) '!Welcome to Hydrodocde!'
WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!'
WRITE(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! Solving for the initial hydrostatic star !
CALL INITIAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We open files related to hydro variables !
CALL OPENFILEHYDRO

! We output the files related to deflagration !
IF(flame_flag == 1) THEN
	CALL OPENFILEFLAME
END IF

! We output the chemical composition !
IF(xisotran_flag == 1) THEN
	CALL OPENXISOFILE
END IF

! Neutrino spectra !	
IF(nuspec_flag == 1) THEN 
	CALL OPENNEUFILE
END IF

! ppt !
filecount_PPT = 0
IF (tracer_flag == 1) THEN
	call outputPPT(filecount_PPT)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Outputing initial star hydro data into files !
CALL OUTPUTHYDRO (global_time)
CALL OUTPUTPROFILEHYDRO (global_time)

! We output the files related to deflagration !
IF(flame_flag == 1) THEN
	CALL OUTPUTFLAME(global_time)
	CALL OUTPUTPROFILEFLAME(global_time)
END IF

! We output the chemical composition !
IF(xisotran_flag == 1) THEN
	CALL OUTPUTXISOPROFILE(global_time)
	CALL OUTPUTXISO(global_time)
END IF

! Neutrino spectra !
IF(nuspec_flag == 1) THEN
	CALL OUTPUTNEU(global_time)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE (*,*) 'Finish preparing initial data. We do the rungekutta loop'
WRITE (*,*) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Output setting !
output_file = .false.

! Solving for hydrodynamics evolution in time !
DO n = 1, total_time_step

	! Find the corresponding time step dt !
	CALL FINDDT

	! Print Out Unusual time step !
	IF (test_model == 0 .AND. (dt > 1.0E4_DP .OR. dt < 1.0E-10_DP)) THEN
		WRITE (*,*) 'unusual time step'
		WRITE (*,*) dt, global_time, 'time'
		WRITE (*,*) dx1, 'DM grid_size'
		WRITE (*,*) dx2, 'NM grid_size'

		! Output everything before such crash !
		! Outputing initial star hydro data into files !
		CALL OUTPUTHYDRO (global_time)
		CALL OUTPUTPROFILEHYDRO (global_time)

		! We output the files related to deflagration !
		IF(flame_flag == 1) THEN
			CALL OUTPUTFLAME(global_time)
			CALL OUTPUTPROFILEFLAME(global_time)
		END IF

		! We output the chemical composition !
		IF(xisotran_flag == 1) THEN
			CALL OUTPUTXISOPROFILE(global_time)
			CALL OUTPUTXISO(global_time)
		END IF

		! Neutrino spectra !
		IF(nuspec_flag == 1) THEN
			CALL OUTPUTNEU(global_time)
		END IF

		! Quit the simulations !
		EXIT
		
	END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Calling the 3rd order rungekutta time solver !
        CALL RUNGEKUTTA (n)

	! Update the global simulation time !
	global_time = global_time + dt
	
	! Write out timestep and global time !
	WRITE (*,*) n, 'time step'
	WRITE (*,*) length_step_part_2, 'length_step_part_2'
	WRITE (*,*) global_time, 'global time'
	WRITE (*,*) dt, 'dt'
	WRITE (*,*) dx2, 'dx2'
	WRITE (*,*) rho2(length_step_2), 'NM ATM Density'

	! Write out atmospheirc density !
	If(RUNDM_flag == 1) THEN
		WRITE (*,*) dx1, 'dx1'
		WRITE (*,*) rho1(length_step_1), 'DM ATM Density'
		WRITE (*,*) length_step_part_1, 'length_step_part_1'
	END IF
	
	WRITE (*,*)
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! We check whether there are unacceptable output !
	CALL CHECKNAN
	If (foundnan > 0) THEN
		WRITE (*,*) 'unusual output'
		WRITE (*,*) global_time, 'global_time for this happens'
	
		! Outputing initial star hydro data into files !
		CALL OUTPUTHYDRO (global_time)
		CALL OUTPUTPROFILEHYDRO (global_time)

		! We output the files related to deflagration !
		IF(flame_flag == 1) THEN
			CALL OUTPUTFLAME(global_time)
			CALL OUTPUTPROFILEFLAME(global_time)
		END IF

		! We output the chemical composition !
		IF(xisotran_flag == 1) THEN
			CALL OUTPUTXISOPROFILE(global_time)
			CALL OUTPUTXISO(global_time)
		END IF

		! Neutrino spectra !
		IF(nuspec_flag == 1) THEN
			CALL OUTPUTNEU(global_time)
		END IF

		! We exit the simulation !
		EXIT
	
	END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! If there is a detonation, output more frequent profile and log !
	If(founddeton_flag == 1) THEN
		IF(global_time - output_ddt_last >= output_ddt .AND. global_time - found_deton_time <= 1.0D5) THEN
			output_file = .TRUE.
			output_ddt_last = global_time
		END IF
	END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Determine whether datas should be ouputed according to output time step !
	IF (global_time - last_outputtime >= output_time .OR. output_file .eqv. .TRUE.) THEN
		WRITE (*,*) global_time, 'global_time for data ouput'
		! Outputing initial star hydro data into files !
		CALL OUTPUTHYDRO (global_time)

		! We output the files related to deflagration !
		IF(flame_flag == 1) THEN
			CALL OUTPUTFLAME (global_time)
		END IF

		! We output the isotope composition !
		IF(xisotran_flag == 1) THEN
			CALL OUTPUTXISO(global_time)
		END IF

		! Spectra !
		IF(nuspec_flag == 1) THEN
			CALL OUTPUTNEU(global_time)
		END IF

		! Adjust time !
		last_outputtime = global_time		

	END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Determine whether profiles should be ouputed according to output time  !
	! We Output here for geometric series for output_time_profile1 !
	IF (global_time >= output_time_profile1 * 1.0E1_DP ** (proindex1) .OR. output_file .eqv. .TRUE.) THEN
	
		! write the global time and dt !
		WRITE (*,*) global_time, 'global_time for profile ouput'

		! Outputing initial star hydro data into files !
		CALL OUTPUTPROFILEHYDRO (global_time)

		! We output the chemical composition !
		IF(xisotran_flag == 1) THEN
			CALL OUTPUTXISOPROFILE(global_time)
		END IF

		! Flame !
		IF (flame_flag == 1) THEN
			CALL OUTPUTPROFILEFLAME(global_time)
		END IF

		! assign the last output profile time !
		last_outputtime_profile1 = global_time
		
		! update the profile index !
		proindex1 = proindex1 + 1.0E0_DP
		
	END IF

	! Determine whether profiles should be ouputed according to output time  !
	! We Output here for geometric series for output_time_profile2 !
	IF (global_time >= output_time_profile2 * 1.0E1_DP ** (proindex2)) THEN

		! write the global time and dt !
		WRITE (*,*) global_time, 'global_time for profile ouput'

		! Outputing initial star hydro data into files !
		CALL OUTPUTPROFILEHYDRO (global_time)

		! We output the chemical composition !
		IF(xisotran_flag == 1) THEN
			CALL OUTPUTXISOPROFILE(global_time)
		END IF

		IF (flame_flag == 1) THEN
			CALL OUTPUTPROFILEFLAME(global_time)
		END IF

		! Assign the last output profile time !
		last_outputtime_profile2 = global_time
	
		! update the profile index !
		proindex2 = proindex2 + 1.0E0_DP
		
	END IF

	   ! Output tracer particle scheme
   	   IF(tracer_flag == 1) THEN
	      IF(ABS(global_time - output_PPTtime_last) >= output_PPTtime .or. output_file .eqv. .true.) then
                 filecount_PPT = filecount_PPT + 1
                 output_PPTtime_last = global_time
	         CALL outputPPT(filecount_PPT)
	      ENDIF
	   ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Reset !
	output_file = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Check if we need to turn on the moving grid !
	IF(found_movinggridnm_flag == 1 .and. movinggridnm_flag == 0) THEN
		If(length_step_2 - length_step_part_2 <= 10) then
			movinggridnm_flag = 1
			checkstepnm_flag = 0
			WRITE (*,*) 'NM Moving Grid Start'
			WRITE (*,*) 'At Global Time = ', global_time
		END IF
	END IF

	! Check if we need to turn on the moving grid !
	IF(found_movinggriddm_flag == 1 .and. movinggriddm_flag == 0) THEN
		If(length_step_1 - length_step_part_1 <= 10) then
			movinggriddm_flag = 1
			checkstepdm_flag = 0
			WRITE (*,*) 'DM Moving Grid Start'
			WRITE (*,*) 'At Global Time = ', global_time
		END IF
	END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Output everything if we reached total_time !
	If (global_time >= total_time) THEN

		WRITE (*,*) 'We reach the end of simulation, output data'
		WRITE (*,*)

		! Outputing initial star hydro data into files !
		CALL OUTPUTHYDRO (global_time)
		CALL OUTPUTPROFILEHYDRO (global_time)

		! We output the files related to deflagration !
		IF(flame_flag == 1) THEN
			CALL OUTPUTFLAME(global_time)
			CALL OUTPUTPROFILEFLAME(global_time)
		END IF

		! We output the chemical composition !
		IF(xisotran_flag == 1) THEN
			CALL OUTPUTXISOPROFILE(global_time)
			CALL OUTPUTXISO(global_time)
		END IF

		! Spectra !
		IF(nuspec_flag == 1) THEN
			CALL OUTPUTNEU(global_time)
		END IF

	END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Determine whether simulation ends !
	If (global_time >= total_time) THEN
		Exit
	END IF
		
END DO

! Output tracer particle scheme
if(tracer_flag == 1) then
   filecount_PPT = filecount_PPT + 1
   call outputPPT(filecount_PPT)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Special section for SNeIa !
If (flame_flag == 1) THEN
	CALL EJECTEDMASS
END IF

! Write out parameters !
CALL PARAMETER(filecount_PPT)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We close files related to hydro variables !
CALL CLOSEFILEHYDRO

! we close the files related to chemical composition !
IF (xisotran_flag == 1) THEN
	CALL CLOSEXISOFILE
END IF

! Neutrino spectra !
IF(nuspec_flag == 1) THEN
	CALL CLOSENEUFILE
END IF

! We close files related to deflagration !
IF (flame_flag == 1) THEN
	CALL CLOSEFILEFLAME
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SSPRK3 variables !
CALL DESTROYWENO

! We deallocate hydrovariables array !
CALL DESTROYHYDRO

! We deallocate nuclear composition array !
IF (helmeos_flag == 1) THEN
	CALL DESTROYNUCLEAR
END IF

! deallocate deflagration !
IF (flame_flag == 1) THEN
	CALL DESTROYFLAME
END IF

! Neutrino variables !
If(nuspec_flag == 1) THEN
	CALL DESTROYNU
END IF

! Tracer !
if(tracer_flag == 1) then
   call destroyPPT
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Write out end message 1
WRITE (*,*) '!!!!!!!!!!!!!!!!'
WRITE (*,*) '!Simulation end!'
WRITE (*,*) '!!!!!!!!!!!!!!!!'
WRITE (*,*)

END PROGRAM
