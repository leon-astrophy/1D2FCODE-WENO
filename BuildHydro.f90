!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine allocate arrays related to hydrodynamics variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDHYDRO
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We allocate some neccesarry arrays for potentials !
ALLOCATE(phinm(-4 : length_step_2 + 5))
ALLOCATE(phinm_1(-4 : length_step_2 + 5))
ALLOCATE(phinm_2(-4 : length_step_2 + 5))
ALLOCATE(phip_nm(-4 : length_step_2 + 5))
ALLOCATE(phipface_2(-4 : length_step_2 + 5))

! GR potentail !
ALLOCATE(phigr(-4 : length_step_2 + 5))

! Distance and volume arrays !
ALLOCATE(r2(-4 : length_step_2 + 5))
ALLOCATE(r2L(-4 : length_step_2 + 5))
ALLOCATE(r2R(-4 : length_step_2 + 5))
ALLOCATE(r2F(-4 : length_step_2 + 5))
ALLOCATE(vol2(-4 : length_step_2 + 5))
ALLOCATE(dvol2(-4 : length_step_2 + 5))

! Mass coordinate
ALLOCATE(m_r(-4 : length_step_2 + 5))
ALLOCATE(m_cell(-4 : length_step_2 + 5))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We determine whether users wants the precense of DM, then allocate array !
If (DM_flag == 1) THEN

	! For DM distance and volume !
	ALLOCATE(r1(-4 : length_step_1 + 5))
	ALLOCATE(r1L(-4 : length_step_1 + 5))
	ALLOCATE(r1R(-4 : length_step_1 + 5))
	ALLOCATE(r1F(-4 : length_step_1 + 5))
	ALLOCATE(vol1(-4 : length_step_1 + 5))
	ALLOCATE(dvol1(-4 : length_step_1 + 5))

	! For DM gravitational potentials !
	ALLOCATE(phidm(-4 : length_step_1 + 5))
	ALLOCATE(phidm_1(-4 : length_step_1 + 5))
	ALLOCATE(phidm_2(-4 : length_step_1 + 5))
	ALLOCATE(phip_dm(-4 : length_step_1 + 5))
	ALLOCATE(phipface_1(-4 : length_step_1 + 5))

	! DM Hydrodynamic Variables !	
	ALLOCATE(rho1(-4 : length_step_1 + 5))
	ALLOCATE(p1(-4 : length_step_1 + 5))
	ALLOCATE(epsilon1(-4 : length_step_1 + 5))
	ALLOCATE(cs1(-4 : length_step_1 + 5))

	! The derivatives !
	ALLOCATE(dpdrho1(-4 : length_step_1 + 5))
	ALLOCATE(dpdepsilon1(-4 : length_step_1 + 5))
	ALLOCATE(dpdr1(-4 : length_step_1 + 5))
		
	! Allocate velocity for movable DM !
	IF (RUNDM_flag == 1) THEN	
		ALLOCATE(vel1(-4 : length_step_1 + 5))
	END IF	
	
	! Build array for fermi eos DM !
	IF (fermieosdm_flag == 1 ) THEN
		ALLOCATE(eostable1(eoslineno,2))
	END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! NM hydrodynamic variables !
ALLOCATE(p2(-4 : length_step_2 + 5))
ALLOCATE(rho2(-4 : length_step_2 + 5))
ALLOCATE(epsilon2(-4 : length_step_2 + 5))
ALLOCATE(vel2(-4 : length_step_2 + 5))
ALLOCATE(temp2(-4:length_step_2 + 5))
ALLOCATE(ye2(-4 : length_step_2 + 5))
ALLOCATE(cs2(-4 : length_step_2 + 5))
ALLOCATE(chempo2(-4 : length_step_2 + 5))

! For dual energy formalism !
ALLOCATE(epsp2(-4 : length_step_2 + 5))
ALLOCATE(rhoe2(-4 : length_step_2 + 5))
ALLOCATE(bige2(-4 : length_step_2 + 5))
ALLOCATE(et2(-4 : length_step_2 + 5))

! Ther derivatives !
ALLOCATE(dpdrho2(-4 : length_step_2 + 5))
ALLOCATE(dpdepsilon2(-4 : length_step_2 + 5))
ALLOCATE(dpdr2(-4 : length_step_2 + 5))

! Build array for fermi eos NM !
IF (fermieosnm_flag == 1) THEN
	ALLOCATE(eostable2(eoslineno,2))
END IF

! Array for CCSN !
IF (ccsneos_flag == 1) THEN
	ALLOCATE(p_p(-4 : length_step_2 + 5))
	ALLOCATE(p_th(-4 : length_step_2 + 5))
	ALLOCATE(eps_p(-4 : length_step_2 + 5))
	ALLOCATE(eps_th(-4 : length_step_2 + 5))
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocate arrays related to hydro variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DESTROYHYDRO
USE DEFINITION
IMPLICIT NONE

! Deallocate NM gravitational potentials !
DEALLOCATE(phinm)
DEALLOCATE(phinm_1)
DEALLOCATE(phinm_2)
DEALLOCATE(phip_nm)
DEALLOCATE(phipface_2)

! GR potentail !
DEALLOCATE(phigr)

! Deallocate NM Distance and volume !
DEALLOCATE(r2)
DEALLOCATE(r2L)
DEALLOCATE(r2R)
DEALLOCATE(r2F)
DEALLOCATE(vol2)
DEALLOCATE(dvol2)

! Deallocate Mass coordinate !
DEALLOCATE(m_r)
DEALLOCATE(m_cell)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Deallocate DM arrays !
If (DM_flag == 1) THEN
	! Deallocate DM Distance and volume !
	DEALLOCATE(r1)
	DEALLOCATE(r1L)
	DEALLOCATE(r1R)
	DEALLOCATE(r1F)
	DEALLOCATE(vol1)
	DEALLOCATE(dvol1)

	! Deallocate DM gravitational potentials !
	DEALLOCATE(phidm)
	DEALLOCATE(phidm_1)
	DEALLOCATE(phidm_2)
	DEALLOCATE(phip_dm)
	DEALLOCATE(phipface_1)

	! Hydro !
	DEALLOCATE(rho1)
	DEALLOCATE(p1)
	DEALLOCATE(epsilon1)	
	DEALLOCATE(cs1)
	DEALLOCATE(dpdrho1)
	DEALLOCATE(dpdepsilon1)
	DEALLOCATE(dpdr1)
	IF (RUNDM_flag == 1) THEN	
		DEALLOCATE(vel1)
	END IF	
	IF (fermieosdm_flag == 1 ) THEN
		DEALLOCATE(eostable1)
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We deallocate NM arrays !	
DEALLOCATE(epsilon2)
DEALLOCATE(rho2)
DEALLOCATE(p2)
DEALLOCATE(vel2)
DEALLOCATE(temp2)
DEALLOCATE(ye2)
DEALLOCATE(chempo2)
DEALLOCATE(cs2)
DEALLOCATE(epsp2)
DEALLOCATE(rhoe2)
DEALLOCATE(bige2)
DEALLOCATE(et2)
DEALLOCATE(dpdrho2)
DEALLOCATE(dpdepsilon2)
DEALLOCATE(dpdr2)
IF (fermieosnm_flag == 2) THEN
	DEALLOCATE(eostable2)
END IF

! Section for CCSN !
IF (ccsneos_flag == 1) THEN
	DEALLOCATE(p_p)
	DEALLOCATE(p_th)
	DEALLOCATE(eps_p)
	DEALLOCATE(eps_th)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE