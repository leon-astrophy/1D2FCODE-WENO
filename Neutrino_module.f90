!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains subroutine that are used to calculate the neutrino !
! emisstivity and neutrino spectrum emitted from the supernovae		  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE NU_MODULE
USE DEFINITION, only : DP
IMPLICIT NONE
INCLUDE 'const.dek'

! theta is sin**2(theta_weinberg) = 0.2223 (2019) !
REAL (DP), PARAMETER :: sin2theta  = 0.2223d0

! xnufam is the number of neutrino flavors OTHER THAN ELECTRON !
REAL (DP), PARAMETER :: xnufam = 2.0d0

! Fermi constant of weak interactions !
REAL (DP), PARAMETER :: gf = 1.435D-49

! Fermi constant of weak interactions in cgs !
REAL (DP), PARAMETER :: gf2 = (1.1663787D-23)*(hbar*clight)**3/(ev2erg)**2

! Electron charge !
REAL (DP), PARAMETER :: charge = 1.60217662D-19

! axial and vector coupling constants !
! for mu and tau, they are cvp and cap !
REAL (DP), PARAMETER :: cvp = -0.5D0 + 2.0D0*sin2theta
REAL (DP), PARAMETER :: cap = -0.5D0
REAL (DP), PARAMETER :: cve = cvp + 1.0D0
REAL (DP), PARAMETER :: cae = cap + 1.0D0
REAL (DP), PARAMETER :: cv2 = cve**2 + xnufam * cvp**2
REAL (DP), PARAMETER :: ca2 = cae**2 + xnufam * cap**2

! Factor in front of pair emisstivity, in cgs !
REAL (DP), PARAMETER :: pair0 = (gf**2*me**8*clight**7)/(1.8D1*pi**5*hbar**10)
REAL (DP), PARAMETER :: pair1 = pair0*me*clight**2

! Emisstivity factor !
REAL (DP), PARAMETER :: plas0 = ((me*clight)**9*gf2**2)/(hbar**10*96.0D0*pi**4*fine)*cv2

! Emisstivity factor for plasma neutrino in cgs (MeV-8cm-3s-1) !
REAL (DP), PARAMETER :: plas1 = 2.115D30
REAL (DP), PARAMETER :: plas2 = plas1/8.0D0
REAL (DP), PARAMETER :: plas3 = plas1/157.5D0

! riemann zeta function !
REAL (DP), PARAMETER :: zeta5 = 1.036927755143369926331
REAL (DP), PARAMETER :: zeta52 = zeta5**2

! Fitting for spectra !
REAL (DP), PARAMETER :: a0 = 1.018192299D0
REAL (DP), PARAMETER :: alpha0 = 3.180657028D0
REAL (DP), PARAMETER :: bigA = 0.1425776426D0

! My patch, liquid metal phase coefficients !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: aiso, iiso
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: eiso, piso
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: biso, jiso
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: fiso, qiso
REAL (DP), ALLOCATABLE, DIMENSION(:) :: ciso, diso, giso, hiso
REAL (DP), ALLOCATABLE, DIMENSION(:) :: kiso, liso, riso, siso

! My patch, liquid metal phase coefficients !
REAL (DP), ALLOCATABLE, DIMENSION(:,:) :: alphaiso, betaiso

! My patch, plasmon neutrino coefficients !
REAL (DP), ALLOCATABLE, DIMENSION(:) :: p_plas
REAL (DP), ALLOCATABLE, DIMENSION(:) :: c_plas
REAL (DP), ALLOCATABLE, DIMENSION(:) :: d_plas
REAL (DP), ALLOCATABLE, DIMENSION(:) :: a_plas
REAL (DP), ALLOCATABLE, DIMENSION(:) :: b_plas
REAL (DP), ALLOCATABLE, DIMENSION(:) :: s_plas
REAL (DP), ALLOCATABLE, DIMENSION(:) :: r_plas
REAL (DP), ALLOCATABLE, DIMENSION(:) :: q_plas

! Energy loss by neutrino !
REAL (DP), allocatable, DIMENSION (:) :: total_nu_qdot
REAL (DP), allocatable, DIMENSION (:,:) :: nu_qdot

! Neutrino number emisstivity !
REAL (DP), allocatable, DIMENSION (:) :: total_nu_pair
REAL (DP), allocatable, DIMENSION (:) :: total_nu_plas
REAL (DP), allocatable, DIMENSION (:,:) :: nu_pair
REAL (DP), allocatable, DIMENSION (:,:) :: nu_plas

! For Back up !
REAL (DP), allocatable, DIMENSION (:) :: nu_qdot_old
REAL (DP), allocatable, DIMENSION (:) :: nu_pair_old
REAL (DP), allocatable, DIMENSION (:) :: nu_plas_old

! Total neutrino number generation !
REAL (DP), allocatable, DIMENSION (:) :: total_loss
REAL (DP), allocatable, DIMENSION (:) :: total_pair
REAL (DP), allocatable, DIMENSION (:) :: total_plas

! Time elapse for neutrino integration !
REAL (DP) :: t_neutrino, t_neutrino_old

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Build the fitting coefficients for liquid metal phase !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDNU
USE NUCLEAR_MODULE, ONLY : total_ion
USE DEFINITION, ONLY : length_step_2
IMPLICIT NONE

! Build liquid metal phase coefficient !
ALLOCATE(aiso(1:total_ion, 0:5))
ALLOCATE(iiso(1:total_ion, 0:5))
ALLOCATE(eiso(1:total_ion, 0:5))
ALLOCATE(piso(1:total_ion, 0:5))
ALLOCATE(biso(1:total_ion, 1:4))
ALLOCATE(jiso(1:total_ion, 1:4))
ALLOCATE(fiso(1:total_ion, 1:4))
ALLOCATE(qiso(1:total_ion, 1:4))

! Build liquid metal phase coefficient !
ALLOCATE(ciso(1:total_ion))
ALLOCATE(diso(1:total_ion))
ALLOCATE(giso(1:total_ion))
ALLOCATE(hiso(1:total_ion))
ALLOCATE(kiso(1:total_ion))
ALLOCATE(liso(1:total_ion))
ALLOCATE(riso(1:total_ion))
ALLOCATE(siso(1:total_ion))

! Build liquid metal phase coefficient !
ALLOCATE(alphaiso(1:total_ion, 0:3))
ALLOCATE(betaiso(1:total_ion, 0:3))

! Build plasmon neutrino coefficients !
ALLOCATE(p_plas(1:10))
ALLOCATE(c_plas(1:2))
ALLOCATE(d_plas(1:2))
ALLOCATE(a_plas(1:2))
ALLOCATE(b_plas(1:2))
ALLOCATE(s_plas(1:3))
ALLOCATE(r_plas(1:3))
ALLOCATE(q_plas(1:5))

! allocate neutrino energy lost !
allocate(total_nu_qdot(1:5))
allocate(nu_qdot(1:5,-4:length_step_2 + 5))

! allocate neutrino energy emisstivity !
allocate(total_nu_pair(1:5))
allocate(total_nu_plas(1:5))
allocate(nu_pair(1:5,-4:length_step_2 + 5))
allocate(nu_plas(1:5,-4:length_step_2 + 5))

! for backup !
allocate(nu_qdot_old(1:5))
allocate(nu_pair_old(1:5))
allocate(nu_plas_old(1:5))

! allocate total neutrino number generated !
allocate(total_loss(1:5))
allocate(total_pair(1:5))
allocate(total_plas(1:5))

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate the fitting coefficients for liquid metal phase !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DESTROYNU
IMPLICIT NONE

! Build liquid metal phase coefficient !
DEALLOCATE(aiso)
DEALLOCATE(iiso)
DEALLOCATE(eiso)
DEALLOCATE(piso)
DEALLOCATE(biso)
DEALLOCATE(jiso)
DEALLOCATE(fiso)
DEALLOCATE(qiso)

! Build liquid metal phase coefficient !
DEALLOCATE(ciso)
DEALLOCATE(diso)
DEALLOCATE(giso)
DEALLOCATE(hiso)
DEALLOCATE(kiso)
DEALLOCATE(liso)
DEALLOCATE(riso)
DEALLOCATE(siso)

! Build liquid metal phase coefficient !
DEALLOCATE(alphaiso)
DEALLOCATE(betaiso)

! Build plasmon neutrino coefficients !
DEALLOCATE(p_plas)
DEALLOCATE(c_plas)
DEALLOCATE(d_plas)
DEALLOCATE(a_plas)
DEALLOCATE(b_plas)
DEALLOCATE(s_plas)
DEALLOCATE(r_plas)
DEALLOCATE(q_plas)

! deallocate neutrino energy lost !
deallocate(total_nu_qdot)
deallocate(nu_qdot)

! deallocate neutrino energy emisstivity !
deallocate(total_nu_pair)
deallocate(total_nu_plas)
deallocate(nu_pair)
deallocate(nu_plas)

! For back up !
deallocate(nu_qdot_old)
deallocate(nu_pair_old)
deallocate(nu_plas_old)

! Deallocate total neutrino number generated !
deallocate(total_loss)
deallocate(total_pair)
deallocate(total_plas)

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine reads in the coefficient from tables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READCOFF
IMPLICIT NONE

! Dummy arrays !
double precision, DIMENSION(1:7,0:23) :: dummy, dummy2
double precision, DIMENSION(1:7,0:7) :: dummy3

! Integer !
INTEGER :: m, n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Open file !
open(unit=601, file='./Library/coff.dat', action='read')
open(unit=602, file='./Library/coff2.dat', action='read')
open(unit=603, file='./Library/coff3.dat', action='read')

! Read in coefficients !
DO m = 0, 23
	read(601,*) (dummy(n,m), n = 1, 7)
END DO
DO m = 0, 23
	read(602,*) (dummy2(n,m), n = 1, 7)
END DO
DO m = 0, 7
	read(603,*) (dummy3(n,m), n = 1, 7)
END DO

! Close !
close(601)
close(602)
close(603)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read in coefficients !
DO m = 0, 23
	If (m <= 5) THEN
		aiso(:, m) = dummy(:, m)
		iiso(:, m) = dummy2(:, m)
	ELSEIF (m > 5 .AND. m <= 9) THEN
		biso(:, m - 5) = dummy(:, m)
		jiso(:, m - 5) = dummy2(:, m)
	ELSEIF (m == 10) THEN
		ciso(:) = dummy(:, m)
		kiso(:) = dummy2(:, m)
	ELSEIF (m == 11) THEN
		diso(:) = dummy(:, m)
		liso(:) = dummy2(:, m)
	ELSEIF (m > 11 .AND. m <= 17) THEN
		eiso(:,m - 12) = dummy(:, m)
		piso(:,m - 12) = dummy2(:, m)
	ELSEIF (m > 17 .AND. m <= 21) THEN
		fiso(:,m - 17) = dummy(:, m)
		qiso(:,m - 17) = dummy2(:, m)
	ELSEIF (m == 22) THEN
		giso(:) = dummy(:, m)
		riso(:) = dummy2(:, m)
	ELSEIF (m == 23) THEN
		hiso(:) = dummy(:, m)
		siso(:) = dummy2(:, m)
	END IF
END DO
DO m = 0, 7
	If (m <= 3) THEN
		alphaiso(:, m) = dummy3(:,m)
	ELSE
		betaiso(:, m - 4) = dummy3(:,m)
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign all fitting constants to plasmon neutrino decay !
p_plas(1) = 1.793D0
p_plas(2) = 0.0645D0
p_plas(3) = 0.433D0
p_plas(4) = 0.01139D0
p_plas(5) = 2.484D6
p_plas(6) = -0.6195D0
p_plas(7) = 0.0009632D0
p_plas(8) = 0.4372D0
p_plas(9) = 1.614D0
p_plas(10) = 8.504D8
s_plas(1) = 9.079D0
s_plas(2) = 1.399D0
s_plas(3) = -0.06592D0
r_plas(1) = 0.3520D0
r_plas(2) = 1.195D0
r_plas(3) = -0.1060D0
q_plas(1) = 0.7886D0
q_plas(2) = 0.2642D0
q_plas(3) = 1.024D0
q_plas(4) = 0.07839D0
q_plas(5) = 0.1784D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialize the time integrated total neutrino emission !
! And total neutrino energy generation rate to zero			 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETNUSPEC
USE DEFINITION 
IMPLICIT NONE

total_loss = 0.0D0
total_pair = 0.0D0
total_plas = 0.0D0

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine open the files for outputing the isotope    !
! composition where we expect it should be changing with time !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENNEUFILE
USE DEFINITION
IMPLICIT NONE

! Integer variables !
INTEGER ::  fileno_len
CHARACTER (len = 256) :: fileno
 
WRITE (fileno, *) n_backup / time_step + 1
fileno = ADJUSTL (fileno)
fileno_len = LEN_TRIM (fileno)

! Energy loss
OPEN (UNIT = 801, FILE = './Outfile/Neutrino/Star_WENO_ThermoNeutrinoEnergyLoss_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 802, FILE = './Outfile/Neutrino/Star_WENO_PairNeutrinoSpectra_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')
OPEN (UNIT = 803, FILE = './Outfile/Neutrino/Star_WENO_PlasNeutrinoSpectra_'//fileno (1 : fileno_len)//'.dat', STATUS = 'REPLACE')

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine close the files for outputing the isotope   !
! composition where we expect it should be changing with time !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CLOSENEUFILE
USE DEFINITION
IMPLICIT NONE

CLOSE(801)
CLOSE(802)
CLOSE(803)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
! The subroutine bridges the GetNuQdot and finds the energy !
! loss rate by neutrino of the hydro grids. 		    !
! Written by Leung Shing Chi in 2016			    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDNEUTRINOLOSS()
USE DEFINITION
USE NUCLEAR_MODULE
IMPLICIT NONE

! dummy variables
integer :: j

! local variables
real (DP) :: temp_in, rho_in, abar_in, zbar_in, ye_in, dummy, eta_dummy
real (DP), DIMENSION(1:total_ion) :: xiso_in   

! initialization
nu_qdot = 0.0D0
nu_pair = 0.0D0
nu_plas = 0.0D0

! Find !
DO j = 1, length_step_part_2

	! Pass the data to local storage
	rho_in = (rho2 (j)/density)
	temp_in = (temp2 (j)/temperature)
	abar_in = abar2 (j)
	zbar_in = min(ye2(j) * abar2(j), zbar2 (j))
	ye_in = ye2(j)
	xiso_in(:) = xiso(:,j)

	! Call GetNuQdot if the matter is sufficiently hot and dense
	IF(temp2 (j)> 1.0D0 .and. rho2 (j) > rho2_a) THEN
		CALL private_sneut5(temp_in, rho_in, abar_in, zbar_in, ye_in, xiso_in, nu_qdot (:, j))
		CALL PAIRNEU(rho_in, ye_in, temp_in, chempo2(j), nu_qdot(1, j), nu_pair(:, j))
		CALL PLASNEU(rho_in, ye_in, temp_in, nu_qdot(2, j), nu_plas(:, j))
		nu_qdot (:,j) = nu_qdot (:,j)*qdot
         ENDIF

ENDDO

! Copy the result to ghost grids
DO j = 1, 5
	CALL BOUNDARY1D_NM(nu_qdot(j, :), even)
END DO

end subroutine findneutrinoloss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the total neutrino loss !
! This one is supplemented to GetNuQdot.f90 	!
! Merged by Leung Shing Chi in 2016		!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDTOTALNEUTRINOLOSS
USE DEFINITION
USE FLAME_MODULE
IMPLICIT NONE

! integer variables !
INTEGER :: j

! Time elapse for neutrino intergration !
REAL (DP) :: dtau

! Initilization !
total_nu_qdot = 0.0D0
total_nu_pair = 0.0D0
total_nu_plas = 0.0D0

! Do the integration !
DO j = 1, length_step_part_2
	total_nu_qdot(:) = total_nu_qdot(:) + vol2 (j) * nu_qdot (:, j) * rho2 (j)
ENDDO

! Find the spectra in cgs (erg-1 s-1)!
DO j = 1, length_step_part_2	
	total_nu_pair(:) = total_nu_pair(:) + (vol2 (j)/length**3) * nu_pair(:, j) 
	total_nu_plas(:) = total_nu_plas(:) + (vol2 (j)/length**3) * nu_plas(:, j)
ENDDO

! Integrate the total spectra !
If(global_time > 0.0E0_DP) THEN
	t_neutrino = global_time
	dtau = t_neutrino - t_neutrino_old
	total_loss(:) = total_loss(:) + 0.5D0*(total_nu_qdot(:)+nu_qdot_old(:))*dtau
	total_pair(:) = total_pair(:) + 0.5D0*(total_nu_pair(:)+nu_pair_old(:))*(dtau/time)
	total_plas(:) = total_plas(:) + 0.5D0*(total_nu_plas(:)+nu_plas_old(:))*(dtau/time)
END IF

! Back up !
nu_qdot_old = total_nu_qdot
nu_pair_old = total_nu_pair
nu_plas_old = total_nu_plas
t_neutrino_old = global_time

END SUBROUTINE FindTotalNeutrinoLoss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do the output of  the isotope               !
! composition where we expect it should be changing with time !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OUTPUTNEU(time)
USE DEFINITION

! Input !
REAL (DP), INTENT(IN) :: time

! Integer !
INTEGER :: n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Done in Main.f90
! Find neutrino loss
!CALL findneutrinoloss
!CALL FindTotalNeutrinoLoss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Output !
WRITE (801, 702) time, (total_nu_qdot(n), n = 1,5)
WRITE (802, 702) time, (total_nu_pair(n), n = 1,5)
WRITE (803, 702) time, (total_nu_plas(n), n = 1,5)
702 FORMAT (F33.15, 5ES33.15)

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find pair neutrino emisstivity and spectrum !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PAIRNEU(rho_in, ye_in, temp_in, beta, emiss, spec)
USE DEFINITION, only : DP
INCLUDE 'implno.dek'
INCLUDE "const.dek"

! Input density, temperature and electron fractions !
REAL (DP), INTENT(IN) :: rho_in, ye_in, temp_in, beta

! Output emisstivity !
REAL (DP), INTENT(OUT), DIMENSION(1:5) :: spec

! Output emisstivity !
REAL (DP), INTENT(OUT) :: emiss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Arrays to store the fermi-dirac functions !
REAL (DP), DIMENSION (-1 : 3) :: gplus, gminus

! Moments funtion !
REAL (DP) :: m10m, m10p, m00

! integration variables for fermi integrals !
REAL (DP) :: alpha, dummy

! Neutrino Energy !
REAL (DP) :: neue

! Integer !
INTEGER :: j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the relativity parameter !
alpha = (me*clight**2)/(kerg*temp_in)

! Find the fermi dirac function !
CALL FERMIDIRAC(alpha, beta, gplus, gminus) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign individual moments !
m00 = (7.0D0*cv2 - 2.0D0*ca2)*gplus(-1)*gminus(-1) + 9.0D0*cv2*gplus(0)*gminus(0) + (cv2 + ca2) * &
	(4.0D0*gplus(1)*gminus(1)-gplus(-1)*gminus(1)-gplus(1)*gminus(-1))

m10m = (7.0D0*cv2 - 2.0D0*ca2)*gplus(0)*gminus(-1) + 9.0D0*cv2*gplus(1)*gminus(0) + (cv2 + ca2) * &
	(4.0D0*gplus(2)*gminus(1)-gplus(0)*gminus(1)-gplus(2)*gminus(-1))

m10p = (7.0D0*cv2 - 2.0D0*ca2)*gplus(-1)*gminus(0) + 9.0D0*cv2*gplus(0)*gminus(1) + (cv2 + ca2) * &
	(4.0D0*gplus(1)*gminus(2)-gplus(-1)*gminus(2)-gplus(1)*gminus(0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign emisstivity in cgs and change to (erg g-1 s-1) !
emiss = (m10m+m10p)*pair1/rho_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the spectra in cgs (erg-1 cm-3 s-1) !
DO j = 1, 5
	! Convert neutrino energy from MeV to erg !
	neue = DBLE(j)*(1.0D6)*ev2erg
	spec(j) = 2.0D0*(pair0*m00)*(bigA/(kerg*temp_in))*(neue/(kerg*temp_in))**(alpha0)*exp(-a0*neue/(kerg*temp_in))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Plasmon Neutrino Emissitivity Using Analytic Fit From Kantor 2007 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PLASNEU(rho_in, ye_in, temp_in, sneup, spec)
USE DEFINITION, only : DP
INCLUDE 'implno.dek'
INCLUDE "const.dek"

! Input density, temperature and electron fractions !
REAL (DP), INTENT(IN) :: rho_in, ye_in, temp_in

! Output emisstivity !
REAL (DP), INTENT(OUT) :: sneup

! Output emisstivity !
REAL (DP), INTENT(OUT), DIMENSION(1:5) :: spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Variables !
Integer :: j

! photon mass !
REAL (DP) :: mt

! Dimensionless fermi momentum !
REAL (DP) :: pfbar, vfbar

! characteristic velocity !
REAL (DP) :: vstar

! dimensionless constant t and f!
REAL (DP) :: t, f

! Plasama frequency !
REAL (DP) :: wp2, wp

! Asymptote !
REAL (DP) :: asy1, asy2, asyt1, asyt2, asyl1, asyl2

! fitting fuction C, D !
REAL (DP) :: funC, funD

! contribution from transverse and longitudinal !
REAL (DP) :: wt, wl

! effective mass density !
REAL (DP) :: rhobar

! beta for fitting !
REAL (DP) :: beta6, beta6s

! emisstivity !
REAL (DP) :: sneut, sneul

! spectra !
REAL (DP) :: spect, specl, integration

! Intermediate !
REAL (DP) :: xin, ain, kappa, din, ein

! Paper input !
REAL (DP) :: kt, neue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First assign density in gram per cc, and temperature !
rhobar = rho_in * ye_in

! Find dimensionless fermi momentum and velocity !
pfbar = (hbar/me/clight) * (3.0D0*pi**2*rhobar/amu) ** (1.0D0/3.0D0)
vfbar = pfbar/sqrt(1.0D0 + pfbar**2)

! dimensionless t !
t = kerg*temp_in/(me*clight**2)

! characteristic velocity !
vstar = ((vfbar**3 + s_plas(1)*(t**s_plas(2))*(rhobar**s_plas(3))) & 
	/(1.0D0 + s_plas(1)*(t**s_plas(2))*(rhobar**s_plas(3)))) ** (1.0D0/3.0D0)

! assign beta !
beta6 = (3.0D0/(2.0D0*vstar**2)*(1.0D0 - (1.0D0 - vstar**2)/(2.0D0*vstar) * log((1.0D0 + vstar)/(1.0D0 - vstar))))**3
beta6s = beta6 + (3.375D0 - beta6)*((t**r_plas(2))*(rhobar**r_plas(3)))/(r_plas(1) + (t**r_plas(2))*(rhobar**r_plas(3)))

! Now assign fitting coefficents !
a_plas(1) = 4.0D0 * beta6s * 1.202056901178332D0 ! Zeta(3) !
a_plas(2) = (8.0D0/105.0D0) + (0.349D0 - (8.0D0/105.0D0))*vstar**10
b_plas(1) = sqrt(2.0D0*pi)*(1.0D0+vstar**2/5.0D0)**(-1.5D0)
b_plas(2) = sqrt(pi/2.0D0)*(3.0D0*vstar**2/5.0D0)**(-1.5D0)
c_plas(1) = p_plas(4)*(1.0D0 + p_plas(5)*(rhobar**p_plas(6)))/(1.0D0 + p_plas(7)*(1.0D0 + p_plas(5)**(rhobar**p_plas(6))))
c_plas(2) = p_plas(8) + p_plas(9)*(rhobar)/(p_plas(10) + rhobar)
d_plas(1) = (6.0D0/pi**2)*(pfbar**2 * (1.0D0 + pfbar**2))/(2.0D0*pfbar**2 + 5.0D0)
d_plas(2) = (pi**2/6.0D0)/(sqrt(1.0D0 + pfbar**2) - 1.0D0)

! Now find the plasma frequency !
funD = t**2/(d_plas(1)*sqrt(1.0D0 + (t*d_plas(2))**2))
funC = 1.0D0 - c_plas(2)*((c_plas(1)*t)**2/(1.0D0 + (c_plas(1)*t)**2))

! asymptote function !
asy1 = (4.0D0*fine)/(3.0D0*pi)*pfbar**3/sqrt(1.0D0 + pfbar**2)
asy2 = (4.0D0*pi*fine)/(9.0D0)*p_plas(2)*(t**2/p_plas(2) + 1.0D0 + p_plas(2)/t**2) & 
	*(1.0D0 + p_plas(3)/(t/sqrt(p_plas(2)))**p_plas(1))**(-10)
wp2 = (me*clight**2/hbar)**2 * sqrt(asy2**2 + (asy1*(1.0D0 - funC*funD))**2)
wp = hbar*sqrt(wp2)

! dimensionless f !
f = hbar*sqrt(wp2)/(kerg*temp_in)

! Now for other asymtpope !
asyt1 = a_plas(1) * f**6
asyt2 = b_plas(1) * f**7.5
asyl1 = a_plas(2) * f**8
asyl2 = b_plas(2) * f**7.5

! Form the contributions from transverse neutrino !
wt = asyt1 + asyt2 * exp(q_plas(3)/((f**q_plas(1))+q_plas(2)))

! From  longitudinal !
wl = asyl2*(asyl1 + q_plas(4)*(1.0D0 + q_plas(5) * vstar**2.5)**3.5 * f**9)& 
	/(asyl2 + (asyl1 + q_plas(4)*(1.0D0 + q_plas(5) * vstar**2.5)**3.5 * f**9))

! Calculate the emissitivity !
sneup = plas0 * t**9 * (wt + wl) * exp(-f)

! Change to erg/g/s !
sneup = sneup/rho_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Photon mass in MeV !
mt = sqrt(1.5D0*(wp/vstar)**2*(1.0D0 - 0.5D0*(1.0D0 - vstar**2)*log((1.0D0 + vstar)/(1.0D0 - vstar))/vstar))
mt = mt/ev2erg*1.0D-6

! kt in MeV !
kt = (kerg*temp_in)/ev2erg*1.0D-6

! Plasma frequency in MeV !
wp = wp/ev2erg*1.0D-6

! Do the loop !
DO j = 1, 5

	! Neutrino energy in MeV !
	neue = DBLE(j)

	! Intemediate !
	xin = neue/mt
	ain = mt/kt
	din = neue/wp
	ein = wp/kt
	kappa = 2.0D0*xin + 0.5D0/xin

	! Assign !
	spect = plas2*mt**7*(kappa - 2.0D0*log(exp(0.5D0*ain*kappa) - 1.0D0)/ain)
	specl = plas3*wp**7*shape(din)/(exp(ein) - 1.0D0)

	! Sum them up !
	spec(j) = spect + specl

	! convert to erg-1 !
	spec(j) = spec(j)/1.60217733E-6

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	REAL (DP) function shape(x)
	REAL (DP), INTENT(IN) :: x
	IF(x > 1.0D0) THEN
		shape = 0.0D0
	ELSEIF(x == 0.5D0) THEN
		shape = 1.0D0
	ELSE
		shape = 4.0D0*x*(x - 1.0D0)*(8.0D0*x**4 - 16.0D0*x**3 + 2.0D0*x**2 + 6.0D0*x - 3.0D0) & 
		+ 3.0D0*(1.0D0 - 2.0D0*x)**2*log((1.0D0 - 2.0D0*x)**2)
	END IF
	shape = shape * 3.28125D0
	END FUNCTION

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the integral to find the spectrum shape !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INTEGRAL (d1, d2, int_out)
IMPLICIT NONE

! Integer !
INTEGER :: j

! Define double precision !
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND (15, 307)

! Input and output !
REAL (DP), INTENT(IN) :: d1, d2
REAL (DP), INTENT(OUT) :: int_out

! number of spatial points !
INTEGER, PARAMETER :: n = 1000

! Spatial seperation !
REAL (DP), PARAMETER :: dh = 2.0D0/DBLE(n)

! Temporal integrand !
REAL (DP) :: xi0, xi1, xi2

! Temporal input !
REAL (DP) :: x_in, y_in
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
x_in = d1/d2
xi0 = integrand(-1.0D0,x_in)/(exp(d1+d2*0.25/x_in)-1.0D0)
xi1 = 0.0D0
xi2 = 0.0D0

! Do the integration !
DO j = 1, n-1
	y_in = -1.0D0 + DBLE(j)*dh
	If(MOD(j,2) == 0.0D0) THEN
		xi2 = xi2 + integrand(y_in,x_in)/(exp(d1+d2*0.5/(1.0D0-y_in)/x_in)-1.0D0)
	ELSE
		xi1 = xi1 + integrand(y_in,x_in)/(exp(d1+d2*0.5/(1.0D0-y_in)/x_in)-1.0D0)
	END IF
END DO

! Assign output !
int_out = dh*(xi0 + 2.0D0*xi2 + 4.0D0*xi1)/3.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	REAL (DP) function integrand(y,x)
	REAL (DP), INTENT(IN) :: x,y
	integrand = (1.0D0+2.0D0*(y-1.0D0)**2*(2*x**2-1.0D0)*x**2)/&
		(x*(y-1.0D0)**2*(1.0D0-2.0D0*y*(y-1.0D0)*x**2+2.0D0*(y-1.0D0)**2*x**4))
	END FUNCTION

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the energ6y loss rate by neutrino    !
! based on the analytic fits by Itoh			     !
! Last Modified by Leung Shing Chi in 2016		     !
! This subroutine is written by Timmes, for more information !
! refer to the cited paper and his webpage		     !
! Muted all unnecssary calculation of derivative 	     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine private_sneut5(temp,den,abar,zbar,ye,xmass,snout)
USE NUCLEAR_MODULE, only : aion, zion, total_ion, xthres
      include 'implno.dek'
      include 'const.dek'

! this routine computes thermal neutrino losses from the analytic fits of
! itoh et al. apjs 102, 411, 1996 and also returns their derivatives.

! input:
! temp = temperature
! den  = density
! abar = mean atomic weight
! zbar = mean charge

! output:
! snu    = total neutrino loss rate in erg/g/sec
! dsnudt = derivative of snu with temperature
! dsnudd = derivative of snu with density
! dsnuda = derivative of snu with abar
! dsnudz = derivative of snu with zbar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! declare the pass
      double precision, INTENT(IN) :: temp,den,abar,zbar,ye

! declare the pass
      double precision, INTENT(IN), DIMENSION(1:total_ion) :: xmass

! declare the pass
      double precision, INTENT(OUT),  DIMENSION(1:5) :: snout

! declare the pass
      double precision snu,dsnudt,dsnudd,dsnuda,dsnudz

! integer 
      integer :: j, k

! local variables
      double precision spair,spairdt,spairdd,spairda,spairdz, & 
                       splas, splasdt,splasdd,splasda,splasdz, &
                       sphot,sphotdt,sphotdd,sphotda,sphotdz, & 
                       sbrem,sbremdt,sbremdd,sbremda,sbremdz, & 
                       sreco,srecodt,srecodd,srecoda,srecodz 

      double precision t9,xl,xldt,xlp5,xl2,xl3,xl4,xl5,xl6,xl7,xl8,xl9, &
                       xlmp5,xlm1,xlm2,xlm3,xlm4,xlnt,cc,den6,tfermi, &
                       a0,a1,a2,a3,b1,b2,c00,c01,c02,c03,c04,c05,c06, &
                       c10,c11,c12,c13,c14,c15,c16,c20,c21,c22,c23,c24, &
                       c25,c26,dd00,dd01,dd02,dd03,dd04,dd05,dd11,dd12, &
                       dd13,dd14,dd15,dd21,dd22,dd23,dd24,dd25,b,c,d,f0, &
                       f1,deni,tempi,abari,zbari,f2,f3,z,xmue, & !,ye, &
                       dum,dumdt,dumdd,dumda,dumdz, &
                       gum,gumdt,gumdd,gumda,gumdz


! pair production
      double precision rm,rmdd,rmda,rmdz,rmi,gl,gldt, &
                       zeta,zetadt,zetadd,zetada,zetadz,zeta2,zeta3, &
                       xnum,xnumdt,xnumdd,xnumda,xnumdz, &
                       xden,xdendt,xdendd,xdenda,xdendz, &
                       fpair,fpairdt,fpairdd,fpairda,fpairdz, &
                       qpair,qpairdt,qpairdd,qpairda,qpairdz

! plasma
      double precision gl2,gl2dt,gl2dd,gl2da,gl2dz,gl12,gl32,gl72,gl6, &
                       ft,ftdt,ftdd,ftda,ftdz,fl,fldt,fldd,flda,fldz, &
                       fxy,fxydt,fxydd,fxyda,fxydz

! photo
      double precision tau,taudt,cos1,cos2,cos3,cos4,cos5,sin1,sin2, &
                       sin3,sin4,sin5,last,xast, &
                       fphot,fphotdt,fphotdd,fphotda,fphotdz, &
                       qphot,qphotdt,qphotdd,qphotda,qphotdz

! brem
      double precision t8,t812,t832,t82,t83,t85,t86,t8m1,t8m2,t8m3,t8m5, &
                       t8m6, &
                       eta,etadt,etadd,etada,etadz,etam1,etam2,etam3, &
                       fbrem,fbremdt,fbremdd,fbremda,fbremdz, &
                       gbrem,gbremdt,gbremdd,gbremda,gbremdz, &
                       u,gm1,gm2,gm13,gm23,gm43,gm53,v,w,fb,gt,gb, &
                       fliq,fliqdt,fliqdd,fliqda,fliqdz, &
                       gliq,gliqdt,gliqdd,gliqda,gliqdz

! recomb
      double precision ifermi12,zfermim12,nu,nudt,nudd,nuda,nudz, &
                       nu2,nu3,bigj,bigjdt,bigjdd,bigjda,bigjdz

! My patch !
double precision, DIMENSION(1:total_ion) :: dumb, zneu

! My patch for brem, the gamma are ion dependents !
double precision, DIMENSION(1:5) :: cosm, sinm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! numerical constants
      double precision fac1,fac2,fac3,oneth,twoth,con1,sixth,iln10
      parameter        (fac1   = 5.0d0 * pi / 3.0d0, &
                        fac2   = 10.0d0 * pi, &
                        fac3   = pi / 5.0d0, &
                        oneth  = 1.0d0/3.0d0, &
                        twoth  = 2.0d0/3.0d0, &
                        con1   = 1.0d0/5.9302d0, &
                        sixth  = 1.0d0/6.0d0, &
                        iln10  = 4.342944819032518d-1)


! theta is sin**2(theta_weinberg) = 0.2223 (2019)
! xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
! change theta and xnufam if need be, and the changes will automatically
! propagate through the routine. cv and ca are the vector and axial currents.

      double precision theta,xnufam,cv,ca,cvp,cap,tfac1,tfac2,tfac3, &
                       tfac4,tfac5,tfac6
      parameter        (theta  = 0.2223d0, &
                        xnufam = 3.0d0, &
                        cv     = 0.5d0 + 2.0d0 * theta, &
                        cvp    = 1.0d0 - cv, &
                        ca     = 0.5d0, &
                        cap    = 1.0d0 - ca, &
                        tfac1  = cv*cv + ca*ca + &
                                 (xnufam-1.0d0) * (cvp*cvp+cap*cap), &
                        tfac2  = cv*cv - ca*ca + &
                                 (xnufam-1.0d0) * (cvp*cvp - cap*cap), &
                        tfac3  = tfac2/tfac1, &
                        tfac4  = 0.5d0 * tfac1, &
                        tfac5  = 0.5d0 * tfac2, &
                        tfac6  = cv*cv + 1.5d0*ca*ca + (xnufam - 1.0d0)* &
                                 (cvp*cvp + 1.5d0*cap*cap))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize
      spair   = 0.0d0
      spairdt = 0.0d0
      spairdd = 0.0d0
      spairda = 0.0d0
      spairdz = 0.0d0

      splas   = 0.0d0
      splasdt = 0.0d0
      splasdd = 0.0d0
      splasda = 0.0d0
      splasdz = 0.0d0

      sphot   = 0.0d0
      sphotdt = 0.0d0
      sphotdd = 0.0d0
      sphotda = 0.0d0
      sphotdz = 0.0d0

      sbrem   = 0.0d0
      sbremdt = 0.0d0
      sbremdd = 0.0d0
      sbremda = 0.0d0
      sbremdz = 0.0d0

      sreco   = 0.0d0
      srecodt = 0.0d0
      srecodd = 0.0d0
      srecoda = 0.0d0
      srecodz = 0.0d0

      snu     = 0.0d0
      dsnudt  = 0.0d0
      dsnudd  = 0.0d0
      dsnuda  = 0.0d0
      dsnudz  = 0.0d0

      !if (temp .lt. 1.0e7) return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! to avoid lots of divisions
      deni  = 1.0d0/den
      tempi = 1.0d0/temp
      abari = 1.0d0/abar
      zbari = 1.0d0/zbar


! some composition variables
      ! ye is given as an input !
      !ye    = zbar*abari
      xmue  = abar*zbari




! some frequent factors
      t9     = temp * 1.0d-9
      xl     = t9 * con1
      xldt   = 1.0d-9 * con1
      xlp5   = sqrt(xl)
      xl2    = xl*xl
      xl3    = xl2*xl
      xl4    = xl3*xl
      xl5    = xl4*xl
      xl6    = xl5*xl
      xl7    = xl6*xl
      xl8    = xl7*xl
      xl9    = xl8*xl
      xlmp5  = 1.0d0/xlp5
      xlm1   = 1.0d0/xl
      xlm2   = xlm1*xlm1
      xlm3   = xlm1*xlm2
      xlm4   = xlm1*xlm3

      rm     = den*ye
      rmdd   = ye
      rmda   = -rm*abari
      rmdz   = den*abari
      rmi    = 1.0d0/rm

      a0     = rm * 1.0d-9
      a1     = a0**oneth
      zeta   = a1 * xlm1
      zetadt = -a1 * xlm2 * xldt
      a2     = oneth * a1*rmi * xlm1
      zetadd = a2 * rmdd
      zetada = a2 * rmda
      zetadz = a2 * rmdz

      zeta2 = zeta * zeta
      zeta3 = zeta2 * zeta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! pair neutrino section
! for reactions like e+ + e- => nu_e + nubar_e
!
! equation 2.8
!      gl   = 1.0d0 - 13.04d0*xl2 +133.5d0*xl4 +1534.0d0*xl6 +918.6d0*xl8
!      !gldt = xldt*(-26.08d0*xl +534.0d0*xl3 +9204.0d0*xl5 +7348.8d0*xl7)
!
! equation 2.7
!
!      a1     = 6.002d19 + 2.084d20*zeta + 1.872d21*zeta2
!      a2     = 2.084d20 + 2.0d0*1.872d21*zeta
!
!      if (t9 .lt. 10.0) then
!       b1     = exp(-5.5924d0*zeta)
!       b2     = -b1*5.5924d0
!      else
!       b1     = exp(-4.9924d0*zeta)
!       b2     = -b1*4.9924d0
!      end if
!
!      xnum   = a1 * b1
!      c      = a2*b1 + a1*b2
!      !xnumdt = c*zetadt
!      !xnumdd = c*zetadd
!      !xnumda = c*zetada
!      !xnumdz = c*zetadz
!
!      if (t9 .lt. 10.0) then
!       a1   = 9.383d-1*xlm1 - 4.141d-1*xlm2 + 5.829d-2*xlm3
!       a2   = -9.383d-1*xlm2 + 2.0d0*4.141d-1*xlm3 - 3.0d0*5.829d-2*xlm4
!      else
!       a1   = 1.2383d0*xlm1 - 8.141d-1*xlm2
!       a2   = -1.2383d0*xlm2 + 2.0d0*8.141d-1*xlm3
!      end if
!
!      b1   = 3.0d0*zeta2
!
!      xden   = zeta3 + a1
!      !xdendt = b1*zetadt + a2*xldt
!      !xdendd = b1*zetadd
!      !xdenda = b1*zetada
!      !xdendz = b1*zetadz
!
!      a1      = 1.0d0/xden
!      fpair   = xnum*a1
!      !fpairdt = (xnumdt - fpair*xdendt)*a1
!      !fpairdd = (xnumdd - fpair*xdendd)*a1
!      !fpairda = (xnumda - fpair*xdenda)*a1
!      !fpairdz = (xnumdz - fpair*xdendz)*a1


! equation 2.6
!      a1     = 10.7480d0*xl2 + 0.3967d0*xlp5 + 1.005d0
!      a2     = xldt*(2.0d0*10.7480d0*xl + 0.5d0*0.3967d0*xlmp5)
!      xnum   = 1.0d0/a1
!      !xnumdt = -xnum*xnum*a2
!
!      a1     = 7.692d7*xl3 + 9.715d6*xlp5
!      a2     = xldt*(3.0d0*7.692d7*xl2 + 0.5d0*9.715d6*xlmp5)
!
!      c      = 1.0d0/a1
!      b1     = 1.0d0 + rm*c
!
!      xden   = b1**(-0.3d0)
!
!      d      = -0.3d0*xden/b1
!      !xdendt = -d*rm*c*c*a2
!      !xdendd = d*rmdd*c
!      !xdenda = d*rmda*c
!      !xdendz = d*rmdz*c
!
!      qpair   = xnum*xden
!      !qpairdt = xnumdt*xden + xnum*xdendt
!      !qpairdd = xnum*xdendd
!      !qpairda = xnum*xdenda
!      !qpairdz = xnum*xdendz



! equation 2.5
!      a1    = exp(-2.0d0*xlm1)
!      a2    = a1*2.0d0*xlm2*xldt
!
!      spair   = a1*fpair
!      !spairdt = a2*fpair + a1*fpairdt
!      !spairdd = a1*fpairdd
!      !spairda = a1*fpairda
!      !spairdz = a1*fpairdz
!
!      a1      = spair
!      spair   = gl*a1
!      !spairdt = gl*spairdt + gldt*a1
!      !spairdd = gl*spairdd
!      !spairda = gl*spairda
!      !spairdz = gl*spairdz
!
!      a1      = tfac4*(1.0d0 + tfac3 * qpair)
!      a2      = tfac4*tfac3
!
!      a3      = spair
!      spair   = a1*a3
!      !spairdt = a1*spairdt + a2*qpairdt*a3
!      !spairdd = a1*spairdd + a2*qpairdd*a3
!      !spairda = a1*spairda + a2*qpairda*a3
!      !spairdz = a1*spairdz + a2*qpairdz*a3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! plasma neutrino section
! for collective reactions like gamma_plasmon => nu_e + nubar_e
! equation 4.6
!
!      a1   = 1.019d-6*rm
!      a2   = a1**twoth
!      a3   = twoth*a2/a1
!
!      b1   =  sqrt(1.0d0 + a2)
!      b2   = 1.0d0/b1
!
!      c00  = 1.0d0/(temp*temp*b1)
!
!      gl2   = 1.1095d11 * rm * c00
!
!      !gl2dt = -2.0d0*gl2*tempi
!      !d     = rm*c00*b2*0.5d0*b2*a3*1.019d-6
!      !gl2dd = 1.1095d11 * (rmdd*c00  - d*rmdd)
!      !gl2da = 1.1095d11 * (rmda*c00  - d*rmda)
!      !gl2dz = 1.1095d11 * (rmdz*c00  - d*rmdz)
!
!
!      gl    = sqrt(gl2)
!      gl12  = sqrt(gl)
!      gl32  = gl * gl12
!      gl72  = gl2 * gl32
!      gl6   = gl2 * gl2 * gl2
!
!
! equation 4.7
!      ft   = 2.4d0 + 0.6d0*gl12 + 0.51d0*gl + 1.25d0*gl32
!      !gum  = 1.0d0/gl2
!      !a1   =(0.25d0*0.6d0*gl12 +0.5d0*0.51d0*gl +0.75d0*1.25d0*gl32)*gum
!      !ftdt = a1*gl2dt
!      !ftdd = a1*gl2dd
!      !ftda = a1*gl2da
!      !ftdz = a1*gl2dz
!
!
! equation 4.8
!      a1   = 8.6d0*gl2 + 1.35d0*gl72
!      a2   = 8.6d0 + 1.75d0*1.35d0*gl72*gum

!      b1   = 225.0d0 - 17.0d0*gl + gl2
!      b2   = -0.5d0*17.0d0*gl*gum + 1.0d0
!
!      c    = 1.0d0/b1
!      fl   = a1*c
!
!      !d    = (a2 - fl*b2)*c
!      !fldt = d*gl2dt
!      !fldd = d*gl2dd
!      !flda = d*gl2da
!      !fldz = d*gl2dz


! equation 4.9 and 4.10
!      cc   = log10(2.0d0*rm)
!      xlnt = log10(temp)

!      xnum   = sixth * (17.5d0 + cc - 3.0d0*xlnt)
!      !xnumdt = -iln10*0.5d0*tempi
!      !a2     = iln10*sixth*rmi
!      !xnumdd = a2*rmdd
!      !xnumda = a2*rmda
!      !xnumdz = a2*rmdz

!      xden   = sixth * (-24.5d0 + cc + 3.0d0*xlnt)
!      !xdendt = iln10*0.5d0*tempi
!      !xdendd = a2*rmdd
!      !xdenda = a2*rmda
!      !xdendz = a2*rmdz


! equation 4.11
!      if (abs(xnum) .gt. 0.7d0  .or.  xden .lt. 0.0d0) then
!       fxy   = 1.0d0
!       !fxydt = 0.0d0
!       !fxydd = 0.0d0
!       !fxydz = 0.0d0
!       !fxyda = 0.0d0
!
!      else
!
!       a1  = 0.39d0 - 1.25d0*xnum - 0.35d0*sin(4.5d0*xnum)
!       a2  = -1.25d0 - 4.5d0*0.35d0*cos(4.5d0*xnum)
!
!       b1  = 0.3d0 * exp(-1.0d0*(4.5d0*xnum + 0.9d0)**2)
!       b2  = -b1*2.0d0*(4.5d0*xnum + 0.9d0)*4.5d0
!
!       c   = min(0.0d0, xden - 1.6d0 + 1.25d0*xnum)
!       !if (c .eq. 0.0) then
!        !dumdt = 0.0d0
!        !dumdd = 0.0d0
!        !dumda = 0.0d0
!        !dumdz = 0.0d0
!       !else
!        !dumdt = xdendt + 1.25d0*xnumdt
!        !dumdd = xdendd + 1.25d0*xnumdd
!        !dumda = xdenda + 1.25d0*xnumda
!        !dumdz = xdendz + 1.25d0*xnumdz
!       !end if
!
!       d   = 0.57d0 - 0.25d0*xnum
!       a3  = c/d
!       c00 = exp(-1.0d0*a3**2)
!
!       f1  = -c00*2.0d0*a3/d
!       !c01 = f1*(dumdt + a3*0.25d0*xnumdt)
!       !c02 = f1*(dumdd + a3*0.25d0*xnumdd)
!       !c03 = f1*(dumda + a3*0.25d0*xnumda)
!       !c04 = f1*(dumdz + a3*0.25d0*xnumdz)
!
!       fxy   = 1.05d0 + (a1 - b1)*c00
!       !fxydt = (a2*xnumdt -  b2*xnumdt)*c00 + (a1-b1)*c01
!       !fxydd = (a2*xnumdd -  b2*xnumdd)*c00 + (a1-b1)*c02
!       !fxyda = (a2*xnumda -  b2*xnumda)*c00 + (a1-b1)*c03
!       !fxydz = (a2*xnumdz -  b2*xnumdz)*c00 + (a1-b1)*c04
!
!      end if



! equation 4.1 and 4.5
!      splas   = (ft + fl) * fxy
!      !splasdt = (ftdt + fldt)*fxy + (ft+fl)*fxydt
!      !splasdd = (ftdd + fldd)*fxy + (ft+fl)*fxydd
!      !splasda = (ftda + flda)*fxy + (ft+fl)*fxyda
!      !splasdz = (ftdz + fldz)*fxy + (ft+fl)*fxydz
!
!      a2      = exp(-gl)
!      a3      = -0.5d0*a2*gl*gum
!
!      a1      = splas
!      splas   = a2*a1
!      !splasdt = a2*splasdt + a3*gl2dt*a1
!      !splasdd = a2*splasdd + a3*gl2dd*a1
!      !splasda = a2*splasda + a3*gl2da*a1
!      !splasdz = a2*splasdz + a3*gl2dz*a1
!
!      a2      = gl6
!      a3      = 3.0d0*gl6*gum
!
!      a1      = splas
!      splas   = a2*a1
!      !splasdt = a2*splasdt + a3*gl2dt*a1
!      !splasdd = a2*splasdd + a3*gl2dd*a1
!      !splasda = a2*splasda + a3*gl2da*a1
!      !splasdz = a2*splasdz + a3*gl2dz*a1
!
!
!      a2      = (cv**2 + xnufam * cvp**2) * 3.0d21 * xl9 !0.93153d0 * 3.0d21 * xl9
!      a3      = (cv**2 + xnufam * cvp**2) * 0.93153d0 * 3.0d21 * 9.0d0*xl8*xldt !0.93153d0 * 3.0d21 * 9.0d0*xl8*xldt
!
!      a1      = splas
!      splas   = a2*a1
!      !splasdt = a2*splasdt + a3*a1
!      !splasdd = a2*splasdd
!      !splasda = a2*splasda
!      !splasdz = a2*splasdz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! photoneutrino process section
! for reactions like e- + gamma => e- + nu_e + nubar_e
!                    e+ + gamma => e+ + nu_e + nubar_e
! equation 3.8 for tau, equation 3.6 for cc,
! and table 2 written out for speed
      if (temp .ge. 1.0d7  .and. temp .lt. 1.0d8) then
       tau  =  log10(temp * 1.0d-7)
       cc   =  0.5654d0 + tau
       c00  =  1.008d11
       c01  =  0.0d0
       c02  =  0.0d0
       c03  =  0.0d0
       c04  =  0.0d0
       c05  =  0.0d0
       c06  =  0.0d0
       c10  =  8.156d10
       c11  =  9.728d8
       c12  = -3.806d9
       c13  = -4.384d9
       c14  = -5.774d9
       c15  = -5.249d9
       c16  = -5.153d9
       c20  =  1.067d11
       c21  = -9.782d9
       c22  = -7.193d9
       c23  = -6.936d9
       c24  = -6.893d9
       c25  = -7.041d9
       c26  = -7.193d9
       dd01 =  0.0d0
       dd02 =  0.0d0
       dd03 =  0.0d0
       dd04 =  0.0d0
       dd05 =  0.0d0
       dd11 = -1.879d10
       dd12 = -9.667d9
       dd13 = -5.602d9
       dd14 = -3.370d9
       dd15 = -1.825d9
       dd21 = -2.919d10
       dd22 = -1.185d10
       dd23 = -7.270d9
       dd24 = -4.222d9
       dd25 = -1.560d9

      else if (temp .ge. 1.0d8  .and. temp .lt. 1.0d9) then
       tau   =  log10(temp * 1.0d-8)
       cc   =  1.5654d0
       c00  =  9.889d10
       c01  = -4.524d8
       c02  = -6.088d6
       c03  =  4.269d7
       c04  =  5.172d7
       c05  =  4.910d7
       c06  =  4.388d7
       c10  =  1.813d11
       c11  = -7.556d9
       c12  = -3.304d9
       c13  = -1.031d9
       c14  = -1.764d9
       c15  = -1.851d9
       c16  = -1.928d9
       c20  =  9.750d10
       c21  =  3.484d10
       c22  =  5.199d9
       c23  = -1.695d9
       c24  = -2.865d9
       c25  = -3.395d9
       c26  = -3.418d9
       dd01 = -1.135d8
       dd02 =  1.256d8
       dd03 =  5.149d7
       dd04 =  3.436d7
       dd05 =  1.005d7
       dd11 =  1.652d9
       dd12 = -3.119d9
       dd13 = -1.839d9
       dd14 = -1.458d9
       dd15 = -8.956d8
       dd21 = -1.548d10
       dd22 = -9.338d9
       dd23 = -5.899d9
       dd24 = -3.035d9
       dd25 = -1.598d9

      else if (temp .ge. 1.0d9) then
       tau  =  log10(t9)
       cc   =  1.5654d0
       c00  =  9.581d10
       c01  =  4.107d8
       c02  =  2.305d8
       c03  =  2.236d8
       c04  =  1.580d8
       c05  =  2.165d8
       c06  =  1.721d8
       c10  =  1.459d12
       c11  =  1.314d11
       c12  = -1.169d11
       c13  = -1.765d11
       c14  = -1.867d11
       c15  = -1.983d11
       c16  = -1.896d11
       c20  =  2.424d11
       c21  = -3.669d9
       c22  = -8.691d9
       c23  = -7.967d9
       c24  = -7.932d9
       c25  = -7.987d9
       c26  = -8.333d9
       dd01 =  4.724d8
       dd02 =  2.976d8
       dd03 =  2.242d8
       dd04 =  7.937d7
       dd05 =  4.859d7
       dd11 = -7.094d11
       dd12 = -3.697d11
       dd13 = -2.189d11
       dd14 = -1.273d11
       dd15 = -5.705d10
       dd21 = -2.254d10
       dd22 = -1.551d10
       dd23 = -7.793d9
       dd24 = -4.489d9
       dd25 = -2.185d9
      end if

      taudt = iln10*tempi


! equation 3.7, compute the expensive trig functions only one time
      cos1 = cos(fac1*tau)
      cos2 = cos(fac1*2.0d0*tau)
      cos3 = cos(fac1*3.0d0*tau)
      cos4 = cos(fac1*4.0d0*tau)
      cos5 = cos(fac1*5.0d0*tau)
      last = cos(fac2*tau)

      sin1 = sin(fac1*tau)
      sin2 = sin(fac1*2.0d0*tau)
      sin3 = sin(fac1*3.0d0*tau)
      sin4 = sin(fac1*4.0d0*tau)
      sin5 = sin(fac1*5.0d0*tau)
      xast = sin(fac2*tau)

      a0 = 0.5d0*c00 &
           + c01*cos1 + dd01*sin1 + c02*cos2 + dd02*sin2 &
           + c03*cos3 + dd03*sin3 + c04*cos4 + dd04*sin4 &
           + c05*cos5 + dd05*sin5 + 0.5d0*c06*last

      f0 =  taudt*fac1*(-c01*sin1 + dd01*cos1 - c02*sin2*2.0d0 &
           + dd02*cos2*2.0d0 - c03*sin3*3.0d0 + dd03*cos3*3.0d0 &
           - c04*sin4*4.0d0 + dd04*cos4*4.0d0 &
           - c05*sin5*5.0d0 + dd05*cos5*5.0d0) &
           - 0.5d0*c06*xast*fac2*taudt

      a1 = 0.5d0*c10 &
           + c11*cos1 + dd11*sin1 + c12*cos2 + dd12*sin2 &
           + c13*cos3 + dd13*sin3 + c14*cos4 + dd14*sin4 &
           + c15*cos5 + dd15*sin5 + 0.5d0*c16*last

      f1 = taudt*fac1*(-c11*sin1 + dd11*cos1 - c12*sin2*2.0d0 &
           + dd12*cos2*2.0d0 - c13*sin3*3.0d0 + dd13*cos3*3.0d0 &
           - c14*sin4*4.0d0 + dd14*cos4*4.0d0 - c15*sin5*5.0d0 &
           + dd15*cos5*5.0d0) - 0.5d0*c16*xast*fac2*taudt

      a2 = 0.5d0*c20 &
           + c21*cos1 + dd21*sin1 + c22*cos2 + dd22*sin2 &
           + c23*cos3 + dd23*sin3 + c24*cos4 + dd24*sin4 &
           + c25*cos5 + dd25*sin5 + 0.5d0*c26*last

      f2 = taudt*fac1*(-c21*sin1 + dd21*cos1 - c22*sin2*2.0d0 &
           + dd22*cos2*2.0d0 - c23*sin3*3.0d0 + dd23*cos3*3.0d0 &
           - c24*sin4*4.0d0 + dd24*cos4*4.0d0 - c25*sin5*5.0d0 &
           + dd25*cos5*5.0d0) - 0.5d0*c26*xast*fac2*taudt

! equation 3.4
      dum   = a0 + a1*zeta + a2*zeta2
      !dumdt = f0 + f1*zeta + a1*zetadt + f2*zeta2 + 2.0d0*a2*zeta*zetadt
      !dumdd = a1*zetadd + 2.0d0*a2*zeta*zetadd
      !dumda = a1*zetada + 2.0d0*a2*zeta*zetada
      !dumdz = a1*zetadz + 2.0d0*a2*zeta*zetadz

      z      = exp(-cc*zeta)

      xnum   = dum*z
      !xnumdt = dumdt*z - dum*z*cc*zetadt
      !xnumdd = dumdd*z - dum*z*cc*zetadd
      !xnumda = dumda*z - dum*z*cc*zetada
      !xnumdz = dumdz*z - dum*z*cc*zetadz

      xden   = zeta3 + 6.290d-3*xlm1 + 7.483d-3*xlm2 + 3.061d-4*xlm3

      dum    = 3.0d0*zeta2
      !xdendt = dum*zetadt - xldt*(6.290d-3*xlm2 &
               !+ 2.0d0*7.483d-3*xlm3 + 3.0d0*3.061d-4*xlm4)
      !xdendd = dum*zetadd
      !xdenda = dum*zetada
      !xdendz = dum*zetadz

      dum      = 1.0d0/xden
      fphot   = xnum*dum
      !fphotdt = (xnumdt - fphot*xdendt)*dum
      !fphotdd = (xnumdd - fphot*xdendd)*dum
      !fphotda = (xnumda - fphot*xdenda)*dum
      !fphotdz = (xnumdz - fphot*xdendz)*dum


! equation 3.3
      a0     = 1.0d0 + 2.045d0 * xl
      xnum   = 0.666d0*a0**(-2.066d0)
      !xnumdt = -2.066d0*xnum/a0 * 2.045d0*xldt

      dum    = 1.875d8*xl + 1.653d8*xl2 + 8.499d8*xl3 - 1.604d8*xl4
      !dumdt  = xldt*(1.875d8 + 2.0d0*1.653d8*xl + 3.0d0*8.449d8*xl2 &
               !- 4.0d0*1.604d8*xl3)

      z      = 1.0d0/dum
      xden   = 1.0d0 + rm*z
      !xdendt =  -rm*z*z*dumdt
      !xdendd =  rmdd*z
      !xdenda =  rmda*z
      !xdendz =  rmdz*z

      z      = 1.0d0/xden
      qphot = xnum*z
      !qphotdt = (xnumdt - qphot*xdendt)*z
      dum      = -qphot*z
      !qphotdd = dum*xdendd
      !qphotda = dum*xdenda
      !qphotdz = dum*xdendz

! equation 3.2
      sphot   = xl5 * fphot
      !sphotdt = 5.0d0*xl4*xldt*fphot + xl5*fphotdt
      !sphotdd = xl5*fphotdd
      !sphotda = xl5*fphotda
      !sphotdz = xl5*fphotdz

      a1      = sphot
      sphot   = rm*a1
      !sphotdt = rm*sphotdt
      !sphotdd = rm*sphotdd + rmdd*a1
      !sphotda = rm*sphotda + rmda*a1
      !sphotdz = rm*sphotdz + rmdz*a1

      a1      = tfac4*(1.0d0 - tfac3 * qphot)
      a2      = -tfac4*tfac3

      a3      = sphot
      sphot   = a1*a3
      !sphotdt = a1*sphotdt + a2*qphotdt*a3
      !sphotdd = a1*sphotdd + a2*qphotdd*a3
      !sphotda = a1*sphotda + a2*qphotda*a3
      !sphotdz = a1*sphotdz + a2*qphotdz*a3

      if (sphot .le. 0.0) then
       sphot   = 0.0d0
       !sphotdt = 0.0d0
       !sphotdd = 0.0d0
       !sphotda = 0.0d0
       !sphotdz = 0.0d0
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! bremsstrahlung neutrino section
! for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
!                    n  + n     => n + n + nu + nubar
!                    n  + p     => n + p + nu + nubar
! equation 4.3

      den6   = den * 1.0d-6
      t8     = temp * 1.0d-8
      t812   = sqrt(t8)
      t832   = t8 * t812
      t82    = t8*t8
      t83    = t82*t8
      t85    = t82*t83
      t86    = t85*t8
      t8m1   = 1.0d0/t8
      t8m2   = t8m1*t8m1
      t8m3   = t8m2*t8m1
      t8m5   = t8m3*t8m2
      t8m6   = t8m5*t8m1


      tfermi = 5.9302d9*(sqrt(1.0d0+1.018d0*(den6*ye)**twoth)-1.0d0)

! "weak" degenerate electrons only
      if (temp .gt. 0.3d0 * tfermi) then

! equation 5.3
       dum   = 7.05d6 * t832 + 5.12d4 * t83
       !dumdt = (1.5d0*7.05d6*t812 + 3.0d0*5.12d4*t82)*1.0d-8

       z     = 1.0d0/dum
       eta   = rm*z
       !etadt = -rm*z*z*dumdt
       !etadd = rmdd*z
       !etada = rmda*z
       !etadz = rmdz*z

       etam1 = 1.0d0/eta
       etam2 = etam1 * etam1
       etam3 = etam2 * etam1


! equation 5.2
       a0    = 23.5d0 + 6.83d4*t8m2 + 7.81d8*t8m5
       f0    = (-2.0d0*6.83d4*t8m3 - 5.0d0*7.81d8*t8m6)*1.0d-8
       xnum  = 1.0d0/a0

       dum   = 1.0d0 + 1.47d0*etam1 + 3.29d-2*etam2
       z     = -1.47d0*etam2 - 2.0d0*3.29d-2*etam3
       !dumdt = z*etadt
       !dumdd = z*etadd
       !dumda = z*etada
       !dumdz = z*etadz

       c00   = 1.26d0 * (1.0d0+etam1)
       z     = -1.26d0*etam2
       c01   = z*etadt
       c02   = z*etadd
       c03   = z*etada
       c04   = z*etadz

       z      = 1.0d0/dum
       xden   = c00*z
       !xdendt = (c01 - xden*dumdt)*z
       !xdendd = (c02 - xden*dumdd)*z
       !xdenda = (c03 - xden*dumda)*z
       !xdendz = (c04 - xden*dumdz)*z

       fbrem   = xnum + xden
       !fbremdt = -xnum*xnum*f0 + xdendt
       !fbremdd = xdendd
       !fbremda = xdenda
       !fbremdz = xdendz


! equation 5.9
       a0    = 230.0d0 + 6.7d5*t8m2 + 7.66d9*t8m5
       f0    = (-2.0d0*6.7d5*t8m3 - 5.0d0*7.66d9*t8m6)*1.0d-8

       z     = 1.0d0 + rm*1.0d-9
       dum   = a0*z
       !dumdt = f0*z
       z     = a0*1.0d-9
       !dumdd = z*rmdd
       !dumda = z*rmda
       !dumdz = z*rmdz

       xnum   = 1.0d0/dum
       z      = -xnum*xnum
       !xnumdt = z*dumdt
       !xnumdd = z*dumdd
       !xnumda = z*dumda
       !xnumdz = z*dumdz

       c00   = 7.75d5*t832 + 247.0d0*t8**(3.85d0)
       dd00  = (1.5d0*7.75d5*t812 + 3.85d0*247.0d0*t8**(2.85d0))*1.0d-8

       c01   = 4.07d0 + 0.0240d0 * t8**(1.4d0)
       dd01  = 1.4d0*0.0240d0*t8**(0.4d0)*1.0d-8

       c02   = 4.59d-5 * t8**(-0.110d0)
       dd02  = -0.11d0*4.59d-5 * t8**(-1.11d0)*1.0d-8

       z     = den**(0.656d0)
       dum   = c00*rmi  + c01  + c02*z
       !dumdt = dd00*rmi + dd01 + dd02*z
       z     = -c00*rmi*rmi
       !dumdd = z*rmdd + 0.656d0*c02*den**(-0.454d0)
       !dumda = z*rmda
       !dumdz = z*rmdz

       xden  = 1.0d0/dum
       z      = -xden*xden
       !xdendt = z*dumdt
       !xdendd = z*dumdd
       !xdenda = z*dumda
       !xdendz = z*dumdz

       gbrem   = xnum + xden
       !gbremdt = xnumdt + xdendt
       !gbremdd = xnumdd + xdendd
       !gbremda = xnumda + xdenda
       !gbremdz = xnumdz + xdendz


! equation 5.1
       dum    = 0.5738d0*zbar*ye*t86*den
       !dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
       !dumdd  = 0.5738d0*zbar*ye*t86
       !dumda  = -dum*abari
       !dumdz  = 0.5738d0*2.0d0*ye*t86*den

       z       = tfac4*fbrem - tfac5*gbrem
       sbrem   = dum * z
       !sbremdt = dumdt*z + dum*(tfac4*fbremdt - tfac5*gbremdt)
       !sbremdd = dumdd*z + dum*(tfac4*fbremdd - tfac5*gbremdd)
       !sbremda = dumda*z + dum*(tfac4*fbremda - tfac5*gbremda)
       !sbremdz = dumdz*z + dum*(tfac4*fbremdz - tfac5*gbremdz)




! liquid metal with c12 parameters (not too different for other elements)
! equation 5.18 and 5.16

      else
       u     = fac3 * (log10(den) - 3.0d0)
       a0    = iln10*fac3*deni

! compute the expensive trig functions of equation 5.21 only once
       !cos1 = cos(u)
       !cos2 = cos(2.0d0*u)
       !cos3 = cos(3.0d0*u)
       !cos4 = cos(4.0d0*u)
       !cos5 = cos(5.0d0*u)

       !sin1 = sin(u)
       !sin2 = sin(2.0d0*u)
       !sin3 = sin(3.0d0*u)
       !sin4 = sin(4.0d0*u)
       !sin5 = sin(5.0d0*u)

! My patch, do loop is better !
DO j = 1, 5 
	cosm(j) = cos(DBLE(j)*u)
END DO
DO j = 1, 5
	sinm(j) = sin(DBLE(j)*u)
END DO

! My patch, calculate the fitting functions for each isotope !
! equation 5.21
DO j = 1, total_ion
	fb =  0.5d0 * aiso(j,0)  + ciso(j) * u + diso(j)
	DO k = 1, 4
		fb = fb + biso(j,k)*sinm(k)
	END DO
	DO k = 1, 5
		fb = fb + aiso(j,k)*cosm(k)
	END DO

	c00 = ciso(j)
	DO k = 1, 4
		c00 = c00 + biso(j,k)*cosm(k)*DBLE(k)
	END DO
		DO k = 1, 5
		c00 = c00 - aiso(j,k)*sinm(k)*DBLE(k)
	END DO
	c00 = c00 * a0

       !fb =  0.5d0 * 0.17946d0  + 0.00945d0*u + 0.34529d0 &
             !- 0.05821d0*cos1 - 0.04969d0*sin1 &
             !- 0.01089d0*cos2 - 0.01584d0*sin2 &
             !- 0.01147d0*cos3 - 0.00504d0*sin3 &
             !- 0.00656d0*cos4 - 0.00281d0*sin4 &
             !- 0.00519d0*cos5

       !c00 =  a0*(0.00945d0 &
             !+ 0.05821d0*sin1       - 0.04969d0*cos1 &
             !+ 0.01089d0*sin2*2.0d0 - 0.01584d0*cos2*2.0d0 &
             !+ 0.01147d0*sin3*3.0d0 - 0.00504d0*cos3*3.0d0 &
             !+ 0.00656d0*sin4*4.0d0 - 0.00281d0*cos4*4.0d0 &
             !+ 0.00519d0*sin5*5.0d0)


! equation 5.22
	ft =  0.5d0 * eiso(j,0)  + giso(j) * u + hiso(j)
	DO k = 1, 4
		ft = ft + fiso(j,k)*sinm(k)
	END DO
	DO k = 1, 5
		ft = ft + eiso(j,k)*cosm(k)
	END DO

	c01 = giso(j)
	DO k = 1, 4
		c01 = c01 + fiso(j,k)*cosm(k)*DBLE(k)
	END DO
		DO k = 1, 5
		c01 = c01 - eiso(j,k)*sinm(k)*DBLE(k)
	END DO
	c01 = c01 * a0

       !ft =  0.5d0 * 0.06781d0 - 0.02342d0*u + 0.24819d0 &
             !- 0.00944d0*cos1 - 0.02213d0*sin1 &
             !- 0.01289d0*cos2 - 0.01136d0*sin2 &
             !- 0.00589d0*cos3 - 0.00467d0*sin3 &
             !- 0.00404d0*cos4 - 0.00131d0*sin4 &
             !- 0.00330d0*cos5

       !c01 = a0*(-0.02342d0 &
             !+ 0.00944d0*sin1       - 0.02213d0*cos1 &
             !+ 0.01289d0*sin2*2.0d0 - 0.01136d0*cos2*2.0d0 &
             !+ 0.00589d0*sin3*3.0d0 - 0.00467d0*cos3*3.0d0 &
             !+ 0.00404d0*sin4*4.0d0 - 0.00131d0*cos4*4.0d0 &
             !+ 0.00330d0*sin5*5.0d0)


! equation 5.23
       	gb =  0.5d0 * iiso(j,0)  + kiso(j) * u + liso(j) !&
	DO k = 1, 4
		gb = gb + jiso(j,k)*sinm(k)
	END DO
	DO k = 1, 5
		gb = gb + iiso(j,k)*cosm(k)
	END DO

     	c02 = kiso(j)
	DO k = 1, 4
		c02 = c02 + jiso(j,k)*cosm(k)*DBLE(k)
	END DO
		DO k = 1, 5
		c02 = c02 - iiso(j,k)*sinm(k)*DBLE(k)
	END DO
	c02 = c02 * a0

       !gb =  0.5d0 * 0.00766d0 - 0.01259d0*u + 0.07917d0 &
             !- 0.00710d0*cos1 + 0.02300d0*sin1 &
             !- 0.00028d0*cos2 - 0.01078d0*sin2 &
             !+ 0.00232d0*cos3 + 0.00118d0*sin3 &
             !+ 0.00044d0*cos4 - 0.00089d0*sin4 &
             !+ 0.00158d0*cos5

       !c02 = a0*(-0.01259d0 &
             !+ 0.00710d0*sin1       + 0.02300d0*cos1 &
             !+ 0.00028d0*sin2*2.0d0 - 0.01078d0*cos2*2.0d0 &
             !- 0.00232d0*sin3*3.0d0 + 0.00118d0*cos3*3.0d0 &
             !- 0.00044d0*sin4*4.0d0 - 0.00089d0*cos4*4.0d0 &
             !- 0.00158d0*sin5*5.0d0)


! equation 5.24
       	gt =  0.5d0 * piso(j,0)  + riso(j) * u + siso(j) !&
	DO k = 1, 4
		gt = gt + qiso(j,k)*sinm(k)
	END DO
	DO k = 1, 5
		gt = gt + piso(j,k)*cosm(k)
	END DO

       	c03 = riso(j)
	DO k = 1, 4
		c03 = c03 + qiso(j,k)*cosm(k)*DBLE(k)
	END DO
		DO k = 1, 5
		c03 = c03 - piso(j,k)*sinm(k)*DBLE(k)
	END DO
	c03 = c03 * a0

       !gt =  -0.5d0 * 0.00769d0  - 0.00829d0*u + 0.05211d0 &
             !+ 0.00356d0*cos1 + 0.01052d0*sin1 &
             !- 0.00184d0*cos2 - 0.00354d0*sin2 &
             !+ 0.00146d0*cos3 - 0.00014d0*sin3 &
             !+ 0.00031d0*cos4 - 0.00018d0*sin4 &
             !+ 0.00069d0*cos5

       !c03 = a0*(-0.00829d0 &
             !- 0.00356d0*sin1       + 0.01052d0*cos1 &
             !+ 0.00184d0*sin2*2.0d0 - 0.00354d0*cos2*2.0d0 &
             !- 0.00146d0*sin3*3.0d0 - 0.00014d0*cos3*3.0d0 &
             !- 0.00031d0*sin4*4.0d0 - 0.00018d0*cos4*4.0d0 &
             !- 0.00069d0*sin5*5.0d0)


       dum   = 2.275d-1 * zbar * zbar*t8m1 * (den6*abari)**oneth
       !dumdt = -dum*tempi
       !dumdd = oneth*dum*deni
       !dumda = -oneth*dum*abari
       !dumdz = 2.0d0*dum*zbari

       gm1   = 1.0d0/dum
       gm2   = gm1*gm1
       gm13  = gm1**oneth
       gm23  = gm13 * gm13
       gm43  = gm13*gm1
       gm53  = gm23*gm1


! equation 5.25 and 5.26
       v  = -0.05483d0 - 0.01946d0*gm13 + 1.86310d0*gm23 - 0.78873d0*gm1
       !a0 = oneth*0.01946d0*gm43 - twoth*1.86310d0*gm53 + 0.78873d0*gm2

       w  = -0.06711d0 + 0.06859d0*gm13 + 1.74360d0*gm23 - 0.74498d0*gm1
       !a1 = -oneth*0.06859d0*gm43 - twoth*1.74360d0*gm53 + 0.74498d0*gm2


! equation 5.19 and 5.20
       fliq   = v*fb + (1.0d0 - v)*ft
       !fliqdt = a0*dumdt*(fb - ft)
       !fliqdd = a0*dumdd*(fb - ft) + v*c00 + (1.0d0 - v)*c01
       !fliqda = a0*dumda*(fb - ft)
       !fliqdz = a0*dumdz*(fb - ft)

       gliq   = w*gb + (1.0d0 - w)*gt
       !gliqdt = a1*dumdt*(gb - gt)
       !gliqdd = a1*dumdd*(gb - gt) + w*c02 + (1.0d0 - w)*c03
       !gliqda = a1*dumda*(gb - gt)
       !gliqdz = a1*dumdz*(gb - gt)


! equation 5.17
	dumb(j) = 0.5738d0*xmass(j)*zion(j)**2/aion(j)*t86*den
       !dum    = 0.5738d0*zbar*ye*t86*den
       !dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
       !dumdd  = 0.5738d0*zbar*ye*t86
       !dumda  = -dum*abari
       !dumdz  = 0.5738d0*2.0d0*ye*t86*den

	zneu(j) = tfac4*fliq - tfac5*gliq
       !z       = tfac4*fliq - tfac5*gliq
       !sbrem   = dum * z
       !sbremdt = dumdt*z + dum*(tfac4*fliqdt - tfac5*gliqdt)
       !sbremdd = dumdd*z + dum*(tfac4*fliqdd - tfac5*gliqdd)
       !sbremda = dumda*z + dum*(tfac4*fliqda - tfac5*gliqda)
       !sbremdz = dumdz*z + dum*(tfac4*fliqdz - tfac5*gliqdz)
END DO

! Sum all the contributions !
DO j = 1, total_ion
	If(xmass(j) <= xthres) THEN
		sbrem = sbrem + 0.0D0
	ELSE
		sbrem   = sbrem + dumb(j) * zneu(j)
	END IF
END DO

      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! recombination neutrino section
! for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
! equation 6.11 solved for nu
      xnum   = 1.10520d8 * den * ye /(temp*sqrt(temp))
      !xnumdt = -1.50d0*xnum*tempi
      !xnumdd = xnum*deni
      !xnumda = -xnum*abari
      !xnumdz = xnum*zbari

! the chemical potential
      nu   = ifermi12(xnum)

! a0 is d(nu)/d(xnum)
      a0 = 1.0d0/(0.5d0*zfermim12(nu))
      !nudt = a0*xnumdt
      !nudd = a0*xnumdd
      !nuda = a0*xnumda
      !nudz = a0*xnumdz

      nu2  = nu * nu
      nu3  = nu2 * nu

! table 12
      if (nu .ge. -20.0  .and. nu .lt. 0.0) then
       a1 = 1.51d-2
       a2 = 2.42d-1
       a3 = 1.21d0
       b  = 3.71d-2
       c  = 9.06e-1
       d  = 9.28d-1
       f1 = 0.0d0
       f2 = 0.0d0
       f3 = 0.0d0
      else if (nu .ge. 0.0  .and. nu .le. 10.0) then
       a1 = 1.23d-2
       a2 = 2.66d-1
       a3 = 1.30d0
       b  = 1.17d-1
       c  = 8.97e-1
       d  = 1.77d-1
       f1 = -1.20d-2
       f2 = 2.29d-2
       f3 = -1.04d-3
      end if


! equation 6.7, 6.13 and 6.14
      if (nu .ge. -20.0  .and.  nu .le. 10.0) then

       zeta   = 1.579d5*zbar*zbar*tempi
       !zetadt = -zeta*tempi
       !zetadd = 0.0d0
       !zetada = 0.0d0
       !zetadz = 2.0d0*zeta*zbari

       c00    = 1.0d0/(1.0d0 + f1*nu + f2*nu2 + f3*nu3)
       c01    = f1 + f2*2.0d0*nu + f3*3.0d0*nu2
       dum    = zeta*c00
       !dumdt  = zetadt*c00 + zeta*c01*nudt
       !dumdd  = zeta*c01*nudd
       !dumda  = zeta*c01*nuda
       !dumdz  = zetadz*c00 + zeta*c01*nudz


       z      = 1.0d0/dum
       dd00   = dum**(-2.25)
       dd01   = dum**(-4.55)
       c00    = a1*z + a2*dd00 + a3*dd01
       c01    = -(a1*z + 2.25*a2*dd00 + 4.55*a3*dd01)*z


       z      = exp(c*nu)
       dd00   = b*z*(1.0d0 + d*dum)
       gum    = 1.0d0 + dd00
       !gumdt  = dd00*c*nudt + b*z*d*dumdt
       !gumdd  = dd00*c*nudd + b*z*d*dumdd
       !gumda  = dd00*c*nuda + b*z*d*dumda
       !gumdz  = dd00*c*nudz + b*z*d*dumdz


       z   = exp(nu)
       a1  = 1.0d0/gum

       bigj   = c00 * z * a1
       !bigjdt = c01*dumdt*z*a1 + c00*z*nudt*a1 - c00*z*a1*a1 * gumdt
       !bigjdd = c01*dumdd*z*a1 + c00*z*nudd*a1 - c00*z*a1*a1 * gumdd
       !bigjda = c01*dumda*z*a1 + c00*z*nuda*a1 - c00*z*a1*a1 * gumda
       !bigjdz = c01*dumdz*z*a1 + c00*z*nudz*a1 - c00*z*a1*a1 * gumdz


! equation 6.5
       z     = exp(zeta + nu)
       dum   = 1.0d0 + z
       a1    = 1.0d0/dum
       a2    = 1.0d0/bigj

       sreco   = tfac6 * 2.649d-18 * ye * zbar**13 * den * bigj*a1
       !srecodt = sreco*(bigjdt*a2 - z*(zetadt + nudt)*a1)
       !srecodd = sreco*(1.0d0*deni + bigjdd*a2 - z*(zetadd + nudd)*a1)
       !srecoda = sreco*(-1.0d0*abari + bigjda*a2 - z*(zetada+nuda)*a1)
       !srecodz = sreco*(14.0d0*zbari + bigjdz*a2 - z*(zetadz+nudz)*a1)

      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! convert from erg/cm^3/s to erg/g/s
! comment these out to duplicate the itoh et al plots
! exclude plasmon neutrino because there are better fits !

      spair   = spair*deni
      !spairdt = spairdt*deni
      !spairdd = spairdd*deni - spair*deni
      !spairda = spairda*deni
      !spairdz = spairdz*deni

      splas   = splas*deni
      !splasdt = splasdt*deni
      !splasdd = splasdd*deni - splas*deni
      !splasda = splasda*deni
      !splasdz = splasdz*deni

      sphot   = sphot*deni
      !sphotdt = sphotdt*deni
      !sphotdd = sphotdd*deni - sphot*deni
      !sphotda = sphotda*deni
      !sphotdz = sphotdz*deni

      sbrem   = sbrem*deni
      !sbremdt = sbremdt*deni
      !sbremdd = sbremdd*deni - sbrem*deni
      !sbremda = sbremda*deni
      !sbremdz = sbremdz*deni

      sreco   = sreco*deni
      !srecodt = srecodt*deni
      !srecodd = srecodd*deni - sreco*deni
      !srecoda = srecoda*deni
      !srecodz = srecodz*deni


! the total neutrino loss rate
! muted because we want individual rates 
      !snu    =  splas + spair + sphot + sbrem + sreco
      !dsnudt =  splasdt + spairdt + sphotdt + sbremdt + srecodt
      !dsnudd =  splasdd + spairdd + sphotdd + sbremdd + srecodd
      !dsnuda =  splasda + spairda + sphotda + sbremda + srecoda
      !dsnudz =  splasdz + spairdz + sphotdz + sbremdz + srecodz

! My patch, pass to the output !
snout (1) = spair
snout (2) = splas
snout (3) = sphot
snout (4) = sbrem
snout (5) = sreco

      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function ifermi12(f)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the inverse
! fermi function of order 1/2 when it is equal to f.
! maximum error is 4.19d-9.   reference: antia apjs 84,101 1993
!
! declare
      integer  :: i,m1,k1,m2,k2
      real*8   :: f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /0.5d0, 4, 3, 6, 5/
      data  (a1(i),i=1,5)/ 1.999266880833d4,   5.702479099336d3, &
                           6.610132843877d2,   3.818838129486d1, &
                           1.0d0/
      data  (b1(i),i=1,4)/ 1.771804140488d4,  -2.014785161019d3, &
                           9.130355392717d1,  -1.670718177489d0/
      data  (a2(i),i=1,7)/-1.277060388085d-2,  7.187946804945d-2, &
                          -4.262314235106d-1,  4.997559426872d-1, &
                          -1.285579118012d0,  -3.930805454272d-1, &
                           1.0d0/
      data  (b2(i),i=1,6)/-9.745794806288d-3,  5.485432756838d-2, &
                          -3.299466243260d-1,  4.077841975923d-1, &
                          -1.145531476975d0,  -6.067091689181d-2/


      if (f .lt. 4.0d0) then
       rn = f + a1(m1)
       do i=m1-1,1,-1
        rn = rn*f + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*f + b1(i)
       enddo
       ifermi12 = log(f * rn/den)

      else
       ff = 1.0d0/f**(1.0d0/(1.0d0 + an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
        rn = rn*ff + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*ff + b2(i)
       enddo
       ifermi12 = rn/(den*ff)
      end if
      return
      end function ifermi12

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8 function zfermim12(x)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the fermi function
! of order -1/2 evaluated at x. maximum error is 1.23d-12.
! reference: antia apjs 84,101 1993
!
! declare
      integer :: i,m1,k1,m2,k2
      real*8  :: x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /-0.5d0, 7, 7, 11, 11/
      data  (a1(i),i=1,8)/ 1.71446374704454d7,    3.88148302324068d7, &
                           3.16743385304962d7,    1.14587609192151d7, &
                           1.83696370756153d6,    1.14980998186874d5, &
                           1.98276889924768d3,    1.0d0/
      data  (b1(i),i=1,8)/ 9.67282587452899d6,    2.87386436731785d7, &
                           3.26070130734158d7,    1.77657027846367d7, &
                           4.81648022267831d6,    6.13709569333207d5, &
                           3.13595854332114d4,    4.35061725080755d2/
      data (a2(i),i=1,12)/-4.46620341924942d-15, -1.58654991146236d-12, &
                          -4.44467627042232d-10, -6.84738791621745d-8, &
                          -6.64932238528105d-6,  -3.69976170193942d-4, &
                          -1.12295393687006d-2,  -1.60926102124442d-1, &
                          -8.52408612877447d-1,  -7.45519953763928d-1, &
                           2.98435207466372d0,    1.0d0/
      data (b2(i),i=1,12)/-2.23310170962369d-15, -7.94193282071464d-13, &
                          -2.22564376956228d-10, -3.43299431079845d-8, &
                          -3.33919612678907d-6,  -1.86432212187088d-4, &
                          -5.69764436880529d-3,  -8.34904593067194d-2, &
                          -4.78770844009440d-1,  -4.99759250374148d-1, &
                           1.86795964993052d0,    4.16485970495288d-1/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermim12 = xx * rn/den
!
      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermim12 = sqrt(x)*rn/den
      end if
      return
      end function zfermim12



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the interface to communicate between the main program and the   !
! Fermi-Dirac function integrator. To be specific, this subroutine assign !
! All the possible order (-1/2, 0, 1/2, 3.2) Fermi-Dirac function	  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FERMIDIRAC(alpha, beta, gplus, gminus) 
IMPLICIT NONE

! Define double precision !
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND (15, 307)

! Input !
REAL (DP), INTENT(IN) :: alpha, beta

! Arrays to store the fermi-dirac functions !
REAL (DP), DIMENSION (-1 : 3), INTENT(OUT) :: gplus, gminus

! Integer !
INTEGER :: j, r

! Temporaily integral !
REAL (DP) :: integral, dummy

! Initialize !
gplus = 0.0D0
gminus = 0.0D0

! Assign values !
DO j = -1, 3
	DO r = 0, j + 1
		CALL dfermi(DBLE(j) - DBLE(r) + 1.5D0, -alpha-beta, 1.0D0/alpha, integral, dummy, dummy, dummy, dummy, dummy)
		gplus (j) = gplus (j) + integral * choose (j+1, r) * alpha**r
		CALL dfermi(DBLE(j) - DBLE(r) + 1.5D0, -alpha+beta, 1.0D0/alpha, integral, dummy, dummy, dummy, dummy, dummy)
		gminus (j)= gminus (j) + integral * choose (j+1, r) * alpha**r
	END DO
	gplus (j) = gplus (j) * sqrt(2.0D0)/alpha**(2.5D0 + DBLE(j))
	gminus (j) = gminus (j) * sqrt(2.0D0)/alpha**(2.5D0 + DBLE(j))
END DO

! Binominal Coefficients !
! Ref : https://rosettacode.org/ !
contains
 
  function factorial (n) result (res)
 
    implicit none
    integer, intent (in) :: n
    integer :: res
    integer :: i
 
    res = product ((/(i, i = 1, n)/))
 
  end function factorial
 
  function choose (n, k) result (res)
 
    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: k
    integer :: res
 
    res = factorial (n) / (factorial (k) * factorial (n - k))
 
  end function choose

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the degeneracy parameter eta = mu/kT where mu !
! is the chemical potential for both electron and positron 	      !
! mutted derivative calculations to speed up the process	      !
! Ref : CoCoCube						      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine dfermi gets the fermi-dirac functions and their derivatives
! routine fdfunc1 forms the integrand of the fermi-dirac functions
! routine fdfunc2 same as fdfunc but with the change of variable z**2=x
! routine dqleg010 does 10 point gauss-legendre integration  9 fig accuracy
! routine dqleg020 does 20 point gauss-legendre integration 14 fig accuracy
! routine dqleg040 does 40 point gauss-legendre integration 18 fig accuracy
! routine dqleg080 does 80 point gauss-legendre integration 32 fig accuracy
! routine dqlag010 does 10 point gauss-laguerre integration  9 fig accuracy
! routine dqlag020 does 20 point gauss-laguerre integration 14 fig accuracy
! routine dqlag040 does 40 point gauss-laguerre integration 18 fig accuracy
! routine dqlag080 does 80 point gauss-laguerre integration 32 fig accuracy


! quad precision 
! routine qfermi gets the fermi-dirac functions and their derivatives
! routine qfdfunc1 forms the integrand of the fermi-dirac functions
! routine qfdfunc2 same as fdfunc but with the change of variable z**2=x
! routine qqleg010 does  10 point gauss-legendre integration  9 fig accuracy
! routine qqleg020 does  20 point gauss-legendre integration 14 fig accuracy
! routine qqleg040 does  40 point gauss-legendre integration 18 fig accuracy
! routine qqleg080 does  80 point gauss-legendre integration 32 fig accuracy
! routine qqleg200 does 200 point gauss-legendre integration 32 fig accuracy
! routine qqlag010 does  10 point gauss-laguerre integration  9 fig accuracy
! routine qqlag020 does  20 point gauss-laguerre integration 14 fig accuracy
! routine qqlag040 does  40 point gauss-laguerre integration 18 fig accuracy
! routine qqlag080 does  80 point gauss-laguerre integration 32 fig accuracy
! routine qqlag200 does 200 point gauss-laguerre integration 32 fig accuracy




      subroutine dfermi(dk,eta,theta,fd,fdeta,fdtheta, &
                        fdeta2,fdtheta2,fdetadtheta)
      include 'implno.dek'

! this routine computes the fermi-dirac integrals f(dk,eta,theta) of
! index dk, with degeneracy parameter eta and relativity parameter theta.

! input is dk the index of the fermi-dirac function,
! eta the degeneracy parameter, and theta the relativity parameter.

! output is fd is computed by applying three 10-point gauss-legendre
! and one 10-point gauss-laguerre rules over four appropriate subintervals.
! fdeta is the derivative of f with respect to eta.
! fdthetha is the derivative of f with theta.
! fdeta2 is the second derivative of f with respect to eta.
! fdtheta2 is the second derivative of f with respect to theta
! fdetadtheta is the mixed second derivative of f with respect to eta and theta.

! reference: j.m. aparicio, apjs 117, 632 1998


! declare
      external :: fdfunc1,fdfunc2
      real*8   :: dk,eta,theta,fd,fdeta,fdtheta, &
                  fdeta2,fdtheta2,fdetadtheta, &
                  d,sg,a1,b1,c1,a2,b2,c2,d2,e2,a3,b3,c3,d3,e3, &
                  eta1,xi,xi2,x1,x2,x3,s1,s2,s3,s12,par(3), &
                  res(4),drde(4),drdt(4),drde2(4),drdt2(4),drdet(4)


!   parameters defining the location of the breakpoints for the
!   subintervals of integration:
      data d   / 3.3609d0 /
      data sg  / 9.1186d-2 /
      data a1  / 6.7774d0 /
      data b1  / 1.1418d0 /
      data c1  / 2.9826d0 /
      data a2  / 3.7601d0 /
      data b2  / 9.3719d-2 /
      data c2  / 2.1063d-2 /
      data d2  / 3.1084d1 /
      data e2  / 1.0056d0 /
      data a3  / 7.5669d0 /
      data b3  / 1.1695d0 /
      data c3  / 7.5416d-1 /
      data d3  / 6.6558d0 /
      data e3  /-1.2819d-1 /


!   integrand parameters:
      par(1)=dk
      par(2)=eta
      par(3)=theta


!   definition of xi:
      eta1=sg*(eta-d)
      if (eta1 .le. 5.0d1) then
        xi=log(1.d0+exp(eta1))/sg
      else
        xi=eta-d
      endif
      xi2=xi*xi

!   definition of the x_i:
      x1=(a1 + b1*xi + c1*xi2)    / (1.d0 +c1*xi)
      x2=(a2 + b2*xi + c2*d2*xi2) / (1.d0 +e2*xi + c2*xi2)
      x3=(a3 + b3*xi + c3*d3*xi2) / (1.d0 +e3*xi + c3*xi2)

!   breakpoints:
      s1=x1-x2
      s2=x1
      s3=x1+x3
      s12=sqrt(s1)

!   quadrature integrations:

! 9 significant figure accuracy
      call dqleg010(fdfunc2, 0.d0,  s12, res(1), drde(1), drdt(1), &
                    drde2(1), drdt2(1), drdet(1), par,3)
      call dqleg010(fdfunc1,   s1,   s2, res(2), drde(2), drdt(2), &
                    drde2(2), drdt2(2), drdet(2), par,3)
      call dqleg010(fdfunc1,   s2,   s3, res(3), drde(3), drdt(3), &
                    drde2(3), drdt2(3), drdet(3), par,3)
      call dqlag010(fdfunc1,   s3, 1.d0, res(4), drde(4), drdt(4), &
                    drde2(4), drdt2(4), drdet(4), par,3)

! 14 significant figure accuracy
!      call dqleg020(fdfunc2, 0.d0,  s12, res(1), drde(1), drdt(1), &
!                    drde2(1), drdt2(1), drdet(1), par,3)
!      call dqleg020(fdfunc1,   s1,   s2, res(2), drde(2), drdt(2), &
!                    drde2(2), drdt2(2), drdet(2), par,3)
!      call dqleg020(fdfunc1,   s2,   s3, res(3), drde(3), drdt(3), &
!                    drde2(3), drdt2(3), drdet(3), par,3)
!      call dqlag020(fdfunc1,   s3, 1.d0, res(4), drde(4), drdt(4), &
!                    drde2(4), drdt2(4), drdet(4), par,3)

! 18 significant figure accuracy
!      call dqleg040(fdfunc2, 0.d0,  s12, res(1), drde(1), drdt(1), &
!                    drde2(1), drdt2(1), drdet(1), par,3)
!      call dqleg040(fdfunc1,   s1,   s2, res(2), drde(2), drdt(2), &
!                    drde2(2), drdt2(2), drdet(2), par,3)
!      call dqleg040(fdfunc1,   s2,   s3, res(3), drde(3), drdt(3), &
!                    drde2(3), drdt2(3), drdet(3), par,3)
!      call dqlag040(fdfunc1,   s3, 1.d0, res(4), drde(4), drdt(4), &
!                    drde2(4), drdt2(4), drdet(4), par,3)

! 32 significant figure accuracy
!      call dqleg080(fdfunc2, 0.d0,  s12, res(1), drde(1), drdt(1), &
!                    drde2(1), drdt2(1), drdet(1), par,3)
!      call dqleg080(fdfunc1,   s1,   s2, res(2), drde(2), drdt(2), &
!                    drde2(2), drdt2(2), drdet(1), par,3)
!      call dqleg080(fdfunc1,   s2,   s3, res(3), drde(3), drdt(3), &
!                    drde2(3), drdt2(3), drdet(1), par,3)
!      call dqlag080(fdfunc1,   s3, 1.d0, res(4), drde(4), drdt(4), &
!                    drde2(4), drdt2(4), drdet(1), par,3)


! sum the contributions
      fd          = res(1)   + res(2)   + res(3)   + res(4)
      fdeta       = drde(1)  + drde(2)  + drde(3)  + drde(4)
      fdtheta     = drdt(1)  + drdt(2)  + drdt(3)  + drdt(4)
      fdeta2      = drde2(1) + drde2(2) + drde2(3) + drde2(4)
      fdtheta2    = drdt2(1) + drdt2(2) + drdt2(3) + drdt2(4)
      fdetadtheta = drdet(1) + drdet(2) + drdet(3) + drdet(4)
      return
      end subroutine dfermi




      subroutine fdfunc1(x,par,n,fd,fdeta,fdtheta, &
                         fdeta2,fdtheta2,fdetadtheta)
      include 'implno.dek'

! forms the fermi-dirac integrand and its first and second
! derivatives with eta and theta.

! input:
! x is the integration variable
! par(1) is the index
! par(2) is the degeneravy parameter
! par(3) is the relativity parameter

! output:
! fd is the integrand
! fdeta is the first derivative with eta
! fdeta2 is the second derivative with eta
! fdtheta is the first derivative with theta
! fdtheta2 is the second derivative with theta
! fdetadtheta is the mixed second derivative

! declare the pass
      integer     :: n
      real*8      :: x,par(n),fd, &
                     fdeta,fdeta2,fdtheta,fdtheta2,fdetadtheta

! local variables
      real*8      :: dk,eta,theta, &
                     factor,xst,dxst,denom,denomi,denom2,xdk


! initialize
      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xdk   = x**dk
      xst   = 1.0d0 + 0.5d0*x*theta
      dxst  = sqrt(xst)

!   avoid overflow in the exponentials at large x
      if ((x-eta) .lt. 1.0d2) then
       factor      = exp(x-eta)
       denom       = factor + 1.0d0
       denomi      = 1.0d0/denom
       fd          = xdk * dxst * denomi
       fdeta       = fd * factor * denomi
       fdeta2      = (2.0d0 * factor * denomi - 1.0d0)*fdeta
       denom2      = 1.0d0/(4.0d0 * xst)
       fdtheta     = fd * x * denom2
       fdtheta2    = -fdtheta * x * denom2
       fdetadtheta = fdtheta * factor * denomi
      else
       factor      = exp(eta-x)
       fd          = xdk * dxst * factor
       fdeta       = fd
       fdeta2      = fd
       denom2      = 1.0d0/(4.0d0 * xst)
       fdtheta     = fd * x * denom2
       fdtheta2    = -fdtheta * x * denom2
       fdetadtheta = fdtheta
      endif

      return
      end subroutine fdfunc1




      subroutine fdfunc2(x,par,n,fd,fdeta,fdtheta, &
                         fdeta2,fdtheta2,fdetadtheta)
      include 'implno.dek'

! forms the fermi-dirac integrand and its first and second
! derivatives with eta and theta when the variable change z**2=x
! has been made.

! input:
! x is the integration variable
! par(1) is the index
! par(2) is the degeneravy parameter
! par(3) is the relativity parameter

! output:
! fd is the integrand
! fdeta is the first derivative with eta
! fdeta2 is the second derivative with eta
! fdtheta is the first derivative with theta
! fdtheta2 is the second derivative with theta
! fdetadtheta is the mixed second derivative

! declare the pass
      integer     ::  n
      real*8      ::  x,par(n),fd, &
                      fdeta,fdeta2,fdtheta,fdtheta2,fdetadtheta

! local variables
      real*8      ::  dk,eta,theta, &
                      factor,xst,dxst,denom,denomi,denom2,xdk,xsq


      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xsq   = x * x
      xdk   = x**(2.0d0 * dk + 1.0d0)
      xst   = 1.0d0 + 0.5d0 * xsq * theta
      dxst  = sqrt(xst)

!   avoid an overflow in the denominator at large x:
      if ((xsq-eta) .lt. 1.0d2) then
       factor      = exp(xsq - eta)
       denom       = factor + 1.0d0
       denomi      = 1.0d0/denom
       fd          = 2.0d0 * xdk * dxst * denomi
       fdeta       = fd * factor * denomi
       fdeta2      = (2.0d0 * factor * denomi - 1.0d0)*fdeta
       denom2      = 1.0d0/(4.0d0 * xst)
       fdtheta     = fd * xsq * denom2
       fdtheta2    = -fdtheta * xsq * denom2
       fdetadtheta = fdtheta * factor * denomi
      else
       factor      = exp(eta - xsq)
       fd          = 2.0d0 * xdk * dxst * factor
       fdeta       = fd
       fdeta2      = fd
       denom2      = 1.0d0/(4.0d0 * xst)
       fdtheta     = fd * xsq * denom2
       fdtheta2    = -fdtheta * xsq * denom2
       fdetadtheta = fdtheta
      endif

      return
      end subroutine fdfunc2





      subroutine dqleg010(f,a,b,result,drdeta,drdtheta, &
                          drdeta2,drdtheta2,drdetadtheta,par,n)
      include 'implno.dek'
!
! 10 point gauss-legendre rule for the fermi-dirac function and
! its derivatives with respect to eta and theta.
! on input f is the name of the subroutine containing the integrand,
! a is the lower end point of the interval, b is the higher end point,
! par is an array of constant parameters to be passed to subroutine f,
! and n is the length of the par array. on output result is the
! approximation from applying the 10-point gauss-legendre rule,
! drdeta is the derivative with respect to eta, and drdtheta is the
! derivative with respect to theta.
!
! note: since the number of nodes is even, zero is not an abscissa.
!
! declare
      external    ::  f
      integer     ::  j,n
      real*8      ::  a,b,result,drdeta,drdtheta, &
                       drdeta2,drdtheta2,drdetadtheta,par(n), &
                       absc1,absc2,center,hlfrun,wg(5),xg(5), &
                       fval(2),deta(2),dtheta(2), &
                       deta2(2),dtheta2(2),detadtheta(2)

! the abscissae and weights are given for the interval (-1,1).
! xg     - abscissae of the 20-point gauss-legendre rule
!          for half of the usual run (-1,1), i.e.
!          the positive nodes of the 20-point rule
! wg     - weights of the 20-point gauss rule.
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.48874338981631210884826001129719984d-1 /
      data xg (  2) /   4.33395394129247190799265943165784162d-1 /
      data xg (  3) /   6.79409568299024406234327365114873575d-1 /
      data xg (  4) /   8.65063366688984510732096688423493048d-1 /
      data xg (  5) /   9.73906528517171720077964012084452053d-1 /

      data wg (  1) /   2.95524224714752870173892994651338329d-1 /
      data wg (  2) /   2.69266719309996355091226921569469352d-1 /
      data wg (  3) /   2.19086362515982043995534934228163192d-1 /
      data wg (  4) /   1.49451349150580593145776339657697332d-1 /
      data wg (  5) /   6.66713443086881375935688098933317928d-2 /


!           list of major variables
!           -----------------------
!
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 10-point gauss formula

      center       = 0.5d0 * (a+b)
      hlfrun       = 0.5d0 * (b-a)
      result       = 0.0d0
      drdeta       = 0.0d0
      drdtheta     = 0.0d0
      drdeta2      = 0.0d0
      drdtheta2    = 0.0d0
      drdetadtheta = 0.0d0
      do j=1,5
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)

        call f(absc1, par, n, fval(1), deta(1), dtheta(1), &
               deta2(1), dtheta2(1), detadtheta(1))
        call f(absc2, par, n, fval(2), deta(2), dtheta(2), &
               deta2(2), dtheta2(2), detadtheta(2))

        result       = result    + (fval(1)    + fval(2))*wg(j)
        drdeta       = drdeta    + (deta(1)    + deta(2))*wg(j)
        drdtheta     = drdtheta  + (dtheta(1)  + dtheta(2))*wg(j)
        drdeta2      = drdeta2   + (deta2(1)   + deta2(2))*wg(j)
        drdtheta2    = drdtheta2 + (dtheta2(1) + dtheta2(2))*wg(j)
        drdetadtheta = drdetadtheta+(detadtheta(1)+detadtheta(2))*wg(j)
      enddo

      result       = result * hlfrun
      drdeta       = drdeta * hlfrun
      drdtheta     = drdtheta * hlfrun
      drdeta2      = drdeta2 * hlfrun
      drdtheta2    = drdtheta2 * hlfrun
      drdetadtheta = drdetadtheta * hlfrun
      return
      end subroutine dqleg010




      subroutine dqleg020(f,a,b,result,drdeta,drdtheta, &
                          drdeta2,drdtheta2,drdetadtheta,par,n)
      include 'implno.dek'
!
! 20 point gauss-legendre rule for the fermi-dirac function and
! its derivatives with respect to eta and theta.
! on input f is the name of the subroutine containing the integrand,
! a is the lower end point of the interval, b is the higher end point,
! par is an array of constant parameters to be passed to subroutine f,
! and n is the length of the par array. on output result is the
! approximation from applying the 20-point gauss-legendre rule,
! drdeta is the derivative with respect to eta, and drdtheta is the
! derivative with respect to theta.
!
! note: since the number of nodes is even, zero is not an abscissa.
!
! declare
      external    ::  f
      integer     ::  j,n
      real*8      ::  a,b,result,drdeta,drdtheta, &
                      drdeta2,drdtheta2,drdetadtheta,par(n), &
                      absc1,absc2,center,hlfrun,wg(10),xg(10), &
                      fval(2),deta(2),dtheta(2), &
                      deta2(2),dtheta2(2),detadtheta(2)


! the abscissae and weights are given for the interval (-1,1).
! xg     - abscissae of the 20-point gauss-legendre rule
!          for half of the usual run (-1,1), i.e.
!          the positive nodes of the 20-point rule
! wg     - weights of the 20-point gauss rule.
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.


      data xg (  1) /   7.65265211334973337546404093988382110d-2 /
      data xg (  2) /   2.27785851141645078080496195368574624d-1 /
      data xg (  3) /   3.73706088715419560672548177024927237d-1 /
      data xg (  4) /   5.10867001950827098004364050955250998d-1 /
      data xg (  5) /   6.36053680726515025452836696226285936d-1 /
      data xg (  6) /   7.46331906460150792614305070355641590d-1 /
      data xg (  7) /   8.39116971822218823394529061701520685d-1 /
      data xg (  8) /   9.12234428251325905867752441203298113d-1 /
      data xg (  9) /   9.63971927277913791267666131197277221d-1 /
      data xg ( 10) /   9.93128599185094924786122388471320278d-1 /

      data wg (  1) /   1.52753387130725850698084331955097593d-1 /
      data wg (  2) /   1.49172986472603746787828737001969436d-1 /
      data wg (  3) /   1.42096109318382051329298325067164933d-1 /
      data wg (  4) /   1.31688638449176626898494499748163134d-1 /
      data wg (  5) /   1.18194531961518417312377377711382287d-1 /
      data wg (  6) /   1.01930119817240435036750135480349876d-1 /
      data wg (  7) /   8.32767415767047487247581432220462061d-2 /
      data wg (  8) /   6.26720483341090635695065351870416063d-2 /
      data wg (  9) /   4.06014298003869413310399522749321098d-2 /
      data wg ( 10) /   1.76140071391521183118619623518528163d-2 /


!           list of major variables
!           -----------------------
!
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 20-point gauss formula


      center       = 0.5d0 * (a+b)
      hlfrun       = 0.5d0 * (b-a)
      result       = 0.0d0
      drdeta       = 0.0d0
      drdtheta     = 0.0d0
      drdeta2      = 0.0d0
      drdtheta2    = 0.0d0
      drdetadtheta = 0.0d0
      do j=1,10
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)

        call f(absc1, par, n, fval(1), deta(1), dtheta(1), &
               deta2(1), dtheta2(1), detadtheta(1))
        call f(absc2, par, n, fval(2), deta(2), dtheta(2), &
               deta2(2), dtheta2(2), detadtheta(2))

        result       = result    + (fval(1)    + fval(2))*wg(j)
        drdeta       = drdeta    + (deta(1)    + deta(2))*wg(j)
        drdtheta     = drdtheta  + (dtheta(1)  + dtheta(2))*wg(j)
        drdeta2      = drdeta2   + (deta2(1)   + deta2(2))*wg(j)
        drdtheta2    = drdtheta2 + (dtheta2(1) + dtheta2(2))*wg(j)
        drdetadtheta = drdetadtheta+(detadtheta(1)+detadtheta(2))*wg(j)
      enddo

      result       = result * hlfrun
      drdeta       = drdeta * hlfrun
      drdtheta     = drdtheta * hlfrun
      drdeta2      = drdeta2 * hlfrun
      drdtheta2    = drdtheta2 * hlfrun
      drdetadtheta = drdetadtheta * hlfrun
      return
      end subroutine dqleg020





      subroutine dqleg040(f,a,b,result,drdeta,drdtheta, &
                          drdeta2,drdtheta2,drdetadtheta,par,n)
      include 'implno.dek'
!
! 40 point gauss-legendre rule for the fermi-dirac function and
! its derivatives with respect to eta and theta.
! on input f is the name of the subroutine containing the integrand,
! a is the lower end point of the interval, b is the higher end point,
! par is an array of constant parameters to be passed to subroutine f,
! and n is the length of the par array. on output result is the
! approximation from applying the 40-point gauss-legendre rule,
! drdeta is the derivative with respect to eta, and drdtheta is the
! derivative with respect to theta.
!
! note: since the number of nodes is even, zero is not an abscissa.
!
! declare
      external    ::  f
      integer     ::  j,n
      real*8      ::  a,b,result,drdeta,drdtheta, &
                      drdeta2,drdtheta2,drdetadtheta,par(n), &
                      absc1,absc2,center,hlfrun,wg(20),xg(20), &
                      fval(2),deta(2),dtheta(2), &
                      deta2(2),dtheta2(2),detadtheta(2)


! the abscissae and weights are given for the interval (-1,1).
! xg     - abscissae of the 40-point gauss-legendre rule
!          for half of the usual run (-1,1), i.e.
!          the positive nodes of the 40-point rule
! wg     - weights of the 40-point gauss rule.
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   3.87724175060508219331934440246232946d-2 /
      data xg (  2) /   1.16084070675255208483451284408024113d-1 /
      data xg (  3) /   1.92697580701371099715516852065149894d-1 /
      data xg (  4) /   2.68152185007253681141184344808596183d-1 /
      data xg (  5) /   3.41994090825758473007492481179194310d-1 /
      data xg (  6) /   4.13779204371605001524879745803713682d-1 /
      data xg (  7) /   4.83075801686178712908566574244823004d-1 /
      data xg (  8) /   5.49467125095128202075931305529517970d-1 /
      data xg (  9) /   6.12553889667980237952612450230694877d-1 /
      data xg ( 10) /   6.71956684614179548379354514961494109d-1 /
      data xg ( 11) /   7.27318255189927103280996451754930548d-1 /
      data xg ( 12) /   7.78305651426519387694971545506494848d-1 /
      data xg ( 13) /   8.24612230833311663196320230666098773d-1 /
      data xg ( 14) /   8.65959503212259503820781808354619963d-1 /
      data xg ( 15) /   9.02098806968874296728253330868493103d-1 /
      data xg ( 16) /   9.32812808278676533360852166845205716d-1 /
      data xg ( 17) /   9.57916819213791655804540999452759285d-1 /
      data xg ( 18) /   9.77259949983774262663370283712903806d-1 /
      data xg ( 19) /   9.90726238699457006453054352221372154d-1 /
      data xg ( 20) /   9.98237709710559200349622702420586492d-1 /

      data wg (  1) /   7.75059479784248112637239629583263269d-2 /
      data wg (  2) /   7.70398181642479655883075342838102485d-2 /
      data wg (  3) /   7.61103619006262423715580759224948230d-2 /
      data wg (  4) /   7.47231690579682642001893362613246731d-2 /
      data wg (  5) /   7.28865823958040590605106834425178358d-2 /
      data wg (  6) /   7.06116473912867796954836308552868323d-2 /
      data wg (  7) /   6.79120458152339038256901082319239859d-2 /
      data wg (  8) /   6.48040134566010380745545295667527300d-2 /
      data wg (  9) /   6.13062424929289391665379964083985959d-2 /
      data wg ( 10) /   5.74397690993915513666177309104259856d-2 /
      data wg ( 11) /   5.32278469839368243549964797722605045d-2 /
      data wg ( 12) /   4.86958076350722320614341604481463880d-2 /
      data wg ( 13) /   4.38709081856732719916746860417154958d-2 /
      data wg ( 14) /   3.87821679744720176399720312904461622d-2 /
      data wg ( 15) /   3.34601952825478473926781830864108489d-2 /
      data wg ( 16) /   2.79370069800234010984891575077210773d-2 /
      data wg ( 17) /   2.22458491941669572615043241842085732d-2 /
      data wg ( 18) /   1.64210583819078887128634848823639272d-2 /
      data wg ( 19) /   1.04982845311528136147421710672796523d-2 /
      data wg ( 20) /   4.52127709853319125847173287818533272d-3 /


!           list of major variables
!           -----------------------
!
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 20-point gauss formula


      center       = 0.5d0 * (a+b)
      hlfrun       = 0.5d0 * (b-a)
      result       = 0.0d0
      drdeta       = 0.0d0
      drdtheta     = 0.0d0
      drdeta2      = 0.0d0
      drdtheta2    = 0.0d0
      drdetadtheta = 0.0d0
      do j=1,20
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)

        call f(absc1, par, n, fval(1), deta(1), dtheta(1), &
               deta2(1), dtheta2(1), detadtheta(1))
        call f(absc2, par, n, fval(2), deta(2), dtheta(2), &
               deta2(2), dtheta2(2), detadtheta(2))

        result       = result    + (fval(1)    + fval(2))*wg(j)
        drdeta       = drdeta    + (deta(1)    + deta(2))*wg(j)
        drdtheta     = drdtheta  + (dtheta(1)  + dtheta(2))*wg(j)
        drdeta2      = drdeta2   + (deta2(1)   + deta2(2))*wg(j)
        drdtheta2    = drdtheta2 + (dtheta2(1) + dtheta2(2))*wg(j)
        drdetadtheta = drdetadtheta+(detadtheta(1)+detadtheta(2))*wg(j)
      enddo

      result       = result * hlfrun
      drdeta       = drdeta * hlfrun
      drdtheta     = drdtheta * hlfrun
      drdeta2      = drdeta2 * hlfrun
      drdtheta2    = drdtheta2 * hlfrun
      drdetadtheta = drdetadtheta * hlfrun
      return
      end subroutine dqleg040





      subroutine dqleg080(f,a,b,result,drdeta,drdtheta, &
                          drdeta2,drdtheta2,drdetadtheta,par,n)
      include 'implno.dek'
!
! 80 point gauss-legendre rule for the fermi-dirac function and
! its derivatives with respect to eta and theta.
! on input f is the name of the subroutine containing the integrand,
! a is the lower end point of the interval, b is the higher end point,
! par is an array of constant parameters to be passed to subroutine f,
! and n is the length of the par array. on output result is the
! approximation from applying the 80-point gauss-legendre rule,
! drdeta is the derivative with respect to eta, and drdtheta is the
! derivative with respect to theta.
!
! note: since the number of nodes is even, zero is not an abscissa.
!
! declare
      external    ::  f
      integer     ::  j,n
      real*8      ::  a,b,result,drdeta,drdtheta, &
                      drdeta2,drdtheta2,drdetadtheta,par(n), &
                      absc1,absc2,center,hlfrun,wg(40),xg(40), &
                      fval(2),deta(2),dtheta(2), &
                      deta2(2),dtheta2(2),detadtheta(2)


! the abscissae and weights are given for the interval (-1,1).
! xg     - abscissae of the 80-point gauss-legendre rule
!          for half of the usual run (-1,1), i.e.
!          the positive nodes of the 80-point rule
! wg     - weights of the 80-point gauss rule.
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.


      data xg (  1) /   1.95113832567939976543512341074545479d-2 /
      data xg (  2) /   5.85044371524206686289933218834177944d-2 /
      data xg (  3) /   9.74083984415845990632784501049369020d-2 /
      data xg (  4) /   1.36164022809143886559241078000717067d-1 /
      data xg (  5) /   1.74712291832646812559339048011286195d-1 /
      data xg (  6) /   2.12994502857666132572388538666321823d-1 /
      data xg (  7) /   2.50952358392272120493158816035004797d-1 /
      data xg (  8) /   2.88528054884511853109139301434713898d-1 /
      data xg (  9) /   3.25664370747701914619112943627358695d-1 /
      data xg ( 10) /   3.62304753499487315619043286358963588d-1 /
      data xg ( 11) /   3.98393405881969227024379642517533757d-1 /
      data xg ( 12) /   4.33875370831756093062386700363181958d-1 /
      data xg ( 13) /   4.68696615170544477036078364935808657d-1 /
      data xg ( 14) /   5.02804111888784987593672750367568003d-1 /
      data xg ( 15) /   5.36145920897131932019857253125400904d-1 /
      data xg ( 16) /   5.68671268122709784725485786624827158d-1 /
      data xg ( 17) /   6.00330622829751743154746299164006848d-1 /
      data xg ( 18) /   6.31075773046871966247928387289336863d-1 /
      data xg ( 19) /   6.60859898986119801735967122844317234d-1 /
      data xg ( 20) /   6.89637644342027600771207612438935266d-1 /
      data xg ( 21) /   7.17365185362099880254068258293815278d-1 /
      data xg ( 22) /   7.44000297583597272316540527930913673d-1 /
      data xg ( 23) /   7.69502420135041373865616068749026083d-1 /
      data xg ( 24) /   7.93832717504605449948639311738454358d-1 /
      data xg ( 25) /   8.16954138681463470371124994012295707d-1 /
      data xg ( 26) /   8.38831473580255275616623043902867064d-1 /
      data xg ( 27) /   8.59431406663111096977192123491656492d-1 /
      data xg ( 28) /   8.78722567678213828703773343639124407d-1 /
      data xg ( 29) /   8.96675579438770683194324071967395986d-1 /
      data xg ( 30) /   9.13263102571757654164733656150947478d-1 /
      data xg ( 31) /   9.28459877172445795953045959075453133d-1 /
      data xg ( 32) /   9.42242761309872674752266004500001735d-1 /
      data xg ( 33) /   9.54590766343634905493481517021029508d-1 /
      data xg ( 34) /   9.65485089043799251452273155671454998d-1 /
      data xg ( 35) /   9.74909140585727793385645230069136276d-1 /
      data xg ( 36) /   9.82848572738629070418288027709116473d-1 /
      data xg ( 37) /   9.89291302499755531026503167136631385d-1 /
      data xg ( 38) /   9.94227540965688277892063503664911698d-1 /
      data xg ( 39) /   9.97649864398237688899494208183122985d-1 /
      data xg ( 40) /   9.99553822651630629880080499094567184d-1 /

      data wg (  1) /   3.90178136563066548112804392527540483d-2 /
      data wg (  2) /   3.89583959627695311986255247722608223d-2 /
      data wg (  3) /   3.88396510590519689317741826687871658d-2 /
      data wg (  4) /   3.86617597740764633270771102671566912d-2 /
      data wg (  5) /   3.84249930069594231852124363294901384d-2 /
      data wg (  6) /   3.81297113144776383442067915657362019d-2 /
      data wg (  7) /   3.77763643620013974897749764263210547d-2 /
      data wg (  8) /   3.73654902387304900267053770578386691d-2 /
      data wg (  9) /   3.68977146382760088391509965734052192d-2 /
      data wg ( 10) /   3.63737499058359780439649910465228136d-2 /
      data wg ( 11) /   3.57943939534160546028615888161544542d-2 /
      data wg ( 12) /   3.51605290447475934955265923886968812d-2 /
      data wg ( 13) /   3.44731204517539287943642267310298320d-2 /
      data wg ( 14) /   3.37332149846115228166751630642387284d-2 /
      data wg ( 15) /   3.29419393976454013828361809019595361d-2 /
      data wg ( 16) /   3.21004986734877731480564902872506960d-2 /
      data wg ( 17) /   3.12101741881147016424428667206035518d-2 /
      data wg ( 18) /   3.02723217595579806612200100909011747d-2 /
      data wg ( 19) /   2.92883695832678476927675860195791396d-2 /
      data wg ( 20) /   2.82598160572768623967531979650145302d-2 /
      data wg ( 21) /   2.71882275004863806744187066805442598d-2 /
      data wg ( 22) /   2.60752357675651179029687436002692871d-2 /
      data wg ( 23) /   2.49225357641154911051178470032198023d-2 /
      data wg ( 24) /   2.37318828659301012931925246135684162d-2 /
      data wg ( 25) /   2.25050902463324619262215896861687390d-2 /
      data wg ( 26) /   2.12440261157820063887107372506131285d-2 /
      data wg ( 27) /   1.99506108781419989288919287151135633d-2 /
      data wg ( 28) /   1.86268142082990314287354141521572090d-2 /
      data wg ( 29) /   1.72746520562693063585842071312909998d-2 /
      data wg ( 30) /   1.58961835837256880449029092291785257d-2 /
      data wg ( 31) /   1.44935080405090761169620745834605500d-2 /
      data wg ( 32) /   1.30687615924013392937868258970563403d-2 /
      data wg ( 33) /   1.16241141207978269164667699954326348d-2 /
      data wg ( 34) /   1.01617660411030645208318503524069436d-2 /
      data wg ( 35) /   8.68394526926085842640945220403428135d-3 /
      data wg ( 36) /   7.19290476811731275267557086795650747d-3 /
      data wg ( 37) /   5.69092245140319864926910711716201847d-3 /
      data wg ( 38) /   4.18031312469489523673930420168135132d-3 /
      data wg ( 39) /   2.66353358951268166929353583166845546d-3 /
      data wg ( 40) /   1.14495000318694153454417194131563611d-3 /


!           list of major variables
!           -----------------------
!
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 20-point gauss formula


      center       = 0.5d0 * (a+b)
      hlfrun       = 0.5d0 * (b-a)
      result       = 0.0d0
      drdeta       = 0.0d0
      drdtheta     = 0.0d0
      drdeta2      = 0.0d0
      drdtheta2    = 0.0d0
      drdetadtheta = 0.0d0
      do j=1,40
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)

        call f(absc1, par, n, fval(1), deta(1), dtheta(1), &
               deta2(1), dtheta2(1), detadtheta(1))
        call f(absc2, par, n, fval(2), deta(2), dtheta(2), &
               deta2(2), dtheta2(2), detadtheta(2))

        result       = result    + (fval(1)    + fval(2))*wg(j)
        drdeta       = drdeta    + (deta(1)    + deta(2))*wg(j)
        drdtheta     = drdtheta  + (dtheta(1)  + dtheta(2))*wg(j)
        drdeta2      = drdeta2   + (deta2(1)   + deta2(2))*wg(j)
        drdtheta2    = drdtheta2 + (dtheta2(1) + dtheta2(2))*wg(j)
        drdetadtheta = drdetadtheta+(detadtheta(1)+detadtheta(2))*wg(j)
      enddo

      result       = result * hlfrun
      drdeta       = drdeta * hlfrun
      drdtheta     = drdtheta * hlfrun
      drdeta2      = drdeta2 * hlfrun
      drdtheta2    = drdtheta2 * hlfrun
      drdetadtheta = drdetadtheta * hlfrun
      return
      end subroutine dqleg080





      subroutine dqlag010(f,a,b,result,drdeta,drdtheta, &
                          drdeta2,drdtheta2,drdetadtheta,par,n)
      include 'implno.dek'
!
! 10 point gauss-laguerre rule for the fermi-dirac function.
! on input f is the external function defining the integrand
! f(x)=g(x)*w(x), where w(x) is the gaussian weight
! w(x)=exp(-(x-a)/b) and g(x) a smooth function,
! a is the lower end point of the interval, b is the higher end point,
! par is an array of constant parameters to be passed to the function f,
! and n is the length of the par array. on output result is the
! approximation from applying the 10-point gauss-laguerre rule.
! since the number of nodes is even, zero is not an abscissa.
!
! declare
      external    ::  f
      integer     ::  j,n
      real*8      ::  a,b,result,drdeta,drdtheta, &
                      drdeta2,drdtheta2,drdetadtheta,par(n), &
                      absc,wg(10),xg(10),fval,deta,dtheta, &
                      deta2,dtheta2,detadtheta


! the abscissae and weights are given for the interval (0,+inf).
! xg     - abscissae of the 10-point gauss-laguerre rule
! wg     - weights of the 10-point gauss rule. since f yet
!          includes the weight function, the values in wg
!          are actually exp(xg) times the standard
!          gauss-laguerre weights
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.37793470540492430830772505652711188d-1 /
      data xg (  2) /   7.29454549503170498160373121676078781d-1 /
      data xg (  3) /   1.80834290174031604823292007575060883d0 /
      data xg (  4) /   3.40143369785489951448253222140839067d0 /
      data xg (  5) /   5.55249614006380363241755848686876285d0 /
      data xg (  6) /   8.33015274676449670023876719727452218d0 /
      data xg (  7) /   1.18437858379000655649185389191416139d1 /
      data xg (  8) /   1.62792578313781020995326539358336223d1 /
      data xg (  9) /   2.19965858119807619512770901955944939d1 /
      data xg ( 10) /   2.99206970122738915599087933407991951d1 /

      data wg (  1) /   3.54009738606996308762226891442067608d-1 /
      data wg (  2) /   8.31902301043580738109829658127849577d-1 /
      data wg (  3) /   1.33028856174932817875279219439399369d0 /
      data wg (  4) /   1.86306390311113098976398873548246693d0 /
      data wg (  5) /   2.45025555808301016607269373165752256d0 /
      data wg (  6) /   3.12276415513518249615081826331455472d0 /
      data wg (  7) /   3.93415269556152109865581245924823077d0 /
      data wg (  8) /   4.99241487219302310201148565243315445d0 /
      data wg (  9) /   6.57220248513080297518766871037611234d0 /
      data wg ( 10) /   9.78469584037463069477008663871859813d0 /


!           list of major variables
!           -----------------------
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 10-point gauss formula

      result       = 0.0d0
      drdeta       = 0.0d0
      drdtheta     = 0.0d0
      drdeta2      = 0.0d0
      drdtheta2    = 0.0d0
      drdetadtheta = 0.0d0

      do j=1,10
       absc = a+b*xg(j)

       call f(absc, par, n, fval, deta, dtheta, &
              deta2,dtheta2,detadtheta)

       result       = result    + fval*wg(j)
       drdeta       = drdeta    + deta*wg(j)
       drdtheta     = drdtheta  + dtheta*wg(j)
       drdeta2      = drdeta2   + deta2*wg(j)
       drdtheta2    = drdtheta2 + dtheta2*wg(j)
       drdetadtheta = drdetadtheta + detadtheta*wg(j)
      enddo

      result       = result*b
      drdeta       = drdeta*b
      drdtheta     = drdtheta*b
      drdeta2      = drdeta2*b
      drdtheta2    = drdtheta2*b
      drdetadtheta = drdetadtheta*b
      return
      end subroutine dqlag010





      subroutine dqlag020(f,a,b,result,drdeta,drdtheta, &
                          drdeta2,drdtheta2,drdetadtheta,par,n)
      include 'implno.dek'
!
! 20 point gauss-laguerre rule for the fermi-dirac function.
! on input f is the external function defining the integrand
! f(x)=g(x)*w(x), where w(x) is the gaussian weight
! w(x)=dexp(-(x-a)/b) and g(x) a smooth function,
! a is the lower end point of the interval, b is the higher end point,
! par is an array of constant parameters to be passed to the function f,
! and n is the length of the par array. on output result is the
! approximation from applying the 20-point gauss-laguerre rule.
! since the number of nodes is even, zero is not an abscissa.
!
! declare
      external    ::  f
      integer     ::  j,n
      real*8      ::  a,b,result,drdeta,drdtheta, &
                      drdeta2,drdtheta2,drdetadtheta,par(n), &
                      absc,wg(20),xg(20),fval,deta,dtheta, &
                      deta2,dtheta2,detadtheta


! the abscissae and weights are given for the interval (0,+inf).
! xg     - abscissae of the 20-point gauss-laguerre rule
! wg     - weights of the 20-point gauss rule. since f yet
!          includes the weight function, the values in wg
!          are actually exp(xg) times the standard
!          gauss-laguerre weights
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   7.05398896919887533666890045842150958d-2 /
      data xg (  2) /   3.72126818001611443794241388761146636d-1 /
      data xg (  3) /   9.16582102483273564667716277074183187d-1 /
      data xg (  4) /   1.70730653102834388068768966741305070d0 /
      data xg (  5) /   2.74919925530943212964503046049481338d0 /
      data xg (  6) /   4.04892531385088692237495336913333219d0 /
      data xg (  7) /   5.61517497086161651410453988565189234d0 /
      data xg (  8) /   7.45901745367106330976886021837181759d0 /
      data xg (  9) /   9.59439286958109677247367273428279837d0 /
      data xg ( 10) /   1.20388025469643163096234092988655158d1 /
      data xg ( 11) /   1.48142934426307399785126797100479756d1 /
      data xg ( 12) /   1.79488955205193760173657909926125096d1 /
      data xg ( 13) /   2.14787882402850109757351703695946692d1 /
      data xg ( 14) /   2.54517027931869055035186774846415418d1 /
      data xg ( 15) /   2.99325546317006120067136561351658232d1 /
      data xg ( 16) /   3.50134342404790000062849359066881395d1 /
      data xg ( 17) /   4.08330570567285710620295677078075526d1 /
      data xg ( 18) /   4.76199940473465021399416271528511211d1 /
      data xg ( 19) /   5.58107957500638988907507734444972356d1 /
      data xg ( 20) /   6.65244165256157538186403187914606659d1 /

      data wg (  1) /   1.81080062418989255451675405913110644d-1 /
      data wg (  2) /   4.22556767878563974520344172566458197d-1 /
      data wg (  3) /   6.66909546701848150373482114992515927d-1 /
      data wg (  4) /   9.15352372783073672670604684771868067d-1 /
      data wg (  5) /   1.16953970719554597380147822239577476d0 /
      data wg (  6) /   1.43135498592820598636844994891514331d0 /
      data wg (  7) /   1.70298113798502272402533261633206720d0 /
      data wg (  8) /   1.98701589079274721410921839275129020d0 /
      data wg (  9) /   2.28663578125343078546222854681495651d0 /
      data wg ( 10) /   2.60583472755383333269498950954033323d0 /
      data wg ( 11) /   2.94978373421395086600235416827285951d0 /
      data wg ( 12) /   3.32539578200931955236951937421751118d0 /
      data wg ( 13) /   3.74225547058981092111707293265377811d0 /
      data wg ( 14) /   4.21423671025188041986808063782478746d0 /
      data wg ( 15) /   4.76251846149020929695292197839096371d0 /
      data wg ( 16) /   5.42172604424557430380308297989981779d0 /
      data wg ( 17) /   6.25401235693242129289518490300707542d0 /
      data wg ( 18) /   7.38731438905443455194030019196464791d0 /
      data wg ( 19) /   9.15132873098747960794348242552950528d0 /
      data wg ( 20) /   1.28933886459399966710262871287485278d1 /


!           list of major variables
!           -----------------------
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 20-point gauss formula

      result       = 0.0d0
      drdeta       = 0.0d0
      drdtheta     = 0.0d0
      drdeta2      = 0.0d0
      drdtheta2    = 0.0d0
      drdetadtheta = 0.0d0

      do j=1,20
       absc = a+b*xg(j)

       call f(absc, par, n, fval, deta, dtheta, &
              deta2,dtheta2,detadtheta)

       result       = result    + fval*wg(j)
       drdeta       = drdeta    + deta*wg(j)
       drdtheta     = drdtheta  + dtheta*wg(j)
       drdeta2      = drdeta2   + deta2*wg(j)
       drdtheta2    = drdtheta2 + dtheta2*wg(j)
       drdetadtheta = drdetadtheta + detadtheta*wg(j)
      enddo

      result       = result*b
      drdeta       = drdeta*b
      drdtheta     = drdtheta*b
      drdeta2      = drdeta2*b
      drdtheta2    = drdtheta2*b
      drdetadtheta = drdetadtheta*b
      return
      end subroutine dqlag020




      subroutine dqlag040(f,a,b,result,drdeta,drdtheta, &
                          drdeta2,drdtheta2,drdetadtheta,par,n)
      include 'implno.dek'
!
! 20 point gauss-laguerre rule for the fermi-dirac function.
! on input f is the external function defining the integrand
! f(x)=g(x)*w(x), where w(x) is the gaussian weight
! w(x)=dexp(-(x-a)/b) and g(x) a smooth function,
! a is the lower end point of the interval, b is the higher end point,
! par is an array of constant parameters to be passed to the function f,
! and n is the length of the par array. on output result is the
! approximation from applying the 20-point gauss-laguerre rule.
! since the number of nodes is even, zero is not an abscissa.
!
! declare
      external    ::  f
      integer     ::  j,n
      real*8      ::  a,b,result,drdeta,drdtheta, &
                      drdeta2,drdtheta2,drdetadtheta,par(n), &
                      absc,wg(40),xg(40),fval,deta,dtheta, &
                      deta2,dtheta2,detadtheta


! the abscissae and weights are given for the interval (0,+inf).
! xg     - abscissae of the 20-point gauss-laguerre rule
! wg     - weights of the 20-point gauss rule. since f yet
!          includes the weight function, the values in wg
!          are actually exp(xg) times the standard
!          gauss-laguerre weights
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   3.57003943088883851220844712866008554d-2 /
      data xg (  2) /   1.88162283158698516003589346219095913d-1 /
      data xg (  3) /   4.62694281314576453564937524561190364d-1 /
      data xg (  4) /   8.59772963972934922257272224688722412d-1 /
      data xg (  5) /   1.38001082052733718649800032959526559d0 /
      data xg (  6) /   2.02420913592282673344206600280013075d0 /
      data xg (  7) /   2.79336935350681645765351448602664039d0 /
      data xg (  8) /   3.68870267790827020959152635190868698d0 /
      data xg (  9) /   4.71164114655497269361872283627747369d0 /
      data xg ( 10) /   5.86385087834371811427316423799582987d0 /
      data xg ( 11) /   7.14724790810228825068569195197942362d0 /
      data xg ( 12) /   8.56401701758616376271852204208813232d0 /
      data xg ( 13) /   1.01166340484519394068496296563952448d1 /
      data xg ( 14) /   1.18078922940045848428415867043606304d1 /
      data xg ( 15) /   1.36409337125370872283716763606501202d1 /
      data xg ( 16) /   1.56192858933390738372019636521880145d1 /
      data xg ( 17) /   1.77469059500956630425738774954243772d1 /
      data xg ( 18) /   2.00282328345748905296126148101751172d1 /
      data xg ( 19) /   2.24682499834984183513717862289945366d1 /
      data xg ( 20) /   2.50725607724262037943960862094009769d1 /
      data xg ( 21) /   2.78474800091688627207517041404557997d1 /
      data xg ( 22) /   3.08001457394454627007543851961911114d1 /
      data xg ( 23) /   3.39386570849137196090988585862819990d1 /
      data xg ( 24) /   3.72722458804760043283207609906074207d1 /
      data xg ( 25) /   4.08114928238869204661556755816006426d1 /
      data xg ( 26) /   4.45686031753344627071230206344983559d1 /
      data xg ( 27) /   4.85577635330599922809620488067067936d1 /
      data xg ( 28) /   5.27956111872169329693520211373917638d1 /
      data xg ( 29) /   5.73018633233936274950337469958921651d1 /
      data xg ( 30) /   6.21001790727751116121681990578989921d1 /
      data xg ( 31) /   6.72193709271269987990802775518887054d1 /
      data xg ( 32) /   7.26951588476124621175219277242619385d1 /
      data xg ( 33) /   7.85728029115713092805438968334812596d1 /
      data xg ( 34) /   8.49112311357049845427015647096663186d1 /
      data xg ( 35) /   9.17898746712363769923371934806273153d1 /
      data xg ( 36) /   9.93208087174468082501090541654868123d1 /
      data xg ( 37) /   1.07672440639388272520796767611322664d2 /
      data xg ( 38) /   1.17122309512690688807650644123550702d2 /
      data xg ( 39) /   1.28201841988255651192541104389631263d2 /
      data xg ( 40) /   1.42280044469159997888348835359541764d2 /

      data wg (  1) /   9.16254711574598973115116980801374830d-2 /
      data wg (  2) /   2.13420584905012080007193367121512341d-1 /
      data wg (  3) /   3.35718116680284673880510701616292191d-1 /
      data wg (  4) /   4.58540935033497560385432380376452497d-1 /
      data wg (  5) /   5.82068165779105168990996365401543283d-1 /
      data wg (  6) /   7.06495216367219392989830015673016682d-1 /
      data wg (  7) /   8.32026903003485238099112947978349523d-1 /
      data wg (  8) /   9.58878198794443111448122679676028906d-1 /
      data wg (  9) /   1.08727616203054971575386933317202661d0 /
      data wg ( 10) /   1.21746232797778097895427785066560948d0 /
      data wg ( 11) /   1.34969549135676530792393859442394519d0 /
      data wg ( 12) /   1.48425492977684671120561178612978719d0 /
      data wg ( 13) /   1.62144416281182197802316884316454527d0 /
      data wg ( 14) /   1.76159537467676961118424220420981598d0 /
      data wg ( 15) /   1.90507466589479967668299320597279371d0 /
      data wg ( 16) /   2.05228834726171671760199582272947454d0 /
      data wg ( 17) /   2.20369055324509588909828344328140570d0 /
      data wg ( 18) /   2.35979253852320332354037375378901497d0 /
      data wg ( 19) /   2.52117414037643299165313690287422820d0 /
      data wg ( 20) /   2.68849805540884226415950544706374659d0 /
      data wg ( 21) /   2.86252781321044881203476395983104311d0 /
      data wg ( 22) /   3.04415066531151710041043967954333670d0 /
      data wg ( 23) /   3.23440709726353194177490239428867111d0 /
      data wg ( 24) /   3.43452939842774809220398481891602464d0 /
      data wg ( 25) /   3.64599282499408907238965646699490434d0 /
      data wg ( 26) /   3.87058459721651656808475320213444338d0 /
      data wg ( 27) /   4.11049868043282265583582247263951577d0 /
      data wg ( 28) /   4.36846872325406347450808338272945025d0 /
      data wg ( 29) /   4.64795898407446688299303399883883991d0 /
      data wg ( 30) /   4.95344611240989326218696150785562721d0 /
      data wg ( 31) /   5.29084840590073657468737365718858968d0 /
      data wg ( 32) /   5.66820460903297677000730529023263795d0 /
      data wg ( 33) /   6.09679641474342030593376010859198806d0 /
      data wg ( 34) /   6.59310886103999953794429664206294899d0 /
      data wg ( 35) /   7.18249599553689315064429801626699574d0 /
      data wg ( 36) /   7.90666631138422877369310742310586595d0 /
      data wg ( 37) /   8.84089249281034652079125595063026792d0 /
      data wg ( 38) /   1.01408992656211694839094600306940468d1 /
      data wg ( 39) /   1.22100212992046038985226485875881108d1 /
      data wg ( 40) /   1.67055206420242974052468774398573553d1 /


!           list of major variables
!           -----------------------
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 20-point gauss formula


      result       = 0.0d0
      drdeta       = 0.0d0
      drdtheta     = 0.0d0
      drdeta2      = 0.0d0
      drdtheta2    = 0.0d0
      drdetadtheta = 0.0d0

      do j=1,40
       absc = a+b*xg(j)

       call f(absc, par, n, fval, deta, dtheta, &
              deta2,dtheta2,detadtheta)

       result       = result    + fval*wg(j)
       drdeta       = drdeta    + deta*wg(j)
       drdtheta     = drdtheta  + dtheta*wg(j)
       drdeta2      = drdeta2   + deta2*wg(j)
       drdtheta2    = drdtheta2 + dtheta2*wg(j)
       drdetadtheta = drdetadtheta + detadtheta*wg(j)
      enddo

      result       = result*b
      drdeta       = drdeta*b
      drdtheta     = drdtheta*b
      drdeta2      = drdeta2*b
      drdtheta2    = drdtheta2*b
      drdetadtheta = drdetadtheta*b
      return
      end subroutine dqlag040





      subroutine dqlag080(f,a,b,result,drdeta,drdtheta, &
                          drdeta2,drdtheta2,drdetadtheta,par,n)
      include 'implno.dek'
!
! 20 point gauss-laguerre rule for the fermi-dirac function.
! on input f is the external function defining the integrand
! f(x)=g(x)*w(x), where w(x) is the gaussian weight
! w(x)=dexp(-(x-a)/b) and g(x) a smooth function,
! a is the lower end point of the interval, b is the higher end point,
! par is an array of constant parameters to be passed to the function f,
! and n is the length of the par array. on output result is the
! approximation from applying the 20-point gauss-laguerre rule.
! since the number of nodes is even, zero is not an abscissa.
!
! declare
      external    ::  f
      integer     ::  j,n
      real*8      ::  a,b,result,drdeta,drdtheta, &
                      drdeta2,drdtheta2,drdetadtheta,par(n), &
                      absc,wg(80),xg(80),fval,deta,dtheta, &
                      deta2,dtheta2,detadtheta


! the abscissae and weights are given for the interval (0,+inf).
! xg     - abscissae of the 20-point gauss-laguerre rule
! wg     - weights of the 20-point gauss rule. since f yet
!          includes the weight function, the values in wg
!          are actually exp(xg) times the standard
!          gauss-laguerre weights
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.79604233006983655540103192474016803d-2 /
      data xg (  2) /   9.46399129943539888113902724652172943d-2 /
      data xg (  3) /   2.32622868125867569207706157216349831d-1 /
      data xg (  4) /   4.31992547802387480255786172497770411d-1 /
      data xg (  5) /   6.92828861352021839905702213635446867d-1 /
      data xg (  6) /   1.01523255618947143744625436859935350d0 /
      data xg (  7) /   1.39932768784287277414419051430978382d0 /
      data xg (  8) /   1.84526230383584513811177117769599966d0 /
      data xg (  9) /   2.35320887160926152447244708016140181d0 /
      data xg ( 10) /   2.92336468655542632483691234259732862d0 /
      data xg ( 11) /   3.55595231404613405944967308324638370d0 /
      data xg ( 12) /   4.25122008230987808316485766448577637d0 /
      data xg ( 13) /   5.00944263362016477243367706818206389d0 /
      data xg ( 14) /   5.83092153860871901982127113295605083d0 /
      data xg ( 15) /   6.71598597785131711156550087635199430d0 /
      data xg ( 16) /   7.66499349489177306073418909047823480d0 /
      data xg ( 17) /   8.67833082516770109543442255542661083d0 /
      data xg ( 18) /   9.75641480574293071316550944617366591d0 /
      data xg ( 19) /   1.08996933712878553774361001021489406d1 /
      data xg ( 20) /   1.21086466423656999007054848698315593d1 /
      data xg ( 21) /   1.33837881127786473701629840603833297d1 /
      data xg ( 22) /   1.47256659435085855393358076261838437d1 /
      data xg ( 23) /   1.61348643716624665791658545428990907d1 /
      data xg ( 24) /   1.76120052438144378598635686943586520d1 /
      data xg ( 25) /   1.91577496842412479221729970205674985d1 /
      data xg ( 26) /   2.07727999097920960924419379010489579d1 /
      data xg ( 27) /   2.24579012045404583114095916950877516d1 /
      data xg ( 28) /   2.42138440689586473771922469392447092d1 /
      data xg ( 29) /   2.60414665601655866929390053565435682d1 /
      data xg ( 30) /   2.79416568418594655558233069293692111d1 /
      data xg ( 31) /   2.99153559649009855011270412115737715d1 /
      data xg ( 32) /   3.19635609022089207107748887542636533d1 /
      data xg ( 33) /   3.40873278647261898749834947342860505d1 /
      data xg ( 34) /   3.62877759287814544588031988436216948d1 /
      data xg ( 35) /   3.85660910092922104582563052172908535d1 /
      data xg ( 36) /   4.09235302180312671999095850595544326d1 /
      data xg ( 37) /   4.33614266517312302957826760468219500d1 /
      data xg ( 38) /   4.58811946612788863456266489974878378d1 /
      data xg ( 39) /   4.84843356608331891358737273353563006d1 /
      data xg ( 40) /   5.11724445446070105959889432334907144d1 /
      data xg ( 41) /   5.39472167895544471206210278787572430d1 /
      data xg ( 42) /   5.68104563346362231341248503244102122d1 /
      data xg ( 43) /   5.97640843421099549427295961277471927d1 /
      data xg ( 44) /   6.28101489639264772036272917590288682d1 /
      data xg ( 45) /   6.59508362574560573434640627160792248d1 /
      data xg ( 46) /   6.91884824202362773741980288648237373d1 /
      data xg ( 47) /   7.25255875442633453588389652616568450d1 /
      data xg ( 48) /   7.59648311278641748269449794974796502d1 /
      data xg ( 49) /   7.95090896290888369620572826259980809d1 /
      data xg ( 50) /   8.31614564010536896630429506875848705d1 /
      data xg ( 51) /   8.69252644196156234481165926040448396d1 /
      data xg ( 52) /   9.08041123009407559518411727820318427d1 /
      data xg ( 53) /   9.48018942159474332072071889138735302d1 /
      data xg ( 54) /   9.89228344469405791648019372738036790d1 /
      data xg ( 55) /   1.03171527508039130233047094167345654d2 /
      data xg ( 56) /   1.07552984977539906327607890798975954d2 /
      data xg ( 57) /   1.12072690484128333623930046166211013d2 /
      data xg ( 58) /   1.16736664673503666318157888130801099d2 /
      data xg ( 59) /   1.21551542490952625566863895752110813d2 /
      data xg ( 60) /   1.26524665796515540341570265431653573d2 /
      data xg ( 61) /   1.31664195252120310870089086308006192d2 /
      data xg ( 62) /   1.36979246686936973947570637289463788d2 /
      data xg ( 63) /   1.42480058912161601930826569200455232d2 /
      data xg ( 64) /   1.48178202455004441818652384836007732d2 /
      data xg ( 65) /   1.54086842281798697859417425265596259d2 /
      data xg ( 66) /   1.60221072870095715935268416893010646d2 /
      data xg ( 67) /   1.66598351934053918744521179733712213d2 /
      data xg ( 68) /   1.73239071334249503830906503775056999d2 /
      data xg ( 69) /   1.80167323049032317982430208997701523d2 /
      data xg ( 70) /   1.87411949676963772390490134588021771d2 /
      data xg ( 71) /   1.95008022441532991450390479600599643d2 /
      data xg ( 72) /   2.02998984195074937824807677823714777d2 /
      data xg ( 73) /   2.11439870494836466691484904695542608d2 /
      data xg ( 74) /   2.20402368151735739654044206677763168d2 /
      data xg ( 75) /   2.29983206075680004348410969675844754d2 /
      data xg ( 76) /   2.40319087055841540417597460479219628d2 /
      data xg ( 77) /   2.51615879330499611167444939310973194d2 /
      data xg ( 78) /   2.64213823883199102097696108691435553d2 /
      data xg ( 79) /   2.78766733046004563652014172530611597d2 /
      data xg ( 80) /   2.96966511995651345758852859155703581d2 /

      data wg (  1) /   4.60931031330609664705251321395510083d-2 /
      data wg (  2) /   1.07313007783932752564150320304398860d-1 /
      data wg (  3) /   1.68664429547948111794220457782702406d-1 /
      data wg (  4) /   2.30088089384940054411257181978193282d-1 /
      data wg (  5) /   2.91601302502437964832169318772943752d-1 /
      data wg (  6) /   3.53226753575408236352723125805647046d-1 /
      data wg (  7) /   4.14988177550940466187197686311280092d-1 /
      data wg (  8) /   4.76909792302936241314777025418505661d-1 /
      data wg (  9) /   5.39016218474955374499507656522327912d-1 /
      data wg ( 10) /   6.01332497447190529086765248840739512d-1 /
      data wg ( 11) /   6.63884136396680571849442240727299214d-1 /
      data wg ( 12) /   7.26697163614156688973567296249140514d-1 /
      data wg ( 13) /   7.89798189428428531349793078398788294d-1 /
      data wg ( 14) /   8.53214471438152298354598162431362968d-1 /
      data wg ( 15) /   9.16973983833892698590342900031553302d-1 /
      data wg ( 16) /   9.81105491004005747195060155984218607d-1 /
      data wg ( 17) /   1.04563862580654218147568445663176029d0 /
      data wg ( 18) /   1.11060397300025890771124763259729371d0 /
      data wg ( 19) /   1.17603315841226175056651076519208666d0 /
      data wg ( 20) /   1.24195894449809359279351761817871338d0 /
      data wg ( 21) /   1.30841533303134064261188542845954645d0 /
      data wg ( 22) /   1.37543767574892843813155917093490796d0 /
      data wg ( 23) /   1.44306279387849270398312417207247308d0 /
      data wg ( 24) /   1.51132910758830693847655020559917703d0 /
      data wg ( 25) /   1.58027677653099415830201878723121659d0 /
      data wg ( 26) /   1.64994785280267874116012042819355036d0 /
      data wg ( 27) /   1.72038644781283277182004281452290770d0 /
      data wg ( 28) /   1.79163891476093832891442620527688915d0 /
      data wg ( 29) /   1.86375404864909708435925709028688162d0 /
      data wg ( 30) /   1.93678330603070923513925434327841646d0 /
      data wg ( 31) /   2.01078104701134222912614988175555546d0 /
      data wg ( 32) /   2.08580480238741046429303978512989079d0 /
      data wg ( 33) /   2.16191556924159897378316344048827763d0 /
      data wg ( 34) /   2.23917813882364652373453997447445645d0 /
      data wg ( 35) /   2.31766146114651854068606048043496370d0 /
      data wg ( 36) /   2.39743905144001430514117238638849980d0 /
      data wg ( 37) /   2.47858944444973417756369164455222527d0 /
      data wg ( 38) /   2.56119670357790455335115509222572643d0 /
      data wg ( 39) /   2.64535099306968892850463441000367534d0 /
      data wg ( 40) /   2.73114922289915138861410287131169260d0 /
      data wg ( 41) /   2.81869577775934171703141873747811157d0 /
      data wg ( 42) /   2.90810334368223018934550276777492687d0 /
      data wg ( 43) /   2.99949384839685626832412451829968724d0 /
      data wg ( 44) /   3.09299953469357468116695108353033660d0 /
      data wg ( 45) /   3.18876418994712376429365271501623466d0 /
      data wg ( 46) /   3.28694455975337531998378107012216956d0 /
      data wg ( 47) /   3.38771197960397652334054908762154571d0 /
      data wg ( 48) /   3.49125426598732012281732423782764895d0 /
      data wg ( 49) /   3.59777791769613046096294730174902943d0 /
      data wg ( 50) /   3.70751069001745708341027155659228179d0 /
      data wg ( 51) /   3.82070461965311695152029959430467622d0 /
      data wg ( 52) /   3.93763959771430720676800540657330923d0 /
      data wg ( 53) /   4.05862761338354481597420116187988679d0 /
      data wg ( 54) /   4.18401782381424031850607692334503121d0 /
      data wg ( 55) /   4.31420264929613425820084573217987912d0 /
      data wg ( 56) /   4.44962515053655906604982820155377774d0 /
      data wg ( 57) /   4.59078802263617511042959849148929810d0 /
      data wg ( 58) /   4.73826464598929537394753873505838770d0 /
      data wg ( 59) /   4.89271277966692168696886936743283567d0 /
      data wg ( 60) /   5.05489168534039512820572507135175938d0 /
      data wg ( 61) /   5.22568375594272391089278010166022467d0 /
      data wg ( 62) /   5.40612213379727909853323512340717863d0 /
      data wg ( 63) /   5.59742640184041404016553694158980053d0 /
      data wg ( 64) /   5.80104932137643943530626162455394841d0 /
      data wg ( 65) /   6.01873893878099108768015151514026344d0 /
      data wg ( 66) /   6.25262247491437403092934213480091928d0 /
      data wg ( 67) /   6.50532173517668675787482719663696133d0 /
      data wg ( 68) /   6.78011521200777294201287347980059368d0 /
      data wg ( 69) /   7.08117122025414518776174311916759402d0 /
      data wg ( 70) /   7.41389244615305421421695606226687752d0 /
      data wg ( 71) /   7.78544154841612700386232740339230532d0 /
      data wg ( 72) /   8.20557347814596472333905086100917119d0 /
      data wg ( 73) /   8.68801383996161871469419958058255237d0 /
      data wg ( 74) /   9.25286973415578523923556506201979918d0 /
      data wg ( 75) /   9.93114471840215736008370986534009772d0 /
      data wg ( 76) /   1.07739736414646829405750843522990655d1 /
      data wg ( 77) /   1.18738912465097447081950887710877400d1 /
      data wg ( 78) /   1.34228858497264236139734940154089734d1 /
      data wg ( 79) /   1.59197801616897924449554252200185978d1 /
      data wg ( 80) /   2.14214542964372259537521036186415127d1 /


!           list of major variables
!           -----------------------
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 20-point gauss formula

      result       = 0.0d0
      drdeta       = 0.0d0
      drdtheta     = 0.0d0
      drdeta2      = 0.0d0
      drdtheta2    = 0.0d0
      drdetadtheta = 0.0d0

      do j=1,80
       absc = a+b*xg(j)

       call f(absc, par, n, fval, deta, dtheta, &
              deta2,dtheta2,detadtheta)

       result       = result    + fval*wg(j)
       drdeta       = drdeta    + deta*wg(j)
       drdtheta     = drdtheta  + dtheta*wg(j)
       drdeta2      = drdeta2   + deta2*wg(j)
       drdtheta2    = drdtheta2 + dtheta2*wg(j)
       drdetadtheta = drdetadtheta + detadtheta*wg(j)
      enddo

      result       = result*b
      drdeta       = drdeta*b
      drdtheta     = drdtheta*b
      drdeta2      = drdeta2*b
      drdtheta2    = drdtheta2*b
      drdetadtheta = drdetadtheta*b
      return
      end subroutine dqlag080
