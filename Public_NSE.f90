!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NSE subroutine taken from Frank Timmes website, with modification including !
! plasma columb corrections, created by Leung et al.			      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getnse(den_input, temp_input, x_output)
	USE DEFINITION, ONLY: density, temperature
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'


! exercises the nse routines


! declare
      character*7      whose
      integer          i,j,indx(abignet),newguess,iprint
      double precision xtemp,xden,xye

      double precision temp_input, den_input, ye_input
      double precision, dimension(ionmax):: x_output

! formats
01    format(1x,a)
02    format(1x,1p2e14.6)
03    format(1x,a4,'=',1pe10.3,'  ',a4,'=',1pe10.3,'  ', &
                a4,'=',1pe10.3,'  ',a4,'=',1pe10.3)
04    format(1x,a,' ',1p2e11.3)





! initialize the network
!      call init_network

! print some statistics of the root find
! and always use new initial guesses

      xtemp = temp_input/temperature
      xden = den_input/density
      xye = 0.5D0 ! Unchange electron fractions !

      iprint   = 0
      newguess = 1 

      call nse_private(xtemp,xden,xye,newguess,xmass_nse,iprint)

! keep coming back to here
!100   write(6,01) 'give temp, den, ye =>'
!      read (5,*)  xtemp,xden,xye
!      if (xtemp .eq. 0.0) stop 'normal termination'

!      if (xtemp .gt. 2.0e9) then
!       newguess = 1
!       call nse_private(xtemp,xden,xye,newguess,xmass_nse,iprint)

!      else
!       write(6,*) 'no nse state'
!       goto 100
!      end if



! nse mass fractions
      !write(6,*) ' '
      !write(6,*) 'abundance vector:'
      !write(6,03) (ionam(i),xmass_nse(i), i=1,ionmax)

      do i = 1, ionmax, 1
         x_output(i) = xmass_nse(i)
         !write(*,*) x_output(i)
      enddo

      !write(6,*) ' '
      !write(6,*) 'top 10 abundances:'
      !call indexx(ionmax,xmass_nse,indx)
      !write(6,03) (ionam(indx(i)), &
      !             xmass_nse(indx(i)), i=ionmax,ionmax-9,-1)
      !write(6,*) ' '


      !goto 100
      end






      subroutine nse_private(tt,dd,yye,newguess,xmass_out,iprint)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

	
! given the temperature tt, density dd, and electron mole number ye
! this routine puts a chosen reaction network into its nse distribution.

! input:
! tt = temperature
! dd = density
! ye = electron mol number
! newguess = 1 = a new initial guess is made
!          = 0 = values from the previous call are used as the initial guess
! iprint = print flag

! output:
! xmass_out  = nse mass fraction


! declare the pass
      integer          newguess,iprint
      double precision tt,dd,yye,xmass_out(1)


! communicate
      double precision temp,den,ye_want,beta
      common /nsec1/   temp,den,ye_want,beta


! locals
      external         nsefunc_private
      logical          check
      integer          ntrial,nfev,ntaken,n,i
      parameter        (ntrial = 200, n = 2)
      double precision x(n),amass,fac1,fac2,tolf,tolx,twopi, &
                       dum,resid(n)
      parameter        (tolf = 1.0d-8, tolx = 1.0d-14, twopi=2.0d0*pi)



! fill the common block
      temp    = tt
      den     = dd
      ye_want = yye
      beta    = 1.0d0/(kerg * temp)

	

! set the partition functions
! these are generally temperature dependent
! here i'll just take the ground state

      !do i=1,ionmax
      ! wpart(i) = 1.0d0
      !enddo
      !wpart(ineut) = 2.0d0
      !wpart(iprot) = 2.0d0

! Use temperature dependent partition functions !
call get_partition_functions(tt/1.0D9)



! here is an initial guess for the neutron and proton chemical potentials,
! (x1) and (x2) respectively. obtained by setting xmass(ini56) = 1,
! setting mup = mun, and inverting the saha equation.
! all nse networks should at least have ni56 and this appears to be a
! decent guess for all temp, rho, ye combinations.

      if (newguess .eq. 1) then
       newguess = 0
       i      = ini56
       amass  = aion(i) * amu
       fac1   = aion(i)/(avo * den) * wpart(i)
       fac2   = (twopi/(beta*h) * amass/h )**1.5d0
       x(1)   = -(log(fac1*fac2)/beta + bion(i)*ev2erg*1.0d6)/aion(i)
       x(2)   = x(1)
      end if



! root find on mass and charge conservation for
! the chemical potentials of protons and neutrons

      call xnewt_nse_private(ntrial,x,n,tolx,tolf,ntaken,check,nfev,nsefunc_private)


! be sure we converged
      if (check .or. ntaken .eq. ntrial) then
       write(6,*)
       write(6,*) 'check convergence of root finder'
       write(6,*)
      end if


! some optional diagnostics
      if (iprint .eq. 1) then
       write(6,*)
       write(6,110) 'iterations taken             =',ntaken
       write(6,110) 'function evals               =',nfev
       write(6,111) 'roots                        =',x(1),x(2)
       call nsefunc_private(dum,x,resid)
       write(6,111) 'mass conservation   residual =',resid(1)
       write(6,111) 'charge conservation residual =',resid(2)

 110   format(1x,a,i4)
 111   format(1x,a,1p2e14.6)
      end if


      if (check .or. ntaken .eq. ntrial) stop



! fill the output array using the converged values
      call nsefunc_private(dum,x,resid)
      do i=1,ionmax
       xmass_out(i) = xmass(i)
      enddo

      return
      end







      subroutine nsefunc_private(x,y,f)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'


! this routine returns the root find functions.
! input is the point x and y a vector of the unknowns.
! output is the vector of root find functions f, which should be the
! zero vector upon convergence.

! y(1) is input as the neutron chemical potential
! y(2) is input as the proton chemical potential


! declare the pass
      double precision x,y(*),f(*)


! locals
      integer          i,indx(abignet),j
      double precision ye,mu,amass,fac1,fac2,fac3,sum,twopi
      parameter        (twopi = 2.0d0 * pi)


! communicate
      double precision temp,den,ye_want,beta
      common /nsec1/   temp,den,ye_want,beta

 !!!!!!!!!!!!!!!!!!!!PATCH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 double precision, parameter :: facA1 = -0.9052D0, &
                     		facA2 = 0.6322D0, &
                     		facA3 = -SQRT(3.0D0) / 2.0D0 - facA1 / SQRT(facA2)

 double precision :: fac0a, fac0g, fac0
 double precision :: facpa, facpg, facp

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


! chemical potential and mass fraction of each isotope
! hartmann et al, 1985 apj 297 837, eq 2
! take the mass of each isotope to be amu * aion, otherwise a mass formula

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PATCH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      facpa  = (4.00 * pi * den * ye_want / 3.0D0 / amu)**(-1.0D0/3.0D0)
      facpg  = 4.803204D-10**2 * beta / facpa  
      facp   = facA1* (SQRT(facpg * (facA2 + facpg)) - &
               facA2 * LOG(SQRT(facpg / facA2) + &
               SQRT(1.0D0 + facpg / facA2))) + &
               2.0D0 * facA3 * (SQRT(facpg) - ATAN(facpg))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,ionmax

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PATCH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       fac0a  = (4.00 * pi * den * ye_want / 3.0D0 / amu)**(-1.0D0/3.0D0)
       fac0g  = 4.803204D-10**2 * beta / fac0a * (zion(i) ** (5.0D0/3.0D0))
       fac0   = facA1* (SQRT(fac0g * (facA2 + fac0g)) - &
                facA2 * LOG(SQRT(fac0g / facA2) + &
                SQRT(1.0D0 + fac0g / facA2))) + &
                2.0D0 * facA3 * (SQRT(fac0g) - ATAN(fac0g))
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       mu       = (aion(i) - zion(i))*y(1) + zion(i)*y(2)
       amass    = aion(i) * amu
       fac1     = aion(i)/(avo * den) * wpart(i)
       fac2     = ( twopi/(beta*h) * amass/h )**1.5d0
       !fac3     = exp( beta * (mu + bion(i)*ev2erg*1.0d6) )
	fac3     = exp( beta * (mu + bion(i)*ev2erg*1.0d6 - fac0 / beta) )
       xmass(i) = fac1 * fac2 * fac3
      enddo


! sum the mass fractions in ascending order to minimize roundoff
      call indexx_private(ionmax,xmass,indx)
      sum   = 0.0d0
      do i=1,ionmax
       sum   = sum + xmass(indx(i))
      enddo

! sum the mass fractions to form ye
      ye = 0.0d0
      do i=1,ionmax
       j = indx(i)
       ye = ye + zion(j)/aion(j) * xmass(j)
      enddo

! mass and charge conservation are the requirements
      f(1) = sum - 1.0d0
      f(2) = ye - ye_want

      return
      end








      subroutine nsejac_private(x,y,f,dfdy,n,np)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! this routine returns the functions and the jacobian to do the root find on
! input is x, and y(n) a vector of the unknowns. output is f(n)
! and its jacobian dfdy(np,np).

! y(1) is the neutron chemical potential
! y(2) is the proton chemical potential


! declare the pass
      integer          n,np
      double precision x,y(*),f(*),dfdy(np,np)


! locals
      integer          indx(ionmax),i,j
      double precision mu,mubn,mubp,amass,fac1,fac2,fac3,fac4,fac5, &
                       xmbn(ionmax),xmbp(ionmax),sum,sumbn,sumbp, &
                       ye,yebn,yebp,twopi
      parameter        (twopi = 2.0d0 * pi)


! communicate
      double precision temp,den,ye_want,beta
      common /nsec1/   temp,den,ye_want,beta


 !!!!!!!!!!!!!!!!!!!!PATCH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          

 double precision, parameter :: facA1 = -0.9052D0, &
                     		facA2 = 0.6322D0, &                 
                     		facA3 = -SQRT(3.0D0) / 2.0D0 - facA1 / SQRT(facA2)

 double precision :: fac0a, fac0g, fac0
 double precision :: facpa, facpg, facp

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! chemical potential and mass fraction of each isotope
! hartmann et al, 1985 apj 297 837, eq 2
! take the mass of each isotope to be amu * aion, otherwise a mass formula

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PATCH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      facpa  = (4.00 * pi * den * ye_want / 3.0D0 / amu)**(-1.0D0/3.0D0)
      facpg  = 4.803204D-10**2 * beta / facpa
      facp   = facA1* (SQRT(facpg * (facA2 + facpg)) - &
               facA2 * LOG(SQRT(facpg / facA2) + &
               SQRT(1.0D0 + facpg / facA2))) + &
               2.0D0 * facA3 * (SQRT(facpg) - ATAN(facpg))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,ionmax

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PATCH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       fac0a  = (4.00 * pi * den * ye_want / 3.0D0 / amu)**(-1.0D0/3.0D0)
       fac0g  = 4.803204D-10**2 * beta / fac0a * (zion(i) ** (5.0D0/3.0D0))
       fac0   = facA1* (SQRT(fac0g * (facA2 + fac0g)) - &
                facA2 * LOG(SQRT(fac0g / facA2) + &
                SQRT(1.0D0 + fac0g / facA2))) + &
                2.0D0 * facA3 * (SQRT(fac0g) - ATAN(fac0g))
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       mu       = (aion(i) - zion(i)) * y(1) + zion(i) * y(2)
       mubn     = aion(i) - zion(i)
       mubp     = zion(i)

       amass    = aion(i) * amu

       fac1     = aion(i)/(avo * den) * wpart(i)
       fac2     = ( twopi/(beta*h) * amass/h )**1.5d0
       !fac3     = exp( beta * (mu + bion(i) * ev2erg * 1.0d6) )
	fac3     = exp( beta * (mu + bion(i) * ev2erg * 1.0d6 - fac0 / beta) )
       fac4     = fac1 * fac2 * fac3

       xmass(i) = fac4
       xmbn(i)  = fac4 * beta * mubn
       xmbp(i)  = fac4 * beta * mubp

      enddo


! sum the mass fractions in ascending order to minimize roundoff
      call indexx_private(ionmax,xmass,indx)
      sum   = 0.0d0
      sumbn = 0.0d0
      sumbp = 0.0d0
      do i=1,ionmax
       j     = indx(i)
       sum   = sum   + xmass(j)
       sumbn = sumbn + xmbn(j)
       sumbp = sumbp + xmbp(j)
      enddo


! sum the mass fractions to form ye
      ye   = 0.0d0
      yebn = 0.0d0
      yebp = 0.0d0
      do i=1,ionmax
       j    = indx(i)
       fac5 = zion(j)/aion(j)
       ye   = ye   + fac5 * xmass(j)
       yebn = yebn + fac5 * xmbn(j)
       yebp = yebp + fac5 * xmbp(j)
      enddo


! mass and charge conservation are the requirements
      f(1) = sum - 1.0d0
      f(2) = ye - ye_want


! jacobian
      dfdy(1,1) = sumbn
      dfdy(1,2) = sumbp
      dfdy(2,1) = yebn
      dfdy(2,2) = yebp

      return
      end







      subroutine xnewt_nse_private(ntrial,x,n,tolx,tolf,ntaken,check,nfev,func)
      include 'implno.dek'


! given an initial guess x(1:n) for the root of n equations, this routine
! finds the root by a globally convergent newtons method. the vector of
! functions to be zeroed, called fvec(1:n) in the routine below, is
! returned by the user supplied routine func. the output quantity check
! is false on nomal return and true if the routine has converged to a
! local minimum of the function xfminx_nse. if so, try restarting from a
! different initial guess.

! np is the maximum number of equations n
! ntrial is the maximum number of iterations to try
! ntaken is the number of iterations done
! tolf sets the convergence on function values
! tolmin sets the criterion for deciding wheather spurious convergence to
!        a false minimum of xfminx_nse has occured
! tolx is the convergence criteria on deltax
! stpmx is the scaled maximum step length allowed in the line searches
! nfev is the number of function evaluations


! declare the pass
      external         func
      logical          check
      integer          ntrial,n,ntaken,nfev
      double precision x(n),tolx,tolf


! common block communicates values from routine xfminx_nse
      integer          nn,np
      parameter        (np = 4)
      double precision fvec(np)
      common /newtnse/ fvec,nn

! locals
      integer          i,its,j,indx(np)
      double precision tolmin,stpmx,d,den,f,fold,stpmax,sum,temp,test, &
                       fjac(np,np),g(np),p(np),xold(np),xfminx_nse_private,dum

      parameter        (tolmin = 1.0d-12, &
                        stpmx = 2.0d0)



! initialize
      if (n .gt. np) stop 'n > np in routine xnewt'
      nn     = n
      f      = xfminx_nse_private(x,func)
      nfev   = 1
      ntaken = 0


!  test for the initial guess being a root, using a more stringent tolf
      test = 0.0d0
      do i=1,n
       if (abs(fvec(i)) .gt. test) test = abs(fvec(i))
      enddo
      if (test .lt. 0.01*tolf) then
       check = .false.
       return
      end if


! get stpmax for the line search
      sum = 0.0d0
      do i=1,n
       sum = sum + x(i)*x(i)
      enddo
      stpmax = stpmx * max(sqrt(sum),dfloat(n))


! start of iteration loop; get the jacobian
      do its = 1, ntrial
       ntaken = its


! second order accurate numerical jacobian
!      call jac_nse_private(dum,x,fjac,n,n,np,np,func)
!      nfev = nfev + 2*n + 1


! analytic jacobian
       call nsejac_private(dum,x,fvec,fjac,n,np)
       nfev = nfev + 1


! compute grad f for the line searches
       do i=1,n
        sum = 0.0d0
        do j=1,n
         sum = sum + fjac(j,i)*fvec(j)
        enddo
        g(i) = sum
       enddo


! store x, and f and form right hand sides
       do i=1,n
        xold(i) = x(i)
       enddo
       fold = f
       do i=1,n
        p(i) = -fvec(i)
       enddo


! solve the linear systems
       call ludcmp_private(fjac,n,np,indx,d)
       call lubksb_private(fjac,n,np,indx,p)


! line search returns new x and f
! it also gets fvec at the new x when it calls xfminx_nse
       call lnsrch_nse_private(n,xold,fold,g,p,x,f,stpmax,check,nfev,func)


! test for convergence on function value
       test = 0.0d0
       do i=1,n
        if (abs(fvec(i)) .gt. test) test = abs(fvec(i))
       enddo
!      write(6,*) its,test,tolf
       if (test .lt. tolf) then
        check = .false.
!       write(6,*) 'returning on converged function value',test,tolf
        return
       end if

! check for zero gradiant, i.e. spurious convergence
       if (check) then
        test = 0.0d0
        den  = max(f, 0.5d0 * n)
        do i=1,n
         temp = abs(g(i)) * max(abs(x(i)),1.0d0)/den
         if (temp .gt. test) test = temp
        enddo
        if (test .lt. tolmin) then
         check = .true.
        else
         check = .false.
        end if
        return
       end if

! test for convergence on deltax
       test = 0.0d0
       do i=1,n
        temp = (abs(x(i)-xold(i)))/max(abs(x(i)),1.0d0)
        if (temp .gt. test) test = temp
       enddo
       if (test .lt. tolx) return

!      write(6,*) its,test,tolx

! back for another iteration
      enddo
      check = .true.
      return
      end





      subroutine lnsrch_nse_private(n,xold,fold,g,p,x,f,stpmax,check,nfev,func)
      include 'implno.dek'

! given an n dimensional point xold(1:n), the value of the function fold
! and the gradient g(1:n) at the point, and a direction p(1:n), this routine
! finds a new point x(1:n) along the direction of p from xold where the
! function xfminx_nse has decreased "sufficiently". the new function value is
! returned in f. stpmax is an input quanity that limits the length of the
! steps so that the function is not evaluated in regions where it is
! undefined or subject to overflow. p is usually the newton direction. the
! output quantity check is false on normal exit, and true when x is too
! close to xold. in a minimization routine, this usually signals
! convergence and can be ignored. however, in a root finding routine, the
! calling routine should check wheather the convergence is spurious.


! declare the pass
      external         func
      logical          check
      integer          n,nfev
      double precision f,fold,stpmax,g(n),p(n),x(n),xold(n)


! locals
      integer          i
      double precision xfminx_nse_private,a,alam,alam2,alamin,b,disc,f2,rhs1, &
                       rhs2,slope,sum,temp,test,tmplam, &
                       alf,tolx
      parameter        (alf  = 1.0d-4, &
                        tolx = 3.0d-15)


! alf ensures sufficient decrease in the function value, tolx is the
! convergence criterion on deltax


! initialize and scale if the attempted step is too big
      check = .false.
      sum   = 0.0d0
      do i=1,n
       sum = sum + p(i)*p(i)
      enddo
      sum = sqrt(sum)
      if (sum .gt. stpmax) then
       do i=1,n
        p(i) = p(i) * stpmax/sum
       enddo
      end if
      slope = 0.0d0
      do i=1,n
       slope = slope + g(i)*p(i)
      enddo
      if (slope .ge. 0.0) stop 'roundoff problem in lnsrch_nse'


! compute lambda_min
      test = 0.0d0
      do i=1,n
       temp = abs(p(i))/max(abs(xold(i)),1.0d0)
       if (temp .gt. test) test = temp
      enddo
      alamin = tolx/test


! always try a full newton step, start of iteration loop
      alam = 1.0d0
1     continue
      do i=1,n
       x(i) = xold(i) + alam*p(i)
      enddo


      f    = xfminx_nse_private(x,func)
      nfev = nfev + 1



! convergence on deltax, for root finding, the calling routine
! should verify the convergence
      if (alam .lt. alamin) then
       do i=1,n
        x(i) = xold(i)
       enddo
       check = .true.
       return

! sufficient function decrease
      else if (f .le. fold + alf*alam*slope) then
       return

! backtrack
      else
       if (alam .eq. 1.0) then
        tmplam = -slope / (2.0d0 * (f-fold-slope))
       else
        rhs1 = f  - fold - alam*slope
        rhs2 = f2 - fold - alam2*slope
        a    = (rhs1/alam**2 - rhs2/alam2**2)/(alam-alam2)
        b    = (-alam2*rhs1/alam**2 + alam*rhs2/alam2**2) / (alam-alam2)
        if (a .eq. 0.0) then
         tmplam = -slope/(2.0d0 * b)
        else
         disc = b*b - 3.0d0 * a * slope
         if (disc .lt. 0.0) then
          tmplam = 0.5d0 * alam
         else if (b .le. 0.0) then
          tmplam = (-b + sqrt(disc)) / (3.0d0 * a)
         else
          tmplam = -slope/(b + sqrt(disc))
         end if
        end if
        if (tmplam .gt. 0.5d0*alam) tmplam = 0.5d0*alam
       end if
      end if

! store for the next trip through
      alam2 = alam
      f2    = f
      alam  = max(tmplam, 0.1d0*alam)
      goto 1
      end






      double precision function xfminx_nse_private(x,func)
      include 'implno.dek'


! returns f = 0.5 f dot f at x. func is a user supplied routine of the
! functions to be root found.

! declare the pass
      external         func
      double precision x(1)


! locals
      integer          i
      double precision sum,dum


! common block communicates values back to routine xnewt
      integer          nn,np
      parameter        (np = 4)
      double precision fvec(np)
      common /newtnse/ fvec,nn



      call func(dum,x,fvec)

      sum = 0.0d0
      do i=1,nn
       sum = sum + fvec(i)*fvec(i)
      enddo

      xfminx_nse_private = 0.5d0 * sum
      return
      end






      subroutine jac_nse_private(x,y,dfdy,mcol,nrow,mmax,nmax,derivs)
      include 'implno.dek'

! this routine computes a second order accurate jacobian matrix
! of the function contained in the routine derivs.
! 
! input is the point x and the the vector y(nrow) at which to compute the
! jacobian dfdy(mcol,nrow).
! 
! uses 2*nrow + 1 function evaluations


! declare the pass
      external         derivs
      integer          mcol,nrow,mmax,nmax
      double precision x,y(nmax),dfdy(mmax,nmax)


! locals
      integer          i,j,imax
      parameter        (imax = 4)
      double precision fminus(imax),fplus(imax),rel,ax,temp,h,hinv
      parameter        (rel = 3.162278d-8, &
                        ax  = 1.0d-16)


! check
       if (nrow .gt. imax) stop 'nrow > imax in jacobian2'


! for each row, get the right stepsize
      do j=1,nrow
       temp = y(j)
       h    = rel * max(abs(y(j)),ax)
       y(j) = temp + h
       h    = y(j) - temp
       call derivs(x,y,fplus)
       y(j) = temp

       temp = y(j)
       y(j) = temp - h
       h    = temp - y(j)
       call derivs(x,y,fminus)
       y(j) = temp

! compute the jth row of the jacobian
        hinv = 1.0d0/(2.0d0 * h)
        do i=1,mcol
         dfdy(i,j) = (fplus(i) - fminus(i)) * hinv
        enddo
       enddo

! restore the original state
      call derivs(x,y,fplus)
      return
      end








! lu decomposition:
! routine ludcmp does a pivoting lower-upper decomposition
! routine lubksb does the backsubstitution from ludcmp



      subroutine ludcmp_private(a,n,np,indx,d)
      implicit none
      save

! given the matrix a(n,n), with physical dimensions a(np,ap) this routine
! replaces a by the lu decompostion of a row-wise permutation of itself.
! input are a,n,np. output is a, indx which records the row
! permutations effected by the partial pivoting, and d which is 1 if
! the number of interchanges is even, -1 if odd.
! use routine lubksb to solve a system of linear equations.
! 
! nmax is the largest expected value of n

! declare
      integer          n,np,indx(np),nmax,i,j,k,imax
      parameter        (nmax=500)
      double precision a(np,np),d,tiny,vv(nmax),aamax,sum,dum
      parameter        (tiny=1.0d-20)


! vv stores the implicit scaling of each row
! loop over the rows to get the scaling information
      d = 1.0d0
      do i=1,n
       aamax = 0.0d0
       do j=1,n
        if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
       enddo
       if (aamax .eq. 0.0) stop 'singular matrix in ludcmp'
       vv(i) = 1.0d0/aamax
      enddo

! for each column apply crouts method; see equation 2.3.12
      do j=1,n
       do i=1,j-1
        sum = a(i,j)
        do k=1,i-1
         sum = sum - a(i,k)*a(k,j)
        enddo
        a(i,j) = sum
       enddo

! find the largest pivot element
       aamax = 0.0d0
       do i=j,n
        sum=a(i,j)
        do k=1,j-1
         sum = sum - a(i,k)*a(k,j)
        enddo
        a(i,j) = sum
        dum = vv(i)*abs(sum)
        if (dum .ge. aamax) then
         imax  = i
         aamax = dum
        end if
       enddo

! if we need to interchange rows
       if (j .ne. imax) then
        do k=1,n
         dum       = a(imax,k)
         a(imax,k) = a(j,k)
         a(j,k)    = dum
        enddo
        d          = -d
        vv(imax)   = vv(j)
       end if

! divide by the pivot element
       indx(j) = imax
       if (a(j,j) .eq. 0.0) a(j,j) = tiny
       if (j .ne. n) then
        dum = 1.0d0/a(j,j)
        do i=j+1,n
         a(i,j) = a(i,j)*dum
        enddo
       end if

! and go back for another column of crouts method
      enddo
      return
      end





      subroutine lubksb_private(a,n,np,indx,b)
      implicit none
      save

! solves a set of n linear equations ax=b. a is input in its lu decomposition
! form, determined by the routine above ludcmp. indx is input as the
! permutation vector also returned by ludcmp. b is input as the right hand
! side vector and returns with the solution vector x.
! a,n ans np are not modified by this routine and thus can be left in place
! for successive calls (i.e matrix inversion)

! declare
      integer           n,np,indx(np),i,ii,j,ll
      double precision  a(np,np),b(np),sum

! when ii is > 0, ii becomes the index of the first nonzero element of b
! this is forward substitution of equation 2.3.6, and unscamble in place
      ii = 0
      do i=1,n
       ll = indx(i)
       sum = b(ll)
       b(ll) = b(i)
       if (ii .ne. 0) then
        do j=ii,i-1
         sum = sum - a(i,j) * b(j)
        enddo

! nonzero element was found, so dos the sums in the loop above
       else if (sum .ne. 0.0) then
        ii  = i
       end if
       b(i) = sum
      enddo

! back substitution equation 2.3.7
      do i = n,1,-1
       sum = b(i)
       if (i .lt. n) then
        do j=i+1,n
         sum = sum - a(i,j) * b(j)
        enddo
       end if
       b(i) = sum/a(i,i)
      enddo
      return
      end







      subroutine indexx_private(n,arr,indx)
      include 'implno.dek'
! 
! indexes an array arr(1:n). that is it outputs the array indx(1:n) such
! that arr(indx(j)) is in ascending order for j=1...n. the input quantities
! are not changed.
! 
! declare
      integer          n,indx(n),m,nstack
      parameter        (m=7, nstack = 50)
      integer          i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
      double precision arr(n),a
! 
! initialize
      do 11 j=1,n
       indx(j) = j
11    continue
      jstack = 0
      l      = 1
      ir     = n
! 
! insertion sort when subbarray small enough
1     if (ir - l .lt. m) then
       do 13 j=l+1,ir
        indxt = indx(j)
        a     = arr(indxt)
        do 12 i=j-1,l,-1
         if (arr(indx(i)) .le. a) go to 2
         indx(i+1) = indx(i)
12      continue
        i = l - 1
2       indx(i+1) = indxt
13     continue
! 
! pop stack and begin a new round of partitioning
       if (jstack .eq. 0) return
       ir     = istack(jstack)
       l      = istack(jstack-1)
       jstack = jstack - 2
! 
! choose median of left, center and right elements as partitioning element
! also rearrange so that a(l+1) < a(l) < a(ir)
      else
       k         = (l + ir)/2
       itemp     = indx(k)
       indx(k)   = indx(l+1)
       indx(l+1) = itemp

       if (arr(indx(l)) .gt. arr(indx(ir))) then
        itemp    = indx(l)
        indx(l)  = indx(ir)
        indx(ir) = itemp
       end if


       if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
       endif
       if(arr(indx(l)).gt.arr(indx(l+1)))then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
       endif

! 
! initialize pointers for partitioning
       i     = l + 1
       j     = ir
       indxt = indx(l+1)
       a     = arr(indxt)
3      continue
       i = i + 1
       if (arr(indx(i)) .lt. a) go to 3
4      continue
       j = j - 1
       if (arr(indx(j)) .gt. a) go to 4
       if (j .lt. i) go to 5
       itemp   = indx(i)
       indx(i) = indx(j)
       indx(j) = itemp
       go to 3
! 
5      indx(l+1) = indx(j)
       indx(j)   = indxt
       jstack    = jstack + 2
! 
! push pointers to larger subarray on stack
       if (jstack .gt. nstack) stop 'jstack > nstack in routine indexx'
       if (ir - i + 1  .ge.  j - l) then
        istack(jstack)   = ir
        istack(jstack-1) = i
        ir               = j - 1
       else
        istack(jstack)   = j-1
        istack(jstack-1) = l
        l                = i
       end if
      end if
      go to 1
      end







  subroutine init_network
      include 'implno.dek'
      include 'network.dek'

! this routine initializes stuff for a network

! declare
      integer          i


! for easy zeroing of the isotope pointers
      integer          isotp(nisotp)
      equivalence      (isotp(1),ih1)


! zero all the isotope pointers
      do i=1,nisotp
       isotp(i)   = 0
      enddo


! set the size of the network and the number of rates
      ionmax  = 7

! for partition function !
      ionbeg = 1
      ionend = ionmax

! open the partition function data file and read it

      call read_partition_functions

! set the id numbers of the elements
      ihe4  = 1
      ic12  = 2
      io16  = 3
      ine20 = 4
      img24 = 5
      isi28 = 6
      ini56 = 7
      !ineut = 8
      !iprot = 9

! set the names of the elements
      ionam(ihe4)  = 'he4 '
      ionam(ic12)  = 'c12 '
      ionam(io16)  = 'o16 '
      ionam(ine20) = 'ne20'
      ionam(img24) = 'mg24'
      ionam(isi28) = 'si28'
      ionam(ini56) = 'ni56'
      !ionam(ineut) = 'neut'
      !ionam(iprot) = 'prot'


! set the number of nucleons in the element
      aion(ihe4)  = 4.0d0
      aion(ic12)  = 12.0d0
      aion(io16)  = 16.0d0
      aion(ine20) = 20.0d0
      aion(img24) = 24.0d0
      aion(isi28) = 28.0d0
      aion(ini56) = 56.0d0
      !aion(ineut) = 1.0d0
      !aion(iprot) = 1.0d0




! set the number of protons in the element
      zion(ihe4)  = 2.0d0
      zion(ic12)  = 6.0d0
      zion(io16)  = 8.0d0
      zion(ine20) = 10.0d0
      zion(img24) = 12.0d0
      zion(isi28) = 14.0d0
      zion(ini56) = 28.0d0
      !zion(ineut) = 0.0d0
      !zion(iprot) = 1.0d0


! set the number of neutrons
       do i=1,ionmax
        nion(i) = aion(i) - zion(i)
       enddo



! set the binding energy of the element
      bion(ihe4)  = 28.29603d0
      bion(ic12)  = 92.16294d0
      bion(io16)  = 127.62093d0
      bion(ine20) = 160.64788d0
      bion(img24) = 198.2579d0
      bion(isi28) = 236.5379d0
      bion(ini56) = 484.003d0
      bion(ini58) = 506.4450d0
      !bion(ineut) = 0.0d0
      !bion(iprot) = 0.0d0

      return
      end

subroutine read_partition_functions
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
!
! this routine reads and initializes the partition function data
!
! declare

      character*132    :: string
      integer          :: i,j,izpf,iapf,inpf
      double precision :: xx


! initialize
      pf(1:pf_tmax,1:pf_zmax,1:pf_nmax) = -101.0d0 
      pf_gspin(1:pf_zmax,1:pf_nmax)     = -101.0d0

! open the partition function data file and read the header
      open(unit=22, file='./Library/partition_function_frdm.dat', status='old')
      do i=1,4
       read(22,11) string
11     format(a)
      enddo

! read the file
      do i=1,40000
       read(22,11,end=30) string
       if (string(1:12) .eq. 'END OF TABLE') goto 30
       read(22,*) izpf,iapf,xx
       inpf = iapf - izpf
       if (izpf .lt. 1) stop 'izpf < 1 during read of partition functions'
       if (inpf .lt. 1) stop 'inpf < 1 during read of partition functions'
       if (izpf .gt. pf_zmax) stop 'izpf > pf_izmax during read of partition functions'
       if (inpf .gt. pf_nmax) stop 'inpf > pf_inmax during read of partition functions'
       pf_gspin(izpf,inpf) = xx
       read(22,*) (pf(j,izpf,inpf), j=1,24)
      end do

! close up
30    close(unit=22)


! set the temperatures on which the partition functions are calculated
      pf_t9 = (/0.01d0, 0.15d0,  0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, &
                 0.8d0,  0.9d0,  1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, 3.5d0, & 
                 4.0d0,  4.5d0,  5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0, 10.0d0/)

      return
      end


subroutine get_partition_functions(t9)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'vector_eos.dek'
      include 'network.dek'


! this routine generates partition function values

! declare the pass
     double precision :: t9


! declare locals
! for the partition function interpolation
      integer            :: i,jat
      integer, parameter :: norder=2
      double precision   :: x,x1,x2,a,b,e,alfa,beta,stat_weight


! waprt(i) is the temperature dependent partition functions
! zwork1(i) is the ratio of wpart(i) to the ground state partition function
! zwork2(i) is the temperature derivative of zwork1
! initialize


      zwork1(ionbeg:ionend) = 1.0d0
      zwork2(ionbeg:ionend) = 0.0d0
      wpart(ionbeg:ionend)  = 1.0d0

      !if (ineut .ne. 0) then
      ! zwork1(ineut) = 1.0d0
      ! zwork2(ineut) = 0.0d0
      ! wpart(ineut)  = 2.0d0
      !end if

      !if (iprot .ne. 0) then
      ! zwork1(iprot) = 1.0d0
      ! zwork2(iprot) = 0.0d0
      ! wpart(iprot)  = 2.0d0
      !end if

      !if (ih1 .ne. 0) then
      ! zwork1(ih1) = 1.0d0
      ! zwork2(ih1) = 0.0d0
      ! wpart(ih1)  = 2.0d0
      !end if

      if (ihe4 .ne. 0) then
       zwork1(ihe4) = 1.0d0
       zwork2(ihe4) = 0.0d0
       wpart(ihe4)  = 1.0d0
      end if



! for every isotope not named n, p or alpha

      do i=ionbeg,ionend
	
! bullet check
       if (int(zion(i)) .lt. 1) stop 'izpf < 1 during read of partition functions'
       if (int(nion(i)) .lt. 1) stop 'izpf < 1 during read of partition functions'
       if (int(zion(i)) .gt. pf_zmax) stop 'izpf > pf_izmax during read of partition functions'
       if (int(nion(i)) .gt. pf_nmax) stop 'inpf > pf_inmax during read of partition functions'


! if isotope is in the network 
       if (pf_gspin(int(zion(i)),int(nion(i))) .lt. 0.0d0) cycle

!       write(6,*) int(zion(i)),int(nion(i)),pf_gspin(int(zion(i)),int(nion(i)))


        stat_weight = 2.0d0 * pf_gspin(int(zion(i)),int(nion(i))) + 1.0d0

! clamp at both ends of the temperature grid
       if (t9 .le. pf_t9(1)) then
        zwork1(i) = pf(1,int(zion(i)),int(nion(i)))
        zwork2(i) = 0.0d0
        wpart(i)  = zwork1(i) * stat_weight

       else if (t9 .ge. pf_t9(pf_tmax)) then
        zwork1(i) = pf(pf_tmax,int(zion(i)),int(nion(i)))
        zwork2(i) = 0.0d0
        wpart(i)  = zwork1(i) * stat_weight

! linear larange interpolation

       else
        call locate(pf_t9,pf_tmax,t9,jat)
        jat = max(1,min(jat - norder/2 + 1,pf_tmax - norder + 1))

!        write(6,*) jat,t9
!        write(6,*) pf_t9(jat),pf_t9(jat+1)
!        write(6,*) pf(jat,int(zion(i)),int(nion(i))),pf(jat+1,int(zion(i)),int(nion(i)))

        x  = t9
        x1 = pf_t9(jat)
        x2 = pf_t9(jat+1)
        a  = x - x1
        b  = x - x2
        e  = x1 - x2
        alfa =  b/e
        beta = -a/e
        zwork1(i) = alfa * pf(jat,int(zion(i)),int(nion(i))) + beta * pf(jat+1,int(zion(i)),int(nion(i)))
        zwork2(i) = 0.0d0
        wpart(i)  = zwork1(i) * stat_weight
       end if

!       write(*,*) i, ionam(i), zwork1(i), wpart(i)
!       write(6,156) i, ionam(i), zwork1(i), wpart(i)
!156    format(1x,i4,' ',a,' ',1p2e12.4)
!       read(5,*)
	
      enddo
      return
      end

 subroutine locate(xx,n,x,j)
      implicit none
      save


! given an array xx of length n, and a value of x, this routine returns
! a value j such that x is between xx(j) and xx(j+1). the array xx must be
! monotonic. j=0 or j=n indicates that x is out of range. bisection is used
! to find the entry

! declare
      integer           n,j,jl,ju,jm
      double precision  xx(n),x

! initialize
      jl = 0
      ju = n+1

! compute a midpoint, and replace either the upper or lower limit
 10   if (ju-jl .gt. 1) then
       jm = (ju+jl)/2
       if ( (xx(n) .ge. xx(1)) .eqv. (x .ge. xx(jm)) ) then
        jl = jm
       else
        ju = jm
       end if
       goto 10
      end if
      if (x .eq. xx(1))then
        j = 1
      else if(x .eq. xx(n))then
        j = n - 1
      else
        j = jl
      end if
      return
      end
!---------------------------------------------------------------------