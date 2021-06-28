!=============================================================================!
!=============================================================================!
!!*ALPS                                                     	             *!!
!!                                                                           !!
!!Kristopher Klein, Daniel Verscharen                                        !!
!!
!!Space Science Center, University of New Hampshire                          !!
!!                                                                           !!
!!*NUMERICAL FUNCTIONS                                                      *!!
!=============================================================================!
!=============================================================================!
!

module alps_fns_rel

  implicit none

  private :: int_T_rel,  int_T_res_rel, integrate_resU_rel
  private :: funct_g_rel, LUDCMP, LUBKSB
  private :: determine_sproc_rel,fact

  public :: derivative_f0_rel, polyharmonic_spline
  public :: integrate_res_rel, landau_integrate_rel, resU_rel
  public :: bessj0, bessj1, bessj


contains


!-=-=-=-=-=-=
subroutine derivative_f0_rel(is,is_rel)
    use alps_var, only : f0, pp, nperp, npar, vA, ms, writeOut, arrayName
    use alps_var, only : f0_rel, df0_rel,gamma_rel, pparbar_rel,nspec_rel,ngamma,npparbar
    use alps_var, only : writeOut, pi
    use alps_io,  only : get_unused_unit
    implicit none

    integer :: is, is_rel, n_coarse, counter, igamma, ipparbar,iperp,ipar
    double precision :: pparbar_max,pparbar_min,gamma_max,gamma_min,gamma_max_use,dgamma,dpparbar
    double precision, allocatable, dimension (:) :: gamma_coarse,pparbar_coarse,grid_coarse
    character (50) :: fmt   !Output format
    character (100) :: writename
    double precision :: integrate, gamma, smoothing
    integer :: unit_f

    ! Determine the minimum and maximum values of Gamma and ppar
	pparbar_min=999999.d0
	pparbar_max=-999999.d0
	gamma_min=999999.d0
	gamma_max=-1.d0

	do iperp=0,nperp
		do ipar=0,npar
			gamma = sqrt((pp(is,iperp,ipar,1)**2 + pp(is,iperp,ipar,2)**2) * vA**2/ms(is)**2 +1.d0)
			if (gamma.GT.gamma_max) gamma_max=gamma
			if (gamma.LT.gamma_min) gamma_min=gamma
			if (pp(is,iperp,ipar,2).LT.pparbar_min) pparbar_min=pp(is,iperp,ipar,2)
			if (pp(is,iperp,ipar,2).GT.pparbar_max) pparbar_max=pp(is,iperp,ipar,2)
		enddo
	enddo

	pparbar_min=pparbar_min*vA/ms(is)
	pparbar_max=pparbar_max*vA/ms(is)

	gamma_max_use = sqrt(1.d0+pp(is,nperp,1,1)**2 *vA**2/ms(is)**2)

	if (writeOut) then
		write (*,'(a,i3,a,1es14.4)') "Maximum Gamma for species",is,":", gamma_max
		write (*,'(a,i3,a,1es14.4)') "Minimum Gamma for species",is,":", gamma_min
		write (*,'(a,i3,a,1es14.4)') "Maximum Pparbar: for species",is,":", pparbar_max
		write (*,'(a,i3,a,1es14.4)') "Minimum Pparbar: for species",is,":", pparbar_min
		write (*,'(a,i3,a,1es14.4)') "Usable Maximum Gamma for species",is,":", gamma_max_use
    endif

    smoothing = 0.d0!0.001d0! make this a user-defined parameter ?

    n_coarse=(nperp+1)*(npar+1)	! This has to be (nperp+2)*(npar+1) if we include a line of zeroes
	allocate (gamma_coarse(n_coarse))
	allocate (pparbar_coarse(n_coarse))
	allocate (grid_coarse(n_coarse))


    counter=0
	do iperp = 0, nperp
		do ipar = 0,npar
			counter=counter+1
			gamma_coarse(counter) = sqrt(1.d0+(pp(is,iperp,ipar,1)**2+pp(is,iperp,ipar,2)**2)*&
					vA*vA/(ms(is)*ms(is)))
			pparbar_coarse(counter) = pp(is,iperp,ipar,2)*vA/ms(is)
  			grid_coarse(counter) = log(f0(is,iperp,ipar))
		enddo
	enddo

    ! add a line of zeroes at the end:
!	do ipar = 0, npar
!		counter=counter+1
!		gamma_coarse(counter)=gamma_max+(gamma_max-gamma_min)/(1.d0*ngamma)
!		pparbar_coarse(counter)=pp(is,nperp,ipar,2)*vA/ms(is)
!		grid_coarse(counter)=0.d0
!	enddo


    do igamma = 0,ngamma
		do ipparbar = 0,npparbar
			gamma_rel(is_rel,igamma,ipparbar) = gamma_min + ((gamma_max_use - gamma_min)*igamma)/(1.d0*ngamma)
			pparbar_rel(is_rel,igamma,ipparbar) = pparbar_min + ((pparbar_max - pparbar_min)*ipparbar)/(1.d0*npparbar)
		enddo
	enddo


	if(writeOut) write (*,'(a)') 'Polyharmonic spline interpolation on relativistic grid...'

	call polyharmonic_spline(grid_coarse,gamma_coarse,pparbar_coarse,n_coarse,gamma_rel,pparbar_rel,&
			ngamma,npparbar,smoothing,f0_rel,is_rel,nspec_rel)


	! Stay within the subluminal cone
	do igamma=0,ngamma
	 do ipparbar=0,npparbar

     f0_rel(is_rel,igamma,ipparbar)=exp(f0_rel(is_rel,igamma,ipparbar))

		if ((gamma_rel(is_rel,igamma,ipparbar)**2-1.d0).LT.(pparbar_rel(is_rel,igamma,ipparbar)**2))&
			 f0_rel(is_rel,igamma,ipparbar)=-1.d0
	 enddo
	enddo


      integrate = 0.d0

      dgamma=gamma_rel(is_rel,2,2)-gamma_rel(is_rel,1,2)
      dpparbar=pparbar_rel(is_rel,2,2)-pparbar_rel(is_rel,2,1)
        do igamma = 0, ngamma
             do ipparbar = 0, npparbar
              if(f0_rel(is_rel,igamma,ipparbar).GT.-1.d0) then
                integrate = integrate + &
                     gamma_rel(is_rel,igamma,ipparbar) * f0_rel(is_rel,igamma,ipparbar) * &
                     2.d0 * pi * dgamma * dpparbar * (ms(is) / vA)**3
                endif
             enddo
        enddo


          write(*,'(a,i3,a, 2es14.4)') 'Integration of species', is,':', integrate


	 	if(writeOut) write (*,'(a)') 'Writing relativistic grid to file...'
        write(fmt,'(a)') '(2es14.4,1es14.4)'
	     write(writeName,'(3a,i0,a)') 'distribution/',trim(arrayName),'_f0_rel.',is,'.array'
	      call get_unused_unit(unit_f)
          open(unit=unit_f,file=trim(writeName),status='replace')
		  do igamma = 0, ngamma
		  	do ipparbar = 0, npparbar
			   if (f0_rel(is_rel,igamma,ipparbar).NE.-1.d0) f0_rel(is_rel,igamma,ipparbar)=f0_rel(is_rel,igamma,ipparbar)/integrate
		  	   write(unit_f,fmt) gamma_rel(is_rel,igamma,ipparbar),pparbar_rel(is_rel,igamma,ipparbar),&
					f0_rel(is_rel,igamma,ipparbar)
		  	enddo
           write(unit_f,*)
		  enddo
	 close(unit_f)


   if(writeOut) write (*,'(a)') 'Determining relativistic derivatives...'

        do igamma = 1, ngamma-1
           do ipparbar = 1, npparbar-1

              !index 1-> gamma derivative
              df0_rel(is_rel,igamma,ipparbar,1) = 0.d0

              if ((f0_rel(is_rel,igamma-1,ipparbar).GT.0.d0).AND.(f0_rel(is_rel,igamma+1,ipparbar).GT.0.d0)) then
                df0_rel(is_rel,igamma,ipparbar,1) = &
                   (f0_rel(is_rel,igamma+1,ipparbar) - f0_rel(is_rel,igamma-1,ipparbar))/&
                   (gamma_rel(is_rel,igamma+1,ipparbar) - gamma_rel(is_rel,igamma-1,ipparbar))
             endif

              !index 2-> pparbar derivative
              df0_rel(is_rel,igamma,ipparbar,2) = 0.d0

              if ((f0_rel(is_rel,igamma,ipparbar+1).GT.0.d0).AND.(f0_rel(is_rel,igamma,ipparbar-1).GT.0.d0)) then
                   df0_rel(is_rel,igamma,ipparbar,2) = &
                     (f0_rel(is_rel,igamma,ipparbar+1) - f0_rel(is_rel,igamma,ipparbar-1))/&
                     (pparbar_rel(is_rel,igamma,ipparbar+1) - pparbar_rel(is_rel,igamma,ipparbar-1))
              endif

              if (f0_rel(is_rel,igamma,ipparbar).GE.0.d0) then
                if ((f0_rel(is_rel,igamma,ipparbar+1).LE.0.d0).AND.(f0_rel(is_rel,igamma,ipparbar-1).GT.0.d0)) then
                  df0_rel(is_rel,igamma,ipparbar,2) = &
                   (f0_rel(is_rel,igamma,ipparbar) - f0_rel(is_rel,igamma,ipparbar-1))/&
                   (pparbar_rel(is_rel,igamma,ipparbar) - pparbar_rel(is_rel,igamma,ipparbar-1))
                 endif
                 if ((f0_rel(is_rel,igamma,ipparbar+1).GT.0.d0).AND.(f0_rel(is_rel,igamma,ipparbar-1).LE.0.d0)) then
                      df0_rel(is_rel,igamma,ipparbar,2) = &
                        (f0_rel(is_rel,igamma,ipparbar+1) - f0_rel(is_rel,igamma,ipparbar))/&
                        (pparbar_rel(is_rel,igamma,ipparbar+1) - pparbar_rel(is_rel,igamma,ipparbar))
                 endif
              endif
           enddo
        enddo

	if(writeOut) write (*,'(a)') 'Writing relativistic derivatives to file...'

	      write(fmt,'(a)') '(2es14.4,2es14.4)'

	     write(writeName,'(3a,i0,a)') 'distribution/',trim(arrayName),'_dfdv_rel.',is,'.array'
        call get_unused_unit(unit_f)
          open(unit=unit_f,file=trim(writeName),status='replace')
          do igamma = 0, ngamma
             do ipparbar = 0, npparbar
               write(unit_f,fmt) gamma_rel(is_rel,igamma,ipparbar),pparbar_rel(is_rel,igamma,ipparbar),&
					df0_rel(is_rel,igamma,ipparbar,1),df0_rel(is_rel,igamma,ipparbar,2)
             enddo
             write(unit_f,*)
          enddo
          close(unit_f)

	if (writeOut) write(*,'(a)') '-=-=-=-=-=-=-=-=-'

end subroutine derivative_f0_rel



! Polyharmonic Spline for relativistic calculation:
subroutine polyharmonic_spline(grid_coarse,gamma_coarse,pparbar_coarse,n_coarse,gamma_rel,pparbar,ngamma,npparbar,&
		smoothing,f0_rel,is_rel,nspec_rel)
	!
	! This soubroutine interpolates the grid with a polyharmonic thin-plate spline.
	!
	! Input:
	! grid_coarse is a vector of length n_coarse that includes the values of f at each point i
	! gamma_coarse is a vector of length n_coarse that includes the values of gamma at each point i
	! pparbar_coarse is a vector of length n_coarse that includes the values of ppar at each point i
	! n_coarse is the total number of points in the coarse grid (ngamma_coarse * npparbar_coarse)
	!
	! gamma_rel is the value of gamma in the fine grid. It is a field of rank (ngamma, npparbar)
	! pparbar is the value of ppar in the fine grid. It is a field of rank (ngamma, npparbar)
	! ngamma is the number of perpendicular data points in the fine grid
	! npparbar is the number of parallel data points in the fine grid
	!
	!
	! Output:
	! f0_rel is the interpolated grid. It is a field of rank (ngamma, npparbar)
	!
	! This subroutine needs the subroutines LUDCMP and LUBKSB -- see below
	!

	! This is the Thin Plate Spline:
	! We use these resources:
	! http://cseweb.ucsd.edu/~sjb/eccv_tps.pdf
	! http://www.univie.ac.at/nuhag-php/bibtex/open_files/po94_M%20J%20D%20Powell%2003%2093.pdf
	! http://vision.ucsd.edu/sites/default/files/fulltext(4).pdf
	implicit none

	integer :: i,j,k,permutation_index(n_coarse+3),odd_even,code
	integer :: ngamma,npparbar,n_coarse,is_rel,nspec_rel
	double precision :: gamma_coarse(n_coarse),pparbar_coarse(n_coarse)
	double precision :: grid_coarse(n_coarse),f0_rel(nspec_rel,0:ngamma,0:npparbar)

	double precision :: fullmatrix(n_coarse+3,n_coarse+3)
	double precision :: grid_vector(n_coarse+3),weight_param(n_coarse+3)
	double precision :: gamma_rel(nspec_rel,0:ngamma,0:npparbar),pparbar(nspec_rel,0:ngamma,0:npparbar)
	double precision :: r,smoothing



	grid_vector=0.d0
	do i=1,n_coarse
		grid_vector(i)=grid_coarse(i)
	enddo
	! grid_vector has three additional entries. The last three entries are all zero.

	fullmatrix=0.d0
	do i=1,n_coarse
		do j=1,n_coarse

			! Do the K-matrix part first:
			r=sqrt((gamma_coarse(i)-gamma_coarse(j))**2+(pparbar_coarse(i)-pparbar_coarse(j))**2)
			if(r.GE.1.d0) then
				fullmatrix(i,j)=r*r*log(r)
			elseif (r.EQ.0.d0) then
				fullmatrix(i,j)=0.d0
			else
				fullmatrix(i,j)=r*log(r**r)
			endif
		enddo

		fullmatrix(i,i)=fullmatrix(i,i)+smoothing

		! Now the P-matrix parts:
		fullmatrix(i,n_coarse+1)=1.d0
		fullmatrix(i,n_coarse+2)=gamma_coarse(i)
		fullmatrix(i,n_coarse+3)=pparbar_coarse(i)

		! and the transposed P-matrix:
		fullmatrix(n_coarse+1,i)=1.d0
		fullmatrix(n_coarse+2,i)=gamma_coarse(i)
		fullmatrix(n_coarse+3,i)=pparbar_coarse(i)
	enddo

	weight_param=grid_vector
	call LUDCMP(fullmatrix,n_coarse+3,permutation_index,odd_even,code)
	call LUBKSB(fullmatrix,n_coarse+3,permutation_index,weight_param)

	f0_rel(is_rel,:,:)=0.d0
	do i=0,ngamma
		do j=0,npparbar

		 do k=1,n_coarse
			r=sqrt((gamma_rel(is_rel,i,j)-gamma_coarse(k))**2+(pparbar(is_rel,i,j)-pparbar_coarse(k))**2)
			if (r.GE.1.d0) then
				f0_rel(is_rel,i,j)=f0_rel(is_rel,i,j)+weight_param(k)*r*r*log(r)
			elseif (r.EQ.0.d0) then
				f0_rel(is_rel,i,j)=f0_rel(is_rel,i,j)
			else
				f0_rel(is_rel,i,j)=f0_rel(is_rel,i,j)+weight_param(k)*r*log(r**r)
			endif
		 enddo


		 f0_rel(is_rel,i,j)=f0_rel(is_rel,i,j)+weight_param(n_coarse+1)+&
		 	weight_param(n_coarse+2)*gamma_rel(is_rel,i,j)+weight_param(n_coarse+3)*pparbar(is_rel,i,j)

		enddo
	enddo

end subroutine


!-=-=-=-=-=-=
! The following two functions are from pre-written sources:
!-=-=-=-=-=-=
!*******************************************************
!*    LU decomposition routines 				       *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     *
!*******************************************************


!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
subroutine LUDCMP(a,n,indx,d,code)
	implicit none
	integer, parameter :: nmax=100000
	double precision, parameter :: tiny=1.5d-16
	double precision :: amax, dum, sum, a(n,n), vv(nmax)
	integer :: n, indx(n), code, d, i, j, k, imax


	d=1; code=0
	do i=1,n
		amax=0.d0
		do j=1,n
			if (dabs(a(i,j)).GT.amax) amax=DABS(a(i,j))
		enddo ! j loop
		if(amax.LT.tiny) then
			code = 1
			return
		endif
		vv(i) = 1.d0 / amax
	enddo ! i loop

	do j=1,n
		do i=1,j-1
			sum = a(i,j)
			do k=1,i-1
				sum = sum - a(i,k)*A(k,j)
			enddo ! k loop
				a(i,j) = sum
		enddo ! i loop
		amax = 0.d0
		do i=j,n
			sum = a(i,j)
			do k=1,j-1
				sum = sum - a(i,k)*a(k,j)
			enddo ! k loop
			a(i,j) = sum
			dum = vv(i)*dabs(sum)
			if(dum.GE.amax) then
				imax = i
				amax = dum
			endif
		enddo ! i loop

		if(j.NE.imax) then
			do k=1,n
				dum = a(imax,k)
				a(imax,k) = a(j,k)
				a(j,k) = dum
			enddo ! k loop
			d = -d
			vv(imax) = vv(j)
		endif

		indx(j) = imax
		if(dabs(a(j,j)) < tiny) a(j,j) = tiny

		if(j.NE.n) then
			dum = 1.d0 / a(j,j)
			do i=j+1,n
				a(i,j) = a(i,j)*dum
			enddo ! i loop
		endif
	enddo ! j loop

	return

 end subroutine LUDCMP


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
subroutine LUBKSB(a,n,indx,b)
	implicit none
	integer :: n, indx(n), ii, i, ll, j
	double precision :: sum, a(n,n),b(n)

	ii = 0

	do i=1,n
		ll = indx(i)
		sum = b(ll)
		b(ll) = b(i)
		if(ii.NE.0) then
			do j=ii,i-1
				sum = sum - a(i,j)*b(j)
			enddo ! j loop
		else if(sum.NE.0.d0) then
			ii = i
		endif
		b(i) = sum
	enddo ! i loop

	do i=n,1,-1
		sum = b(i)
		if(i < n) then
			do j=i+1,n
				sum = sum - a(i,j)*b(j)
			enddo ! j loop
		endif
		b(i) = sum / a(i,i)
	enddo ! i loop

	return

end subroutine LUBKSB

!-=-=-=-=-=-=
!
!-=-=-=-=-=-=


! This subroutine determines sproc_rel for the given process
subroutine determine_sproc_rel(sproc_rel)
	use alps_var, only : relativistic, sproc, nspec
	integer :: is_rel, is, sproc_rel

	is_rel=0

	do is=1,nspec
		if (relativistic(is)) then
			is_rel=is_rel+1
			if (is.EQ.sproc) sproc_rel=is_rel
		endif
	enddo

end subroutine



! This function does the relativistic integration around resonances if necessary:
double complex function integrate_res_rel(om,nn,mode)
	use alps_var, only : ngamma, gamma_rel, pparbar_rel
	implicit none
	integer :: igamma
	integer :: nn,mode,sproc_rel
	double precision :: dgamma,dpparbar
	double complex :: ii,om

	integrate_res_rel=cmplx(0.d0,0.d0,kind(1.d0))
	ii=cmplx(0.d0,1.d0,kind(1.d0))

	call determine_sproc_rel(sproc_rel)

	dgamma = gamma_rel(sproc_rel,2,2)-gamma_rel(sproc_rel,1,2)
	dpparbar = pparbar_rel(sproc_rel,2,2)-pparbar_rel(sproc_rel,2,1)


	! At igamma=1, we are already missing the part from igamma=0, where we should actually start.
  ! Therefore, we use 2 instead of 1 in the trapezoid integration:

	do igamma = 1,ngamma - 2
		integrate_res_rel = integrate_res_rel + 2.d0 * integrate_resU_rel(sproc_rel,om,nn,mode,igamma)
	enddo

	igamma = ngamma - 1
	integrate_res_rel = integrate_res_rel + integrate_resU_rel(sproc_rel,om,nn,mode,igamma)

	integrate_res_rel = integrate_res_rel * dgamma * 0.25d0

	return

end function integrate_res_rel




! This function performs the integration of the resU term in (gamma,pparbar) space at a constant igamma
double complex function integrate_resU_rel(sproc_rel,om,nn,mode,igamma)
	use alps_var, only : npparbar,sproc,positions_principal,f0_rel,kpar,vA,ms,qs
	use alps_var, only : gamma_rel,pparbar_rel
	use alps_io, only : alps_error
	implicit none
	integer :: igamma,ipparbar,int_start,int_end
	integer :: nn,mode,ipparbar_res,sproc_rel,lowerlimit,upperlimit,ipparbar_lower,ipparbar_upper
	double precision :: dgamma,dpparbar
	double complex :: ii,om,pparbar_res
	logical :: found_res,found_lower,found_upper

	integrate_resU_rel=cmplx(0.d0,0.d0,kind(1.d0))
	ii=cmplx(0.d0,1.d0,kind(1.d0))


	dgamma = gamma_rel(sproc_rel,2,2)-gamma_rel(sproc_rel,1,2)
	dpparbar = pparbar_rel(sproc_rel,2,2)-pparbar_rel(sproc_rel,2,1)



	! determine the position of the resonance (i.e., the step LEFT of it):
	ipparbar = 0
	ipparbar_res=0
	found_res = .FALSE.

	pparbar_res = (gamma_rel(sproc_rel,igamma,1)*om-(1.d0*nn)*qs(sproc)/ms(sproc))*vA/kpar

  if ((real(pparbar_res)**2).LE.(gamma_rel(sproc_rel,igamma,1)**2-1.d0)) then
	   do while ((ipparbar.LT.(npparbar-2)).AND.(.NOT.found_res))
		     ipparbar = ipparbar + 1

		       if ((pparbar_rel(sproc_rel,2,ipparbar+1).GT.real(pparbar_res)).and.&
			        (pparbar_rel(sproc_rel,2,ipparbar).LE.real(pparbar_res))) then
				          ipparbar_res=ipparbar
				          found_res = .TRUE.
		       endif
	    enddo
  endif

	! Handle resonances that are right outside the integration domain:

	do ipparbar=0,positions_principal
		if ((real(pparbar_res).GE.(pparbar_rel(sproc_rel,2,0)-dpparbar*ipparbar)).AND.&
		  (real(pparbar_res).LT.(pparbar_rel(sproc_rel,2,0)-dpparbar*(ipparbar-1)))) then
				ipparbar_res = -ipparbar
				found_res = .TRUE.
			endif

		if ((real(pparbar_res).GE.(pparbar_rel(sproc_rel,2,npparbar-1)+dpparbar*ipparbar)).AND.&
		  (real(pparbar_res).LT.(pparbar_rel(sproc_rel,2,npparbar-1)+dpparbar*(ipparbar+1)))) then
				 ipparbar_res = npparbar-1+ipparbar
				 found_res = .TRUE.
			endif
	enddo



	! determine the limits for the integration:

	!! What are the relevant ranges in pparbar:
	found_lower=.FALSE.
	found_upper=.FALSE.
	ipparbar_lower = 1
	ipparbar_upper = npparbar - 1
	do ipparbar=1,npparbar-1
		if((.NOT.found_lower).AND.(f0_rel(sproc_rel,igamma,ipparbar-1).LE.-1.d0).AND.(f0_rel(sproc_rel,igamma,ipparbar).GT.-1.d0)) then
			ipparbar_lower = ipparbar
			found_lower=.TRUE.
		endif
		if((.NOT.found_upper).AND.(f0_rel(sproc_rel,igamma,ipparbar).GT.-1.d0).AND.(f0_rel(sproc_rel,igamma,ipparbar+1).LE.-1.d0)) then
			ipparbar_upper = ipparbar
			found_upper=.TRUE.
		endif
	enddo


	! define int_start, int_end, lowerlimit, and upperlimit:
	if (found_res) then
		! These four are in general correct:
		int_start = ipparbar_lower
		int_end = ipparbar_upper
		lowerlimit = ipparbar_res - positions_principal
		upperlimit = ipparbar_res + positions_principal + 1

    if ((ipparbar_res.GE.0).and.(ipparbar_res.LE.npparbar)) then
			if(abs(real(pparbar_res)-pparbar_rel(sproc_rel,2,ipparbar_res)).GT.(0.5d0*dpparbar)) upperlimit = upperlimit + 1
		endif

    ! But there are special circumstances:
		if ((lowerlimit.LT.ipparbar_lower).AND.(upperlimit.GT.ipparbar_upper)) then
      call alps_error(8)
		elseif (lowerlimit.LE.ipparbar_lower) then	! resonance outside or near the left end of the subluminal cone
				int_start = 1
				lowerlimit = 0
				!upperlimit = max(upperlimit, ipparbar_lower)
        upperlimit = ipparbar_lower
			elseif (upperlimit.GE.ipparbar_upper) then ! resonance outside or near the right end of the subluminal cone
				!lowerlimit = min(lowerlimit, ipparbar_upper)
        lowerlimit = ipparbar_upper
				upperlimit = npparbar
				int_end = npparbar - 1
			endif



	else ! no resonance (only integrate from int_start to lowerlimit)
		int_start = ipparbar_lower
		lowerlimit = ipparbar_upper
		int_end = npparbar-1
		upperlimit = npparbar
	endif


	! Now comes the normal integration:

	if (int_start.LE.lowerlimit) then
		ipparbar = int_start
		integrate_resU_rel = integrate_resU_rel + resU_rel(sproc_rel,om,nn,igamma,ipparbar) &
			*  int_T_rel(sproc_rel,nn,igamma,ipparbar,mode)
	endif

	do ipparbar = int_start+1, lowerlimit-1
	  integrate_resU_rel = integrate_resU_rel + 2.d0 * resU_rel(sproc_rel,om,nn,igamma,ipparbar) &
		*  int_T_rel(sproc_rel,nn,igamma,ipparbar,mode)
	enddo

	if (int_start.LT.lowerlimit) then
		ipparbar = lowerlimit
		integrate_resU_rel = integrate_resU_rel + resU_rel(sproc_rel,om,nn,igamma,ipparbar) &
			*  int_T_rel(sproc_rel,nn,igamma,ipparbar,mode)
	endif

	if (upperlimit.LE.int_end) then
		ipparbar = upperlimit
		integrate_resU_rel = integrate_resU_rel + resU_rel(sproc_rel,om,nn,igamma,ipparbar) &
			*  int_T_rel(sproc_rel,nn,igamma,ipparbar,mode)
	endif

	do ipparbar = upperlimit+1, int_end-1
		integrate_resU_rel = integrate_resU_rel + 2.d0 * resU_rel(sproc_rel,om,nn,igamma,ipparbar) &
			* int_T_rel(sproc_rel,nn,igamma,ipparbar,mode)
	enddo

	if (upperlimit.LT.int_end) then
		ipparbar = int_end
		integrate_resU_rel = integrate_resU_rel + resU_rel(sproc_rel,om,nn,igamma,ipparbar) &
			* int_T_rel(sproc_rel,nn,igamma,ipparbar,mode)
	endif

	integrate_resU_rel = integrate_resU_rel * dpparbar


	! And this is the resonance integration:
	if (found_res.AND.(lowerlimit.GE.int_start).AND.(upperlimit.LE.int_end)) &
		integrate_resU_rel = integrate_resU_rel + principal_integral_rel(sproc_rel,om,nn,mode,igamma,ipparbar_res,upperlimit)

	return

end function integrate_resU_rel






double complex function principal_integral_rel(sproc_rel,om,nn,mode,igamma,ipparbar_res,upperlimit)
	use alps_var, only : positions_principal,n_resonance_interval, sproc
	use alps_var, only : gamma_rel, pparbar_rel,vA,kpar,qs,ms, Tlim, pi
	implicit none
	double complex :: om, ii, pparbar_res, gprimetr
	double precision :: denomR,denomI,capDelta,smdelta,correction,pparbar, dpparbar
	integer :: nn, mode, ipparbar_res,ipparbar,igamma,upperlimit, sproc_rel, ntiny


	ii = cmplx(0.d0,1.d0,kind(1.d0))

	principal_integral_rel=cmplx(0.d0,0.d0,kind(1.d0))

    dpparbar = pparbar_rel(sproc_rel,2,2)-pparbar_rel(sproc_rel,2,1)
	! split the range between the resonance and the upper limit into n_resonance_interval steps:

	! Now comes the resonance part:
	! We call the function that needs to be integrated WITHOUT the resonance part funct_g.
	! We linearize this function. Now we can calculate the even part of the integration.
	! We set Delta so that it starts at ipparbar_res-positions_principal. In that way, there is only
	! a tiny rest left on the right side that needs to be integrated.
	! split the range between the resonance and the upper limit into n_resonance_interval steps:
	! the denominator is:
	! (pparbar - gamma * om/kpar + (1.d0*nn) * qs(sproc)/ms(sproc)) * vA /kpar
	! we define
	! denomR=real(gamma *om/kpar - (1.d0*nn) * qs(sproc)/ms(sproc))vA /kpar )
	! denomI=aimag(gamma * ms(sproc)*om/kpar)
	! so that the denominator is
	! (ppar-denomR-ii*denomI)


	denomR=real( gamma_rel(sproc_rel,igamma,ipparbar_res) * om * vA / kpar - (1.d0*nn) * (qs(sproc)/ms(sproc)) * vA / kpar )
	denomI=aimag(gamma_rel(sproc_rel,igamma,ipparbar_res) * om * vA / kpar)

	pparbar_res = denomR+denomI*ii

	capDelta=real(pparbar_res)-pparbar_rel(sproc_rel,1,ipparbar_res-positions_principal)
	smdelta=capDelta/(1.d0*n_resonance_interval)


	if (abs(denomI).GT.Tlim) then ! regular integration:
		! Integrate the boundaries:
		pparbar=real(pparbar_res)

		principal_integral_rel=principal_integral_rel + 1.d0 * &
			funct_g_rel(sproc_rel,pparbar,igamma,om,nn,mode)/(pparbar-denomR-ii*denomI)

		principal_integral_rel=principal_integral_rel - 1.d0 * &
			funct_g_rel(sproc_rel,2.d0*denomR-pparbar,igamma,om,nn,mode)/(pparbar-denomR+ii*denomI)

		pparbar=real(pparbar_res)+capDelta

		principal_integral_rel=principal_integral_rel + 1.d0 * &
			funct_g_rel(sproc_rel,pparbar,igamma,om,nn,mode)/(pparbar-denomR-ii*denomI)

		principal_integral_rel=principal_integral_rel - 1.d0 * &
			funct_g_rel(sproc_rel,2.d0*denomR-pparbar,igamma,om,nn,mode)/(pparbar-denomR+ii*denomI)

		do ipparbar = 1, n_resonance_interval-1
			pparbar=real(pparbar_res)+smdelta*ipparbar

			principal_integral_rel=principal_integral_rel + 2.d0 * &
				funct_g_rel(sproc_rel,pparbar,igamma,om,nn,mode)/(pparbar-denomR-ii*denomI)

			principal_integral_rel=principal_integral_rel - 2.d0 * &
				funct_g_rel(sproc_rel,2.d0*denomR-pparbar,igamma,om,nn,mode)/(pparbar-denomR+ii*denomI)

		enddo


	else ! analytical approximation

		gprimetr = (funct_g_rel(sproc_rel,denomR+dpparbar,igamma,om,nn,mode)&
			-funct_g_rel(sproc_rel,denomR-dpparbar,igamma,om,nn,mode))/(2.d0*dpparbar)

		! Integrate the edges:
		pparbar=real(pparbar_res)
		if (denomI.NE.0.d0) then
			principal_integral_rel=principal_integral_rel &
				+ 2.d0 * gprimetr*(pparbar-denomR)**2 / ((pparbar-denomR)**2+denomI**2)
		else
			principal_integral_rel=principal_integral_rel + 2.d0 * gprimetr
		endif


		pparbar=real(pparbar_res)+capDelta
		principal_integral_rel=principal_integral_rel &
			+ 2.d0 * gprimetr*(pparbar-denomR)**2 / ((pparbar-denomR)**2+denomI**2)

		if (denomI.NE.0.d0) then ! the factor 2 is for normalization reasons in the trapezoidal rule
			principal_integral_rel=principal_integral_rel &
				+ 2.d0 * ii * pi * funct_g_rel(sproc_rel,denomR,igamma,om,nn,mode) * denomI/(abs(denomI)*smdelta)
		endif
		! end of edges.

		do ipparbar = 1, n_resonance_interval-1
			pparbar=real(pparbar_res)+smdelta*ipparbar

			principal_integral_rel=principal_integral_rel &
				+ 2.d0 * 2.d0 * gprimetr*(pparbar-denomR)**2 / ((pparbar-denomR)**2+denomI**2)

		enddo

	endif







	! There is a tiny rest left between the point real(pparbar_res)+capDelta and the position
	! pparbar_rel(sproc_rel,iperp,upperlimit). We split this interval into steps of roughly size smdelta:

	ntiny=int((pparbar_rel(sproc_rel,igamma,upperlimit)-real(pparbar_res)-capDelta)/smdelta)


	if (ntiny.GT.0) then

		! Correct for the fact that smdelta is not exactly the step width in the tiny-rest integration:
		correction=((pparbar_rel(sproc_rel,igamma,upperlimit)-real(pparbar_res)-capDelta)/(1.d0*ntiny))/smdelta


		pparbar=real(pparbar_res)+capDelta

		principal_integral_rel=principal_integral_rel + 1.d0*&
		correction*(funct_g_rel(sproc_rel,pparbar,igamma,om,nn,mode)/(pparbar-denomR-ii*denomI))


		pparbar=real(pparbar_res)+capDelta+correction*smdelta*ntiny

		principal_integral_rel=principal_integral_rel + 1.d0*&
		correction*(funct_g_rel(sproc_rel,pparbar,igamma,om,nn,mode)/(pparbar-denomR-ii*denomI))

		do ipparbar=1,ntiny-1
			pparbar=real(pparbar_res)+capDelta+correction*smdelta*ipparbar

			principal_integral_rel=principal_integral_rel + 2.d0*&
			correction*(funct_g_rel(sproc_rel,pparbar,igamma,om,nn,mode)/(pparbar-denomR-ii*denomI))

		enddo
	endif


	principal_integral_rel = principal_integral_rel * smdelta

	return

end function



! Linearized integrand WITHOUT the resonance part.
! It is  resU * int_T * Omega but without the denominator.
! It can be evaluated at any real value of pp (ppar_real) between the grid around the resonance.
double complex function funct_g_rel(sproc_rel,pparbar,igamma,om,nn,mode)
	use alps_var,only : ms,qs,kpar,df0_rel,sproc,vA,npparbar,pparbar_rel,pi,f0_rel
	implicit none
	integer :: nn, mode, ipparbar,ipparbar_close,sproc_rel,igamma
	double complex :: om,integrandplus,integrandminus,integrand
	double precision :: pparbar,dpparbar

	dpparbar = pparbar_rel(sproc_rel,2,2)-pparbar_rel(sproc_rel,2,1)

	ipparbar_close=-2
	! determine the closest ipar (on the left) to this p_res_real:
	do ipparbar=0,npparbar-1
		if ((pparbar_rel(sproc_rel,igamma,ipparbar+1).GT.pparbar).AND.(pparbar_rel(sproc_rel,igamma,ipparbar).LE.pparbar)) then
			ipparbar_close=ipparbar
		endif
	enddo


	if (f0_rel(sproc_rel,igamma,ipparbar_close+1).LE.-1.d0) ipparbar_close = ipparbar_close - 1
	if (f0_rel(sproc_rel,igamma,ipparbar_close-1).LE.-1.d0) ipparbar_close = ipparbar_close + 1

	if (pparbar.EQ.pparbar_rel(sproc_rel,igamma,npparbar)) ipparbar_close = npparbar-2
	if (ipparbar_close.GE.(npparbar-1)) ipparbar_close=npparbar-2
	if (ipparbar_close.LE.1) ipparbar_close=2

	! calculate the function on the grid (left and right of pparbar):

	integrandplus = -2.d0*pi*(ms(sproc)/vA)**3 * (qs(sproc)*vA/(kpar*ms(sproc))) * &
		(om*df0_rel(sproc_rel,igamma,ipparbar_close+1,1) + &
			 (kpar/(vA))*df0_rel(sproc_rel,igamma,ipparbar_close+1,2))&
			 *int_T_rel(sproc_rel,nn, igamma, ipparbar_close+1, mode)

	integrand = -2.d0*pi*(ms(sproc)/vA)**3 * (qs(sproc)*vA/(kpar*ms(sproc))) * &
		(om*df0_rel(sproc_rel,igamma,ipparbar_close,1) + &
			 (kpar/(vA))*df0_rel(sproc_rel,igamma,ipparbar_close,2))&
			 *int_T_rel(sproc_rel,nn, igamma, ipparbar_close, mode)

	integrandminus = -2.d0*pi*(ms(sproc)/vA)**3 * (qs(sproc)*vA/(kpar*ms(sproc))) * &
		(om*df0_rel(sproc_rel,igamma,ipparbar_close-1,1) + &
			 (kpar/(vA))*df0_rel(sproc_rel,igamma,ipparbar_close-1,2))&
			 *int_T_rel(sproc_rel,nn, igamma, ipparbar_close-1, mode)

	funct_g_rel = integrand+ &
		0.5d0*((integrandplus-integrandminus)/dpparbar) * (pparbar - pparbar_rel(sproc_rel,igamma,ipparbar_close))

	return

end function funct_g_rel





double complex function landau_integrate_rel(om, nn, mode)
	use alps_var, only : ngamma, gamma_rel,pparbar_rel, pi, ms, qs, kpar, sproc, vA
	use alps_analyt, only: eval_fit
	implicit none
	!Passed
	double complex :: om   !complex frequency
	integer :: nn          !Bessel N
	integer :: mode        !index in T tensor
	!Local
	integer :: igamma, sproc_rel
	double precision :: dgamma,dpparbar
	double precision :: h
	double complex :: ii
	double complex :: dfgamma_C,dfpparbar_C,pparbar_res

	ii = cmplx(0.d0,1.d0,kind(1.d0))

	landau_integrate_rel = cmplx(0.d0,0.d0,kind(1.d0))


	call determine_sproc_rel(sproc_rel)


	dgamma = gamma_rel(sproc_rel,2,2)-gamma_rel(sproc_rel,1,2)
	dpparbar = pparbar_rel(sproc_rel,2,2)-pparbar_rel(sproc_rel,2,1)

	! Landau contour integral:
	do igamma = 1, ngamma-1

		pparbar_res = gamma_rel(sproc_rel,igamma,1)*om*vA/kpar - (1.d0*nn)*qs(sproc)*vA/(kpar*ms(sproc))


		if ((real(pparbar_res)**2).LE.(gamma_rel(sproc_rel,igamma,1)**2-1.d0)) then

			h=1.d0

      ! At igamma=1, we are already missing the part from igamma=0, where we should actually start.
      ! Therefore, we use 4 instead of 2 in the trapezoid integration:
      if (igamma.EQ.(ngamma-1)) h = 0.5d0

      if (igamma.EQ.1) then
        dfgamma_C=(eval_fit(sproc,igamma+1,pparbar_res)-eval_fit(sproc,igamma,pparbar_res))/dgamma
      else
          dfgamma_C=(eval_fit(sproc,igamma+1,pparbar_res)-eval_fit(sproc,igamma-1,pparbar_res))/(2.d0*dgamma)
      endif

      dfpparbar_C=(eval_fit(sproc,igamma,pparbar_res+dpparbar)-eval_fit(sproc,igamma,pparbar_res-dpparbar))/(2.d0*dpparbar)

			landau_integrate_rel = landau_integrate_rel - h * ( &
				om*dfgamma_C+(kpar/(vA))*dfpparbar_C) * int_T_res_rel(sproc_rel,nn, igamma,pparbar_res, mode)

      endif

	enddo

	landau_integrate_rel = landau_integrate_rel * ii * dgamma * pi * 2.d0 * pi &
		* (qs(sproc)*vA/(kpar*ms(sproc))) * (ms(sproc)/vA)**3

	return

end function landau_integrate_rel




double complex function int_ee_rel(om)
	use alps_var, only : qs, ms,pi, f0_rel,df0_rel, vA, sproc, gamma_rel, pparbar_rel
	use alps_var, only : ngamma, npparbar
	implicit none
	double complex :: om   !complex frequency
	integer :: sproc_rel, igamma, ipparbar
	integer :: ipparbar_lower, ipparbar_upper
	double precision :: dgamma, dpparbar
	logical :: found_lower, found_upper

	call determine_sproc_rel(sproc_rel)

	dgamma = gamma_rel(sproc_rel,2,2)-gamma_rel(sproc_rel,1,2)
	dpparbar = pparbar_rel(sproc_rel,2,2)-pparbar_rel(sproc_rel,2,1)


	int_ee_rel = cmplx (0.d0, 0.d0,kind(1.d0))

	!! What are the relevant ranges in pparbar:
	igamma = ngamma-1

	found_lower=.FALSE.
	found_upper=.FALSE.
	ipparbar_lower = 1
	ipparbar_upper = npparbar - 1
	do ipparbar=1,npparbar-1
		if((.NOT.found_lower).AND.(f0_rel(sproc_rel,igamma,ipparbar-1).LE.-1.d0).AND.(f0_rel(sproc_rel,igamma,ipparbar).GT.-1.d0)) then
			ipparbar_lower = ipparbar
			found_lower=.TRUE.
		endif
		if((.NOT.found_upper).AND.(f0_rel(sproc_rel,igamma,ipparbar).GT.-1.d0).AND.(f0_rel(sproc_rel,igamma,ipparbar+1).LE.-1.d0)) then
			ipparbar_upper = ipparbar
			found_upper=.TRUE.
		endif
	enddo


	int_ee_rel = int_ee_rel + &
		pparbar_rel(sproc_rel,igamma,ipparbar_lower)*df0_rel(sproc_rel,igamma,ipparbar_lower,2)

	int_ee_rel = int_ee_rel + &
		pparbar_rel(sproc_rel,igamma,ipparbar_upper)*df0_rel(sproc_rel,igamma,ipparbar_upper,2)



	do ipparbar = ipparbar_lower+1,ipparbar_upper-1
			int_ee_rel = int_ee_rel + &
				2.d0 * pparbar_rel(sproc_rel,igamma,ipparbar)*df0_rel(sproc_rel,igamma,ipparbar,2)
	enddo


  ! At igamma=1, we are already missing the part from igamma=0, where we should actually start.
  ! Therefore, we use 4 instead of 2 in the trapezoid integration:
	do igamma = 1, ngamma-2

	found_lower=.FALSE.
	found_upper=.FALSE.
	ipparbar_lower = 1
	ipparbar_upper = npparbar - 1
		do ipparbar=1,npparbar-1
			if((.NOT.found_lower).AND.(f0_rel(sproc_rel,igamma,ipparbar-1).LE.-1.d0).AND.&
				(f0_rel(sproc_rel,igamma,ipparbar).GT.-1.d0)) then
					ipparbar_lower = ipparbar
					found_lower=.TRUE.
			endif
			if((.NOT.found_upper).AND.(f0_rel(sproc_rel,igamma,ipparbar).GT.-1.d0).AND.&
				(f0_rel(sproc_rel,igamma,ipparbar+1).LE.-1.d0)) then
					ipparbar_upper = ipparbar
					found_upper=.TRUE.
			endif
		enddo


	   do ipparbar = ipparbar_lower+1, ipparbar_upper-1
		int_ee_rel = int_ee_rel + &
			4.d0 * pparbar_rel(sproc_rel,igamma,ipparbar)*df0_rel(sproc_rel,igamma,ipparbar,2)
	   enddo

	   int_ee_rel = int_ee_rel + &
			2.d0 * pparbar_rel(sproc_rel,igamma,ipparbar_lower)*df0_rel(sproc_rel,igamma,ipparbar_lower,2)
	   int_ee_rel = int_ee_rel + &
			2.d0 * pparbar_rel(sproc_rel,igamma,ipparbar_upper)*df0_rel(sproc_rel,igamma,ipparbar_upper,2)
	enddo


	int_ee_rel = int_ee_rel * 2.d0 * pi * qs(sproc) / (ms(sproc))
	int_ee_rel = int_ee_rel * dgamma * dpparbar * 0.25d0 * (ms(sproc)/vA)**3


	return

end function int_ee_rel



!-=-=-=-=-=-=
!Functions for resonant term in integral
!-=-=-=-=-=-=

!-=-=-=-=-=-=
!Functions for resonant term in integral for the relativistic calculation:
!-=-=-=-=-=-=
double complex function resU_rel(sproc_rel,om, nn, igamma, ipparbar)
	use ALPS_var, only : kpar, ms, qs, df0_rel, vA, sproc, gamma_rel
	use ALPS_var, only : pi, pparbar_rel
	implicit none
	!Passed
	integer :: nn          !Bessel N
	integer :: sproc_rel,igamma,ipparbar
	double complex :: om   !complex frequency
	!Local

	resU_rel = -2.d0 * pi * (ms(sproc)/vA)**3 * (qs(sproc)*vA/(kpar*ms(sproc))) * &
		(om*df0_rel(sproc_rel,igamma,ipparbar,1) + (kpar/(vA))*df0_rel(sproc_rel,igamma,ipparbar,2)) / &
		(pparbar_rel(sproc_rel,igamma,ipparbar)-gamma_rel(sproc_rel,igamma,ipparbar)*om*vA/kpar &
		+ (1.d0*nn)*qs(sproc)*vA/(kpar*ms(sproc)))

	return

end function resU_rel

!-=-=-=-=-=-=
!Functions for Tij-
!-=-=-=-=-=-=

!Function to pass T_ij into integrator
double complex function int_T_rel(sproc_rel,nn, igamma, ipparbar, mode)
	use ALPS_var, only : kperp, qs, sproc, pparbar_rel, gamma_rel, vA, ms
	implicit none
	!Passed
	integer :: nn, igamma, ipparbar , sproc_rel
	integer :: mode        !index in T tensor
	!Local
	double precision :: z,zbar  !Bessel Argument
	double precision :: pperpbar
	double precision :: bessel 	! Bessel function for nn and z
	double precision :: besselP	! first derivative of Bessel function for nn and z
	double complex :: ii = cmplx(0.d0,1.d0,kind(1.d0))



	!Bessel Fn Argument
	pperpbar = sqrt(gamma_rel(sproc_rel,igamma,ipparbar)**2-1.d0-pparbar_rel(sproc_rel,igamma,ipparbar)**2)
	z= (kperp*ms(sproc)/(vA*qs(sproc)))*pperpbar
	zbar = kperp*ms(sproc)/(vA*qs(sproc))

	! Look up array of Bessel functions:
	if (nn.LT.0) then
	   bessel=((-1.d0)**nn)*BESSJ(-nn,z)
	else
	   bessel=BESSJ(nn,z)
	endif

	! determine derivative of Bessel function:
	if (nn.GE.1) then
		besselP = 0.5d0 * (BESSJ(nn-1,z) - BESSJ(nn+1,z))
	else if (nn.LT.-1) then
		besselP = 0.5d0 * ((((-1.d0)**(nn-1))*BESSJ(-(nn-1),z))&
			-(((-1.d0)**(nn+1))*BESSJ(-(nn+1),z)))
	else if (nn.EQ.0) then
	   besselP = -BESSJ(1,z)
	else if (nn.EQ.-1) then
		besselP = 0.5d0 * (BESSJ(2,z) - BESSJ(0,z))
	endif


	select case(mode)

	  case(1) !T xx
		 int_T_rel = 1.d0 * (nn * nn) * bessel * bessel / (zbar * zbar)

	  case(2) !T yy
		 int_T_rel = besselP * besselP * pperpbar * pperpbar

	  case(3) !T zz
		 int_T_rel = bessel * bessel * pparbar_rel( sproc_rel, igamma, ipparbar)**2

	  case(4) !T xy
		 int_T_rel = ii*(1.d0 * (nn)) * bessel * besselP * pperpbar / zbar

	  case(5) !T xz
		 int_T_rel = (1.d0 * nn) * bessel * bessel* pparbar_rel(sproc_rel,igamma,ipparbar)  / zbar

	  case(6) !T yz
		 int_T_rel = (-1.d0 * ii) * bessel * besselP * pparbar_rel(sproc_rel,igamma,ipparbar) * pperpbar

	end select

	return

end function int_T_rel



!Function to pass T_ij into integrator
double complex function int_T_res_rel(sproc_rel,nn, igamma, pparbar, mode)
	use ALPS_var, only : kperp, qs, sproc, gamma_rel, vA, ms
	implicit none
	!Passed
	integer :: nn, igamma,sproc_rel
	integer :: mode        !index in T tensor
	!Local
	double precision :: zbar
	double complex :: z
	double complex :: pperpbar
	double complex :: bessel 	! Bessel function for nn and z
	double complex :: besselP	! first derivative of Bessel function for nn and z
	double complex :: besselH ! Help variable
	double complex :: ii = cmplx(0.d0,1.d0,kind(1.d0))
	double complex :: pparbar


	!Bessel Fn Argument
	pperpbar = sqrt(gamma_rel(sproc_rel,igamma,1)**2-1.d0-pparbar**2)
	z= (kperp*ms(sproc)/(vA*qs(sproc)))*pperpbar
	zbar = kperp*ms(sproc)/(vA*qs(sproc))

	! Look up array of Bessel functions:
	if (nn.LT.0) then
	   call CBESSJ(z,-nn,bessel)
	   bessel=bessel*(-1.d0)**nn
	else
	   call CBESSJ(z, nn, bessel)
	endif

	! determine derivative of Bessel function:
	if (nn.GE.1) then
		call CBESSJ(z,nn-1,besselP)
		call CBESSJ(z,nn+1,besselH)
		besselP= 0.5d0*(besselP - besselH)
	elseif (nn.LT.-1) then
		call CBESSJ(z,-(nn-1),besselP)
		call CBESSJ(z,-(nn+1),besselH)
		besselP = 0.5d0 * ((((-1.d0)**(nn-1))*besselP)&
			-(((-1.d0)**(nn+1))*besselH))
	elseif (nn.EQ.0) then
	  	call CBESSJ(z,1,besselP)
		besselP = -besselP
	elseif (nn.EQ.-1) then
		call CBESSJ(z,2,besselP)
		call CBESSJ(z,0,besselH)
		besselP = 0.5d0 * (besselP - besselH)
	endif


	select case(mode)

	  case(1) !T xx
		 int_T_res_rel = 1.d0 * (nn * nn) * bessel * bessel / (zbar * zbar)

	  case(2) !T yy
		 int_T_res_rel = besselP * besselP * pperpbar * pperpbar

	  case(3) !T zz
		 int_T_res_rel = bessel * bessel * pparbar**2

	  case(4) !T xy
		 int_T_res_rel = ii*(1.d0 * (nn)) * bessel * besselP * pperpbar / zbar

	  case(5) !T xz
		 int_T_res_rel = (1.d0 * nn) * bessel * bessel* pparbar  / zbar

	  case(6) !T yz
		 int_T_res_rel = (-1.d0 * ii) * bessel * besselP * pparbar * pperpbar
	end select
	return

end function int_T_res_rel





! Calculate the complex Bessel function
subroutine CBESSJ(z, nu, z1)
!---------------------------------------------------
!                       inf.     (-z^2/4)^k
!   Jnu(z) = (z/2)^nu x Sum  ------------------
!                       k=0  k! x Gamma(nu+k+1)
!  (nu must be >= 0).
!---------------------------------------------------
implicit none

  double complex :: z,z1
  integer :: k, MAXK, nu
  double complex :: sum,tmp
  double precision :: ZERO
  parameter(MAXK=20,ZERO=0.d0)
  sum = cmplx(0.d0,0.d0,kind(1.d0))
  do k=0, MAXK
    !calculate (-z**2/4)**k
	tmp = (-z*z/4.d0)**k
    !divide by k!
	tmp = tmp / Fact(k)
    !divide by Gamma(nu+k+1)
    tmp = tmp / Gamma(1.d0*(nu+k+1))
    !actualize sum
	sum = sum + tmp
  end do
  !calculate (z/2)**nu
  tmp = (z/2.d0)**nu
  !multiply (z/2)**nu by sum
  z1 = tmp*sum
return
end subroutine


double precision function Fact(K)
implicit none
  integer :: i,k
  double precision :: f
    f=1.d0
    do i=2, k
	  f=f*(1.d0*i)
    end do
    Fact=f
  return
end function





DOUBLE PRECISION FUNCTION BESSJ (N,X)
!double precision FUNCTION BESSJ (N,X)

!     This subroutine calculates the first kind Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows.
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.

      IMPLICIT NONE
      INTEGER, PARAMETER :: IACC = 40
      REAL*8, PARAMETER :: BIGNO = 1.D10, BIGNI = 1.D-10
      INTEGER M, N, J, JSUM
      !REAL *8 X,BESSJ,TOX,BJM,BJ,BJP,SUM
      REAL *8 X,TOX,BJM,BJ,BJP,SUM
      !REAL *8 X,BESSJ0,BESSJ1,TOX,BJM,BJ,BJP,SUM

!      REAL *8 X,BESSJ,BESSJ0,BESSJ1,TOX,BJM,BJ,BJP,SUM
      !write(*,*)'z: ', x
      IF (N.EQ.0) THEN
      BESSJ = BESSJ0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSJ = BESSJ1(X)
      RETURN
      ENDIF
      IF (X.EQ.0.) THEN
      BESSJ = 0.
      RETURN
      ENDIF
      TOX = 2./X
      IF (X.GT.FLOAT(N)) THEN
      BJM = BESSJ0(X)
      BJ  = BESSJ1(X)
      DO 11 J = 1,N-1
      BJP = J*TOX*BJ-BJM
      BJM = BJ
      BJ  = BJP
   11 CONTINUE
      BESSJ = BJ
      ELSE
      M = 2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
      BESSJ = 0.
      JSUM = 0
      SUM = 0.
      BJP = 0.
      BJ  = 1.
      DO 12 J = M,1,-1
      BJM = J*TOX*BJ-BJP
      BJP = BJ
      BJ  = BJM
      IF (ABS(BJ).GT.BIGNO) THEN
      BJ  = BJ*BIGNI
      BJP = BJP*BIGNI
      BESSJ = BESSJ*BIGNI
      SUM = SUM*BIGNI
      ENDIF
      IF (JSUM.NE.0) SUM = SUM+BJ
      JSUM = 1-JSUM
      IF (J.EQ.N) BESSJ = BJP
   12 CONTINUE
      SUM = 2.*SUM-BJ
      BESSJ = BESSJ/SUM
      ENDIF
      RETURN
      END FUNCTION



      !double precision FUNCTION BESSJ0 (X)
      DOUBLE PRECISION FUNCTION BESSJ0 (X)
      IMPLICIT NONE
      REAL *8 X,AX,FR,FS,Z,FP,FQ,XX
      !REAL *8 X,BESSJ0,AX,FR,FS,Z,FP,FQ,XX

!     This subroutine calculates the First Kind Bessel Function of
!     order 0, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.

      REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
      -.2073370639D-5,.2093887211D-6 /
      DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
      DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
      651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
      DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
      9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /

      IF(X.EQ.0.D0) GO TO 1
      AX = ABS (X)
      IF (AX.LT.8.D0) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ0 = FR/FS
      ELSE
      Z = 8.D0/AX
      Y = Z*Z
      XX = AX-.785398164D0
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ0 = SQRT(.636619772D0/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
      ENDIF
      RETURN
1     BESSJ0 = 1.D0
      RETURN
    END FUNCTION BESSJ0
! ---------------------------------------------------------------------------
    !double precision FUNCTION BESSJ1 (X)
     double precision FUNCTION BESSJ1 (X)
      IMPLICIT NONE
      !REAL *8 X,BESSJ1,AX,FR,FS,Z,FP,FQ,XX
      REAL *8 X,AX,FR,FS,Z,FP,FQ,XX
!     This subroutine calculates the First Kind Bessel Function of
!     order 1, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
      REAL *8 Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
      .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
      DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
      .8449199096D-5,-.88228987D-6,.105787412D-6 /
      DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, &
      242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
      DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
      18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /

      AX = ABS(X)
      IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ1 = X*(FR/FS)
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-2.35619491
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
      ENDIF
      RETURN
    END FUNCTION BESSJ1


!*******************************************
!*           FUNCTION  GAMMA(X)            *
!* --------------------------------------- *
!* Returns the value of Gamma(x) in double *
!* precision as EXP(LN(GAMMA(X))) for X>0. *
!*******************************************
double precision Function Gamma(xx)
implicit none

  double precision :: xx, ONE, FPF, HALF
  double precision :: cof(6)
  double precision :: stp,x,tmp,ser
  integer :: j

parameter(ONE=1.d0,FPF=5.5d0,HALF=0.5d0)
  cof(1)=76.18009173d0
  cof(2)=-86.50532033d0
  cof(3)=24.01409822d0
  cof(4)=-1.231739516d0
  cof(5)=0.120858003d-2
  cof(6)=-0.536382d-5
  stp=2.50662827465d0

  x=xx-ONE
  tmp=x+FPF
  tmp=(x+HALF)*LOG(tmp)-tmp
  ser=ONE
  do j=1, 6
    x=x+ONE
    ser=ser+cof(j)/x
  end do
  Gamma = EXP(tmp+LOG(stp*ser))
return
end function



end module alps_fns_rel
