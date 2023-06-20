!Copyright (c) 2023, Kristopher G. Klein and Daniel Verscharen
!All rights reserved.
!
!This source code is licensed under the BSD-style license found in the
!LICENSE file in the root directory of this source tree.
!
!===============================================================================
!I                                                                             I
!I                              A  L  P  S                                     I
!I                     Arbitrary Linear Plasma Solver                          I
!I                                                                             I
!I                              Version 1.0                                    I
!I                                                                             I
!I  Kristopher Klein   (kgklein@arizona.edu)                                   I
!I  Daniel Verscharen  (d.verscharen@ucl.ac.uk)                                I
!I                                                                             I
!===============================================================================

module alps_analyt
	!! This module contains functions and subroutines for the hybrid analytical continuation.
	implicit none

	public  :: eval_fit, determine_param_fit
	private :: fit_function, output_fit, determine_JT, LM_nonlinear_fit

contains



double complex function eval_fit(is,iperp,ppar_valC)
	!! This function evaluates the fit to f0 at and the complex parallel
	!! momentum [[ppar_valC(variable):eval_fit(function)]]. It requires the fit parameters that will be determined
	!! by the subroutine [[determine_param_fit(subroutine)]].

	use alps_var, only : fit_type, pp, param_fit, n_fits, gamma_rel,nspec,relativistic
	use alps_distribution_analyt, only : distribution_analyt
	implicit none
	integer :: is,iperp,par_ind,ifit,n_params,is_rel,sproc_rel,is_run
	double precision :: pperp_val
	double precision,allocatable, dimension(:) :: params
	double complex :: ppar_valC



	! Use the pre-coded distribution from distribution/distribution_analyt.f90
	if (n_fits(is).EQ.0) then
		eval_fit=distribution_analyt(is,pp(is,iperp,1,1),ppar_valC)
		return
	endif



	n_params=0	! total number of fit_parameters
	do ifit=1,n_fits(is)
		if (fit_type(is,ifit).EQ.1) n_params=n_params+3	! Maxwell
		if (fit_type(is,ifit).EQ.2) n_params=n_params+5	! kappa
		if (fit_type(is,ifit).EQ.3) n_params=n_params+3	! Juettner with pperp and ppar
		if (fit_type(is,ifit).EQ.4) n_params=n_params+1	! Juettner with gamma and pparbar, constant in pparbar
		if (fit_type(is,ifit).EQ.5) n_params=n_params+3	! Juettner with gamma and pparbar with pparbar-dependence
		if (fit_type(is,ifit).EQ.6) n_params=n_params+4	! Bi-Moyal distribution
	enddo

	allocate (params(n_params))


	par_ind=0
	do ifit=1,n_fits(is)
		if (fit_type(is,ifit).EQ.1) then	! Maxwell
			pperp_val=pp(is,iperp,1,1)
			params(ifit+par_ind+0)=param_fit(is,iperp,1,ifit)
			params(ifit+par_ind+1)=param_fit(is,iperp,2,ifit)
			params(ifit+par_ind+2)=param_fit(is,iperp,3,ifit)
			par_ind=par_ind+2
		endif

		if (fit_type(is,ifit).EQ.2) then	! kappa
			pperp_val=pp(is,iperp,1,1)
			params(ifit+par_ind+0)=param_fit(is,iperp,1,ifit)
			params(ifit+par_ind+1)=param_fit(is,iperp,2,ifit)
			params(ifit+par_ind+2)=param_fit(is,iperp,3,ifit)
			params(ifit+par_ind+3)=param_fit(is,iperp,4,ifit)
			params(ifit+par_ind+4)=param_fit(is,iperp,5,ifit)
			par_ind=par_ind+4
		endif

		if (fit_type(is,ifit).EQ.3) then	! Juettner in pperp and ppar
			pperp_val=pp(is,iperp,1,1)
			params(ifit+par_ind+0)=param_fit(is,iperp,1,ifit)
			params(ifit+par_ind+1)=param_fit(is,iperp,2,ifit)
			params(ifit+par_ind+2)=param_fit(is,iperp,3,ifit)
			par_ind=par_ind+2
		endif

		if (fit_type(is,ifit).EQ.4) then	! Juettner in gamma and pparbar, constant in pparbar
			! Determine your sproc_rel:
 			 is_rel=0
			 do is_run=1,nspec
			 	if (relativistic(is_run)) then
		 			is_rel=is_rel+1
	 				if (is_run.EQ.is) sproc_rel=is_rel
		 		endif
			 enddo
			pperp_val=gamma_rel(sproc_rel,iperp,1)
			params(ifit+par_ind+0)=param_fit(is,iperp,1,ifit)
			par_ind=par_ind+0
		endif

		if (fit_type(is,ifit).EQ.5) then	! Juettner in gamma and pparbar with pparbar-dependence
			! Determine your sproc_rel:
 			 is_rel=0
			 do is_run=1,nspec
			 	if (relativistic(is_run)) then
		 			is_rel=is_rel+1
	 				if (is_run.EQ.is) sproc_rel=is_rel
		 		endif
			 enddo
			pperp_val=gamma_rel(sproc_rel,iperp,1)
			params(ifit+par_ind+0)=param_fit(is,iperp,1,ifit)
			params(ifit+par_ind+1)=param_fit(is,iperp,2,ifit)
			params(ifit+par_ind+2)=param_fit(is,iperp,3,ifit)
			par_ind=par_ind+2
		endif

		if (fit_type(is,ifit).EQ.6) then	! Bi-Moyal
			pperp_val=pp(is,iperp,1,1)
			params(ifit+par_ind+0)=param_fit(is,iperp,1,ifit)
			params(ifit+par_ind+1)=param_fit(is,iperp,2,ifit)
			params(ifit+par_ind+2)=param_fit(is,iperp,3,ifit)
			params(ifit+par_ind+3)=param_fit(is,iperp,4,ifit)
			par_ind=par_ind+3
		endif
	enddo

	eval_fit=fit_function(is,n_params,params,pperp_val,ppar_valC)

	deallocate (params)

	return
end function




!-=-=-=-=-=-=
! This function evaluates the fit to f0 at real pperp_val and complex ppar_val
! provided that the one-dimensional fit-parameter array params is fed into the
! function. This is only used during the fitting. For the evaluation in ALPS,
! use eval_fit.
!-=-=-=-=-=-=
double complex function fit_function(is,n_params,params,pperp_val,ppar_val)
	use alps_var, only : fit_type, n_fits, ms, vA, perp_correction
	implicit none
	integer :: n_params,ifit,is,par_ind
	double precision :: params(n_params),pperp_val
	double complex :: ppar_val,sqrtpart,kappapart

	fit_function=cmplx(0.d0,0.d0,kind(1.d0))

	par_ind=0
	do ifit=1,n_fits(is)
		if (fit_type(is,ifit).EQ.1) then	! Maxwell
			fit_function=fit_function+params(ifit+par_ind+0)*exp(-perp_correction(is,ifit)*pperp_val**2)&
				*exp(-params(ifit+par_ind+1)*(ppar_val-params(ifit+par_ind+2))**2)
			par_ind=par_ind+2
		endif

		if (fit_type(is,ifit).EQ.2) then	! kappa
			kappapart=1.d0+params(ifit+par_ind+1)*(ppar_val-params(ifit+par_ind+2))**2&
				+ perp_correction(is,ifit)*params(ifit+par_ind+4)*pperp_val**2
			fit_function=fit_function+params(ifit+par_ind+0)*kappapart**params(ifit+par_ind+3)
			par_ind=par_ind+4
		endif

		if (fit_type(is,ifit).EQ.3) then	! Juettner in pperp and ppar
			sqrtpart=sqrt(1.d0+(pperp_val**2+(ppar_val-params(ifit+par_ind+2))**2)*vA*vA/(ms(is)*ms(is)))
			fit_function=fit_function+params(ifit+par_ind+0)*exp(-params(ifit+par_ind+1)*sqrtpart)
			par_ind=par_ind+2
		endif

		if (fit_type(is,ifit).EQ.4) then ! Juettner in gamma only, constant in pparbar
			fit_function=fit_function+params(ifit+par_ind+0)*exp(-perp_correction(is,ifit)*pperp_val)
			par_ind=par_ind+0
		endif

		if (fit_type(is,ifit).EQ.5) then ! Juettner in gamma and pparbar with pparbar-dependence
			fit_function=fit_function+params(ifit+par_ind+0)*exp(-perp_correction(is,ifit)*pperp_val)* &
				exp(-params(ifit+par_ind+1)*(ppar_val-params(ifit+par_ind+2))**2)
			par_ind=par_ind+2
		endif

		if (fit_type(is,ifit).EQ.6) then	! Bi-Moyal
			fit_function=fit_function+params(ifit+par_ind+0)*exp(0.5d0*(params(ifit+par_ind+3)* &
					perp_correction(is,ifit)*pperp_val**2 + params(ifit+par_ind+1)* &
					(ppar_val-params(ifit+par_ind+2))**2 - &
					exp(params(ifit+par_ind+3)*perp_correction(is,ifit)*pperp_val**2 + &
					params(ifit+par_ind+1)*(ppar_val-params(ifit+par_ind+2))**2) ))
			par_ind=par_ind+3
		endif

	enddo
	return
end function



!-=-=-=-=-=-=
! This is the fitting routine for the hybrid analytic continuation. It determines
! the full field param_fit.
!-=-=-=-=-=-=
subroutine determine_param_fit
	use alps_var, only : writeOut, fit_type, param_fit, n_fits, nspec, f0, nperp, npar, logfit, runname
	use alps_var, only : relativistic,npparbar,f0_rel,ngamma, perp_correction, gamma_rel, usebM
	implicit none
	integer :: ifit,n_params,par_ind,iperp,is,is_run,is_rel,sproc_rel,nJT
	integer :: ipparbar,ipparbar_lower,ipparbar_upper,upperlimit,unit_spec
	logical :: found_lower, found_upper
	double precision :: quality,qualitytotal
	logical, allocatable, dimension (:) :: param_mask
	double precision,allocatable,dimension (:) :: g
	double precision,allocatable,dimension(:) :: params
	character (10) :: specwrite

	if (writeOut) then
		write (*,'(a)') 'Determine fit parameters for hybrid analytic continuation...'
	endif

	qualitytotal=0.d0

	do is=1,nspec

		if (usebM(is)) then
							write (*,'(a,i2)') ' Bi-Maxwellian calcuation: no fits necessary for species ',is
							param_fit(is,:,:,:)=0.d0
							fit_type(is,:)=1
							perp_correction(is,:)=1.d0

		 elseif (n_fits(is).EQ.0) then

			 	write (*,'(a,i2)') ' Using analytical function from distribution/distribution_analyt.f90: no fits necessary for species ',is

				param_fit(is,:,:,:)=0.d0
				fit_type(is,:)=1
				perp_correction(is,:)=1.d0


		 else
		! For all fit types that include a fit parameter for the perpendicular momentum (kappa and Moyal),
		! we must not fit this parameter when pperp=0. Otherwise, the LM matrix is singular:

					n_params=0	! total number of fit_parameters
					do ifit=1,n_fits(is)
						if (fit_type(is,ifit).EQ.1) n_params=n_params+3	! Maxwell
						if (fit_type(is,ifit).EQ.2) n_params=n_params+5	! kappa
						if (fit_type(is,ifit).EQ.3) n_params=n_params+3	! Juettner with pperp and ppar
						if (fit_type(is,ifit).EQ.4) n_params=n_params+1	! Juettner with gamma and pparbar, constant in pparbar
						if (fit_type(is,ifit).EQ.5) n_params=n_params+3	! Juettner with gamma and pparbar with pparbar-dependence
						if (fit_type(is,ifit).EQ.6) n_params=n_params+4	! Bi-Moyal
					enddo

			allocate (params(n_params))
			allocate (param_mask(n_params))


		if (relativistic(is)) then
			upperlimit=ngamma
		else
			upperlimit=nperp
		endif


	  unit_spec=2500+is
	  write(specwrite,'(i0)') is
	  open(unit = unit_spec,file = 'distribution/'//trim(runname)//'.fit_parameters.'//trim(specwrite)//'.out', status = 'replace')


		do iperp=0,upperlimit

			! Every step that is not iperp = 0 should use the previous result as a start value:
			if (iperp.NE.0) param_fit(is,iperp,:,:)=param_fit(is,iperp-1,:,:)

			par_ind=0
			nJT=0
			param_mask = .TRUE.

			do ifit=1,n_fits(is)
				if (fit_type(is,ifit).EQ.1) then	! Maxwell
					params(ifit+par_ind+0)=param_fit(is,iperp,1,ifit)
					params(ifit+par_ind+1)=param_fit(is,iperp,2,ifit)
					params(ifit+par_ind+2)=param_fit(is,iperp,3,ifit)
					par_ind=par_ind+2
				endif

				if (fit_type(is,ifit).EQ.2) then	! kappa
					params(ifit+par_ind+0)=param_fit(is,iperp,1,ifit)
					params(ifit+par_ind+1)=param_fit(is,iperp,2,ifit)
					params(ifit+par_ind+2)=param_fit(is,iperp,3,ifit)
					params(ifit+par_ind+3)=param_fit(is,iperp,4,ifit)
					params(ifit+par_ind+4)=param_fit(is,iperp,5,ifit)
					if (iperp.EQ.0) then
							nJT=nJT-1
							param_mask(ifit+par_ind+4)=.FALSE.
					endif
					par_ind=par_ind+4
				endif

				if (fit_type(is,ifit).EQ.3) then	! Juettner in pperp and ppar
					params(ifit+par_ind+0)=param_fit(is,iperp,1,ifit)
					params(ifit+par_ind+1)=param_fit(is,iperp,2,ifit)
					params(ifit+par_ind+2)=param_fit(is,iperp,3,ifit)
					par_ind=par_ind+2
				endif

				if (fit_type(is,ifit).EQ.4) then	! Juettner in gamma and pparbar, constant in pparbar
					params(ifit+par_ind+0)=param_fit(is,iperp,1,ifit)
					par_ind=par_ind+0
				endif

				if (fit_type(is,ifit).EQ.5) then	! Juettner in gamma and pparbar with pparbar-dependence
					params(ifit+par_ind+0)=param_fit(is,iperp,1,ifit)
					params(ifit+par_ind+1)=param_fit(is,iperp,2,ifit)
					params(ifit+par_ind+2)=param_fit(is,iperp,3,ifit)
					par_ind=par_ind+2
				endif

				if (fit_type(is,ifit).EQ.6) then	! Bi-Moyal
					params(ifit+par_ind+0)=param_fit(is,iperp,1,ifit)
					params(ifit+par_ind+1)=param_fit(is,iperp,2,ifit)
					params(ifit+par_ind+2)=param_fit(is,iperp,3,ifit)
					params(ifit+par_ind+3)=param_fit(is,iperp,4,ifit)
					if (iperp.EQ.0) then
							nJT=nJT-1
							param_mask(ifit+par_ind+3)=.FALSE.
					endif
					par_ind=par_ind+3
				endif

			enddo


			nJT=nJT+n_params


			! Fit and return everything in one array "params":
			if (relativistic(is)) then
				 ! Determine your sproc_rel:
				 is_rel=0
				 do is_run=1,nspec
					if (relativistic(is_run)) then
						is_rel=is_rel+1
						if (is_run.EQ.is) sproc_rel=is_rel
					endif
				 enddo

				 ! What are the relevant ranges in pparbar:
				found_lower=.FALSE.
				found_upper=.FALSE.
				ipparbar_lower = 1
				ipparbar_upper = npparbar - 1
				do ipparbar=1,npparbar-1
					if((.NOT.found_lower).AND.(f0_rel(sproc_rel,iperp,ipparbar-1).LE.-1.d0).AND.(f0_rel(sproc_rel,iperp,ipparbar).GT.-1.d0)) then
						ipparbar_lower = ipparbar
						found_lower=.TRUE.
					endif
					if((.NOT.found_upper).AND.(f0_rel(sproc_rel,iperp,ipparbar).GT.-1.d0).AND.(f0_rel(sproc_rel,iperp,ipparbar+1).LE.-1.d0)) then
						ipparbar_upper = ipparbar
						found_upper=.TRUE.
					endif
				enddo

				if ((ipparbar_upper-ipparbar_lower).GT.2) then

					allocate(g(0:ipparbar_upper-ipparbar_lower))


					if (logfit(is)) then
						g=log(f0_rel(sproc_rel,iperp,ipparbar_lower:ipparbar_upper))
					else
						g=f0_rel(sproc_rel,iperp,ipparbar_lower:ipparbar_upper)
					endif


					call LM_nonlinear_fit(is,g,n_params,nJT,params,param_mask,iperp,&
											(ipparbar_upper-ipparbar_lower),ipparbar_lower,quality)
					deallocate(g)

				else


					par_ind=0
					do ifit=1,n_fits(is)
						params(ifit+par_ind+0)=f0_rel(sproc_rel,iperp,(ipparbar_upper+ipparbar_lower)/2) /&
							exp(-perp_correction(is,ifit)*gamma_rel(sproc_rel,iperp,1))
						if (fit_type(is,ifit).EQ.5) then
							params(ifit+par_ind+1)=1.d-12
							params(ifit+par_ind+2)=0.d0
							par_ind=par_ind+2
						endif
					enddo

				endif

			else	! non-relativistic


				allocate(g(0:npar))

				if (logfit(is)) then
					g=log(f0(is,iperp,:))
				else
					g=f0(is,iperp,:)
				endif


				call LM_nonlinear_fit(is,g,n_params,nJT,params,param_mask,iperp,npar,0,quality)
				deallocate(g)

			endif

			qualitytotal=qualitytotal+quality

			! Write  Fit parameters to output files
			write (unit_spec,*) iperp,params


			! Fill it back into the param_fit field:
			par_ind=0
			do ifit=1,n_fits(is)
				if (fit_type(is,ifit).EQ.1) then	! Maxwell
					param_fit(is,iperp,1,ifit)=params(ifit+par_ind+0)
					param_fit(is,iperp,2,ifit)=params(ifit+par_ind+1)
					param_fit(is,iperp,3,ifit)=params(ifit+par_ind+2)
					par_ind=par_ind+2
				endif

				if (fit_type(is,ifit).EQ.2) then	! kappa
					param_fit(is,iperp,1,ifit)=params(ifit+par_ind+0)
					param_fit(is,iperp,2,ifit)=params(ifit+par_ind+1)
					param_fit(is,iperp,3,ifit)=params(ifit+par_ind+2)
					param_fit(is,iperp,4,ifit)=params(ifit+par_ind+3)
					param_fit(is,iperp,5,ifit)=params(ifit+par_ind+4)
					par_ind=par_ind+4
				endif

				if (fit_type(is,ifit).EQ.3) then	! Juettner in pperp and ppar
					param_fit(is,iperp,1,ifit)=params(ifit+par_ind+0)
					param_fit(is,iperp,2,ifit)=params(ifit+par_ind+1)
					param_fit(is,iperp,3,ifit)=params(ifit+par_ind+2)
					par_ind=par_ind+2
				endif

				if (fit_type(is,ifit).EQ.4) then	! Juettner in gamma and pparbar, constant in pparbar
					param_fit(is,iperp,1,ifit)=params(ifit+par_ind+0)
					par_ind=par_ind+0
				endif

				if (fit_type(is,ifit).EQ.5) then	! Juettner in gamma and pparbar with pparbar-dependence
					param_fit(is,iperp,1,ifit)=params(ifit+par_ind+0)
					param_fit(is,iperp,2,ifit)=params(ifit+par_ind+1)
					param_fit(is,iperp,3,ifit)=params(ifit+par_ind+2)
					par_ind=par_ind+2
				endif

				if (fit_type(is,ifit).EQ.6) then	! Bi-Moyal
					param_fit(is,iperp,1,ifit)=params(ifit+par_ind+0)
					param_fit(is,iperp,2,ifit)=params(ifit+par_ind+1)
					param_fit(is,iperp,3,ifit)=params(ifit+par_ind+2)
					param_fit(is,iperp,4,ifit)=params(ifit+par_ind+3)
					par_ind=par_ind+3
				endif

			enddo
		enddo	! End loop over iperp

		close (unit_spec)

		deallocate (params)
		deallocate (param_mask)

		endif
	enddo	! End loop over is





	if(writeOut) then
		call output_fit(qualitytotal)
		   write(*,'(a)') '-=-=-=-=-=-=-=-=-'
	endif

end subroutine



!-=-=-=-=-=-=
! This subroutine outputs the fit parameters for iperp=0 to stdout to monitor the fit
!-=-=-=-=-=-=
subroutine output_fit(qualitytotal)
	use alps_io, only : isnancheck, alps_error
	use alps_var, only : fit_type, param_fit, n_fits, nspec, nperp, npar, pp, f0, pi, vA, runname
	use alps_var, only : relativistic,gamma_rel,pparbar_rel,ngamma,npparbar,f0_rel, ms, usebM


	implicit none
	integer :: is,ifit,n_params,iparam,unit_spec,ipar,iperp, is_rel,sproc_rel,igamma,ipparbar,is_run
	double precision :: qualitytotal, integrate, dpperp, dppar, dgamma, dpparbar
	double complex :: ppar_comp
	character (10) :: specwrite

	write (*,'(a)') ' Results of the fit for hybrid analytic continuation at iperp = 1:'
	do is=1,nspec
		do ifit=1,n_fits(is)
			if(fit_type(is,ifit).EQ.1) n_params=3
			if(fit_type(is,ifit).EQ.2) n_params=5
			if(fit_type(is,ifit).EQ.3) n_params=3
			if(fit_type(is,ifit).EQ.4) n_params=1
			if(fit_type(is,ifit).EQ.5) n_params=3
			if(fit_type(is,ifit).EQ.6) n_params=4

			if (.not.usebM(is)) then
				do iparam=1,n_params
					write (*,'(a,i2,a,i2,a,i2,a,2es14.4)') '  param_fit(',is,', 1,',iparam,',',ifit,') = ',param_fit(is,1,iparam,ifit)
				enddo
			write (*,'(a)') ' '
			endif
		enddo
	enddo

	write (*,'(a,es14.4)') ' Sum of all least-squares: ', qualitytotal
	write (*,'(a,es14.4)') ' Standard error of the estimate: ', sqrt(qualitytotal/(1.d0*(nspec*nperp*npar)))


	if (isnancheck(qualitytotal)) call alps_error(7)

	write (*,'(a)') ' Writing fit result to files'
	do is=1,nspec

		if (.not.usebM(is)) then

	  unit_spec=2000+is
	  write(specwrite,'(i0)') is
	  open(unit = unit_spec,file = 'distribution/'//trim(runname)//'.fit_result.'//trim(specwrite)//'.out', status = 'replace')


	  integrate = 0.d0

	  if (relativistic(is)) then

		 ! Determine your sproc_rel:
		 is_rel=0
		 do is_run=1,nspec
				if (relativistic(is_run)) then
					is_rel=is_rel+1
					if (is_run.EQ.is) sproc_rel=is_rel
				endif
		 enddo

		dgamma = gamma_rel(sproc_rel,2,2)-gamma_rel(sproc_rel,1,2)
		dpparbar = pparbar_rel(sproc_rel,2,2)-pparbar_rel(sproc_rel,2,1)

		do igamma=0,ngamma
		 do ipparbar=0,npparbar
			ppar_comp=pparbar_rel(sproc_rel,igamma,ipparbar)
			if (f0_rel(sproc_rel,igamma,ipparbar).EQ.-1.d0) then
				write (unit_spec,*) gamma_rel(sproc_rel,igamma,ipparbar),pparbar_rel(sproc_rel,igamma,ipparbar),&
					"-1.0",	abs(-1.d0-f0_rel(sproc_rel,igamma,ipparbar))
			else

				write (unit_spec,*) gamma_rel(sproc_rel,igamma,ipparbar),pparbar_rel(sproc_rel,igamma,ipparbar),&
					real(eval_fit(is,igamma,ppar_comp)),&
						abs(real(eval_fit(is,igamma,ppar_comp))-f0_rel(sproc_rel,igamma,ipparbar))/f0_rel(sproc_rel,igamma,ipparbar)

						  integrate = integrate + &
                     gamma_rel(is_rel,igamma,ipparbar) * real(eval_fit(is,igamma,ppar_comp)) * &
                     2.d0 * pi * dgamma * dpparbar * (ms(is) / vA)**3
			endif
		 enddo
		enddo

	  else 		! non-relativistic:

	    dpperp = pp(is,2,2,1)-pp(is,1,2,1)
	    dppar = pp(is,2,2,2)-pp(is,2,1,2)

		do iperp=0,nperp
			do ipar=0,npar
				ppar_comp=pp(is,iperp,ipar,2)
				write (unit_spec,*) pp(is,iperp,ipar,1), pp(is,iperp,ipar,2), real(eval_fit(is,iperp,ppar_comp)), &
						abs(real(eval_fit(is,iperp,ppar_comp))-f0(is,iperp,ipar))/f0(is,iperp,ipar)
			     integrate = integrate + pp(is,iperp,ipar,1) * real(eval_fit(is,iperp,ppar_comp)) * &
                     2.d0 * pi * dpperp * dppar
			enddo
		enddo
	  endif
	  close(unit_spec)

       write(*,'(a,i3,a, 2es14.4)') ' Integration of fit/analytical function for species', is,':', integrate
	endif
	enddo

end subroutine



!-=-=-=-=-=-=
! This subroutine calculates the transposed Jacobian matrix of the fit function w.r.t. the
! fit parameter array.
!-=-=-=-=-=-=!
subroutine determine_JT(is,n_params,nJT,JT,params,iperp,upper_limit,ipparbar_lower)
	use alps_var, only : fit_type, n_fits, pp, ms, vA, perp_correction
	use alps_var, only : gamma_rel,pparbar_rel,nspec,relativistic
	implicit none
	integer :: n_params,ifit,par_ind,iperp,ipar,is,upper_limit,is_rel,sproc_rel,is_run,ipparbar_lower,nJT,JT_ind
	double precision :: ppar_val,pperp_val,params(n_params)
	double precision :: sqrtpart,expterm,kappapart, JT(nJT,0:upper_limit)


	if (relativistic(is)) then
		 ! Determine your sproc_rel:
		 is_rel=0
		 do is_run=1,nspec
				if (relativistic(is_run)) then
					is_rel=is_rel+1
					if (is_run.EQ.is) sproc_rel=is_rel
				endif
		 enddo
	endif

	do ipar=0,upper_limit

		if (relativistic(is)) then
			pperp_val=gamma_rel(sproc_rel,iperp,ipar+ipparbar_lower)
			ppar_val=pparbar_rel(sproc_rel,iperp,ipar+ipparbar_lower)
		else
			pperp_val=pp(is,iperp,ipar,1)
			ppar_val=pp(is,iperp,ipar,2)
		endif


		par_ind=0
		JT_ind=0

		do ifit=1,n_fits(is)

			if (fit_type(is,ifit).EQ.1) then	! Maxwell
				expterm=exp(-params(ifit+par_ind+1)*(ppar_val-params(ifit+par_ind+2))**2 &
							-perp_correction(is,ifit)*pperp_val**2)

				JT(ifit+JT_ind+0,ipar)=expterm

				JT(ifit+JT_ind+1,ipar)=-((ppar_val-params(ifit+par_ind+2))**2)*params(ifit+par_ind+0)*expterm

				JT(ifit+JT_ind+2,ipar)=2.d0*params(ifit+par_ind+1)*(ppar_val-params(ifit+par_ind+2))*&
							params(ifit+par_ind+0)*expterm

				par_ind=par_ind+2
				JT_ind=JT_ind+2
			endif

			if (fit_type(is,ifit).EQ.2) then	! kappa
				kappapart=1.d0+params(ifit+par_ind+1)*(ppar_val-params(ifit+par_ind+2))**2 &
						+ perp_correction(is,ifit)*params(ifit+par_ind+4)*pperp_val**2

				JT(ifit+JT_ind+0,ipar)=kappapart**params(ifit+par_ind+3)

				JT(ifit+JT_ind+1,ipar)=params(ifit+par_ind+0)*params(ifit+par_ind+3)*kappapart**(params(ifit+par_ind+3)-1.d0)*&
							(ppar_val-params(ifit+par_ind+2))**2

				JT(ifit+JT_ind+2,ipar)=params(ifit+par_ind+0)*params(ifit+par_ind+3)*kappapart**(params(ifit+par_ind+3)-1.d0)*&
							2.d0*params(ifit+par_ind+1)*(params(ifit+par_ind+2)-ppar_val)

				JT(ifit+JT_ind+3,ipar)=log(kappapart)*params(ifit+par_ind+0)*kappapart**params(ifit+par_ind+3)

				if (iperp.EQ.0) then
					JT_ind=JT_ind+3
				else
					JT(ifit+JT_ind+4,ipar)=params(ifit+par_ind+0)*params(ifit+par_ind+3)*&
									kappapart**(params(ifit+par_ind+3)-1.d0) *perp_correction(is,ifit)*pperp_val**2

					JT_ind=JT_ind+4
				endif

				par_ind=par_ind+4

			endif

			if (fit_type(is,ifit).EQ.3) then	! Juettner with pperp and ppar
				sqrtpart=sqrt(1.d0+(pperp_val**2+(ppar_val-params(ifit+par_ind+2))**2)*vA*vA/(ms(is)*ms(is)))

				JT(ifit+JT_ind+0,ipar)=exp(-params(ifit+par_ind+1)*sqrtpart)

				JT(ifit+JT_ind+1,ipar)=-params(ifit+par_ind+0)*exp(-params(ifit+par_ind+1)*sqrtpart)*sqrtpart

				JT(ifit+JT_ind+2,ipar)=params(ifit+par_ind+0)*exp(-params(ifit+par_ind+1)*sqrtpart)*&
							(params(ifit+par_ind+1)/sqrtpart)*(ppar_val-params(ifit+par_ind+2))*vA*vA/(ms(is)*ms(is))

				par_ind=par_ind+2
				JT_ind=JT_ind+2
			endif

			if (fit_type(is,ifit).EQ.4) then ! Juettner with gamma and pparbar, constant in pparbar
				JT(ifit+JT_ind+0,ipar)=exp(-perp_correction(is,ifit)*pperp_val)

				par_ind=par_ind+0
				JT_ind=JT_ind+0
			endif

			if (fit_type(is,ifit).EQ.5) then ! Juettner with gamma and pparbar with pparbar-dependence
				expterm=exp(-params(ifit+par_ind+1)*(ppar_val-params(ifit+par_ind+2))**2)* &
				exp(-perp_correction(is,ifit)*pperp_val)

				JT(ifit+JT_ind+0,ipar)=expterm

				JT(ifit+JT_ind+1,ipar)=-params(ifit+par_ind+0)*expterm*(ppar_val-params(ifit+par_ind+2))**2

				JT(ifit+JT_ind+2,ipar)=2.d0*params(ifit+par_ind+0)*(ppar_val-params(ifit+par_ind+2))*&
								params(ifit+par_ind+1)*expterm

				par_ind=par_ind+2
				JT_ind=JT_ind+2
			endif


			if (fit_type(is,ifit).EQ.6) then	! bi-Moyal
			expterm=exp(0.5d0*(params(ifit+par_ind+3)*perp_correction(is,ifit)*pperp_val**2 + &
					params(ifit+par_ind+1)*(ppar_val-params(ifit+par_ind+2))**2 - &
					exp(params(ifit+par_ind+3)*perp_correction(is,ifit)*pperp_val**2 + &
					params(ifit+par_ind+1)*(ppar_val-params(ifit+par_ind+2))**2) ))

				JT(ifit+JT_ind+0,ipar)=expterm

				JT(ifit+JT_ind+1,ipar)=params(ifit+par_ind+0)*expterm*0.5d0*(ppar_val-params(ifit+par_ind+2))**2 * &
							(1.d0 - exp(params(ifit+par_ind+3)*perp_correction(is,ifit)*pperp_val**2 + &
							params(ifit+par_ind+1) * (ppar_val-params(ifit+par_ind+2))**2) )

				JT(ifit+JT_ind+2,ipar)=params(ifit+par_ind+0)*expterm*params(ifit+par_ind+1)*&
							(params(ifit+par_ind+2)-ppar_val)* &
							(1.d0-exp(params(ifit+par_ind+3)*perp_correction(is,ifit)*pperp_val**2 + &
				 			params(ifit+par_ind+1) * (ppar_val-params(ifit+par_ind+2))**2) )

			 if (iperp.EQ.0) then
						JT_ind=JT_ind+2
					else
				 		JT(ifit+JT_ind+3,ipar)=params(ifit+par_ind+0)*expterm*0.5d0*perp_correction(is,ifit)*pperp_val**2 * &
 								(1.d0 - exp(params(ifit+par_ind+3)*perp_correction(is,ifit)*pperp_val**2 + &
									params(ifit+par_ind+1) * (ppar_val-params(ifit+par_ind+2))**2)  )

					JT_ind=JT_ind+3
				endif

				par_ind=par_ind+3
			endif

		enddo
	enddo

end subroutine





!-=-=-=-=-=-=
! This subroutine processes the nonlinear Levenberg-Marquart algorithm and returns
! the one-dimensional array params at a given iperp.
! quality is the sum of the squares of all residuals.
!-=-=-=-=-=-=
subroutine LM_nonlinear_fit(is,g,n_params,nJT,params,param_mask,iperp,npar,ipparbar_lower,quality)
	use alps_var, only : lambda_initial_fit, pp, lambdafac_fit, logfit
	use alps_var, only : epsilon_fit, maxsteps_fit
	use alps_var, only : gamma_rel,pparbar_rel,relativistic,nspec
	use alps_io, only : alps_error
	use mpi
	implicit none

	logical :: converged,param_mask(n_params)
	integer :: ipar,k,counter,is,n_params,nJT,iperp,npar,is_rel,is_run,sproc_rel,ipparbar_lower,l
	double precision :: LSQ,LSQnew,lambda_fit,quality,pperp_val,ppar_val
	double precision :: g(0:npar),params(n_params)
	double precision :: residuals(0:npar),deltaparam_fit(nJT)
	double precision :: JTJ(nJT,nJT),JT(nJT,0:npar)
	double precision :: diagmat(nJT,nJT), Amat(nJT,nJT)
	double precision :: work_array(nJT)
	integer :: ipiv(nJT), info


	converged=.FALSE.
	counter=0
	lambda_fit=lambda_initial_fit

	if ((1+npar).LT.n_params) call alps_error(2)

	if (relativistic(is)) then
		! Determine your sproc_rel:
		 is_rel=0
		 do is_run=1,nspec
				if (relativistic(is_run)) then
					is_rel=is_rel+1
					if (is_run.EQ.is) sproc_rel=is_rel
				endif
		 enddo
	endif

	do while(.NOT.converged)
	counter=counter+1

	LSQ=0.d0

	! Determine the transposed Jacobian and the residuals:

	call determine_JT(is,n_params,nJT,JT,params,iperp,npar,ipparbar_lower)



	do ipar=0,npar

		if (relativistic(is)) then
			pperp_val=gamma_rel(sproc_rel,iperp,ipar+ipparbar_lower)
			ppar_val=pparbar_rel(sproc_rel,iperp,ipar+ipparbar_lower)
		else
			pperp_val=pp(is,iperp,ipar,1)
			ppar_val=pp(is,iperp,ipar,2)
		endif



		if (logfit(is)) then

			do k=1,nJT
					JT(k,ipar)=JT(k,ipar)/fit_function(is,n_params,params,pperp_val,cmplx(ppar_val,0.d0,kind(1.d0)))
			enddo

			residuals(ipar)=g(ipar)-log(real(fit_function(is,n_params,params,pperp_val,cmplx(ppar_val,0.d0,kind(1.d0)))))
		else
			residuals(ipar)=g(ipar)-real(fit_function(is,n_params,params,pperp_val,cmplx(ppar_val,0.d0,kind(1.d0))))
		endif

		! Least squares:
		LSQ=LSQ+residuals(ipar)*residuals(ipar)
	enddo


	JTJ=matmul(JT,transpose(JT))


	diagmat=0.d0
	do k=1,nJT
		diagmat(k,k)=JTJ(k,k)
	enddo


	Amat=JTJ+lambda_fit*diagmat

	! The following routine is from LAPACK to invert the matrix:

	call dgetrf(nJT,nJT,Amat,nJT,ipiv,info)
	if (info.NE.0) stop "Fit matrix is numerically singular."

	call dgetri(nJT,Amat,nJT,ipiv,work_array,nJT,info)
	if (info.NE.0) stop "Fit matrix inversion failed."

	! Now Amat is the inverse of JTJ+lambda_fit*diagmat

	deltaparam_fit=matmul(Amat,matmul(JT,residuals))

	l=0
	do k=1,n_params
		if (param_mask(k)) then
			l=l+1
			params(k)=params(k)+deltaparam_fit(l)
		endif
	enddo

	! With the new param_fit, what is the new mean square error:
	LSQnew=0.d0
	do ipar=0,npar

		if (relativistic(is)) then
			pperp_val=gamma_rel(sproc_rel,iperp,ipar+ipparbar_lower)
			ppar_val=pparbar_rel(sproc_rel,iperp,ipar+ipparbar_lower)
		else
			pperp_val=pp(is,iperp,ipar,1)
			ppar_val=pp(is,iperp,ipar,2)
		endif


		if (logfit(is)) then
			residuals(ipar)=g(ipar)-log(real(fit_function(is,n_params,params,pperp_val,cmplx(ppar_val,0.d0,kind(1.d0)))))
		else
			residuals(ipar)=g(ipar)-real(fit_function(is,n_params,params,pperp_val,cmplx(ppar_val,0.d0,kind(1.d0))))
		endif

		! Least squares:
		LSQnew=LSQnew+residuals(ipar)*residuals(ipar)
	enddo


	if (LSQnew.GT.LSQ) then

		l=0
		do k=1,n_params
			if (param_mask(k)) then
				l=l+1
				params(k)=params(k)-deltaparam_fit(l)
			endif
		enddo

		lambda_fit=lambda_fit*lambdafac_fit
	else
		lambda_fit=lambda_fit/lambdafac_fit
	endif


	! Check if it converged (we can add further break criteria):
	if (abs(LSQnew-LSQ).LT.epsilon_fit) converged=.TRUE.
	if (counter.EQ.maxsteps_fit) converged=.TRUE.
	enddo

	quality=LSQ

end subroutine


end module alps_analyt
