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

module alps_fns

  implicit none

  private :: int_T, int_T_res, integrate
  private :: integrate_res, landau_integrate, resU

  public :: derivative_f0, disp, determine_bessel_array
  public :: determine_nmax, split_processes, secant, map_search

  public :: om_scan, om_double_scan

contains

!-=-=-=-=-
!Calculate parallel and perpendicular derivative
!of background velocity distribution functions
!-=-=-=-=-
subroutine derivative_f0
    use alps_var, only : f0, pp, df0, nperp, npar, nspec, arrayName
    use alps_var, only : f0_rel, gamma_rel, pparbar_rel,nspec_rel,df0_rel,ngamma,npparbar
    use alps_var, only : writeOut, pi, relativistic, usebM
    use alps_io,  only : get_unused_unit
    use alps_fns_rel, only : derivative_f0_rel
    implicit none
    !Local
    integer :: iperp, ipar, is, is_rel
    logical :: OutDF        !T-> output df0/dv to file
    logical :: any_relativistic
    character (50) :: fmt   !Output format
    character (100) :: writename
    double precision :: integrate, dpperp, dppar
    integer :: unit_f

    if (writeOut) then
       write(*,'(a)')&
            '-=-=-=-=-=-=-=-=-'
       write(*,'(a)')&
            'Calculating df0/dpperp, df0/ppar if necessary...'
    endif


    do is = 1, nspec

      if (usebM(is)) then
                write (*,'(a,i2)') ' Bi-Maxwellian calcuation: no derivatives necessary for species ',is
                df0(is,:,:,:)=0.d0
       else

       do iperp = 1, nperp-1
          do ipar = 1, npar-1
             !index 1-> vperp derivative
             df0(is,iperp,ipar,1) = &
                  (f0(is,iperp+1,ipar) - f0(is,iperp-1,ipar))/&
                  (pp(is,iperp+1,ipar,1) - pp(is,iperp-1,ipar,1))

             !index 2-> vpar derivative
             df0(is,iperp,ipar,2) = &
                  (f0(is,iperp,ipar+1) - f0(is,iperp,ipar-1))/&
                  (pp(is,iperp,ipar+1,2) - pp(is,iperp,ipar-1,2))
          enddo
       enddo
     endif

    enddo
    write(*,*)'Derivatives calculated'




    !Output df/dv to file
    OutDF = .false.
    if (OutDF) then

       if (writeOut) &
            write(*,'(a)') 'Outputing df0/dpperp, df0/dppar'
       write(fmt,'(a)') '(2es14.4,2es14.4)'
       do is = 1, nspec


         if (.not.usebM(is)) then
          dpperp = pp(is, 2, 2, 1) - pp(is, 1, 2, 1)
          dppar  = abs(pp(is, 2, 2, 2) - pp(is, 2, 1, 2))

          write(writeName,'(3a,i0,a)')&
               'distribution/',trim(arrayName),'_dfdv.',is,'.array'


          call get_unused_unit(unit_f)
          open(unit=unit_f,file=trim(writeName),status='replace')

          do iperp = 1, nperp-1
             do ipar = 1, npar-1

                   write(unit_f,fmt) pp(is,iperp,ipar,1),pp(is,iperp,ipar,2),&
                    df0(is, iperp,ipar,1),df0(is, iperp,ipar,2)
             enddo
             write(unit_f,*)
          enddo
          close(unit_f)
        endif
       enddo
    endif
    !End output



    do is = 1, nspec
       integrate = 0.d0
       dpperp = pp(is, 2, 2, 1) - pp(is, 1, 2, 1)
       dppar  = abs(pp(is, 2, 2, 2) - pp(is, 2, 1, 2))

       do iperp = 0, nperp
          do ipar = 0, npar
             integrate = integrate + &
                  pp(is,iperp,ipar,1) * f0(is,iperp,ipar) * &
                  2.d0 * pi * dpperp * dppar
          enddo
       enddo
       write(*,'(a,i3,a, 2es14.4)') 'Integration of species', is,':', integrate
     enddo




    if (writeOut) write(*,'(a)') '-=-=-=-=-=-=-=-=-'

	if (writeOut) write (*,'(a)') 'Determine if relativistic calculation necessary'
	any_relativistic = .FALSE.

	is_rel = 0
	do is = 1, nspec
		if (relativistic(is)) then
			any_relativistic=.TRUE.
			is_rel=is_rel+1
			write (*,'(a,i3,a,i3)') 'Species',is, ' requires relativistic calculation. is_rel=',is_rel
		endif
	enddo


	if (any_relativistic) then

		if(writeOut) write (*,'(a,i4)') "ngamma   = ",ngamma
		if(writeOut) write (*,'(a,i4)') "npparbar = ",npparbar

		 nspec_rel=is_rel

		 ! Allocate the relativistic fields (only on proc0 for now):
		 allocate(gamma_rel(nspec_rel,0:ngamma,0:npparbar))
		 allocate(pparbar_rel(nspec_rel,0:ngamma,0:npparbar))
		 allocate(f0_rel(nspec_rel,0:ngamma,0:npparbar))
		 allocate(df0_rel(nspec_rel,0:ngamma,0:npparbar,2))

		 is_rel=0
		 do is=1,nspec
		 	if(relativistic(is)) then
		 		is_rel=is_rel+1
		 		call derivative_f0_rel(is,is_rel)
		 	endif
		 enddo
	else
		 nspec_rel=0
		 write (*,'(a)') 'No relativistic calculation necessary.'
		 if (writeOut) write(*,'(a)') '-=-=-=-=-=-=-=-=-'
	endif

end subroutine derivative_f0
!-=-=-=-=-=-=
!
!-=-=-=-=-=-=
 double complex function disp(om)
    use alps_var, only : nlim, proc0, nspec, ierror, sproc, relativistic, iproc
    use alps_var, only : wave, kperp, kpar, ns, qs, vA, chi0, usebM
    use alps_nhds, only: calc_chi
    use alps_fns_rel, only : int_ee_rel
    use mpi
    implicit none

    !Passed
    double complex :: om     !frequency for dispersion solution
    !Local
    double complex :: chi_NHDS(3,3)
    double complex, dimension(1:nspec,1:3,1:3) :: schi,chi
    double complex, dimension(1:3,1:3) :: eps
    double complex :: enx2, enz2, enxnz            !Indices of refraction
    integer :: is, nn        !Indices for species, bessel n
    double complex, dimension(1:nspec) :: norm !normalization for tensors
    logical :: found_res_plus,found_res_minus

    chi=cmplx(0.d0,0.d0,kind(1.d0))
    if (proc0) chi0=cmplx(0.d0,0.d0,kind(1.d0))
    schi=cmplx(0.d0,0.d0,kind(1.d0))

    if (proc0)  then
       !-=-=-=-=-=-
       !proc0 doesn't concern itself with the task
       !of integration, or the calculation of chi_ns
       !-=-=-=-=-=-

       !Indices of refraction for the dispersion relation

       ! old NHDS normalization
       !enx2=((kperp/(om))**2)
       !enz2=((kpar/(om))**2)
       !enxnz=(kperp*kpar/(om**2))

       !new NHDS normalization
       enx2=kperp**2
       enz2=kpar**2
       enxnz=kperp*kpar

       !PLUME normalization
       !enx2=((kperp/(om*vA))**2)
       !enz2=((kpar/(om*vA))**2)
       !enxnz=(kperp*kpar/(om*vA)**2)

    else
       !integrate
       !c.f. Stix Equation 10.48; pg 255

       !Calculate only for assigned species sproc

       !-=-=-=-=-=-
       !The processor runs over the n indicies defined by split_processes
       !-=-=-=-=-=


       ! Split into NHDS or ALPS routines:

       ! Only run the NHDS routine if useBM is on for the species
       !   and if you are handling n=0 according to split_processes
       if (usebM(sproc).and.(nlim(2).GE.0).and.(nlim(1).EQ.0)) then

          ! This is the case to use NHDS for the calculation of chi:

          call calc_chi(chi_NHDS,sproc,kpar,kperp,om)

          ! Account for norm below, which is already included in NHDS:
          schi(sproc,1,1)=chi_NHDS(1,1)/(ns(sproc) * qs(sproc))
          schi(sproc,2,2)=chi_NHDS(2,2)/(ns(sproc) * qs(sproc))
          schi(sproc,3,3)=chi_NHDS(3,3)/(ns(sproc) * qs(sproc))
          schi(sproc,1,2)=chi_NHDS(1,2)/(ns(sproc) * qs(sproc))
          schi(sproc,1,3)=chi_NHDS(1,3)/(ns(sproc) * qs(sproc))
          schi(sproc,2,3)=chi_NHDS(2,3)/(ns(sproc) * qs(sproc))

       else

       do nn = nlim(1),nlim(2)

          call determine_resonances(om,nn,found_res_plus,found_res_minus)

          !CHIij(nn) function calls

          if (nn == 0) then

             !xx term is zero
             !xy term is zero
             !xz term is zero

             !yy term
             schi(sproc,2,2) = schi(sproc,2,2) + &
                  full_integrate(om,nn,2,found_res_plus)

             !zz term
             schi(sproc,3,3) = schi(sproc,3,3) + &
                  full_integrate(om,nn,3,found_res_plus)

             !yz term
             schi(sproc,2,3) = schi(sproc,2,3) + &
                  full_integrate(om,nn,6,found_res_plus)

             !if (sproc==2) then
             !   write(*,'(i4,6es14.4)') nn,schi(sproc,2,2),schi(sproc,2,3),schi(sproc,3,3)
             !endif

          else
             !n != 0 loop

             !xx term
             schi(sproc,1,1) = schi(sproc,1,1) + &
                  full_integrate(om,nn,1,found_res_plus) + &
                  full_integrate(om,-nn,1,found_res_minus)

             !yy term
             schi(sproc,2,2) = schi(sproc,2,2) + &
                  full_integrate(om,nn,2,found_res_plus) + &
                  full_integrate(om,-nn,2,found_res_minus)

             !zz term
             schi(sproc,3,3) = schi(sproc,3,3) + &
                  full_integrate(om,nn,3,found_res_plus) + &
                  full_integrate(om,-nn,3,found_res_minus)

             !xy term
             schi(sproc,1,2) = schi(sproc,1,2) + &
                  full_integrate(om,nn,4,found_res_plus) + &
                  full_integrate(om,-nn,4,found_res_minus)

             !xz term
             schi(sproc,1,3) = schi(sproc,1,3) + &
                  full_integrate(om,nn,5,found_res_plus) + &
                  full_integrate(om,-nn,5,found_res_minus)

             !yz term
             schi(sproc,2,3) = schi(sproc,2,3) + &
                  full_integrate(om,nn,6,found_res_plus) + &
                  full_integrate(om,-nn,6,found_res_minus)

             !End n != 0 loop
          endif

          !if (sproc==2) then
          !   write(*,'(i4,12es14.4)') nn,schi(sproc,1,1),schi(sproc,2,2),schi(sproc,3,3),&
          !        schi(sproc,1,2),schi(sproc,1,3),schi(sproc,2,3)
          !endif

       enddo !End n index loop


       !-=-=-=-=-=-
       !add in non-T zz term
       !(if processor is responsible for the n=0 term)
       !EDIT: Add in ee term
       !-=-=-=-=-=-
       if (nlim(1)==0) then
	     if(relativistic(sproc)) then
	     	    schi(sproc,3,3)=schi(sproc,3,3) + int_ee_rel(om)
	     else
	     		schi(sproc,3,3)=schi(sproc,3,3) + int_ee(om)
	     endif
       endif
       !-=-=-=-=-=-

     endif

       ! NOTE ON NORMALISATION:
       ! NOW WE ARE WORKING WITH A NEW NORMALISATION:
       ! WE HAVE ALREADY MULTIPLIED THE schi TERMS AND THE int_ee TERMS INTERNALLY
       ! BY om. The following norm factors account for this:

       ! old NHDS normalization
       !norm(sproc) = ns(sproc) * qs(sproc) / (om*om)

       !new NHDS normalization
       norm(sproc) = ns(sproc) * qs(sproc)

       !PLUME normalization
       !norm(sproc) = ns(sproc) * qs(sproc)/(om*om*vA**2)

       !-=-=-=-=-=-
       !Multiply chi_s by the desired normalization
       !-=-=-=-=-=-
       schi(sproc,1,1) = schi(sproc,1,1) * norm(sproc)

       schi(sproc,2,2) = schi(sproc,2,2) * norm(sproc)

       schi(sproc,3,3) = schi(sproc,3,3) * norm(sproc)

       schi(sproc,1,2) = schi(sproc,1,2) * norm(sproc)

       schi(sproc,1,3) = schi(sproc,1,3) * norm(sproc)

       schi(sproc,2,3) = schi(sproc,2,3) * norm(sproc)
       !-=-=-=-=-=-

    endif

    ! Return the schi to proc0

    call MPI_REDUCE (schi, chi, size(chi),&
        MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierror)


    if (proc0) then
       !Calculate eps
       !eps  = (1 0 0)         ( chi_xx  chi_xy  chi_xz )
       !       (0 1 0) + SUM_s ( chi_yx  chi_yy  chi_yz )
       !       (0 0 1)         ( chi_zx  chi_zy  chi_zz )|s

       !The global variable 'chi0' is used
       !for heating & eigenfunction calculation


       !old NHDS normalization
       !chi0=chi/(vA*vA)

       !new NHDS normalization
       chi0=chi/(om*om*vA*vA)

       !PLUME normalization
       !chi0=chi

       chi0(:,2,1)=-chi0(:,1,2)
       chi0(:,3,1)=-chi0(:,1,3)
       chi0(:,3,2)=-chi0(:,2,3)

    	wave=cmplx(0.d0,0.d0,kind(1.d0))
      eps=cmplx(0.d0,0.d0,kind(1.d0))


       !Sum over species!
       do is = 1, nspec
          eps(1,1) = eps(1,1) + chi(is,1,1)
          eps(2,2) = eps(2,2) + chi(is,2,2)
          eps(3,3) = eps(3,3) + chi(is,3,3)

          eps(1,2) = eps(1,2) + chi(is,1,2)
          eps(1,3) = eps(1,3) + chi(is,1,3)
          eps(2,3) = eps(2,3) + chi(is,2,3)

       enddo

       !-=-=-=-=-=-
       !The susceptibility tensor has
       !    a nice symmetry
       !-=-=-=-=-=-
       eps(2,1) = -eps(1,2)
       eps(3,1) =  eps(1,3)
       eps(3,2) = -eps(2,3)
       !-=-=-=-=-=-

       !-=-=-=-=-=-=-=-=-=
       !Add the unit tensor-
       !-=-=-=-=-=-=-=-=-=

       !NHDS normalization
       !eps(1,1) = eps(1,1) + vA**2
       !eps(2,2) = eps(2,2) + vA**2
       !eps(3,3) = eps(3,3) + vA**2

       !NHDS normalization multiplied with om**2
       eps(1,1) = eps(1,1) + (om*vA)**2
       eps(2,2) = eps(2,2) + (om*vA)**2
       eps(3,3) = eps(3,3) + (om*vA)**2

       !PLUME normalization
       !eps(1,1) = eps(1,1) + 1.d0
       !eps(2,2) = eps(2,2) + 1.d0
       !eps(3,3) = eps(3,3) + 1.d0

       !-=-=-=-=-=-=-=-=-=

       !Calculate wave
       !wave = ( eps_xx - nz^2  eps_xy              eps_xz + nxnz )
       !       ( eps_yx         eps_yy -nz^2 -nx^2  eps_yz        )
       !       ( eps_zx + nxnz  eps_zy              eps_zz - nx^2 )


       wave(1,1) = eps(1,1) - enz2
       wave(2,2) = eps(2,2) - enz2 - enx2
       wave(3,3) = eps(3,3) - enx2
       wave(1,3) = eps(1,3) + enxnz

       wave(1,2) = eps(1,2)
       wave(2,3) = eps(2,3)

       !-=-=-=-=-=-
       !The wave tensor has
       !    a nice symmetry
       !-=-=-=-=-=-
       wave(2,1) = -wave(1,2)
       wave(3,1) =  wave(1,3)
       wave(3,2) = -wave(2,3)
       !-=-=-=-=-=-

       !Calculate D(k,omega)
       !The below relies on the symmetries of the T_n tensor:
       !Again, c.f. Stix Equation 10.48; pg 255
       !---------------------------------------------------------------------
       disp = wave(1,1)*( wave(2,2)*wave(3,3) + wave(2,3)**2 ) + &
            2.d0*wave(1,2)*wave(2,3)*wave(1,3) - wave(1,3)**2*wave(2,2) + &
            wave(1,2)**2*wave(3,3)
       !---------------------------------------------------------------------
       !write(*,'(a,2es14.4)')'disp: ',disp

    endif

    call mpi_bcast(disp, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)

    !Make sure all processors have completed calculation to
    !cross contamination with map and root searches
    call mpi_barrier(mpi_comm_world,ierror)

    return

end function disp




!-=-=-=-=-=-=
! Determine if there are resonances
!-=-=-=-=-=-=
subroutine determine_resonances(om,nn,found_res_plus,found_res_minus)
	use alps_var, only : npar, pp, vA, ms, qs, sproc, kpar, nperp
	use alps_var, only : positions_principal, relativistic
	use alps_io, only  : alps_error
	implicit none
	integer :: nn, ipar, iperp
	logical :: found_res_plus,found_res_minus
	double precision :: dppar, gamma
	double complex :: p_res, om

	 dppar = pp(sproc,2,2,2)-pp(sproc,2,1,2)
	 found_res_plus = .FALSE.
	 found_res_minus = .FALSE.

	if (relativistic(sproc)) then

  !call determine_sproc_rel(sproc_rel)

  !do igamma=0,ngamma
	  ! Note that this is pparbar in reality, but it does not make a difference for this section:
    ! positive n:
    !p_res = gamma_rel(sproc_rel,igamma,1)*om*vA/kpar - (1.d0*nn)*qs(sproc)*vA/(kpar*ms(sproc))
    !if ((real(p_res)**2).LE.(gamma_rel(sproc_rel,igamma,1)**2-1.d0)) found_res_plus = .TRUE.

    ! negative n:
	  !p_res = gamma_rel(sproc_rel,igamma,1)*om*vA/kpar + (1.d0*nn)*qs(sproc)*vA/(kpar*ms(sproc))
	  !if ((real(p_res)**2).LE.(gamma_rel(sproc_rel,igamma,1)**2-1.d0)) found_res_minus = .TRUE.
!	enddo

! The following checks for resonances in (pperp,ppar)-space rather than in (Gamma,pparbar)-space:
  do iperp=0,nperp
    do ipar=0,npar-1

     gamma = sqrt((pp(sproc, iperp, ipar, 1)**2 + &
        pp(sproc, iperp, ipar, 2)**2) * vA**2/ms(sproc)**2 +1.d0)

     ! positive n:
     p_res=(gamma*ms(sproc) * om - 1.d0 * nn * qs(sproc))/kpar

     if ((pp(sproc,2,ipar,2).LE.real(p_res)).and.&
			  (pp(sproc,2,ipar+1,2).GT.real(p_res))) found_res_plus = .TRUE.


    ! negative n:
    p_res=(gamma*ms(sproc) * om + 1.d0 * nn * qs(sproc))/kpar

        if ((pp(sproc,2,ipar,2).LE.real(p_res)).and.&
   			  (pp(sproc,2,ipar+1,2).GT.real(p_res))) found_res_minus = .TRUE.
	  enddo
  enddo


	else ! non-relativistic case:

	! positive n
	  ipar = 0
	  p_res = (ms(sproc) * om - 1.d0 * nn * qs(sproc))/kpar

	  do while ((ipar.LE.(npar-2)).AND.(.NOT.found_res_plus))
		 ipar = ipar + 1
		 if ((pp(sproc,2,ipar,2).LE.real(p_res)).and.&
			  (pp(sproc,2,ipar+1,2).GT.real(p_res))) found_res_plus = .TRUE.

	  enddo

	! negative n:
	  ipar = 0
	  p_res = (ms(sproc) * om + 1.d0 * nn * qs(sproc))/kpar

	  do while ((ipar.LE.(npar-2)).AND.(.NOT.found_res_minus))
		 ipar = ipar + 1
		 if ((pp(sproc,2,ipar,2).LE.real(p_res)).and.&
			  (pp(sproc,2,ipar+1,2).GT.real(p_res))) found_res_minus = .TRUE.

	  enddo


	  ! Check if there is a resonance right outside the integration domain:
	  p_res = (ms(sproc) * om - 1.d0 * nn * qs(sproc))/kpar

	  if ((real(p_res).LT.pp(sproc,2,1,2)).AND.&
		(real(p_res).GE.(pp(sproc,2,1,2)-(1.d0*positions_principal)*dppar))) found_res_plus=.TRUE.
	  if ((real(p_res).GE.pp(sproc,2,npar-1,2)).AND.&
		(real(p_res).LT.(pp(sproc,2,npar-1,2)+(1.d0*positions_principal)*dppar))) found_res_plus=.TRUE.

	   p_res = (ms(sproc) * om + 1.d0 * nn * qs(sproc))/kpar

	  if ((real(p_res).LT.pp(sproc,2,1,2)).AND.&
		(real(p_res).GE.(pp(sproc,2,1,2)-(1.d0*positions_principal)*dppar))) found_res_minus=.TRUE.
	  if ((real(p_res).GE.pp(sproc,2,npar-1,2)).AND.&
		(real(p_res).LT.(pp(sproc,2,npar-1,2)+(1.d0*positions_principal)*dppar))) found_res_minus=.TRUE.
	endif

end subroutine determine_resonances




!-=-=-=-=-=-=
!THE INTEGRATOR
!-=-=-=-=-=-=
double complex function full_integrate(om, nn, mode, found_res)
  use alps_var, only : npar, relativistic, sproc
  use alps_fns_rel, only : integrate_res_rel,landau_integrate_rel

  implicit none
  !Passed
  double complex :: om   !complex frequency
  integer :: nn          !Bessel N
  integer :: mode        !index in T tensor
  logical :: found_res   !use brute force integral or landau & principal value

  full_integrate = cmplx(0.d0,0.d0,kind(1.d0))

  if (.not. found_res) then
     !Brute force integrate
     full_integrate = integrate(om, nn, mode, 1, npar-1)
	 elseif (found_res.and.relativistic(sproc)) then
       	 	if (aimag(om).GT.0.d0) then
	    		full_integrate = integrate_res_rel(om,nn,mode)
       	 	elseif (aimag(om).LT.0.d0) then
      	 		full_integrate = integrate_res_rel(om,nn,mode) + 2.d0 * landau_integrate_rel(om, nn, mode)
       	 	elseif (aimag(om).EQ.0.d0) then
	       	 	full_integrate = integrate_res_rel(om,nn,mode) + landau_integrate_rel(om, nn, mode)
      	 	endif
  elseif ((found_res).and.(aimag(om).GT.0.d0)) then
     !Brute force integrate
	 full_integrate = integrate_res(om,nn,mode)
  elseif ((found_res).and.(aimag(om).LT.0.d0)) then
     !Landau Integral
   	 full_integrate = integrate_res(om,nn,mode)+2.d0 * landau_integrate(om, nn, mode)
  elseif ((found_res).and.(aimag(om).EQ.0.d0)) then
   	 full_integrate = integrate_res(om,nn,mode)+landau_integrate(om, nn, mode)
  endif
  return
end function full_integrate






double complex function integrate(om, nn, mode, iparmin, iparmax)
  use alps_var, only : nperp, pp, pi, sproc
  implicit none
  !Passed
  double complex :: om   !complex frequency
  integer :: nn          !Bessel N
  integer :: mode        !index in T tensor
  integer :: iparmin,iparmax
  !Local
  integer :: iperp, ipar !p_perp, p_par index
  double precision :: dpperp, dppar !delta p_perp, delta p_par

  !choose numerical integration method...

  integrate = cmplx (0.d0, 0.d0,kind(1.d0))

  dpperp = pp(sproc, 2, 2, 1) - pp(sproc, 1, 2, 1)
  dppar  = abs(pp(sproc, 2, 2, 2) - pp(sproc, 2, 1, 2))

  integrate = integrate + &
       2.d0 * resU(om, nn, 1, iparmin) * int_T(nn, 1, iparmin, mode) + &
       2.d0 * resU(om, nn, 1, iparmax) * int_T(nn, 1, iparmax, mode) + &
       resU(om, nn, nperp -1, iparmin) * int_T(nn, nperp -1, iparmin, mode) + &
       resU(om, nn, nperp -1, iparmax) * int_T(nn, nperp -1, iparmax, mode)

  do iperp = 2, nperp-2
     do ipar = iparmin+1, iparmax-1
        integrate = integrate + 4.d0 * resU(om, nn, iperp, ipar) * int_T(nn, iperp, ipar, mode)
     enddo

     integrate = integrate + &
      2.d0 * ( resU(om, nn, iperp, iparmin) * int_T(nn, iperp, iparmin, mode) + &
      resU(om, nn, iperp, iparmax) * int_T(nn, iperp, iparmax, mode) )
  enddo


  do ipar = iparmin+1, iparmax-1
        integrate = integrate + &
             2.d0 * (2.d0 * resU(om, nn, 1, ipar) * int_T(nn, 1, ipar, mode) + &
             resU(om, nn, nperp -1, ipar) * int_T(nn, nperp -1, ipar, mode) )
  enddo

  integrate = 2.d0 * pi * integrate * dpperp * dppar * 0.25d0

  return

end function integrate






! This function does the integration around resonances if necessary:
double complex function integrate_res(om,nn,mode)
	use alps_var, only : nperp,npar,pp,ms,qs,kpar,pi,sproc
	use alps_var, only : positions_principal,n_resonance_interval, Tlim
	implicit none
	integer :: ipar_res, ipar,iperp, ntiny
	integer :: lowerlimit,upperlimit,nn,mode
	double precision :: dpperp,dppar,capDelta,smdelta,denomR,denomI,ppar,correction
	double complex :: p_res,ii,om,integrate_norm, gprimetr
	logical :: found_res

	dpperp = pp(sproc, 2, 2, 1) - pp(sproc, 1, 2, 1)
	dppar  = pp(sproc, 2, 2, 2) - pp(sproc, 2, 1, 2)

	integrate_res=cmplx(0.d0,0.d0,kind(1.d0))
	integrate_norm=cmplx(0.d0,0.d0,kind(1.d0))

	ii=cmplx(0.d0,1.d0,kind(1.d0))

	! determine the position of the resonance (i.e., the step LEFT of it):
	ipar = 0
	ipar_res=0
	found_res = .FALSE.

	p_res = (ms(sproc) * om - 1.d0*nn * qs(sproc))/kpar

	do while ((ipar.LT.(npar-2)).AND.(.NOT.found_res))
		ipar = ipar + 1

		if ((pp(sproc,2,ipar+1,2).GT.real(p_res)).and.&
			  (pp(sproc,2,ipar,2).LE.real(p_res))) then
				ipar_res=ipar
				found_res = .TRUE.
		endif
	enddo


	! Handle resonances that are right outside the integration domain:
	p_res = (ms(sproc) * om - (1.d0*nn) * qs(sproc))/kpar

	do ipar=0,positions_principal

		if ((real(p_res).GE.(pp(sproc,2,0,2)-dppar*ipar)).AND.&
		  (real(p_res).LT.(pp(sproc,2,0,2)-dppar*(ipar-1)))) ipar_res = -ipar

		if ((real(p_res).GE.(pp(sproc,2,npar-1,2)+dppar*ipar)).AND.&
		  (real(p_res).LT.(pp(sproc,2,npar-1,2)+dppar*(ipar+1)))) ipar_res = npar-1+ipar

	enddo


	! If the resonance is close to the edge, just do the normal integration:
	if ((ipar_res-positions_principal).LE.2) then
		integrate_res=integrate(om, nn, mode, ipar_res+positions_principal,npar-1)
		return
	endif


	if ((ipar_res+positions_principal).GE.(npar-2)) then
		integrate_res=integrate(om, nn, mode, 1, ipar_res-positions_principal)
		return
	endif

	! positions_principal defines how close we can go to ipar_res with the "normal" integration.
	! the following part is basically the normal function "integrate" on the left and on the right of ipar_res:
	! left:
	lowerlimit=ipar_res-positions_principal
	integrate_norm = integrate(om, nn, mode, 1, lowerlimit)


	! right:
	if(abs(real(p_res)-pp(sproc,2,ipar_res,2)).LT.(0.5d0*dppar)) then
		upperlimit=ipar_res+positions_principal+1
	else
		upperlimit=ipar_res+positions_principal+2
	endif
	integrate_norm = integrate_norm + integrate(om, nn, mode, upperlimit, npar-1)



	! The following part includes the analytic switch described in the paper!
	! Now comes the resonance part:
	! We call the function that needs to be integrated WITHOUT the resonance part funct_g.
	! We linearize this function. Now we can calculate the even part of the integration.
	! We set Delta so that it starts at ipar_res-positions_principal. In that way, there is only
	! a tiny rest left on the right side that needs to be integrated.
	! split the range between the resonance and the upper limit into n_resonance_interval steps:
	! the denominator is:
	! (ppar - gamma * ms(sproc)*om/kpar + (1.d0*nn) * qs(sproc) /kpar )
	! we define
	! denomR=real(gamma * ms(sproc)*om/kpar - (1.d0*nn) * qs(sproc) /kpar )
	! denomI=aimag(gamma * ms(sproc)*om/kpar)
	! so that the denominator is
	! (ppar-denomR-ii*denomI)

	denomR=real(ms(sproc)*om/kpar - (1.d0*nn) * qs(sproc) /kpar )
	denomI=aimag(ms(sproc)*om/kpar)

	capDelta=real(p_res)-pp(sproc,1,ipar_res-positions_principal,2)
	smdelta=capDelta/(1.d0*n_resonance_interval)

	if (abs(denomI).GT.Tlim) then ! regular integration:
		! Integrate the boundaries:
		ppar=real(p_res) ! left end of integration interval

    ! At iperp=1, we are already missing the part from iperp=0, where we should actually start. Therefore, we use 4 instead of 2 in the trapezoid integration:
		integrate_res = integrate_res + 2.d0 * funct_g(ppar,1,om,nn,mode)/(ppar-denomR-ii*denomI)
		integrate_res = integrate_res - 2.d0 * funct_g(2.d0*denomR-ppar,1,om,nn,mode)/(ppar-denomR+ii*denomI)

		integrate_res = integrate_res + funct_g(ppar,nperp-1,om,nn,mode)/(ppar-denomR-ii*denomI)
		integrate_res = integrate_res - funct_g(2.d0*denomR-ppar,nperp-1,om,nn,mode)/(ppar-denomR+ii*denomI)

		ppar=capDelta+real(p_res)	! right end of integration interval

		integrate_res = integrate_res + 2.d0 * funct_g(ppar,1,om,nn,mode)/(ppar-denomR-ii*denomI)
		integrate_res = integrate_res - 2.d0 * funct_g(2.d0*denomR-ppar,1,om,nn,mode)/(ppar-denomR+ii*denomI)

		integrate_res = integrate_res + funct_g(ppar,nperp-1,om,nn,mode)/(ppar-denomR-ii*denomI)
		integrate_res = integrate_res - funct_g(2.d0*denomR-ppar,nperp-1,om,nn,mode)/(ppar-denomR+ii*denomI)


		! end of edges.

		do iperp = 2, nperp-2
			do ipar = 1, n_resonance_interval-1
				ppar=real(p_res)+smdelta*ipar

				integrate_res = integrate_res + 4.d0 * funct_g(ppar,iperp,om,nn,mode)/(ppar-denomR-ii*denomI)
				integrate_res = integrate_res - 4.d0 * funct_g(2.d0*denomR-ppar,iperp,om,nn,mode)/(ppar-denomR+ii*denomI)

			enddo

			ppar=real(p_res)+smdelta
			ppar=real(p_res)

			integrate_res = integrate_res + 2.d0 * funct_g(ppar,iperp,om,nn,mode)/(ppar-denomR-ii*denomI)
			integrate_res = integrate_res - 2.d0 * funct_g(2.d0*denomR-ppar,iperp,om,nn,mode)/(ppar-denomR+ii*denomI)

			ppar=real(p_res)+capDelta

			integrate_res = integrate_res + 2.d0 * funct_g(ppar,iperp,om,nn,mode)/(ppar-denomR-ii*denomI)
			integrate_res = integrate_res - 2.d0 * funct_g(2.d0*denomR-ppar,iperp,om,nn,mode)/(ppar-denomR+ii*denomI)

		enddo


		do ipar = 1, n_resonance_interval-1
			ppar=real(p_res)+smdelta*ipar

      ! At iperp=1, we are already missing the part from iperp=0, where we should actually start. Therefore, we use 4 instead of 2 in the trapezoid integration:
			integrate_res = integrate_res + 4.d0 * funct_g(ppar,1,om,nn,mode)/(ppar-denomR-ii*denomI)
			integrate_res = integrate_res - 4.d0 * funct_g(2.d0*denomR-ppar,1,om,nn,mode)/(ppar-denomR+ii*denomI)

			integrate_res = integrate_res + 2.d0 * funct_g(ppar,nperp-1,om,nn,mode)/(ppar-denomR-ii*denomI)
			integrate_res = integrate_res - 2.d0 * funct_g(2.d0*denomR-ppar,nperp-1,om,nn,mode)/(ppar-denomR+ii*denomI)

		  enddo


	else ! analytic approximation

		! Integrate the edges:
		ppar=real(p_res) ! left end of integration interval
    ! At iperp=1, we are already missing the part from iperp=0, where we should actually start. Therefore, we use 4 instead of 2 in the trapezoid integration:
		gprimetr = (funct_g(denomR+dppar,1,om,nn,mode)-funct_g(denomR-dppar,1,om,nn,mode))/(2.d0*dppar)
		if (denomI.NE.0.d0) then
			integrate_res = integrate_res + 2.d0 * 2.d0 * gprimetr * (ppar-denomR)**2 / ((ppar-denomR)**2+denomI**2)
		else
			integrate_res = integrate_res + 2.d0 * 2.d0 * gprimetr
		endif

		gprimetr = (funct_g(denomR+dppar,nperp-1,om,nn,mode)-funct_g(denomR-dppar,nperp-1,om,nn,mode))/(2.d0*dppar)
		if (denomI.NE.0.d0) then
			integrate_res = integrate_res + 2.d0 * gprimetr * (ppar-denomR)**2 / ((ppar-denomR)**2+denomI**2)
		else
			integrate_res = integrate_res + 2.d0 * gprimetr
		endif


		ppar=capDelta+real(p_res)	! right end of integration interval
		gprimetr = (funct_g(denomR+dppar,1,om,nn,mode)-funct_g(denomR-dppar,1,om,nn,mode))/(2.d0*dppar)
		integrate_res = integrate_res + 2.d0 * 2.d0*gprimetr*(ppar-denomR)**2 / ((ppar-denomR)**2+denomI**2)

		gprimetr = (funct_g(denomR+dppar,nperp-1,om,nn,mode)-funct_g(denomR-dppar,nperp-1,om,nn,mode))/(2.d0*dppar)
		integrate_res = integrate_res + 2.d0*gprimetr*(ppar-denomR)**2 / ((ppar-denomR)**2+denomI**2)

    ! The following lines account for Eq. (3.7) in the paper:
    if (denomI.GT.0.d0) then
      integrate_res = integrate_res + 2.d0 * 2.d0 * ii * pi * funct_g(denomR,1,om,nn,mode)/smdelta
			integrate_res = integrate_res + 2.d0 * ii * pi * funct_g(denomR,nperp-1,om,nn,mode)/smdelta
    else if (denomI.LT.0.d0) then
      integrate_res = integrate_res - 2.d0 * 2.d0 * ii * pi * funct_g(denomR,1,om,nn,mode)/smdelta
      integrate_res = integrate_res - 2.d0 * ii * pi * funct_g(denomR,nperp-1,om,nn,mode)/smdelta
    else if (denomI.EQ.0.d0) then
      integrate_res = integrate_res+0.d0
    endif
		! end of edges.

		do iperp = 2, nperp-2
			do ipar = 1, n_resonance_interval-1
				ppar=real(p_res)+smdelta*ipar

				gprimetr = (funct_g(denomR+dppar,iperp,om,nn,mode)-&
					funct_g(denomR-dppar,iperp,om,nn,mode))/(2.d0*dppar)
				integrate_res = integrate_res + 4.d0 * 2.d0 * gprimetr * (ppar-denomR)**2 / ((ppar-denomR)**2+denomI**2)
			enddo

			ppar=real(p_res)
      ! In this case, ppar is equal to denomR, so: no integration needed

	!		gprimetr = (funct_g(denomR+dppar,iperp,om,nn,mode)-funct_g(denomR-dppar,iperp,om,nn,mode))/(2.d0*dppar)
	!		if (denomI.NE.0.d0) then
	!			integrate_res = integrate_res + 2.d0 * 2.d0 * gprimetr * (ppar-denomR)**2 / ((ppar-denomR)**2+denomI**2)
	!		else
	!			integrate_res = integrate_res + 2.d0 * 2.d0 * gprimetr
	!		endif

			ppar=real(p_res)+capDelta
			gprimetr = (funct_g(denomR+dppar,iperp,om,nn,mode)-funct_g(denomR-dppar,iperp,om,nn,mode))/(2.d0*dppar)
			integrate_res = integrate_res + 2.d0*2.d0*gprimetr*(ppar-denomR)**2 / ((ppar-denomR)**2+denomI**2)


      ! The following lines account for Eq. (3.7) in the paper:
      if (denomI.GT.0.d0) then
        integrate_res = integrate_res + 4.d0 * ii * pi * funct_g(denomR,iperp,om,nn,mode)/smdelta
      else if (denomI.LT.0.d0) then
        integrate_res = integrate_res - 4.d0 * ii * pi * funct_g(denomR,iperp,om,nn,mode)/smdelta
      else if (denomI.EQ.0.d0) then
        integrate_res = integrate_res+0.d0
      endif

		enddo



		do ipar = 1, n_resonance_interval-1
			ppar=real(p_res)+smdelta*ipar

			gprimetr = (funct_g(denomR+dppar,1,om,nn,mode)-funct_g(denomR-dppar,1,om,nn,mode))/(2.d0*dppar)
			integrate_res = integrate_res + 4.d0*2.d0*gprimetr*(ppar-denomR)**2 / ((ppar-denomR)**2+denomI**2)

			gprimetr = (funct_g(denomR+dppar,nperp-1,om,nn,mode)-funct_g(denomR-dppar,nperp-1,om,nn,mode))/(2.d0*dppar)
			integrate_res = integrate_res + 2.d0*2.d0*gprimetr*(ppar-denomR)**2 / ((ppar-denomR)**2+denomI**2)

		  enddo


	endif





	! There is a tiny rest left between the point real(p_res)+capDelta and the position
	! pp(sproc,2,upperlimit,2). We split this interval into steps of roughly size smdelta:
	ntiny=int((pp(sproc,2,upperlimit,2)-real(p_res)-capDelta)/smdelta)


	if (ntiny.GT.0) then

		! Correct for the fact that smdelta is not exactly the step width in the tiny-rest integration:
		correction=((pp(sproc,2,upperlimit,2)-real(p_res)-capDelta)/(1.d0*ntiny))/smdelta

		ppar=real(p_res)+capDelta

		integrate_res=integrate_res + &
		 2.d0 * correction*(funct_g(ppar,1,om,nn,mode)/(ppar-denomR-ii*denomI))

		integrate_res=integrate_res + &
		correction*(funct_g(ppar,nperp-1,om,nn,mode)/(ppar-denomR-ii*denomI))

		ppar=real(p_res)+capDelta+correction*smdelta*ntiny
		! this should be the same as pp(upperlimit)

		integrate_res=integrate_res + &
		2.d0 * correction*(funct_g(ppar,1,om,nn,mode)/(ppar-denomR-ii*denomI))

		integrate_res=integrate_res + &
		correction*(funct_g(ppar,nperp-1,om,nn,mode)/(ppar-denomR-ii*denomI))


		do iperp=2,nperp-2
			do ipar=1,ntiny-1
				ppar=real(p_res)+capDelta+correction*smdelta*ipar

				integrate_res=integrate_res + 4.d0*&
				correction*(funct_g(ppar,iperp,om,nn,mode)/(ppar-denomR-ii*denomI))
			enddo

			ppar=real(p_res)+capDelta

			integrate_res=integrate_res + 2.d0*&
			correction*(funct_g(ppar,iperp,om,nn,mode)/(ppar-denomR-ii*denomI))

			ppar=real(p_res)+capDelta+correction*smdelta*ntiny

			integrate_res=integrate_res + 2.d0*&
			correction*(funct_g(ppar,iperp,om,nn,mode)/(ppar-denomR-ii*denomI))

		enddo

		do ipar=1,ntiny-1
			ppar=real(p_res)+capDelta+correction*smdelta*ipar

			integrate_res=integrate_res + 2.d0 * 2.d0*&
			correction*(funct_g(ppar,1,om,nn,mode)/(ppar-denomR-ii*denomI))

			integrate_res=integrate_res + 2.d0*&
			correction*(funct_g(ppar,nperp-1,om,nn,mode)/(ppar-denomR-ii*denomI))
		enddo

	endif


	integrate_res = 2.d0 * pi * integrate_res * smdelta * dpperp * 0.25d0
	integrate_res = integrate_res + integrate_norm

	return

end function integrate_res






! Linearized integrand WITHOUT the resonance part.
! It is - resU * int_T / kpar but without the denominator.
! It can be evaluated at any real value of pp (ppar_real) between the grid around the resonance.
double complex function funct_g(ppar_real,iperp,om,nn,mode)
	use alps_var,only : npar,pp,ms,qs,kpar,df0,sproc
	implicit none
	integer :: nn, mode,ipar,iperp,ipar_close
	double complex :: om,integrandplus,integrandminus,integrand
	double precision :: ppar_real,dppar

	dppar  = abs(pp(sproc, 2, 2, 2) - pp(sproc, 2, 1, 2))

	ipar_close=0
	! determine the closest ipar (on the left) to this p_res_real:
	do ipar=1,npar-1
		if ((pp(sproc,iperp,ipar+1,2).GT.ppar_real).AND.(pp(sproc,iperp,ipar,2).LE.ppar_real)) then
			ipar_close=ipar
		endif
	enddo

	if (ipar_close.GE.(npar-1)) ipar_close=npar-2
	if (ipar_close.LE.1) ipar_close=2


	! calculate the function on the grid (left and right of ppar_real):

	integrandplus=-qs(sproc) * &
		   (om*df0(sproc, iperp, ipar_close+1, 1) + (kpar / ( ms(sproc)) ) * &
		   (pp(sproc, iperp, ipar_close+1, 1) * df0(sproc, iperp, ipar_close+1, 2) -&
		   pp(sproc, iperp, ipar_close+1, 2) * df0(sproc, iperp, ipar_close+1, 1) ) ) * &
		   int_T(nn,iperp,ipar_close+1,mode)/kpar


	integrand=-qs(sproc) * &
		   (om*df0(sproc, iperp, ipar_close, 1) + (kpar / ( ms(sproc)) ) * &
		   (pp(sproc, iperp, ipar_close, 1) * df0(sproc, iperp, ipar_close, 2) -&
		   pp(sproc, iperp, ipar_close, 2) * df0(sproc, iperp, ipar_close, 1) ) ) * &
		   int_T(nn,iperp,ipar_close,mode)/kpar



	integrandminus=-qs(sproc) * &
		   (om*df0(sproc, iperp, ipar_close-1, 1) + (kpar / ( ms(sproc)) ) * &
		   (pp(sproc, iperp, ipar_close-1, 1) * df0(sproc, iperp, ipar_close-1, 2) -&
		   pp(sproc, iperp, ipar_close-1, 2) * df0(sproc, iperp, ipar_close-1, 1) ) ) * &
		   int_T(nn,iperp,ipar_close-1,mode)/kpar


	funct_g = integrand+ &
		0.5d0*((integrandplus-integrandminus)/dppar) * (ppar_real - pp(sproc,iperp,ipar_close,2))

	return

end function funct_g





double complex function landau_integrate(om, nn, mode)
	use alps_var, only : nperp, pp, pi, ms, qs, kpar, sproc
	use alps_analyt, only: eval_fit
	implicit none
	!Passed
	double complex :: om   !complex frequency
	integer :: nn          !Bessel N
	integer :: mode        !index in T tensor
	!Local
	integer :: iperp
	double precision :: dpperp ,dppar
	double precision :: h !delta p_perp
	double complex :: ii
	double complex :: p_res,dfperp_C,dfpar_C


	ii = cmplx(0.d0,1.d0,kind(1.d0))

	!choose numerical integration method...

	landau_integrate = cmplx(0.d0,0.d0,kind(1.d0))

	dpperp = pp(sproc, 2, 2, 1) - pp(sproc, 1, 2, 1)
	dppar = abs(pp(sproc, 2, 2, 2) - pp(sproc, 2, 1, 2))
	! Landau contour integral:

  ! At iperp=1, we are already missing the part from iperp=0, where we should actually start. Therefore, we use 4 instead of 2 in the trapezoid integration:
	do iperp = 1, nperp-1
		if ((iperp .EQ. 0).or.(iperp .EQ. (nperp -1))) then
		   h = 0.5d0
		else
		   h = 1.d0
		endif

		p_res = (ms(sproc) * (om) - 1.d0*nn * qs(sproc))/kpar

		! Calculate the derivatives of f0 at the complex p_res:

		dfperp_C=(eval_fit(sproc,iperp+1,p_res)-eval_fit(sproc,iperp-1,p_res))/(2.d0*dpperp)
		dfpar_C=(eval_fit(sproc,iperp,p_res+dppar)-eval_fit(sproc,iperp,p_res-dppar))/(2.d0*dppar)


    landau_integrate = landau_integrate - h * int_T_res(nn, iperp, p_res, mode)*&
        (qs(sproc) /abs(kpar)) *( (pp(sproc, iperp, 1, 1) * dfpar_C - &
        p_res * dfperp_C) * kpar  /  ( ms(sproc)) + om*dfperp_C )

	enddo

	landau_integrate = landau_integrate * ii * dpperp * pi * 2.d0 * pi

	return

end function landau_integrate





double complex function int_ee(om)
	use alps_var, only : qs, ms, nperp, npar, pp, pi, df0, sproc
	implicit none
	!Passed
	double complex :: om   !complex frequency
	!Local
	integer :: iperp, ipar !p_perp, p_par index
	double precision :: dpperp, dppar !delta p_perp, delta p_par


	!choose numerical integration method...

	int_ee = cmplx (0.d0, 0.d0,kind(1.d0))

	dpperp = pp(sproc, 2, 2, 1) - pp(sproc, 1, 2, 1)
	dppar  = abs(pp(sproc, 2, 2, 2) - pp(sproc, 2, 1, 2))

  ! At iperp=1, we are already missing the part from iperp=0, where we should actually start. Therefore, we use 4 instead of 2 in the trapezoid integration:

	int_ee = int_ee + &
		   2.d0 * pp(sproc, 1, 1, 2) *&
		   (df0(sproc, 1, 1, 2)*pp(sproc, 1, 1, 1) - &
		   pp(sproc, 1, 1, 1)*df0(sproc, 1, 1, 1))


	int_ee = int_ee + &
		   2.d0 * pp(sproc, 1, npar -1, 2) *&
		   (df0(sproc, 1, npar -1, 2)*pp(sproc, 1, npar -1, 1) - &
		   pp(sproc, 1, npar -1, 2)*df0(sproc, 1, npar -1, 1))


	int_ee = int_ee + &
		   pp(sproc, nperp -1, 1, 2) *&
		   (df0(sproc, nperp -1, 1, 2)*pp(sproc, nperp -1, 1, 1) - &
		   pp(sproc, nperp -1, 1, 2)*df0(sproc, nperp -1, 1, 1))


	int_ee = int_ee + &
		   pp(sproc, nperp -1, npar -1, 2) *&
		   (df0(sproc, nperp -1, npar -1, 2)*pp(sproc, nperp -1, npar -1, 1) - &
		   pp(sproc, nperp -1, npar -1, 2)*df0(sproc, nperp -1, npar -1, 1))

	do iperp = 2, nperp-2
		do ipar = 2, npar-2
			  int_ee = int_ee + &
				   4.d0*(&
				   pp(sproc, iperp, ipar, 2) *&
				   (df0(sproc, iperp, ipar, 2)*pp(sproc, iperp, ipar, 1) - &
				   pp(sproc, iperp, ipar, 2)*df0(sproc, iperp, ipar, 1)) )
		enddo
	enddo


	do ipar = 2, npar-2

    ! At iperp=1, we are already missing the part from iperp=0, where we should actually start. Therefore, we use 4 instead of 2 in the trapezoid integration:
		   int_ee = int_ee + &
				2.d0 * 2.d0*(&
				pp(sproc, 1, ipar, 2) *&
				(df0(sproc, 1, ipar, 2)*pp(sproc, 1, ipar, 1) - &
				pp(sproc, 1, ipar, 2)*df0(sproc, 1, ipar, 1)) )


			int_ee = int_ee + &
				2.d0*(&
				pp(sproc, nperp -1, ipar, 2) *&
				(df0(sproc, nperp -1, ipar, 2)*pp(sproc, nperp -1, ipar, 1) - &
				pp(sproc, nperp -1, ipar, 2)*df0(sproc, nperp -1, ipar, 1)) )
	enddo

	do iperp = 2, nperp-2

		  int_ee = int_ee + &
				2.d0*(&
				pp(sproc, iperp, 1, 2) *&
				(df0(sproc, iperp, 1, 2)*pp(sproc, iperp, 1, 1) - &
				pp(sproc, iperp, 1, 2)*df0(sproc, iperp, 1, 1)) )

			int_ee = int_ee + &
				2.d0*(&
				pp(sproc, iperp, npar -1, 2) *&
				(df0(sproc, iperp, npar -1, 2)*pp(sproc, iperp, npar -1, 1) - &
				pp(sproc, iperp, npar -1, 2)*df0(sproc, iperp, npar -1, 1)))
	enddo



	int_ee = int_ee * 2.d0 * pi * qs(sproc) / ms(sproc)
	int_ee = int_ee * dpperp * dppar * 0.25d0


	return

end function int_ee



!-=-=-=-=-=-=
!Functions for resonant term in integral
!-=-=-=-=-=-=

double complex function resU(om, nn, iperp, ipar)
	use ALPS_var, only : pp, kpar, ms, qs, df0, vA, sproc, relativistic
	implicit none
	!Passed
	integer :: nn          !Bessel N
	integer :: iperp, ipar !p_perp, p_par index
	double complex :: om   !complex frequency
	!Local
	double precision :: gamma !Relativistic Factor


	gamma = 1.d0
	if (relativistic(sproc)) gamma = sqrt((pp(sproc, iperp, ipar, 1)**2 + &
		pp(sproc, iperp, ipar, 2)**2) * vA**2/ms(sproc)**2 +1.d0)

	resU = qs(sproc) * &
	   (om*df0(sproc, iperp, ipar, 1) + (kpar / (gamma * ms(sproc)) ) * &
	   (pp(sproc, iperp, ipar, 1) * df0(sproc, iperp, ipar, 2) -&
	   pp(sproc, iperp, ipar, 2) * df0(sproc, iperp, ipar, 1) ) )/&
	   (gamma * ms(sproc)*om - kpar * pp(sproc, iperp, ipar, 2) - &
	   (1.d0*nn) * qs(sproc)  )

	return

end function resU


!-=-=-=-=-=-=
!Functions for Tij-
!-=-=-=-=-=-=
!Function to pass T_ij into integrator
double complex function int_T(nn, iperp, ipar, mode)
	use ALPS_var, only : pp, kperp, qs, bessel_array, sproc
	implicit none
	!Passed
	integer :: nn          !Bessel N
	integer :: iperp, ipar !p_perp, p_par index
	integer :: mode        !index in T tensor
	!Local
	double precision :: z  !Bessel Argument
	double precision :: bessel 	! Bessel function for nn and z
	double precision :: besselP	! first derivative of Bessel function for nn and z
	double complex :: ii = cmplx(0.d0,1.d0,kind(1.d0))


	!Bessel Fn Argument
	z= kperp/qs(sproc)
	! Look up array of Bessel functions:
	if (nn.LT.0) then
	   bessel=((-1.d0)**nn)*bessel_array(-nn,iperp)
	else
	   bessel=bessel_array(nn,iperp)
	endif

	! determine derivative of Bessel function:
	if (nn.GE.1) then
		besselP = 0.5d0 * (bessel_array(nn-1,iperp)-bessel_array(nn+1,iperp))
	else if (nn.LT.-1) then
		besselP = 0.5d0 * ((((-1.d0)**(nn-1))*bessel_array(-(nn-1),iperp))&
			-(((-1.d0)**(nn+1))*bessel_array(-(nn+1),iperp)))
	else if (nn.EQ.0) then
	   besselP = -bessel_array(1,iperp) !NIST 10.6.2- yields same result.
	else if (nn.EQ.-1) then
		besselP = 0.5d0 * (bessel_array(2,iperp)-bessel_array(0,iperp))
	endif


	select case(mode)

	  case(1) !T xx
		 int_T = 1.d0 * (nn * nn) * bessel * bessel / (z * z)

	  case(2) !T yy
		 int_T = (pp(sproc, iperp, ipar, 1)**2)*besselP * besselP

	  case(3) !T zz
		 int_T = bessel * bessel * pp( sproc, iperp, ipar, 2)**2

	  case(4) !T xy
		 int_T = (pp(sproc, iperp, ipar, 1))*ii*(1.d0 * (nn)) * bessel * besselP / z

	  case(5) !T xz
		 int_T = (1.d0 * nn) * bessel * bessel* pp(sproc, iperp, ipar, 2)  / z

	  case(6) !T yz
		 int_T = (-1.d0 * ii) * bessel * besselP * pp(sproc, iperp, ipar, 2) *pp(sproc, iperp, ipar,1)

	end select

	return

end function int_T
!++++

!-=-=-=-=-=-=
!Functions for Tij-
!-=-=-=-=-=-=
!Function to pass T_ij into integrator
double complex function int_T_res(nn, iperp, p_res, mode)
	use ALPS_var, only : pp, kperp, qs, bessel_array,sproc
	implicit none
	!Passed
	integer :: nn          !Bessel N
	integer :: iperp !p_perp, p_par index
	integer :: mode        !index in T tensor
	double complex :: p_res  !Bessel Argument
	!Local
	double precision :: z  !Bessel Argument
	double precision :: bessel 	! Bessel function for nn and z
	double precision :: besselP	! first derivative of Bessel function for nn and z
	double complex :: ii = cmplx(0.d0,1.d0,kind(1.d0))

	!Bessel Fn Argument
	z = kperp / qs(sproc)

	! Look up array of Bessel functions:
	if (nn.LT.0) then
	   bessel=((-1.d0)**(nn))*bessel_array(-nn,iperp)
	else
	   bessel=bessel_array(nn,iperp)
	endif


	! determine derivative of Bessel function:
	if (nn.GE.1) then
		besselP = 0.5d0 * (bessel_array(nn-1,iperp)-bessel_array(nn+1,iperp))
	else if (nn.LT.-1) then
		besselP = 0.5d0 * ((((-1.d0)**(nn-1))*bessel_array(-(nn-1),iperp))&
			-(((-1.d0)**(nn+1))*bessel_array(-(nn+1),iperp)))
	else if (nn.EQ.0) then
	   besselP = -bessel_array(1,iperp) !NIST 10.6.2- yields same result.
	else if (nn.EQ.-1) then
		besselP = 0.5d0 * (bessel_array(2,iperp)-bessel_array(0,iperp))
	endif

	select case(mode)

	  case(1) !T xx
		 int_T_res = 1.d0 * (nn * nn) * bessel * bessel / (z * z)
	  case(2) !T yy
		 int_T_res = (pp(sproc, iperp, 1, 1)**2)*besselP * besselP
	  case(3) !T zz
		 int_T_res = bessel * bessel * p_res**2
	  case(4) !T xy
		 int_T_res = (pp(sproc, iperp, 1, 1))*ii*(1.d0 * (nn)) * bessel * besselP  /z
	  case(5) !T xz
		 int_T_res = (1.d0 * nn) * bessel * bessel* p_res /z
	  case(6) !T yz
		 int_T_res = (-1.d0 * ii) * bessel * besselP * p_res*pp(sproc, iperp, 1, 1)
	end select

	return

end function int_T_res
!++++
!++++

!-=-=-=-=-=-=
!Secant method
!-=-=-=-=-=-=
subroutine secant(om)
	use ALPS_var, only : numiter, D_threshold, ierror, proc0, writeOut, D_prec
   use mpi
	implicit none

	double complex :: om, prevom, ii, D, Dprev, jump, minom, minD
	integer :: iter
	logical :: go_for_secant

	ii=(0.d0,1.d0)

	!prevom = om - 1.d-5 - 1.d-5*ii
	prevom=om*(1.d0-D_prec)
	Dprev = disp(prevom)
	minD=Dprev
	minom=prevom
	call mpi_barrier(mpi_comm_world,ierror)

	iter = 0
	go_for_secant = .TRUE.

	do while ((iter.LE.(numiter-1)).AND.(go_for_secant))
		iter = iter + 1
		D = disp(om)

	 if ((abs(D-Dprev).LT.1.d-80)) then
		prevom = prevom + 1.d-8
		Dprev = disp(prevom)
	 endif

		if ((abs(D).LT.D_threshold)) then
			jump = 0.d0
			go_for_secant = .FALSE.
			if (proc0.AND.writeOut) then
				write(*,'(a,i4)') ' Converged after iteration ',iter
				write(*,'(a,2es14.4,a,2es14.4)') ' D(',real(om),aimag(om),')= ',D
			endif

		else
			jump = D*(om-prevom)/(D-Dprev)
		endif
		prevom = om
		om = om - jump
		Dprev = D
		if (abs(D).LT.abs(minD)) then
			minom=om
			minD=D
		endif
	enddo

	if (proc0.AND.writeOut.AND.(iter.GE.numiter)) then
		write(*,'(a,i4,a)') ' Maximum iteration ',iter,' reached.'
		om=minom
		write(*,'(a,2es14.4,a,2es14.4)') ' D(',real(om),aimag(om),')= ',D
	endif

end subroutine secant


!-=-=-=-=-=-=
!Routine for scanning along a single perscribed path in wavevector space
!-=-=-=-=-=-=
!KGK:
subroutine om_scan(ik)
  use ALPS_var, only : proc0, nroots, runname, ierror, wroots, scan, sproc
  use ALPS_var, only : kperp,kpar,kperp_last,kpar_last, D_prec, D_gap
  use ALPS_var, only : use_secant, nspec
  use ALPS_io,  only : get_unused_unit, isnancheck, alps_error
  use mpi
  implicit none

  !Passed
  integer :: ik !Scan Number

  !Local
  integer :: it, nt !step of scan, number of scans
  integer :: in !number of roots
  character(100), dimension(:), allocatable :: scanName  !Output file name
  character(100), dimension(:), allocatable :: heatName  !Output file name
  character(100), dimension(:), allocatable :: eigenName  !Output file name
  character(6) :: scan_ID
  double precision :: theta_0, theta_1, k_0


  double complex :: omega    !Complex Frequency
  double complex :: om1, om2 !Bounding frequencies for secant search
  integer :: iflag                           !Flag for Root search
  integer, dimension(:), allocatable :: scan_unit
  integer, dimension(:), allocatable :: heat_unit  ! file for outputting damping rates
  integer, dimension(:), allocatable :: eigen_unit ! file for outputting eigenfn values
  integer :: imm
  logical, dimension(:),allocatable :: jump
  logical :: alljump
  double complex :: tmp

      !Eigenfunctions
  double complex, dimension(1:3)       :: ef, bf !E, B
  double complex, dimension(1:nspec)     :: ds     !density
  double complex, dimension(1:3,1:nspec) :: Us     !Velocity
  !Heating
  double precision, dimension(1:nspec) :: Ps !Power into/out of species
  character (50) :: fmt_eigen, fmt_heat   !Output format

  allocate(jump(1:nroots));jump=.true.

  if (proc0) then
     allocate(scan_unit(nroots))
     allocate(scanName(nroots))
     select case(scan(ik)%type_s)
     case (0) !k_0 to k_1
        write(scan_ID,'(a)')'k1_k2_'
     case (1) !theta_0 to theta_1
        write(scan_ID,'(a)')'theta_'
     case (2) !|k_0| to |k_1| @ constant theta
        write(scan_ID,'(a)')'kcstq_'
     case (3) !kperp scan
        write(scan_ID,'(a)')'kperp_'
     case (4) !kpar scan
        write(scan_ID,'(a)')'kpara_'
     end select
          if (scan(ik)%eigen_s) then
        write(fmt_eigen,'(a,i0,a)') '(4es14.4,12es14.4,',nspec*8,'es14.4)'
        allocate(eigen_unit(nroots))
        allocate(eigenName(nroots))
     endif
     if (scan(ik)%heat_s) then
        write(fmt_heat,'(a,i0,a)') '(4es14.4,',nspec,'es14.4)'
        allocate(heat_unit(nroots))
        allocate(heatName(nroots))
     endif
     do in=1,nroots
        write(scanName(in),'(4a,i0,a,i0)')&
             'solution/',trim(runname),'.scan_',scan_ID,ik,'.root_',in
        write(*,'(2a)')' => ',trim(scanName(in))
        call get_unused_unit(scan_unit(in))
        open(unit=scan_unit(in),file=trim(scanName(in)),status='replace')
        write(scan_unit(in),'(4es14.4)') &
             kperp,kpar,wroots(in)
        close(scan_unit(in))
     enddo
  endif

  if ((scan(ik)%eigen_s).or.(scan(ik)%heat_s)) then

     do in=1,nroots
        omega=wroots(in)
        tmp = disp(omega)

        call calc_eigen(omega,ef,bf,Us,ds,Ps,scan(ik)%eigen_s,scan(ik)%heat_s)

        !reassign omega
        omega=wroots(in)
        !this is necessary, as calc_eigen evaluates
        !the wave tensor at several different values of omega

        call mpi_barrier(mpi_comm_world,ierror)

        if (proc0) then
           if (scan(ik)%eigen_s) then
              write(eigenName(in),'(4a,i0,a,i0)')&
                   'solution/',trim(runname),'.eigen_',scan_ID,ik,'.root_',in
              write(*,'(2a)')' => ',trim(eigenName(in))
              call get_unused_unit(eigen_unit(in))
              open(unit=eigen_unit(in),file=trim(eigenName(in)),status='replace')
              write(eigen_unit(in),trim(fmt_eigen)) &
                   kperp,kpar,wroots(in),ef,bf,Us,ds
              close(eigen_unit(in))
           endif

           if (scan(ik)%heat_s) then
              write(heatName(in),'(4a,i0,a,i0)')&
                   'solution/',trim(runname),'.heat_',scan_ID,ik,'.root_',in
              write(*,'(2a)')' => ',trim(heatName(in))
              call get_unused_unit(heat_unit(in))
              open(unit=heat_unit(in),file=trim(heatName(in)),status='replace')
              write(heat_unit(in),trim(fmt_heat)) &
                   kperp,kpar,wroots(in),Ps
              close(heat_unit(in))
           endif

        endif
     enddo
  endif

  nt = scan(ik)%n_out*scan(ik)%n_res

  kperp_last=kperp;kpar_last=kpar

  theta_0=atan(kperp_last/kpar_last)
  k_0=sqrt(kperp_last**2+kpar_last**2)

  do it = 1, nt
     !Scan through wavevector
     select case(scan(ik)%type_s)
     case (0) !k_0 to k_1
        if (scan(ik)%log_scan) then
           kperp=10.d0**(log10(kperp_last)+scan(ik)%diff*it)
           kpar=10.d0**(log10(kpar_last)+scan(ik)%diff2*it)
        else
           kperp=kperp_last+scan(ik)%diff*it
           kpar= kpar_last +scan(ik)%diff2*it
        endif
     case (1) !theta_0 to theta_1
        if (scan(ik)%log_scan) then
           theta_1=10.d0**(log10(theta_0)+scan(ik)%diff*it)
        else
           theta_1=theta_0+scan(ik)%diff*it
        endif
        kperp=k_0*sin(theta_1)
        kpar=k_0*cos(theta_1)
     case (2) !|k_0| to |k_1| @ constant theta
        if (scan(ik)%log_scan) then
           kperp=10.d0**(log10(kperp_last)+scan(ik)%diff*it)
           kpar=10.d0**(log10(kpar_last)+scan(ik)%diff2*it)
        else
           kperp=kperp_last+scan(ik)%diff*it
           kpar= kpar_last +scan(ik)%diff2*it
        endif
     case (3) !kperp scan
        if (scan(ik)%log_scan) then
           kperp=10.d0**(log10(kperp_last)+scan(ik)%diff*it)
        else
           kperp=kperp_last+scan(ik)%diff*it
        endif
     case (4) !kpar scan
        if (scan(ik)%log_scan) then
           kpar=10.d0**(log10(kpar_last)+scan(ik)%diff*it)
        else
           kpar=kpar_last+scan(ik)%diff*it
        endif
     end select

     if (scan(ik)%type_s.ne.4) then
        ! Once we know kperp, we can determine nmax and split the processes:
        ! These three routines must be called when kperp changes
        call determine_nmax
        call split_processes

        ! All processes determine their Bessel function array:
        if(.NOT.(sproc.EQ.0)) call determine_bessel_array
     endif

     call mpi_barrier(mpi_comm_world,ierror)

     if (proc0) write(*,'(a,es14.4,a,es14.4)')'kperp: ',kperp,' kpar: ',kpar

	! Check if all jumps are set to .false.:
	alljump=.FALSE.
	do in = 1,nroots
		alljump=alljump.OR.jump(in)
	enddo
	if (alljump.EQV..FALSE.) call alps_error(9)


     do in = 1,nroots
        !Search for new roots

        if (jump(in)) then

           omega=wroots(in)

           if (use_secant) then
              call secant(omega)
           else
               om1=omega*(1.d0-D_prec)
        	   om2=omega*(1.d0+D_prec)
               omega=rtsec(disp,om1,om2,iflag)
           endif

           wroots(in)=omega

           call mpi_bcast(wroots(in), 1, &
                MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)

           tmp = disp(omega)

           call mpi_barrier(mpi_comm_world,ierror)

                      !Calculate eigenfunctions and heating rates
           !KGK: 200526

           !only call on wavevector steps that will be output
           if (mod(it,scan(ik)%n_res)==0) then

              if ((scan(ik)%eigen_s).or.((scan(ik)%heat_s))) then
                 !the susceptability tensor, chi0
                 !and wave tensor, wave
                 !are set by the
                 !tmp = disp(omega)
                 !call above mpi_barrier

                 call calc_eigen(omega,ef,bf,Us,ds,Ps,scan(ik)%eigen_s,scan(ik)%heat_s)

                 !reassign omega
                 omega=wroots(in)
                 !this is necessary, as calc_eigen evaluates
                 !the wave tensor at several different values of omega

              endif
           endif
           call mpi_barrier(mpi_comm_world,ierror)

           !Output and check for root jumps and NaNs
           if (proc0) then
              !NaN Check
!              tmp=real(omega)
              if(isnancheck(real(omega))) then
              	  omega=cmplx(0.d0,0.d0,kind(1.d0));jump(in)=.false.
              endif
!              if (.not.(tmp .ne. omega)) then
!                 omega=cmplx(0.d0,0.d0);jump(in)=.false.
!              endif
!              !infty Check
!              if (abs(tmp) .gt. 1.d100) then
!                 omega=cmplx(0.d0,0.d0);jump(in)=.false.
!              endif
              !compare to previous roots
	      !KGK: 200811: Updated to separately compare real and imaginary components of roots
	      !previous version would reject all roots except the first solution for scans of
	      !relatvely small wavevectors
              do imm=1,in-1
                 if (abs(wroots(in)-wroots(imm)).lt.D_gap) then
                    write(*,'(a,6es14.4)')'Root too close!',&
                         wroots(in),wroots(imm),&
                         real(wroots(in))-real(wroots(imm)), &
                         aimag(wroots(in))-aimag(wroots(imm))
                    wroots(in)=cmplx(0.d0,0.d0);jump(in)=.false.
                 endif
              enddo

              if (mod(it,scan(ik)%n_res)==0) then
                 open(unit=scan_unit(in),file=trim(scanName(in)),status='old',position='append')
                 write(scan_unit(in),'(4es14.4)') &
                      kperp,kpar,wroots(in)
                 close(scan_unit(in))

                 if (scan(ik)%eigen_s) then
                    open(unit=eigen_unit(in),file=trim(eigenName(in)),status='old',position='append')
                    write(eigen_unit(in),trim(fmt_eigen)) &
                         kperp,kpar,wroots(in),ef,bf,Us,ds
                    close(eigen_unit(in))
                 endif

                 if (scan(ik)%heat_s) then
                    open(unit=heat_unit(in),file=trim(heatName(in)),status='old',position='append')
                    write(heat_unit(in),trim(fmt_heat)) &
                         kperp,kpar,wroots(in),Ps
                    close(heat_unit(in))
                 endif

              endif
           endif
           call mpi_bcast(jump(in),1,MPI_LOGICAL,0,MPI_COMM_WORLD, ierror)

        end if
        call mpi_barrier(mpi_comm_world,ierror)

     enddo
  enddo

  if (proc0) then
     deallocate(scan_unit)
     deallocate(scanName)
  endif

end subroutine om_scan

!-=-=-=-=-=-=
!  Calculates the electric and magnetic fields as well as species
!     velocities and density fluctuations for (omega,gamma)
!     and particle heating/cooling from a given wave
!-=-=-=-=-=-=
!KGK: 200522
!Based upon calc_eigen routine by Greg Howes and Kris Klein
!found in the PLUME linear dispersion solver
subroutine calc_eigen(omega,electric,magnetic,vmean,ds,Ps,eigen_L,heat_L)
  use ALPS_var, only : proc0, nspec, ns, qs, wave, chi0, kperp, kpar, vA

  implicit none
  !Passed
  !Frequency
  double complex :: omega
  !Electromagnetic Eigenfns
  double complex, dimension(1:3), intent(out)       :: electric, magnetic !E, B
  double complex, dimension(1:nspec), intent(out)     :: ds     !density fluctuation
  double complex, dimension(1:3,1:nspec), intent(out) :: vmean     !Velocity fluctiation

  !Heating
  double precision, dimension(1:nspec), intent(out) :: Ps   !Power into/out of species
  logical :: eigen_L, heat_L !Logical for calculating eigenvalues, heating

  !Local
  integer :: ii,j,jj
  double complex :: temp1
  double complex, dimension(nspec,3,3) :: chia
  double complex, dimension(3,3) :: chih,chihold,dchih
  double complex, dimension(nspec,3) :: term
  double complex, dimension(3) :: term1
  double precision :: ewave

  if (proc0) then

     !Note that the electric and magnetic fields are needed for the heating
     !rate calculation; thus, we calculate them for both the eigen and heating
     !logical flags
        !CALCULATE FIELDS FLUCTUATIONS==========================================
        !Calculate Electric Fields, normalized to E_x
        electric(1) = cmplx(1.d0,0.d0,kind(1.d0))
        electric(3)=-electric(1)*(wave(2,1)*wave(3,2)-wave(3,1)*wave(2,2))
        electric(3)= electric(3)/(wave(2,3)*wave(3,2)-wave(3,3)*wave(2,2))
        electric(2) = -electric(3)*wave(3,3) - electric(1)*wave(3,1)
        electric(2) = electric(2)/wave(3,2)

        !Calculate Magnetic Fields, normalized to E_x
        magnetic(1) = -1.d0* kpar*electric(2)/(omega*vA)
        magnetic(2) = -1.d0* (kperp*electric(3) - kpar*electric(1))/(omega*vA)
        magnetic(3) = kperp*electric(2)/(omega*vA)

        !KGK: The magnetic field normalization factors are different from PLUME,
        !as spatial scales are normalized to d_ref, rather than rho_ref.
        !thus, w_perp, ref/c in PLUME becomes v_A/c here.

        if (eigen_L) then

        !CALCULATE VELOCITY FLUCTUATIONS========================================
        !vmean is the velocity perturbutation due to the wave for each species
        !vmean = (delta V_s/v_A)(B_0/E_x)(v_A/c)
        vmean(:,:)=0.d0
        do j=1,3!x,y,z
           do jj = 1,nspec !Species velocity fluctuations
              vmean(j,jj) = -(vA**2.d0/(qs(jj)*ns(jj)))*&
                   cmplx(0.d0,1.d0,kind(1.d0))*&
                   omega*sum(electric(:)*chi0(jj,j,:))

           enddo
        enddo

        !CALCULATE DENSITY FLUCTUATIONS========================================
        ! This is (ns/ns0)(B_0/E_x)(v_A/c)
        do jj=1,nspec
           !ds(jj) = (vmean(1,jj)*kperp+vmean(3,jj)*kpar)/&
           !     (omega-kpar * spec(jj)%vv_s)

           !TO-DO:
           !add in correct effects of drift on density fluctuation
           !e.g. omega-kpar V_0
           !with a physically meaningful V_0 used.
           ds(jj) = 0.d0
           !ds(jj) = (vmean(1,jj)*kperp+vmean(3,jj)*kpar)/&
           !     (omega)
        enddo

        !EndIf (scan(is)%eigen_s) loop
     endif
  endif

!If (scan(is)%heat_s) loop
!Greg Howes, 2006; Kristopher Klein, 2015
if (heat_L) then
   !CALCULATE COMPONENT HEATING======================================
   !evaulate at omega_r=real(omega), gamma=0
   temp1 = cmplx(real(omega),0.d0,kind(1.d0))
   !temp1 = omega
   temp1 = disp(temp1)

   if (proc0) then

      do ii = 1, 3 !tensor index
         do j = 1, 3 !tensor index
            do jj = 1, nspec !species index
               chia(jj,ii,j) = -0.5d0*cmplx(0.d0,1.d0)* &
                    (chi0(jj,ii,j) - conjg(chi0(jj,j,ii)))
            enddo
            chihold(ii,j) = 0.5*(sum(chi0(:,ii,j)) + &
                 sum(conjg(chi0(:,j,ii))))
         enddo
      enddo

      term(:,:)=0.d0
      term1(:)=0.d0
      do ii = 1, 3
         do jj = 1, nspec
            term(jj,ii) = sum(conjg(electric(:))*chia(jj,:,ii))
         enddo
      enddo

      Ps = 0.d0
      do jj = 1, nspec
         Ps(jj) = sum(term(jj,:)*electric(:))
      enddo

   endif

   !recall that disp requires /all/ processors
   temp1 = disp(cmplx(real(omega*1.000001d0),0.d0,kind(1.d0)))
   !but only proc0 'knows' the correct value of chi0

   if (proc0) then
      do ii = 1, 3
         do j = 1, 3
            chih(ii,j) = 0.5d0*(sum(chi0(:,ii,j)) + &
                 sum(conjg(chi0(:,j,ii))))
            dchih(ii,j)=(1.000001d0*chih(ii,j)-chihold(ii,j))/0.000001d0
         enddo
      enddo

      ewave = 0.d0
      do ii = 1, 3
         term1(ii) = sum(conjg(electric(:))*dchih(:,ii))
      enddo

      ewave = sum(term1(:)*electric(:)) + sum(magnetic(:)*conjg(magnetic(:)))

      Ps = Ps/ewave
   endif

endif    !EndIf (scan(is)%heat_s) loop

end subroutine calc_eigen

!-=-=-=-=-=-=
!Routine for scanning along a perscribed plane in wavevector space
!Must have n_scan=2
!-=-=-=-=-=-=
!KGK:
subroutine om_double_scan
  use ALPS_var, only : proc0, nroots, runname, ierror, wroots, scan, sproc
  use ALPS_var, only : kperp,kpar,kperp_last,kpar_last, D_prec, D_gap
  use ALPS_var, only : use_secant, ierror
  use ALPS_io,  only : get_unused_unit, alps_error
  use mpi
  implicit none

  !Local
  integer :: it, it2, nt, nt2 !step of scan, number of scans
  integer :: in !number of roots
  character(100), dimension(:), allocatable :: scanName  !Output file name
  character(6) :: scan_ID
  character(5) :: scan_ID2
  double precision :: theta_0, theta_1, k_0, theta_i, k_i
  double precision :: kperpi,kpari
  double precision, dimension(:), allocatable :: om_tmp


  double complex :: omega    !Complex Frequency
  double complex :: om1, om2 !Bounding frequencies for secant search
  integer :: iflag                           !Flag for Root search
  integer, dimension(:), allocatable :: scan_unit

  double precision :: tmp
  integer :: imm
  logical, dimension(:),allocatable :: jump

  allocate(jump(1:nroots));jump=.true.

  if (scan(1)%type_s==scan(2)%type_s) then
     call alps_error(5)
  endif

  if ((scan(1)%type_s==0).or.(scan(2)%type_s==0)) then
     call alps_error(6)
  endif

  allocate(om_tmp(nroots))

  if (proc0) then
     allocate(scan_unit(nroots))
     allocate(scanName(nroots))
     !Name of first scan
     select case(scan(1)%type_s)
     case (0) !k_0 to k_1
        write(scan_ID,'(a)')'k1_k2_'
     case (1) !theta_0 to theta_1
        write(scan_ID,'(a)')'theta_'
     case (2) !|k_0| to |k_1| @ constant theta
        write(scan_ID,'(a)')'kcstq_'
     case (3) !kperp scan
        write(scan_ID,'(a)')'kperp_'
     case (4) !kpar scan
        write(scan_ID,'(a)')'kpara_'
     end select

     !Name of second scan
     select case(scan(2)%type_s)
     case (0) !k_0 to k_1
        write(scan_ID2,'(a)')'k1_k2'
     case (1) !theta_0 to theta_1
        write(scan_ID2,'(a)')'theta'
     case (2) !|k_0| to |k_1| @ constant theta
        write(scan_ID2,'(a)')'kcstq'
     case (3) !kperp scan
        write(scan_ID2,'(a)')'kperp'
     case (4) !kpar scan
        write(scan_ID2,'(a)')'kpara'
     end select

     do in=1,nroots
        write(scanName(in),'(6a,i0)')&
             'solution/',trim(runname),'.scan_',scan_ID,scan_ID2,'.root_',in
        write(*,'(2a)')' => ',trim(scanName(in))
        call get_unused_unit(scan_unit(in))
        open(unit=scan_unit(in),file=trim(scanName(in)),status='replace')
        write(scan_unit(in),'(4es14.4)') &
             kperp,kpar,wroots(in)
        close(scan_unit(in))
     enddo
  endif

  nt = scan(1)%n_out*scan(1)%n_res
  nt2 = scan(2)%n_out*scan(2)%n_res

  kperp_last=kperp;kpar_last=kpar

  theta_0=atan(kperp_last/kpar_last)
  k_0=sqrt(kperp_last**2+kpar_last**2)

  do it = 1, nt
     !Scan through wavevector
     select case(scan(1)%type_s)
     case (0) !k_0 to k_1
        if (scan(1)%log_scan) then
           kperp=10.d0**(log10(kperp_last)+scan(1)%diff*it)
           kpar=10.d0**(log10(kpar_last)+scan(1)%diff2*it)
        else
           kperp=kperp_last+scan(1)%diff*it
           kpar= kpar_last +scan(1)%diff2*it
        endif
     case (1) !theta_0 to theta_1
        if (scan(1)%log_scan) then
           theta_1=10.d0**(log10(theta_0)+scan(1)%diff*it)
        else
           theta_1=theta_0+scan(1)%diff*it
        endif
        kperp=k_0*sin(theta_1)
        kpar=k_0*cos(theta_1)
     case (2) !|k_0| to |k_1| @ constant theta
        if (scan(1)%log_scan) then
           kperp=10.d0**(log10(kperp_last)+scan(1)%diff*it)
           kpar=10.d0**(log10(kpar_last)+scan(1)%diff2*it)
        else
           kperp=kperp_last+scan(1)%diff*it
           kpar= kpar_last +scan(1)%diff2*it
        endif
     case (3) !kperp scan
        if (scan(1)%log_scan) then
           kperp=10.d0**(log10(kperp_last)+scan(1)%diff*it)
        else
           kperp=kperp_last+scan(1)%diff*it
        endif
     case (4) !kpar scan
        if (scan(1)%log_scan) then
           kpar=10.d0**(log10(kpar_last)+scan(1)%diff*it)
        else
           kpar=kpar_last+scan(1)%diff*it
        endif
     end select

     if (.true.) then

     if (scan(1)%type_s.ne.4) then
        ! Once we know kperp, we can determine nmax and split the processes:
        ! These three routines must be called when kperp changes
        call determine_nmax
        call split_processes

        ! All processes determine their Bessel function array:
        if(.NOT.(sproc.EQ.0)) call determine_bessel_array
     endif

     call mpi_barrier(mpi_comm_world,ierror)

     endif

     if (proc0) write(*,'(a,es14.4,a,es14.4)')'kperp: ',kperp,' kpar: ',kpar
     om_tmp=wroots
     !SECOND SCAN
     do it2 = 1, nt2
        !Scan through wavevector
        if (it2==1) then
           kperpi=kperp; kpari=kpar
           kperpi=kperp; kpari=kpar
           theta_i=atan(kperpi/kpari)
           k_i=sqrt(kperpi**2+kpari**2)
           wroots=om_tmp
        endif
        select case(scan(2)%type_s)
        case (0) !k_0 to k_1
           if (scan(2)%log_scan) then
              kperp=10.d0**(log10(kperpi)+scan(2)%diff*it2)
              kpar=10.d0**(log10(kpari)+scan(2)%diff2*it2)
           else
              kperp=kperpi+scan(2)%diff*it2
              kpar= kpari +scan(2)%diff2*it2
           endif
        case (1) !theta_i to theta_1
           if (scan(2)%log_scan) then
              theta_1=10.d0**(log10(theta_i)+scan(2)%diff*it2)
           else
              theta_1=theta_i+scan(2)%diff*it2
           endif
           kperp=k_i*sin(theta_1)
           kpar=k_i*cos(theta_1)
        case (2) !|k_0| to |k_1| @ constant theta
           if (scan(2)%log_scan) then
              kperp=10.d0**(log10(kperpi)+scan(2)%diff*it2)
              kpar=10.d0**(log10(kpari)+scan(2)%diff2*it2)
           else
              kperp=kperpi+scan(2)%diff*it2
              kpar= kpari +scan(2)%diff2*it2
           endif
        case (3) !kperp scan
           if (scan(2)%log_scan) then
              kperp=10.d0**(log10(kperpi)+scan(2)%diff*it2)
           else
              kperp=kperpi+scan(2)%diff*it2
           endif
        case (4) !kpar scan
           if (scan(2)%log_scan) then
              kpar=10.d0**(log10(kpari)+scan(2)%diff*it2)
           else
              kpar=kpari+scan(2)%diff*it2
           endif
        end select

        if (.true.) then

        if (scan(2)%type_s.ne.4) then
           ! Once we know kperp, we can determine nmax and split the processes:
           ! These three routines must be called when kperp changes
           call determine_nmax
           call split_processes

           ! All processes determine their Bessel function array:
           if(.NOT.(sproc.EQ.0)) call determine_bessel_array
        endif

        call mpi_barrier(mpi_comm_world,ierror)

     endif

          if (proc0) write(*,'(a,es14.4,a,es14.4)')'kperp: ',kperp,' kpar: ',kpar

          do in = 1,nroots
             !Search for new roots

             if (jump(in)) then


                omega=wroots(in)

                om1=omega*(1.d0-D_prec)
                om2=omega*(1.d0+D_prec)

                if (use_secant) then
                   call secant(omega)
                else
                   omega=rtsec(disp,om1,om2,iflag)
                endif


                wroots(in)=omega

                call mpi_bcast(wroots(in), 1, &
                     MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)

                tmp = disp(omega)

                call mpi_barrier(mpi_comm_world,ierror)

           !Output...
           if (proc0) then
              !NaN Check
              tmp=real(omega)
              if (.not.(tmp .ne. omega)) then
                 omega=cmplx(0.d0,0.d0);jump(in)=.false.
              endif
              !infty Check
              if (abs(tmp) .gt. 1.d100) then
                 omega=cmplx(0.d0,0.d0);jump(in)=.false.
              endif


              do imm=1,in-1
                 if (abs(wroots(in)-wroots(imm)).lt.D_gap) then
                    write(*,'(a,6es14.4)')'Root too close!',&
                         wroots(in),wroots(imm),&
                         real(wroots(in))-real(wroots(imm)), &
                         aimag(wroots(in))-aimag(wroots(imm))
                    wroots(in)=cmplx(0.d0,0.d0);jump(in)=.false.
                 endif
              enddo
              !Aplus~=Aminus
              if (mod(it2,scan(2)%n_res)==0) then
                 open(unit=scan_unit(in),file=trim(scanName(in)),&
                      status='old',position='append')
                 write(scan_unit(in),'(4es14.4)') &
                      kperp,kpar,wroots(in)
                 close(scan_unit(in))
              endif
           endif
           call mpi_bcast(jump(in),1,MPI_LOGICAL,0,MPI_COMM_WORLD, ierror)

        endif
           call mpi_barrier(mpi_comm_world,ierror)

        enddo !roots

     enddo !second scan
     if (proc0) then
        do in = 1,nroots
           open(unit=scan_unit(in),file=trim(scanName(in)),&
                status='old',position='append')
           write(scan_unit(in),*)
           close(scan_unit(in))
        enddo
     endif
     kperp=kperp_last;kpar=kpar_last
  enddo !first scan

  if (proc0) then
     deallocate(scan_unit)
     deallocate(scanName)
  endif

  deallocate(om_tmp)

end subroutine om_double_scan

!-=-=-=-=-=-=
!Routine to find solutions over a region of complex frequency space-
!-=-=-=-=-=-=
!KGK:
subroutine map_search
  use ALPS_var, only : ierror
  use ALPS_var, only : omi, omf, gami, gamf, loggridw, loggridg, determine_minima
  use ALPS_var, only : ni, nr, proc0, ms, ns, qs, runname, nspec
  use ALPS_var, only : writeOut, kperp, kpar, wroots, numroots, nroots, nroots_max
  use ALPS_io,  only : get_unused_unit
  use mpi
  implicit none

  double precision :: dr, di                         !Spacing
  double precision :: wr, wi                         !Real,imaginary omega
  double precision, dimension(:,:), pointer :: val   !Value of Dispersion relation
  double complex, dimension(:,:), allocatable :: cal   ! (complex) Value of Dispersion relation
  double complex, dimension(:,:), allocatable :: om     !Complex Frequency
  double complex :: omega    !Complex Frequency
  integer :: ii, ir          !indices for frequency grids
  integer :: is              !species index
  integer :: iw              !root index
  character(100) :: mapName                  !Output file names
  integer, dimension(1:2,1:numroots) :: iroots !Indices of roots in local min. search
  integer :: unit_map
  double precision :: tmp

  if (writeOut .and. proc0.and. .true.) then
     write(*,'(a)')'-=-=-=-=-=-=-=-=-=-'
     write(*,'(a)')      'Global Plasma Parameters:'
     write(*,'(a,es12.3)')' k_perp d_p   = ',kperp
     write(*,'(a,es12.3)')' k_par  d_p   = ',kpar
     do is = 1, nspec
        write(*,'(a)')'-=-=-=-=-=-=-=-=-=-'
        write(*,'(a,i3)')      'Parameters for species',is
        write(*,'(a,es12.3)')' m_s/m_m =        ',ms(is)
        write(*,'(a,es12.3)')' q_s/q_p =        ',qs(is)
        write(*,'(a,es12.3)')' n_s/n_p =        ',ns(is)
     enddo
     write(*,'(a)')'-=-=-=-=-=-=-=-=-=-'
     write(*,'(a)')'Searching over:'
     write(*,'(a,es12.3,a,es12.3,a)')' om  in [',omi,',',omf,']'
     write(*,'(a,es12.3,a,es12.3,a)')' gam in [',gami,',',gamf,']'
     write(*,'(a)')'-=-=-=-=-=-=-=-=-=-'
  endif

  !Allocate array for map values
  !Value of dispersion relation on frequency grid
  allocate(cal(nr,ni)); cal(:,:)=cmplx(0.d0,0.d0,kind(1.d0))
  !magnitude of cal
  allocate(val(nr,ni)); val(:,:)=0.d0
  !Array of complex frequencies
  allocate(om(nr,ni)); om(:,:)=cmplx(0.d0,0.d0,kind(1.d0))

  !Determine spacing in complex omega space (Normal or log)
  dr=omf-omi
  di=gamf-gami
  if (nr.GT.1) dr=(omf-omi)/(1.d0*(nr-1))
  if (ni.GT.1) di=(gamf-gami)/(1.d0*(ni-1))

  if (proc0) then
     write(mapName,'(3a)') 'solution/',trim(runname),'.map'
     call get_unused_unit(unit_map)
     open(unit=unit_map,file=trim(mapName),status='replace')
     close(unit_map)
  endif

  !Scan over complex frequency space and calculate dispersion relation
  do ir=1,nr
     if (loggridw) then
        wr=omi
        if (nr.GT.1) wr=omi*((omf/omi)**((1.d0*(ir-1))/(1.d0*(nr-1))))
     else
        wr=omi+dr*(1.d0*(ir-1))
     endif
     if (proc0.and.writeOut)&
          write(*,'(a,es11.4)')' omega_real = ',wr
     do ii=1,ni
        if (loggridg) then
		   wi=gami
           if(ni.GT.1) wi=gami*((gamf/gami)**((1.d0*(ii-1))/(1.d0*(ni-1))))
        else
           wi=gami+di*(1.d0*(ii-1))
        endif
        !if (proc0.and.writeOut)&
        !     write(*,'(a,es11.4)')' omega_aimag = ',wi
        omega=cmplx(wr,wi,kind(1.d0))
        om(ir,ii)=omega
        cal(ir,ii)=disp(omega)
        val(ir,ii)=abs(cal(ir,ii))
        if (proc0) then
           tmp=cal(ir,ii)
           val(ir,ii)=log10(val(ir,ii))
           !NaN Check
           if (aimag(cal(ir,ii)).ne.0.d0) then
              if (.not.(tmp .ne. cal(ir,ii)) ) then
                 !write(*,*) 'fail check one: ',tmp, cal(ir,ii)
                 cal(ir,ii)=999999.d0;val(ir,ii)=999999.d0
              endif
           else
              if ((tmp .ne. cal(ir,ii)) ) then
                 !write(*,*) 'fail check one: ',tmp, cal(ir,ii)
                 cal(ir,ii)=999999.d0;val(ir,ii)=999999.d0
              endif
           endif
           !infty Check
           if (abs(tmp) .gt. 1.d100) then
              cal(ir,ii)=899999.d0;val(ir,ii)=899999.d0
           endif
           open(unit=unit_map,file=trim(mapName),status='old',position='append')
           write(unit_map,'(2i6,5es14.6)') &
                ir,ii,om(ir,ii),val(ir,ii),cal(ir,ii)
           close(unit_map)
        endif
     enddo
     if (proc0) then
        open(unit=unit_map,file=trim(mapName),status='old',position='append')
        write(unit_map,*)
        close(unit_map)
     endif
  enddo

  if(determine_minima) then
    !Search for Local Minima in Dispersion Surface
    if (proc0.and.(nr.gt.1).and.(ni.gt.1)) then
      write(*,*)'finding minima'
      call find_minima(val,numroots,iroots,nroots_max)
      if (writeOut) &
           write(*,'(i2,a)')nroots_max,'  possible local minimum found'
      do iw=1,nroots_max
         wroots(iw)=om(iroots(1,iw),iroots(2,iw))
         if (writeOut) then
            write(*,'(a,i4,a,i4)')'ir = ',iroots(1,iw),'    ii = ',iroots(2,iw)
            write(*,'(4es14.4)') wroots(iw) ,cal(iroots(1,iw),iroots(2,iw))
         endif
      enddo

      nroots=min(nroots,nroots_max)

   endif

    if (proc0) write(*,*)'testing!'

    call mpi_bcast(nroots, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call mpi_bcast(wroots(:), size(wroots(:)), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)

    call refine_guess

    call mpi_bcast(nroots, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call mpi_bcast(wroots(:), size(wroots(:)), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)

  endif

end subroutine map_search
!-=-=-=-=-=-=
!
!-=-=-=-=-=-=
subroutine refine_guess
  use alps_var, only : wroots, nroots, writeOut, proc0
  use alps_var, only : ierror,runname, D_prec
  use ALPS_var, only : use_secant
  use alps_io,  only : get_unused_unit
  use mpi
  implicit none

  double complex :: omega    !Complex Frequency
  double complex :: om1, om2 !Bounding frequencies for secant search
  integer :: iflag                           !Flag for Root search
  character(100) :: mapName                  !Output file names
  double complex :: tmpDisp  !holding variable for dispersion relation solution
  integer :: iw   !loop index
  integer :: unit_refine

  !Refine

  if (proc0) then
     if (writeOut) write(*,'(a)')' Refining Roots:'
     write(mapName,'(3a)') 'solution/',trim(runname),'.roots'
     call get_unused_unit(unit_refine)
     open(unit=unit_refine,file=trim(mapName),status='replace')
  endif


  do iw=1,nroots
     if (proc0.and.writeOut) write(*,'(a,i0)')'Root ',iw
     call mpi_barrier(mpi_comm_world,ierror)
     omega=wroots(iw)

     om1=omega*(1.d0-D_prec)
     om2=omega*(1.d0+D_prec)


     if (use_secant) then
           call secant(omega)
        else
           omega=rtsec(disp,om1,om2,iflag)
        endif

     wroots(iw)=omega

     tmpDisp=disp(wroots(iw))
     if (proc0.and.(abs(tmpDisp).NE.0.d0)) then
        write(unit_refine,'(i4,5es14.4)') iw,wroots(iw),log10(abs(tmpDisp)),tmpDisp
        write(*,'(i4,5es14.4)') iw,wroots(iw),log10(abs(tmpDisp)),tmpDisp
 !       if (writeOut) write(*,'(a,2es14.4,a,2es14.4)')'D(',wroots(iw),')= ',tmpDisp
     endif

  enddo
  if (proc0) close(unit_refine)



end subroutine refine_guess

!-=-=-=-=-=-=
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!     Based upon routine from Greg Howes, 2005
!     which was adapted from f77 routine by Eliot Quataert
   double complex function rtsec(func,x1,x2,iflag)
     use alps_var, only : proc0,writeOut,D_threshold,numiter
     implicit none
     double complex :: func, x1, xl, x2
     double complex :: fl, f, swap, dx
     integer :: iflag,j

     fl=func(x1)
     f=func(x2)

     if(abs(fl).lt.abs(f))then
        rtsec=x1
        xl=x2
        swap=fl
        fl=f
        f=swap
     else
        xl=x1
        rtsec=x2
     endif
     do  j=1,numiter-1
        iflag = j
        if (abs(f-fl) .GT. 1.d-80) then
           dx=(xl-rtsec)*f/(f-fl)
				else
           dx = (x2-x1)/25.d0
        end if
        xl=rtsec
        fl=f

        rtsec=rtsec+dx/2.d0

        f=func(rtsec)
        if (proc0) write (*,*) j,f,rtsec
        if((abs(dx).LT.D_threshold).OR.(abs(f).EQ.0.d0)) then
	     	if (proc0.AND.writeOut) write(*,'(a,i4)') 'Converged after iteration ',j
	        return
        endif
     enddo

     if (proc0.AND.writeOut) write(*,'(a,i4,a)') 'Maximum iteration ',j,' reached.'

     return

   end function rtsec
!------------------------------------------------------------------------------


!
!-=-=-=-=-=-=
!                          Based upon routine by
!                          Greg Howes, 2006
!-=-=-=-=-=-=
   subroutine find_minima(val,numroots,iroots,nroots)
     use ALPS_var, only : ni,nr
     implicit none
     !Passed
     double precision, dimension(:,:), pointer :: val       !Value of Dispersion relation
     integer :: numroots                        !Number of roots
     integer, dimension(1:2,1:numroots) :: iroots     !Indices of roots
     integer, intent(out) :: nroots             !Number of roots found
     !Local
     integer :: ir,ii                           !Counters

     !Find local minima in map
     iroots=0
     nroots=0
     do ii=ni,1,-1
        do ir=1,nr
           !ir = 0
           if (ir ==1) then
              if (val(ir,ii) .lt. val(ir+1,ii)) then
                 if (ii==1) then
                    if (val(ir,ii) .lt. val(ir,ii+1)) then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 elseif (ii==ni) then
                    if (val(ir,ii) .lt. val(ir,ii-1)) then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 else
                    if (val(ir,ii) .lt. val(ir,ii-1) .and.  &
                         val(ir,ii) .lt. val(ir,ii+1))then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 endif
              endif
           elseif (ir == nr) then
              if (val(ir,ii) .lt. val(ir-1,ii)) then
                 if (ii==1) then
                    if (val(ir,ii) .lt. val(ir,ii+1)) then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 elseif (ii==ni) then
                    if (val(ir,ii) .lt. val(ir,ii-1)) then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 else
                    if (val(ir,ii) .lt. val(ir,ii-1) .and.  &
                         val(ir,ii) .lt. val(ir,ii+1))then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 endif
              endif
           else
              if (val(ir,ii) .lt. val(ir-1,ii) .and.  &
                   val(ir,ii) .lt. val(ir+1,ii))then
                 if (ii==1) then
                    if (val(ir,ii) .lt. val(ir,ii+1)) then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 elseif (ii==ni) then
                    if (val(ir,ii) .lt. val(ir,ii-1)) then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 else
                    if (val(ir,ii) .lt. val(ir,ii-1) .and.  &
                         val(ir,ii) .lt. val(ir,ii+1))then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 endif
              endif
           endif
        enddo
     enddo

     write(*,*)nroots

   end subroutine find_minima

!-=-=-=-=-=-=
!Bessel Functions
!-=-=-=-=-=-=
! determine nmax:
subroutine determine_nmax()
use ALPS_var, only : pp, kperp, qs, Bessel_zero, nmax, ierror
use ALPS_var, only : proc0, nperp, nspec, writeOut, nproc, usebM
use ALPS_fns_rel, only : BESSJ
use mpi
implicit none

integer :: is 				! Species index
integer :: nn				! Bessel n
double precision :: bessel, besselprev, besselmax
double precision :: z		!Bessel Argument
integer :: iperp, ipar 		!p_perp, p_par index
integer :: max_procs

logical :: maximum			! Did we find a maximum?
logical :: modified_nmax

ipar = 1

max_procs=nspec
do is = 1, nspec

  if (usebM(is)) then

         nmax(is)=1

      else


    nn = 0
    besselmax = 10.d0

    do while (besselmax.GT.Bessel_zero)
       nn = nn + 1

       iperp = 0
       maximum = .FALSE.
       bessel = 0.d0
       if (nn.GT.1000) write (*,*) 'Bessel-function n is greater than 1000.'

       do while ((iperp .LT. nperp).AND.(.NOT.maximum))
          iperp=iperp+1
          z = kperp * pp(is, iperp, ipar, 1) / abs(qs(is))

          besselprev = bessel
          bessel = BESSJ(nn,z)
          if (bessel .LT. besselprev) maximum = .TRUE.
       enddo
       besselmax = bessel
    enddo

    nmax(is) = nn

    endif


    if (writeOut .and. proc0) &
         write(*,'(a,i0,a,i0,a)') 'Required nMax for species ',is,' : ',nmax(is)

    max_procs=max_procs+nmax(is)
enddo

! If we have too many processes, then allow for additional ns:
! The condition is that max_procs has to be greater or equal (nproc-1):
is=1
modified_nmax=.FALSE.
do while (max_procs.LT.(nproc-1))
  if (usebM(is)) then

  else
    modified_nmax=.TRUE.
	  nmax(is)=nmax(is)+1
  endif
	  is=is+1
	  max_procs=max_procs+1
	if (is.GT.nspec) is = 1
enddo

if (modified_nmax.AND.proc0.AND.writeOut) then
	write (*,'(a)') "More processes than required. nMax is adjusted."
	do is=1,nspec
	 	write(*,'(a,i0,a,i0,a)') ' Adjusted nMax for Species ',is,' : ',nmax(is)
	enddo
endif

call mpi_barrier(mpi_comm_world,ierror)

end subroutine 	determine_nmax

!-=-=-=-=-=-=
!split the processes on different n for parallelization
!-=-=-=-=-=-=
subroutine split_processes()
use alps_var, only : nproc, iproc, nmax, nlim
use alps_var, only : nspec, sproc, writeOut, proc0
implicit none
integer :: is   		! Species index
integer :: max_procs ! maximum number of processes (all ns for all species)
integer :: ideal_ns_per_proc
integer :: proc_per_spec(nspec) ! number of processes for each species
integer :: ideal_splitting(nspec) ! ideal splitting for all processes associated with spec 1
integer :: splitting_rest(nspec)	! Rest of ideal splitting
integer :: largest_rest			! Largest rest of splitting
integer :: largest_spec			! This species has the largest rest
integer :: used_procs			! Number of used procs
integer :: proc_count,prev_proc_count		! Process Counter
integer :: local_iproc			! Local process number (in species field)
integer :: rest_sum

max_procs = nspec	! to include the Bessel functions with n = 0 for all species
do is = 1,nspec
	max_procs = max_procs + nmax(is)
enddo

ideal_ns_per_proc = ceiling((1.*max_procs)/(1.*nproc-1.))

used_procs = 0
rest_sum = 0
! how many processes does each species get?
do is = 1, nspec
	if ((nmax(is)+1).LE.ideal_ns_per_proc) then
		proc_per_spec(is) = 1
	else
  	  proc_per_spec(is) = (nmax(is)+1)/ideal_ns_per_proc ! COULD LEAVE A REST
  	endif
	  ideal_splitting(is) = (nmax(is)+1)/proc_per_spec(is)  ! is the ideal splitting of species is
	  splitting_rest(is) = modulo((nmax(is)+1),proc_per_spec(is))
	  ! Every process for species is should get as close as possible to this number!
	  ! This is a little bit better than just using ideal_ns_per_proc
	  ! The last one will get the rest.
  	  !  Do we have processes left?
	used_procs = used_procs + proc_per_spec(is)
  rest_sum = rest_sum + splitting_rest(is)
enddo


if (proc0.AND.writeOut) then
	if(modulo((max_procs + rest_sum + 1),2).EQ.0) then
	 write(*,'(a,i0)')&
     'Ideal Number of Processors: ',max_procs + rest_sum + 1
     else
		write(*,'(a,i0)')&
     'Ideal Number of Processors: ',max_procs + rest_sum
     endif
     write(*,'(a,i0)') '-=-=-=-=-=-=-=-=-=-'
endif

! Find out which species has the largest rest of n's:
largest_spec=1
largest_rest=0
do is = 1,nspec
	if (splitting_rest(is).GT.largest_rest) then
		largest_spec = is
		largest_rest = splitting_rest(is)
	endif
enddo

! the rest of our processes will go to species largest_spec
 proc_per_spec(largest_spec) = proc_per_spec(largest_spec) + ((nproc-1)-used_procs)
 ideal_splitting(largest_spec) = nint((1.*nmax(largest_spec)+1.)/(1.*proc_per_spec(largest_spec)))


proc_count = 0
prev_proc_count = 0
do is = 1,nspec
	proc_count = proc_count + proc_per_spec(is)

		if ((iproc.LE.proc_count).AND.(iproc.GT.prev_proc_count)) then
			sproc = is
			local_iproc = (iproc - prev_proc_count)

			nlim(1) = (local_iproc - 1) * (ideal_splitting(is))
      nlim(2) = nlim(1) + ideal_splitting(is)-1

      if ((local_iproc.EQ.proc_per_spec(is)).AND.(nlim(1).LE.nmax(is))) nlim(2)=nmax(is)

		endif
	prev_proc_count = proc_count
enddo

  if (writeOut) &
       write(*,'(a,i4,a,i4,a,i4,a,2i4,a)') &
       'Processor ',iproc,' of ',nproc,' ready. Species ',sproc,': n in [', nlim(1:2),']'

end subroutine split_processes





! determine bessel_array (this process only works on one species
subroutine determine_bessel_array()
  use ALPS_var, only : pp, kperp, qs, sproc,bessel_array, nperp, nlim, iproc
  use ALPS_io, only : get_unused_unit
  use ALPS_fns_rel, only : BESSJ
  implicit none
  integer :: nn				! Bessel n
  !double precision :: BESSJ	! Bessel function
  double precision :: z		!Bessel Argument
  character (10) :: procwrite
  integer :: iperp, ipar !p_perp, p_par index
  !integer :: unit_bessel


  ipar = 1

  write(procwrite,'(i0)') iproc

  ! Allocate bessel_array:
  if (allocated(bessel_array)) deallocate(bessel_array)
  allocate(bessel_array(nlim(1)-1:nlim(2)+1,0:nperp)); bessel_array = 0.d0
  !unit_bessel=1001+iproc

  !open(unit = unit_bessel,file = 'solution/besselArray.'//trim(procwrite)//'.out', status = 'replace')

  !running into an issue for low n, larger kperp, with electrons

  ! Fill this array with values:
  do nn = nlim(1)-1, nlim(2)+1
     do iperp = 0, nperp
        z = kperp * pp(sproc, iperp, ipar, 1) / qs(sproc)
        bessel_array(nn,iperp) = BESSJ(nn,z)
        if (nn.EQ.-1) bessel_array(nn,iperp)=-BESSJ(1,z)
        !write(unit_bessel,'(i4,i4,2es17.9)')&
        !     nn, iperp, z, bessel_array(nn,iperp)
     enddo
     !write(unit_bessel,*); write(unit_bessel,*)
  enddo
  !close(unit_bessel)

end subroutine determine_bessel_array





end module alps_fns
