!=============================================================================!
!=============================================================================!
!!*ALPS                                                                  *!!
!!                                                                           !!
!!Kristopher Klein, Daniel Verscharen                                        !!
!!
!!Space Science Center, University of New Hampshire                          !!
!!                                                                           !!
!!*MPI COMMUNICATIONS MODULE                                                 *!!
!=============================================================================!
!=============================================================================!
!

module alps_com
  implicit none
  
  public :: pass_instructions
  public :: pass_distribution  

contains

!-=-=-=-=-=-=-=-=-=
!Routine for passing information
!-=-=-=-=-=-=-=-=-=
subroutine pass_instructions
  use alps_var, only : proc0, ierror, nroots, n_fits, param_fit, fit_type, perp_correction
  use alps_var, only : writeOut, nperp, npar, nmax, nlim, nspec, numroots, ngamma,npparbar
  use alps_var, only : ns, qs, ms, wroots, kperp, kpar, bessel_zero, D_prec
  use alps_var, only : wave, chi0, pp, df0, vA, pi, numiter, D_threshold, D_gap
  use alps_var, only : use_map
  use alps_var, only : ni, nr, omi, omf, gami, gamf, loggridg, loggridw
  use alps_var, only : determine_minima, n_resonance_interval, positions_principal, Tlim
  use alps_var, only : n_scan, scan, scan_option, use_secant, relativistic
  implicit none
  include 'mpif.h' !Include MPI library variables
  !Local
  integer :: is


  !Broadcast Global Variables needed for code execution
  call mpi_bcast(nperp,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(npar,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(ngamma,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(npparbar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(nspec,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(nroots,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(numiter,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(n_resonance_interval, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(positions_principal, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(writeOut, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(kperp,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(kpar,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(vA,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(bessel_zero,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(D_threshold,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(D_prec,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(D_gap,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(Tlim,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(use_map, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(n_scan, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(scan_option, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(use_secant,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  
  if (use_map) then
     call mpi_bcast(omi,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
     call mpi_bcast(omf,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
     call mpi_bcast(gami,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
     call mpi_bcast(gamf,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

     call mpi_bcast(loggridw, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
     call mpi_bcast(loggridg, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)

     call mpi_bcast(ni,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
     call mpi_bcast(nr,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
     
     call mpi_bcast(determine_minima, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)

  endif

  pi = 4.d0*atan(1.d0)

  allocate(nmax(1:nspec)); nmax = 0
  nlim = 0

  !Point to Point Communication:
  if (proc0) then
     !+ Send instructions

     !wroots allocated for proc0 in ALPS_io
     !to allow for reading in of guesses for the dispersion solution

  else
     !+ Receive instructions
     allocate(ns(1:nspec)); ns = 0.d0
     allocate(qs(1:nspec)); qs = 0.d0
     allocate(ms(1:nspec)); ms = 0.d0
     allocate(n_fits(1:nspec))
     allocate(relativistic(1:nspec)); relativistic=.FALSE.

     !Allocate solution space for nroots dispersion solutions
     allocate(wroots(1:numroots));wroots=cmplx(0.d0,0.d0,kind(1.d0))

     if (n_scan.gt.0) allocate(scan(n_scan))
  endif
  call mpi_bcast(ns(:), size(ns(:)),&
       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(qs(:), size(qs(:)),&
       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(ms(:), size(ms(:)),&
       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(n_fits(:), size(n_fits(:)),&
       MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call mpi_bcast(relativistic(:), size(relativistic(:)),&
       MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)


  call mpi_bcast(wroots(:), size(wroots(:)),&
       MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)

  !There is a 'better' way to do this with derived data types and mpi, but this way works,
  !and only has to be called once.
  if (n_scan.gt.0) then
     do is=1,n_scan
        call mpi_bcast(scan(is)%range_i,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)        
        call mpi_bcast(scan(is)%range_f,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)        
        call mpi_bcast(scan(is)%diff,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)        
        call mpi_bcast(scan(is)%diff2,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)        
        call mpi_bcast(scan(is)%log_scan, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)        
        call mpi_bcast(scan(is)%eigen_s,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)        
        call mpi_bcast(scan(is)%heat_s,   1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)        
        call mpi_bcast(scan(is)%type_s,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
        call mpi_bcast(scan(is)%n_out,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
        call mpi_bcast(scan(is)%n_res,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
     enddo     
  endif

  !Allocation of Arrays:

  if (.true..neqv.proc0) then
     allocate(param_fit(1:nspec,0:(max(nperp,ngamma)),4,maxval(n_fits)))
     allocate(fit_type(1:nspec,maxval(n_fits))) 
     allocate(perp_correction(1:nspec,maxval(n_fits)))
  endif
  
  !Allocate Wave Equation Tensor
  !  Only needed on proc0
  if (proc0) then
     allocate(wave(1:3,1:3))
     !KGK: added chi0 as global character
     allocate(chi0(nspec,1:3,1:3))
     wave=cmplx(0.d0,0.d0,kind(1.d0))
  endif
  
  !Allocate derivative of f0
  allocate(df0(1:nspec,1:nperp-1,1:npar-1,1:2)); df0=0.d0
  
  !Allocate fit parameters
  
   
  !Allocate velocity grid
  allocate(pp(1:nspec,0:nperp,0:npar,1:2)); pp=0.d0
    
end subroutine pass_instructions
!-=-=-=-=-=-=-=-=-=
!Routine for passing distribution function parameters
!-=-=-=-=-=-=-=-=-=
subroutine pass_distribution
	use alps_var,    only : df0, pp, param_fit, fit_type, perp_correction,proc0, writeOut, ierror
	use alps_var,    only : df0_rel, gamma_rel, pparbar_rel, f0_rel
	use alps_var,    only : relativistic, nspec, ngamma, npparbar
	implicit none
	include 'mpif.h' !Include MPI library variables
	integer :: is_rel, is, nspec_rel
	logical :: any_relativistic


	!+ Broadcast derivative array
	if (writeOut.and.proc0)&
	   write(*,'(a)')'Broadcasting df0/dp...'



	any_relativistic = .FALSE.
	is_rel = 0
	do is = 1, nspec
		if (relativistic(is)) then
			any_relativistic=.TRUE.
			is_rel=is_rel+1
		endif
	enddo

	if (any_relativistic) then
		 nspec_rel=is_rel
		 if (.not.(proc0)) then
			 ! Allocate the relativistic fields:
			 allocate(gamma_rel(nspec_rel,0:ngamma,0:npparbar))
			 allocate(pparbar_rel(nspec_rel,0:ngamma,0:npparbar))
			 allocate(df0_rel(nspec_rel,0:ngamma,0:npparbar,2))
			 allocate(f0_rel(nspec_rel,0:ngamma,0:npparbar))
		endif
	
		call mpi_bcast(f0_rel(:,:,:), size(f0_rel(:,:,:)),&
		   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
		call mpi_bcast(df0_rel(:,:,:,:), size(df0_rel(:,:,:,:)),&
		   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
		call mpi_bcast(gamma_rel(:,:,:),  size(gamma_rel(:,:,:)),&
		   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
		call mpi_bcast(pparbar_rel(:,:,:),  size(pparbar_rel(:,:,:)),&
		   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)       
	endif

	  call mpi_bcast(df0(:,:,:,:), size(df0(:,:,:,:)),&
		   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
	  call mpi_bcast(pp(:,:,:,:),  size(pp(:,:,:,:)),&
		   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
	  call mpi_bcast(param_fit(:,:,:,:),  size(param_fit(:,:,:,:)),&
		   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
	  call mpi_bcast(fit_type(:,:),  size(fit_type(:,:)),&
		   MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
	  call mpi_bcast(perp_correction(:,:),  size(perp_correction(:,:)),&
		   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
	
	  if (writeOut.and.proc0)&
		   write(*,'(a)')' df0/dp received'
       
end subroutine pass_distribution

end module alps_com
