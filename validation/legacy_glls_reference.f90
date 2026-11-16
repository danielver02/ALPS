! Historical validation only. Never linked into the production ALPS target.
! Arithmetic copied from the pre-B01 GLLS implementation; full-weight c0.
module legacy_glls_reference
 implicit none
 private
 public :: legacy_determine_param_fit,reference_coefficients,legacy_forward_evaluate,historical_log_max
 double precision,allocatable :: historical_log_max(:)
 public :: historical_coefficients,initialize_historical_nr,ensure_historical_storage,select_nr_test_evaluation
 public :: distribute_historical_coefficients,history_enabled,history_transfer_enabled
 logical :: history_enabled=.true.,history_transfer_enabled=.true.
 double precision, allocatable :: historical_coefficients(:,:,:)
 double precision, allocatable :: polynomials(:,:,:),reference_coefficients(:,:,:)
contains
subroutine ensure_historical_cutoff
 use alps_var,only:nspec
 if(allocated(historical_log_max))then
 if(size(historical_log_max)/=nspec)deallocate(historical_log_max)
 endif
 if(.not.allocated(historical_log_max))then
 allocate(historical_log_max(nspec));historical_log_max=18.d0
 endif
end subroutine
subroutine ensure_historical_storage
 use alps_var,only:nspec,nperp,poly_order
 call ensure_historical_cutoff
 if(.not.allocated(historical_coefficients))then
 allocate(historical_coefficients(nspec,0:nperp,0:maxval(poly_order)))
 historical_coefficients=0.d0
 endif
end subroutine
subroutine select_nr_test_evaluation(enabled)
 use alps_analyt,only:set_nr_cheb_evaluation
 logical,intent(in)::enabled
 history_enabled=.not.enabled
 call set_nr_cheb_evaluation(enabled)
end subroutine
subroutine distribute_historical_coefficients
 use alps_var,only:proc0,nspec,nperp,ACmethod,relativistic,useBM,ierror
 use mpi
 integer::degree
 character(32)::backend
 if(.not.history_transfer_enabled)return
 call get_environment_variable("ALPS_NR_SVD_BACKEND",backend)
 if(trim(backend)=="svd_native".or.trim(backend)=="svd_compat")return
 if(.not.any(ACmethod==2.and..not.relativistic.and..not.useBM))return
 degree=-1
 if(proc0.and.allocated(historical_coefficients))degree=ubound(historical_coefficients,3)
 call mpi_bcast(degree,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
 if(degree<0)return
 if(.not.proc0)then
 if(allocated(historical_coefficients))deallocate(historical_coefficients)
 allocate(historical_coefficients(nspec,0:nperp,0:degree))
 endif
 call mpi_bcast(historical_coefficients,size(historical_coefficients),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
 call ensure_historical_cutoff
 call mpi_bcast(historical_log_max,size(historical_log_max),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
end subroutine
subroutine initialize_historical_nr(status)
 use alps_var,only:nspec,nperp,npar,pp,poly_order,logfit
 use alps_analyt,only:cheb_model,cheb_from_coefficients,nr_cheb_import,nr_cheb_species_eligible, &
 CHEB_SUCCESS,CHEB_INVALID_MODEL,CHEB_INVALID_DEGREE,CHEB_DIRECT,CHEB_LOG10,NR_FIT_IMPORTED,NR_ROW_FITTED, &
 nr_cheb_evaluation_active,set_nr_cheb_evaluation,release_nr_cheb
 integer,intent(out)::status
 integer,allocatable::metadata(:,:)
 double precision,allocatable::mapping(:,:),coefficients(:,:)
 type(cheb_model)::model
 integer::s,r,i,degree,rep
 logical::was_enabled
 was_enabled=.false.
 do s=1,nspec
 was_enabled=was_enabled.or.nr_cheb_evaluation_active(s)
 enddo
 call release_nr_cheb
 status=CHEB_INVALID_MODEL
 if(.not.allocated(historical_coefficients))return
 allocate(metadata(17,nspec*(nperp+1)),mapping(6,nspec*(nperp+1)))
 allocate(coefficients(size(historical_coefficients,3),nspec*(nperp+1)))
 metadata=0;mapping=0.d0;coefficients=0.d0;i=0
 do r=0,nperp
 do s=1,nspec
 i=i+1;metadata(3,i)=CHEB_INVALID_MODEL;metadata(5,i)=-1
 if(.not.nr_cheb_species_eligible(s))cycle
 degree=poly_order(s);status=CHEB_INVALID_DEGREE
 if(degree<0.or.degree>ubound(historical_coefficients,3))return
 rep=CHEB_DIRECT;if(logfit(s))rep=CHEB_LOG10
 call cheb_from_coefficients(historical_coefficients(s,r,0:degree),degree, &
 .5d0*pp(s,r,npar,2)+.5d0*pp(s,r,0,2),.5d0*pp(s,r,npar,2)-.5d0*pp(s,r,0,2),rep,model,status)
 if(status/=CHEB_SUCCESS)return
 metadata(:,i)=[NR_FIT_IMPORTED,NR_ROW_FITTED,status,degree,degree,rep,model%rank,model%initial_rank, &
 model%original_sample_count,model%retained_sample_count,model%sample_count,0,0,0,0,0,0]
 mapping(:,i)=[model%center,model%scale,pp(s,r,0,2),pp(s,r,npar,2),model%retained_min,model%retained_max]
 coefficients(1:degree+1,i)=model%coefficients
 enddo
 enddo
 call nr_cheb_import(nspec,nperp,metadata,mapping,coefficients,status)
 if(status==CHEB_SUCCESS)call set_nr_cheb_evaluation(was_enabled)
end subroutine
subroutine legacy_determine_param_fit
 use alps_var, only: nspec,nperp,npar,poly_order,ACmethod,relativistic,useBM
 use alps_analyt, only: determine_param_fit
 integer :: species
 call ensure_historical_cutoff
 call determine_param_fit
 if(allocated(polynomials))deallocate(polynomials,reference_coefficients)
 allocate(polynomials(nspec,0:npar,0:maxval(poly_order)))
 allocate(reference_coefficients(nspec,0:nperp,0:maxval(poly_order)))
 polynomials=0.d0;reference_coefficients=0.d0
 do species=1,nspec
 if(ACmethod(species)/=2.or.relativistic(species).or.useBM(species))cycle
 call legacy_set_basis(species)
 call legacy_fit_species(species)
 enddo
 ! Validation-owned coefficients; no production buffer or fit backend.
 if(allocated(historical_coefficients))deallocate(historical_coefficients)
 historical_coefficients=reference_coefficients
end subroutine
subroutine legacy_fit_species(is)
  !! This subroutine evaluates the General Linear Least Squares fit
  !! to the distribution function for component 'is' using the selected
  !! polynomial basis functions.
  use alps_var, only : nperp, npar, f0, poly_order
  use alps_var, only : logfit
  implicit none

  integer :: is
  !! Index of species.

  integer :: iperp
  !! Perpendicular index.

  integer :: ipar
  !! parallel index

  double precision, dimension(:), allocatable :: f0_fit

  double precision :: min_val

  write(*,'(a,i2)')'Determining GLLS Coefficients for Component: ',is
  
  allocate(f0_fit(0:npar));f0_fit=0.d0
  do iperp=0,nperp
     if (logfit(is)) then

        if (minval(f0(is,iperp,:)) .gt. 0.d0) then
           f0_fit(:)=log10(f0(is,iperp,:))
        else
           ! Get min nonzero value from f0(is,iperp,:) before log10 conversion:
           if (any(f0(is,iperp,:) > 0.d0)) then
              min_val = minval(f0(is,iperp,:), mask=(f0(is,iperp,:) > 0.d0))
           else
              ! Fallback value if alrel entries are zero; adjust as needed.
              min_val = 1.0d-40
           endif
           
           ! Copy the slice and replace zeros with 0.01*min_val:
           f0_fit(:) = f0(is,iperp,:)
           do ipar = 0, npar
              if (f0_fit(ipar) == 0.d0) then
                 f0_fit(ipar) = 0.01d0 * min_val
              endif
           enddo
           
           do ipar = 0, npar
              if (f0_fit(ipar) <= 0.0d0) then
                 print*, 'Non-positive value at index', ipar, f0_fit(ipar)
              endif
           enddo
           f0_fit(:)=log10(f0_fit(:))
        !else
        !   f0_fit(:)=log10(f0(is,iperp,:))
        endif
     else
        f0_fit(:)=f0(is,iperp,:)
     endif
     call legacy_least_squares_fit(polynomials(is,:,:),f0_fit,reference_coefficients(is,iperp,:),poly_order(is))
  enddo
  deallocate(f0_fit)
  
end subroutine legacy_fit_species

subroutine legacy_least_squares_fit(AA,BB,coeffs,npoly)
  !!Solves General Linear Least Squares Normal Equation
  use alps_var, only : npar
  implicit none

  integer :: npoly
  !! Order for polyhedral representation
  
  double precision, intent(in) :: BB(0:npar)
  !! f0(is,iperp,:)

    double precision, intent(out) :: coeffs(0:npoly)
    !! fit coefficients

    double precision, intent(in) :: AA(0:npar, 0:npoly)
    !! Polynomial Basis
   
    integer, dimension(0:npoly) :: ipiv
    !!Pivot Indices

    integer :: info
    !!LAPACK error diagnostics

    double precision :: alpha(0:npoly,0:npoly)
    !! A^T.A component of General Linear Least Square calculation
    !! e.g. Eqn 15.4.7 from Numerical Recipies
    
    double precision :: beta(0:npoly)
    !! A^T.b component of General Linear Least Square calculation
    !! e.g. Eqn 15.4.7 from Numerical Recipies

    alpha=0.d0
    beta=0.d0
    
    ! Calculate (A^T * A) for the normal equation solution
    call dgemm('T', 'N', npoly+1, npoly+1, npar+1, 1.0d0, AA(0:npar, 0:npoly), &
         npar+1, AA(0:npar, 0:npoly), npar+1, 0.0d0, alpha(0:npoly, 0:npoly), npoly+1)

    ! Calculate (A^T * B) for the normal equation solution
    call dgemv('T', npar+1, npoly+1, 1.0d0, AA(0:npar, 0:npoly), &
         npar+1, BB(0:npar), 1, 0.0d0, beta(0:npoly), 1)
    
    ! Solve for the coefficients using the LAPACK routine
    call dgesv(npoly+1, 1, alpha(0:npoly, 0:npoly), npoly+1, &
         ipiv(0:npoly), beta(0:npoly), npoly+1, info)

    coeffs(0:npoly)=beta(0:npoly)

    if (info /= 0) then
        print *, "Error in dgesv:", info
        stop
    end if
  
end subroutine legacy_least_squares_fit

subroutine legacy_set_basis(is)
  !! This subroutine evaluates the General Linear Least Squares fit
  !! to the distribution function for component 'is' using the selected
  !! polynomial basis functions.
  use alps_var, only : poly_kind, poly_order
  use alps_var, only : writeOut, npar
  use alps_io, only : alps_error
  implicit none

  integer :: is
  !! Index of species.

  integer :: ipar
  !! Index of parallel momentum.

  integer :: n
  !! Polynomial Order

  double precision :: yy
  !! Argument of Polynomial

  select case (poly_kind(is))
  case (1) !Chebyshev Polynomial Basis
     if (writeOut) & 
          write(*,'(a,i2)')'Constructing Chebyshev Basis for Component ',is
     
     polynomials(is,:,0) = 1.0
     do ipar = 0, npar
        yy=-1.d0+ipar*(2.d0/npar)
        polynomials(is,ipar,1) = yy
        do n = 2, poly_order(is)
           polynomials(is,ipar,n) = &
                2.0 * yy * polynomials(is,ipar,n-1) - polynomials(is,ipar,n-2)
        end do
     enddo
  case default
     call alps_error(10)
  end select
  
end subroutine legacy_set_basis

double complex function legacy_forward_evaluate(is,iperp,ppar_val,n_poly,fit_coeffs,dfpar)
  !! This function evaluates the orthogonal polynomical fit to f0
  !! at real pperp_val and complex ppar_val,
  !! with the one-dimensional polynomical coefficient array fit_coeffs is fed into the
  !! function.
  !! For the evaluation in ALPS, use [[eval_fit(function)]].
  use alps_var, only : pp, poly_kind, npar, logfit
 implicit none

 integer :: n_poly
 !! Maximum Polynomial Order
 
 integer, intent(in) :: is
 !! Index of species for which [[eval_fit(function)]] is executed.

 double precision, intent(in) :: fit_coeffs(0:n_poly)
 !! Array of polynomial coefficients

 double complex :: poly_basis(0:n_poly)
 !! Array of polynomial coefficients

  integer :: n
 !! Polynomial Order Index

 integer, intent(in) :: iperp
 !! Index of perpendicular momentum at which [[legacy_forward_evaluate(function)]] is executed.

 double complex :: ppar_val
 !! Complex parallel momentum.

 double complex :: ppar_val_tmp
 !! Complex parallel momentum.

 double precision :: norm_1, norm_2
 !! Range of p_parallel for rescaling polynomials

 double complex, intent(out), optional :: dfpar
 !! Derivative of polynomial representations.

 double complex :: poly_basis_deriv(0:n_poly)
 !! Derivatives of the polynomial basis.
 
 double complex :: poly_deriv
 !! Total Derivative

 call ensure_historical_cutoff
 legacy_forward_evaluate=cmplx(0.d0,0.d0,kind(1.d0))
 if (present(dfpar)) dfpar = (0.d0,0.d0)
 poly_deriv = (0.d0,0.d0)
 
 select case (poly_kind(is))
 case (1)
    !Chebyshev polynomials range from [-1,1]
    norm_1=5.d-1*(pp(is,iperp,npar,2)+pp(is,iperp,0,2))
    norm_2=5.d-1*(pp(is,iperp,npar,2)-pp(is,iperp,0,2))
    ppar_val_tmp=(ppar_val-norm_1)/norm_2

    if (abs(ppar_val_tmp).gt.1.d0) then !!kgk: should this be .gt. or .ge. ?
       legacy_forward_evaluate=0.d0
    else
	!Set n=0
    n=0
    poly_basis(n)=cmplx(1.d0,0.d0,kind(1.d0))
	if (present(dfpar)) &
		poly_basis_deriv(n) = (0.d0,0.d0)
    legacy_forward_evaluate=legacy_forward_evaluate+fit_coeffs(n)*poly_basis(n)
    
	!Set n=1
	n=1
    poly_basis(n)=ppar_val_tmp
	if (present(dfpar)) &
		poly_basis_deriv(n) = (1.d0,0.d0)
    legacy_forward_evaluate=legacy_forward_evaluate+fit_coeffs(n)*poly_basis(n)

	!Start Summing Derivatives
	if (present(dfpar)) &
		poly_deriv = fit_coeffs(n)

	!Start recursion
    do n=2,n_poly
       poly_basis(n) = cmplx(2.d0,0.d0,kind(1.d0)) * &
			ppar_val_tmp * poly_basis(n-1) - poly_basis(n-2)

       legacy_forward_evaluate=legacy_forward_evaluate+&
	   		fit_coeffs(n)*poly_basis(n)
		if (present(dfpar)) then
			poly_basis_deriv(n) = cmplx(2.d0,0.d0,kind(1.d0))*poly_basis(n-1) &
            + cmplx(2.d0,0.d0,kind(1.d0))*ppar_val_tmp*poly_basis_deriv(n-1) &
            - poly_basis_deriv(n-2)
	   		poly_deriv = poly_deriv &
              + fit_coeffs(n)*poly_basis_deriv(n)
		endif
    enddo
    if (logfit(is)) then

       !if (abs(legacy_forward_evaluate).gt.20.) then
       !KGK: how does this impact the 'noise'
       !if ((real(legacy_forward_evaluate).gt.1.).or.(aimag(legacy_forward_evaluate).gt.1.)&
       !     .or.(real(legacy_forward_evaluate).lt.-18.).or.(aimag(legacy_forward_evaluate).lt.-18.)) then

       !if ((real(legacy_forward_evaluate).gt.10.).or.(aimag(legacy_forward_evaluate).gt.10.)&
       !     .or.(real(legacy_forward_evaluate).lt.-22.).or.(aimag(legacy_forward_evaluate).lt.-22.)) then

       !KGK: Need to investigate the appropriate threshold here for when this functional
       !representation ceases to work.
       !if ((real(legacy_forward_evaluate).gt.17.).or.(aimag(legacy_forward_evaluate).gt.17.)&
         !if ((real(legacy_forward_evaluate).lt.-historical_log_max(is)).or.& !??
              !(aimag(legacy_forward_evaluate).lt.-historical_log_max(is))) then

       !if ((real(legacy_forward_evaluate).lt.-historical_log_max(is)).or.& !??
       !     (aimag(legacy_forward_evaluate).lt.-historical_log_max(is)).or.&
       !     (real(legacy_forward_evaluate).gt. historical_log_max(is)).or.& !??
       !     (aimag(legacy_forward_evaluate).gt.historical_log_max(is))) then

       !KGK: 2026-10-12
       if ((real(legacy_forward_evaluate).gt. historical_log_max(is))) then
       !if ((abs(legacy_forward_evaluate).lt.-100.d0).or.(abs(legacy_forward_evaluate).gt.100.d0)) then
       !if ((abs(legacy_forward_evaluate).lt.-50.d0).or.(abs(legacy_forward_evaluate).gt.50.d0)) then
       !if (abs(legacy_forward_evaluate).gt.historical_log_max(is)) then
          legacy_forward_evaluate=cmplx(0.d0,0.d0,kind(1.d0))
       else
          !write(*,*)iproc,is,iperp,ppar_val,ppar_val_tmp,legacy_forward_evaluate
          legacy_forward_evaluate=10.d0**(legacy_forward_evaluate)
          !write(*,*)iproc,is,iperp,ppar_val,ppar_val_tmp,legacy_forward_evaluate
       endif
    endif
	!Calculate the Actual Derivative
	if (present(dfpar)) then
   		if (logfit(is)) then
      		dfpar = log(10.d0)*legacy_forward_evaluate*poly_deriv/norm_2
   		else
      		dfpar = poly_deriv/norm_2
   		endif
	endif
 endif
end select

 return
end function legacy_forward_evaluate


end module

! External hooks resolve only in instrumented validation executables. No module
! dependency from the production numerical kernel back into validation code.
subroutine historical_distribution
 use legacy_glls_reference,only:distribute_historical_coefficients
 call distribute_historical_coefficients
end subroutine
logical function historical_active(species)
 use legacy_glls_reference,only:historical_coefficients,history_enabled
 use alps_analyt,only:nr_cheb_evaluation_active,nr_cheb_species_eligible
 integer,intent(in)::species
 historical_active=.false.
 if(.not.allocated(historical_coefficients).or..not.history_enabled)return
 historical_active=nr_cheb_species_eligible(species).and..not.nr_cheb_evaluation_active(species)
end function
double complex function historical_value(species,row,u,want_derivative,derivative)
 use legacy_glls_reference,only:historical_coefficients,legacy_forward_evaluate
 use alps_var,only:poly_order
 integer,intent(in)::species,row
 double complex,intent(in)::u
 logical,intent(in)::want_derivative
 double complex,intent(out)::derivative
 derivative=(0.d0,0.d0)
 if(want_derivative)then
 historical_value=legacy_forward_evaluate(species,row,u,poly_order(species), &
 historical_coefficients(species,row,0:poly_order(species)),derivative)
 else
 historical_value=legacy_forward_evaluate(species,row,u,poly_order(species), &
 historical_coefficients(species,row,0:poly_order(species)))
 endif
end function
