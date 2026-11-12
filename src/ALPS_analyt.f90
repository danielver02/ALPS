!Copyright (c) 2026, Kristopher G. Klein and Daniel Verscharen
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
!I                              Version 1.2                                    I
!I                                                                             I
!I  Kristopher Klein   (kgklein@arizona.edu)                                   I
!I  Daniel Verscharen  (daniel.verscharen@rub.de)                              I
!I                                                                             I
!===============================================================================

module alps_analyt
	!! This module contains functions and subroutines for the hybrid analytical continuation.
	implicit none

        public  :: eval_fit, determine_param_fit
        private :: set_polynomial_basis
        private :: fit_function, output_fit, determine_JT, LM_nonlinear_fit
        private :: determine_GLLS, factorial

  integer, parameter, private :: cheb_dp=kind(1.d0)
  ! Status values belong to the independent numerical kernel only.
  integer, parameter, public :: CHEB_SUCCESS=0, CHEB_REDUCED=1
  integer, parameter, public :: CHEB_NO_SAMPLES=2, CHEB_INVALID_COORDINATES=3
  integer, parameter, public :: CHEB_INVALID_NORMALIZATION=4, CHEB_INVALID_DEGREE=5
  integer, parameter, public :: CHEB_RANK_FAILURE=6, CHEB_LAPACK_FAILURE=7
  integer, parameter, public :: CHEB_NONFINITE_COEFFICIENTS=8, CHEB_INVALID_MODEL=9
  integer, parameter, public :: CHEB_EVALUATION_FAILURE=10, CHEB_INVALID_VALUES=11
  integer, parameter, public :: CHEB_INVALID_TOLERANCE=12, CHEB_ALLOCATION_FAILURE=13
  integer, parameter, public :: CHEB_INVALID_REPRESENTATION=14, CHEB_NO_POSITIVE_SAMPLES=15
  integer, parameter, public :: CHEB_NEGATIVE_SAMPLES=16, CHEB_INTERIOR_ZERO=17
  integer, parameter, public :: CHEB_LOG_OVERFLOW=18, CHEB_LOG_UNDERFLOW=19
  integer, parameter, public :: CHEB_LOG_VALUE_FAILURE=20, CHEB_LOG_DERIVATIVE_FAILURE=21
  integer, parameter, public :: CHEB_PHASE_FAILURE=22
  integer, parameter, public :: CHEB_NOT_REQUESTED=23
  integer, parameter, public :: CHEB_DIRECT=0, CHEB_LOG10=1
  ! Diagnostic availability is independent of fitting success and rank.
  integer, parameter, public :: CHEB_DIAG_UNAVAILABLE=0, CHEB_DIAG_VALID=1
  integer, parameter, public :: CHEB_DIAG_RANGE=2, CHEB_DIAG_FAILED=3

  type, public :: cheb_model
     !! First-kind series, with full-weight c(0); no species/global state.
     integer :: representation=CHEB_DIRECT
     real(cheb_dp), allocatable :: coefficients(:)       ! 0:effective_degree
     real(cheb_dp), allocatable :: derivative_coefficients(:) ! dP/dx, 0:max(0,degree-1)
     real(cheb_dp) :: center=0.d0, scale=1.d0
     integer :: requested_degree=-1, effective_degree=-1, sample_count=0, rank=0
     ! sample_count is the retained solver row size; diagnostics for LOG10
     ! measure residuals of log10(f), not of the reconstructed distribution.
     integer :: original_sample_count=0, retained_sample_count=0
     integer :: lower_tail_zeros=0, upper_tail_zeros=0
     real(cheb_dp) :: retained_min=0.d0, retained_max=0.d0
     integer :: fit_status=CHEB_INVALID_MODEL, lapack_info=0, fit_attempts=0, initial_rank=0
     logical :: sample_limited=.false., rank_reduced=.false., degenerate_interval=.false.
     real(cheb_dp) :: svd_tolerance=0.d0, singular_value_min=0.d0, singular_value_max=0.d0
     real(cheb_dp) :: residual_norm=0.d0, relative_residual=0.d0, rms_residual=0.d0
     ! Read a diagnostic only when its status is VALID. Fitted unavailable
     ! diagnostics are quiet NaNs; an unfitted/default model is UNAVAILABLE.
     integer :: residual_norm_status=CHEB_DIAG_UNAVAILABLE
     integer :: relative_residual_status=CHEB_DIAG_UNAVAILABLE, rms_residual_status=CHEB_DIAG_UNAVAILABLE
  end type cheb_model

  ! Row availability is distinct from the kernel fit status. No unavailable
  ! row is represented by fabricated coefficients or an implicit zero model.
  integer, parameter, public :: REL_CHEB_NOT_BUILT=-1, REL_CHEB_NOT_SELECTED=0
  integer, parameter, public :: REL_CHEB_FITTED=1, REL_CHEB_REDUCED=2
  integer, parameter, public :: REL_CHEB_NO_GEOMETRY=3, REL_CHEB_ZERO=4
  integer, parameter, public :: REL_CHEB_INVALID_SAMPLES=5, REL_CHEB_FIT_FAILURE=6
  type, public :: rel_cheb_row_info
     integer :: state=REL_CHEB_NOT_BUILT, kernel_status=CHEB_INVALID_MODEL
     integer :: species=0, packed_species=0, gamma_index=-1
     integer :: input_samples=0, geometric_samples=0, masked_samples=0, positive_samples=0
     real(cheb_dp) :: gamma=0.d0
  end type rel_cheb_row_info
  type(cheb_model), allocatable, private :: rel_cheb(:,:)
  type(rel_cheb_row_info), allocatable, private :: rel_cheb_info(:,:)
  integer, private :: rel_cheb_fit_calls=0
  integer, private :: rel_cheb_build_calls=0, rel_cheb_setup_status=CHEB_SUCCESS
  integer(kind=selected_int_kind(18)), private :: rel_cheb_evaluations=0
  public :: build_rel_cheb_rows, get_rel_cheb_row, release_rel_cheb_rows, rel_cheb_collection_info
  public :: initialize_rel_cheb, eval_rel_cheb, rel_cheb_initialization_status
  ! Narrow collection transfer helpers, used only by ALPS_com initialization.
  public :: rel_cheb_export, rel_cheb_import

  public :: cheb_fit, cheb_evaluate, cheb_evaluate_polynomial
  private :: cheb_fit_direct, cheb_polynomial_real, cheb_polynomial_complex
  private :: cheb_reconstruct, cheb_exp_component
  private :: cheb_evaluate_real, cheb_evaluate_complex, cheb_svd, cheb_model_valid
  private :: cheb_clenshaw, cheb_add, cheb_multiply, cheb_divide, cheb_complex_add, cheb_complex_multiply
  private :: cheb_residual_diagnostics, cheb_diagnostic_product
  private :: dgelsd, ilaenv

  interface cheb_evaluate
     module procedure cheb_evaluate_real, cheb_evaluate_complex
  end interface

  interface cheb_evaluate_polynomial
     module procedure cheb_polynomial_real, cheb_polynomial_complex
  end interface

  interface
     subroutine dgelsd(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,iwork,info)
       import cheb_dp
       integer, intent(in) :: m,n,nrhs,lda,ldb,lwork
       real(cheb_dp), intent(inout) :: a(lda,*),b(ldb,*)
       real(cheb_dp), intent(out) :: s(*),work(*)
       real(cheb_dp), intent(in) :: rcond
       integer, intent(out) :: rank,iwork(*),info
     end subroutine
     integer function ilaenv(ispec,name,opts,n1,n2,n3,n4)
       integer, intent(in) :: ispec,n1,n2,n3,n4
       character(*), intent(in) :: name,opts
     end function
  end interface

contains


subroutine initialize_rel_cheb(status)
  !! Rank-0 hook after derivative_f0, common to fresh and cached VDFs.
  !! Normal input still rejects relativistic method 2. Tests explicitly set
  !! internal selections before invoking this same hook; no guard is bypassed.
  use alps_var, only : proc0,relativistic,ACmethod,usebM,poly_order,gamma_rel,pparbar_rel,f0_rel
  implicit none
  integer, intent(out) :: status
  integer, allocatable :: species(:)
  logical, allocatable :: selected(:)
  integer :: is,s,ns,istat
  call release_rel_cheb_rows
  status=CHEB_SUCCESS
  setup: block
     if(.not.proc0) exit setup
     status=CHEB_INVALID_VALUES
     if(.not.allocated(relativistic).or..not.allocated(ACmethod).or..not.allocated(usebM)) exit setup
     ns=size(relativistic)
     if(size(ACmethod)/=ns.or.size(usebM)/=ns) exit setup
     status=CHEB_SUCCESS
     if(.not.any(relativistic.and.ACmethod==2.and..not.usebM)) exit setup
     status=CHEB_INVALID_VALUES
     if(.not.allocated(poly_order)) exit setup
     if(size(poly_order)/=ns) exit setup
     if(.not.allocated(gamma_rel).or..not.allocated(pparbar_rel).or..not.allocated(f0_rel)) exit setup
     allocate(species(count(relativistic)),selected(count(relativistic)),stat=istat)
     if(istat/=0) then
        status=CHEB_ALLOCATION_FAILURE; exit setup
     endif
     s=0
     do is=1,ns
        if(.not.relativistic(is)) cycle
        s=s+1; species(s)=is; selected(s)=ACmethod(is)==2.and..not.usebM(is)
     enddo
     call build_rel_cheb_rows(gamma_rel,pparbar_rel,f0_rel,species,poly_order(species),status,selected)
  end block setup
  rel_cheb_setup_status=status
end subroutine initialize_rel_cheb


integer function rel_cheb_initialization_status() result(status)
  implicit none
  status=rel_cheb_setup_status
end function rel_cheb_initialization_status


subroutine eval_rel_cheb(packed_species,igamma,u,value,status,derivative)
  !! Hot path: direct reference to stored coefficients, no model copy, fitting,
  !! allocation, transformed-table access, or communication. Diagnostics use
  !! get_rel_cheb_row instead. The counter is for lifecycle validation only.
  use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan
  implicit none
  integer, intent(in) :: packed_species,igamma
  complex(cheb_dp), intent(in) :: u
  complex(cheb_dp), intent(out) :: value
  integer, intent(out) :: status
  complex(cheb_dp), intent(out), optional :: derivative
  real(cheb_dp) :: nan
  rel_cheb_evaluations=rel_cheb_evaluations+1
  nan=ieee_value(0.d0,ieee_quiet_nan); value=cmplx(nan,nan,cheb_dp)
  if(present(derivative)) derivative=value
  status=CHEB_INVALID_MODEL
  if(.not.allocated(rel_cheb)) return
  if(packed_species<1.or.packed_species>size(rel_cheb,1)) return
  if(igamma<0.or.igamma>ubound(rel_cheb,2)) return
  status=rel_cheb_info(packed_species,igamma)%kernel_status
  if(rel_cheb_info(packed_species,igamma)%state/=REL_CHEB_FITTED.and. &
       rel_cheb_info(packed_species,igamma)%state/=REL_CHEB_REDUCED) return
  call cheb_evaluate(rel_cheb(packed_species,igamma),u,value,status,derivative)
end subroutine eval_rel_cheb


subroutine rel_cheb_export(metadata,mapping,coefficients,status)
  !! One contiguous column per row, Gamma-major then packed-species order.
  !! Minimal evaluator metadata: species/state/status/requested/effective/rank/
  !! sample_count/representation. Mapping: center/scale/Gamma. Two padded
  !! coefficient blocks hold c(0:Nmax) and d(0:Nmax). No fit diagnostics/workspaces.
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  integer, allocatable, intent(out) :: metadata(:,:)
  real(cheb_dp), allocatable, intent(out) :: mapping(:,:),coefficients(:,:)
  integer, intent(out) :: status
  integer :: s,g,row,nmax,rows,istat,n
  integer(int64) :: nr,nc
  status=CHEB_INVALID_MODEL
  if(.not.allocated(rel_cheb)) return
  nr=int(size(rel_cheb,1),int64)*size(rel_cheb,2)
  nmax=0
  do g=0,ubound(rel_cheb,2)
     do s=1,size(rel_cheb,1)
        nmax=max(nmax,rel_cheb(s,g)%effective_degree)
     enddo
  enddo
  nc=2_int64*(int(nmax,int64)+1)
  status=CHEB_ALLOCATION_FAILURE
  if(nr>huge(0)/8.or.nc>huge(0)) return
  if(nr*nc>huge(0)) return ! MPI count uses the project's default INTEGER ABI.
  rows=int(nr)
  allocate(metadata(8,rows),mapping(3,rows),coefficients(int(nc),rows),stat=istat)
  if(istat/=0) return
  coefficients=0.d0; row=0
  do g=0,ubound(rel_cheb,2)
     do s=1,size(rel_cheb,1)
        row=row+1
        metadata(:,row)=[rel_cheb_info(s,g)%species,rel_cheb_info(s,g)%state, &
             rel_cheb_info(s,g)%kernel_status,rel_cheb(s,g)%requested_degree, &
             rel_cheb(s,g)%effective_degree,rel_cheb(s,g)%rank,rel_cheb(s,g)%sample_count, &
             rel_cheb(s,g)%representation]
        mapping(:,row)=[rel_cheb(s,g)%center,rel_cheb(s,g)%scale,rel_cheb_info(s,g)%gamma]
        if(.not.allocated(rel_cheb(s,g)%coefficients)) cycle
        n=size(rel_cheb(s,g)%coefficients)
        coefficients(1:n,row)=rel_cheb(s,g)%coefficients
        n=size(rel_cheb(s,g)%derivative_coefficients)
        coefficients(nmax+2:nmax+1+n,row)=rel_cheb(s,g)%derivative_coefficients
     enddo
  enddo
  status=CHEB_SUCCESS
end subroutine rel_cheb_export


subroutine rel_cheb_import(nspecies,last_gamma,metadata,mapping,coefficients,status)
  !! Rebuild local allocatables without fitting. Omitted fit diagnostics are
  !! explicitly UNAVAILABLE/NaN; omitted sample-selection counts are -1.
  use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  integer, intent(in) :: nspecies,last_gamma,metadata(:,:)
  real(cheb_dp), intent(in) :: mapping(:,:),coefficients(:,:)
  integer, intent(out) :: status
  integer :: s,g,row,nmax,n,d,istat
  real(cheb_dp) :: nan
  call release_rel_cheb_rows
  status=CHEB_INVALID_VALUES
  if(nspecies<1.or.last_gamma<0) return
  if(size(metadata,1)/=8.or.size(mapping,1)/=3) return
  if(int(nspecies,int64)*(int(last_gamma,int64)+1)/=size(metadata,2)) return
  if(size(mapping,2)/=size(metadata,2).or.size(coefficients,2)/=size(metadata,2)) return
  if(mod(size(coefficients,1),2)/=0.or.size(coefficients,1)<2) return
  nmax=size(coefficients,1)/2-1
  allocate(rel_cheb(nspecies,0:last_gamma),rel_cheb_info(nspecies,0:last_gamma),stat=istat)
  status=CHEB_ALLOCATION_FAILURE
  if(istat/=0) then
     call release_rel_cheb_rows; return
  endif
  nan=ieee_value(0.d0,ieee_quiet_nan); row=0; status=CHEB_INVALID_MODEL
  do g=0,last_gamma
     do s=1,nspecies
        row=row+1
        rel_cheb_info(s,g)%species=metadata(1,row)
        rel_cheb_info(s,g)%state=metadata(2,row); rel_cheb_info(s,g)%kernel_status=metadata(3,row)
        rel_cheb_info(s,g)%packed_species=s; rel_cheb_info(s,g)%gamma_index=g
        rel_cheb_info(s,g)%gamma=mapping(3,row)
        rel_cheb_info(s,g)%input_samples=-1; rel_cheb_info(s,g)%geometric_samples=-1
        rel_cheb_info(s,g)%masked_samples=-1; rel_cheb_info(s,g)%positive_samples=-1
        rel_cheb(s,g)%requested_degree=metadata(4,row); rel_cheb(s,g)%effective_degree=metadata(5,row)
        rel_cheb(s,g)%rank=metadata(6,row); rel_cheb(s,g)%sample_count=metadata(7,row)
        rel_cheb(s,g)%representation=metadata(8,row); rel_cheb(s,g)%fit_status=metadata(3,row)
        rel_cheb(s,g)%center=mapping(1,row); rel_cheb(s,g)%scale=mapping(2,row)
        rel_cheb(s,g)%residual_norm=nan; rel_cheb(s,g)%relative_residual=nan; rel_cheb(s,g)%rms_residual=nan
        rel_cheb(s,g)%svd_tolerance=nan; rel_cheb(s,g)%singular_value_min=nan; rel_cheb(s,g)%singular_value_max=nan
        if(metadata(1,row)<=0) goto 900
        if(metadata(1,row)/=metadata(1,s)) goto 900
        if(s>1) then
           if(metadata(1,row)<=metadata(1,s-1)) goto 900
        endif
        if(metadata(8,row)/=CHEB_LOG10) goto 900
        select case(metadata(2,row))
        case(REL_CHEB_FITTED,REL_CHEB_REDUCED)
           if(metadata(2,row)==REL_CHEB_FITTED.and.metadata(3,row)/=CHEB_SUCCESS) goto 900
           if(metadata(2,row)==REL_CHEB_REDUCED.and.metadata(3,row)/=CHEB_REDUCED) goto 900
           d=metadata(5,row)
           if(d<0.or.d>nmax) goto 900
           allocate(rel_cheb(s,g)%coefficients(0:d), &
                rel_cheb(s,g)%derivative_coefficients(0:max(0,d-1)),stat=istat)
           if(istat/=0) then
              status=CHEB_ALLOCATION_FAILURE; goto 900
           endif
           rel_cheb(s,g)%coefficients=coefficients(1:d+1,row)
           n=max(1,d)
           rel_cheb(s,g)%derivative_coefficients=coefficients(nmax+2:nmax+1+n,row)
           if(.not.cheb_model_valid(rel_cheb(s,g))) goto 900
        case(REL_CHEB_NO_GEOMETRY)
           if(metadata(3,row)/=CHEB_NO_SAMPLES) goto 900
        case(REL_CHEB_ZERO)
           if(metadata(3,row)/=CHEB_NO_POSITIVE_SAMPLES) goto 900
        case(REL_CHEB_NOT_SELECTED)
           if(metadata(3,row)/=CHEB_INVALID_MODEL) goto 900
        case default
           goto 900 ! Unexpected root failures must be communicated before payload transfer.
        end select
     enddo
  enddo
  status=CHEB_SUCCESS
  return
900 call release_rel_cheb_rows
end subroutine rel_cheb_import


subroutine build_rel_cheb_rows(gamma,parallel,distribution,species,degrees,status,selected,diagnostic_unit)
  !! Explicit initialization entry: invoke on rank 0 after derivative_f0 has
  !! loaded/constructed and normalized the transformed tables. The conditional
  !! initialization hook remains unreachable through guarded relativistic input.
  !! species maps packed slots to original species; degrees=poly_order(species).
  !! selected optionally restricts applicable slots, without changing packing.
  !! Rebuilding first releases this collection only. Inputs are never modified.
  !! SUCCESS means construction completed without unexpected row failures;
  !! expected empty/zero rows still have unavailable models and explicit states.
  !! Otherwise return the first unexpected kernel error, retain classifications
  !! for every row, and report each failure. Diagnostics remain in log10 space.
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite,ieee_value,ieee_quiet_nan
  use, intrinsic :: iso_fortran_env, only : error_unit
  implicit none
  real(cheb_dp), intent(in) :: gamma(:,0:,0:),parallel(:,0:,0:),distribution(:,0:,0:)
  integer, intent(in) :: species(:),degrees(:)
  integer, intent(out) :: status
  logical, intent(in), optional :: selected(:)
  integer, intent(in), optional :: diagnostic_unit
  real(cheb_dp), allocatable :: nodes(:),values(:)
  integer :: ns,ng,np,s,g,j,n,code,istat,unit

  call release_rel_cheb_rows
  rel_cheb_build_calls=1
  status=CHEB_INVALID_VALUES; unit=error_unit
  if(present(diagnostic_unit)) unit=diagnostic_unit
  ns=size(distribution,1); ng=ubound(distribution,2); np=ubound(distribution,3)
  if(any(shape(gamma)/=shape(distribution)).or.any(shape(parallel)/=shape(distribution))) return
  if(size(species)/=ns.or.size(degrees)/=ns) return
  if(present(selected)) then
     if(size(selected)/=ns) return
  endif
  if(any(species<=0)) return
  do s=2,ns
     if(species(s)<=species(s-1)) return
  enddo
  if(ng<0.or.np<0) return
  allocate(rel_cheb(ns,0:ng),rel_cheb_info(ns,0:ng),nodes(0:np),values(0:np),stat=istat)
  if(istat/=0) then
     call release_rel_cheb_rows
     status=CHEB_ALLOCATION_FAILURE; return
  endif
  status=CHEB_SUCCESS
  do s=1,ns
     do g=0,ng
        rel_cheb_info(s,g)%species=species(s)
        rel_cheb_info(s,g)%packed_species=s; rel_cheb_info(s,g)%gamma_index=g
        rel_cheb_info(s,g)%input_samples=np+1
        rel_cheb_info(s,g)%gamma=ieee_value(0.d0,ieee_quiet_nan)
        rel_cheb(s,g)%representation=CHEB_LOG10
        rel_cheb(s,g)%requested_degree=degrees(s)
        if(present(selected)) then
           if(.not.selected(s)) then
              rel_cheb_info(s,g)%state=REL_CHEB_NOT_SELECTED
              cycle
           endif
        endif
        n=0
        do j=0,np
           ! Test finiteness before comparing: quiet NaNs must reach the
           ! kernel's validation, not raise invalid in a sentinel comparison.
           if(ieee_is_finite(distribution(s,g,j))) then
              if(distribution(s,g,j)==-1.d0) then
                 rel_cheb_info(s,g)%masked_samples=rel_cheb_info(s,g)%masked_samples+1
                 cycle
              endif
              if(distribution(s,g,j)>0.d0) &
                   rel_cheb_info(s,g)%positive_samples=rel_cheb_info(s,g)%positive_samples+1
           endif
           nodes(n)=parallel(s,g,j); values(n)=distribution(s,g,j); n=n+1
        enddo
        rel_cheb_info(s,g)%geometric_samples=n
        if(degrees(s)<0) then
           code=CHEB_INVALID_DEGREE
        elseif(.not.all(ieee_is_finite(gamma(s,g,:)))) then
           code=CHEB_INVALID_COORDINATES
        elseif(any(gamma(s,g,:)<1.d0)) then
           code=CHEB_INVALID_COORDINATES
        elseif(any(gamma(s,g,:)/=gamma(s,g,0))) then
           code=CHEB_INVALID_COORDINATES
        else
           rel_cheb_info(s,g)%gamma=gamma(s,g,0)
           if(n==0) then
              code=CHEB_NO_SAMPLES
           else
              rel_cheb_fit_calls=rel_cheb_fit_calls+1
              call cheb_fit(nodes(:n-1),values(:n-1),degrees(s),rel_cheb(s,g),code,representation=CHEB_LOG10)
           endif
        endif
        rel_cheb(s,g)%fit_status=code
        rel_cheb_info(s,g)%kernel_status=code
        select case(code)
        case(CHEB_SUCCESS)
           rel_cheb_info(s,g)%state=REL_CHEB_FITTED
        case(CHEB_REDUCED)
           rel_cheb_info(s,g)%state=REL_CHEB_REDUCED
        case(CHEB_NO_SAMPLES)
           rel_cheb_info(s,g)%state=REL_CHEB_NO_GEOMETRY
        case(CHEB_NO_POSITIVE_SAMPLES)
           rel_cheb_info(s,g)%state=REL_CHEB_ZERO
        case(CHEB_INVALID_COORDINATES,CHEB_INVALID_VALUES,CHEB_NEGATIVE_SAMPLES,CHEB_INTERIOR_ZERO)
           rel_cheb_info(s,g)%state=REL_CHEB_INVALID_SAMPLES
        case default
           rel_cheb_info(s,g)%state=REL_CHEB_FIT_FAILURE
        end select
        if(rel_cheb_info(s,g)%state==REL_CHEB_INVALID_SAMPLES.or. &
             rel_cheb_info(s,g)%state==REL_CHEB_FIT_FAILURE) then
           if(status==CHEB_SUCCESS) status=code
           write(unit,'(a,7(a,i0))') 'Relativistic Chebyshev row failure:', &
                ' species=',species(s),' packed_species=',s,' igamma=',g,' kernel_status=',code, &
                ' input=',np+1,' unmasked=',n,' positive=',rel_cheb_info(s,g)%positive_samples
        endif
     enddo
  enddo
end subroutine build_rel_cheb_rows


subroutine get_rel_cheb_row(packed_species,igamma,model,info,status)
  !! Read-only access by deep copy. Changes to the returned model cannot change
  !! stored coefficients. Unavailable rows return their kernel status and state;
  !! absent collections/invalid indices return INVALID_MODEL, NOT_BUILT.
  implicit none
  integer, intent(in) :: packed_species,igamma
  type(cheb_model), intent(out) :: model
  type(rel_cheb_row_info), intent(out) :: info
  integer, intent(out) :: status
  status=CHEB_INVALID_MODEL
  if(.not.allocated(rel_cheb)) return
  if(packed_species<1.or.packed_species>size(rel_cheb,1)) return
  if(igamma<0.or.igamma>ubound(rel_cheb,2)) return
  model=rel_cheb(packed_species,igamma); info=rel_cheb_info(packed_species,igamma)
  status=info%kernel_status
end subroutine get_rel_cheb_row


subroutine rel_cheb_collection_info(nspecies,last_gamma,fit_calls,coefficient_bytes,build_calls,evaluations)
  !! Introspection only; coefficient_bytes excludes derived-type/descriptor
  !! overhead and counts both value and derivative payloads. No fitting occurs.
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  integer, intent(out) :: nspecies,last_gamma,fit_calls
  integer(int64), intent(out) :: coefficient_bytes
  integer, intent(out), optional :: build_calls
  integer(int64), intent(out), optional :: evaluations
  integer :: s,g
  nspecies=0; last_gamma=-1; fit_calls=rel_cheb_fit_calls; coefficient_bytes=0_int64
  if(present(build_calls)) build_calls=rel_cheb_build_calls
  if(present(evaluations)) evaluations=rel_cheb_evaluations
  if(.not.allocated(rel_cheb)) return
  nspecies=size(rel_cheb,1); last_gamma=ubound(rel_cheb,2)
  do g=0,last_gamma
     do s=1,nspecies
        if(allocated(rel_cheb(s,g)%coefficients)) coefficient_bytes=coefficient_bytes+ &
             int(size(rel_cheb(s,g)%coefficients),int64)*(storage_size(0.d0)/8)
        if(allocated(rel_cheb(s,g)%derivative_coefficients)) coefficient_bytes=coefficient_bytes+ &
             int(size(rel_cheb(s,g)%derivative_coefficients),int64)*(storage_size(0.d0)/8)
     enddo
  enddo
end subroutine rel_cheb_collection_info


subroutine release_rel_cheb_rows
  !! Releasing parent arrays also releases all allocatable model components.
  implicit none
  if(allocated(rel_cheb)) deallocate(rel_cheb)
  if(allocated(rel_cheb_info)) deallocate(rel_cheb_info)
  rel_cheb_fit_calls=0
  rel_cheb_build_calls=0; rel_cheb_evaluations=0; rel_cheb_setup_status=CHEB_SUCCESS
end subroutine release_rel_cheb_rows


subroutine cheb_fit(u,values,requested_degree,model,status,rcond,representation)
  !! Default DIRECT fits values; LOG10 fits log10 of strictly positive samples.
  !! Only exact zeros strictly outside the positive coordinate interval are
  !! excluded. A zero coincident with a positive sample is an interior zero.
  !! No floors, physical support selection, or extrapolation cutoff is applied.
  !! SUCCESS/REDUCED apply to either representation; inspect representation.
  !! Finite data are supported. Quiet NaNs/infinities return errors; signaling
  !! NaNs can raise IEEE invalid before status return. FP modes are unchanged.
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  real(cheb_dp), intent(in) :: u(:),values(:)
  integer, intent(in) :: requested_degree
  type(cheb_model), intent(out) :: model
  integer, intent(out) :: status
  real(cheb_dp), intent(in), optional :: rcond
  integer, intent(in), optional :: representation
  real(cheb_dp), allocatable :: nodes(:),logs(:)
  real(cheb_dp) :: lower,upper
  integer :: rep,m,n,i,j,istat,nlower,nupper

  rep=CHEB_DIRECT
  if(present(representation)) rep=representation
  m=size(u)
  model%requested_degree=requested_degree; model%original_sample_count=m
  model%representation=rep
  model%residual_norm=ieee_value(0.d0,ieee_quiet_nan)
  model%relative_residual=model%residual_norm; model%rms_residual=model%residual_norm
  status=CHEB_INVALID_REPRESENTATION
  selection: block
     if(rep/=CHEB_DIRECT.and.rep/=CHEB_LOG10) exit selection
     if(rep==CHEB_DIRECT) then
        call cheb_fit_direct(u,values,requested_degree,model,status,rcond)
        exit selection
     endif
     status=CHEB_INVALID_DEGREE
     if(requested_degree<0) exit selection
     status=CHEB_INVALID_VALUES
     if(size(values)/=m) exit selection
     status=CHEB_NO_POSITIVE_SAMPLES
     if(m==0) exit selection
     status=CHEB_INVALID_COORDINATES
     if(.not.all(ieee_is_finite(u))) exit selection
     status=CHEB_INVALID_VALUES
     if(.not.all(ieee_is_finite(values))) exit selection
     status=CHEB_NEGATIVE_SAMPLES
     if(any(values<0.d0)) exit selection
     n=count(values>0.d0)
     status=CHEB_NO_POSITIVE_SAMPLES
     if(n==0) exit selection
     lower=minval(u,mask=values>0.d0); upper=maxval(u,mask=values>0.d0)
     model%retained_min=lower; model%retained_max=upper
     status=CHEB_INTERIOR_ZERO
     if(any(values==0.d0.and.u>=lower.and.u<=upper)) exit selection
     nlower=count(values==0.d0.and.u<lower); nupper=count(values==0.d0.and.u>upper)
     allocate(nodes(n),logs(n),stat=istat)
     status=CHEB_ALLOCATION_FAILURE
     if(istat/=0) exit selection
     j=0
     do i=1,m
        if(values(i)<=0.d0) cycle
        j=j+1; nodes(j)=u(i); logs(j)=log10(values(i))
     enddo
     call cheb_fit_direct(nodes,logs,requested_degree,model,status,rcond)
     model%lower_tail_zeros=nlower; model%upper_tail_zeros=nupper
  end block selection
  model%representation=rep; model%original_sample_count=m
  model%retained_sample_count=model%sample_count
  if(status==CHEB_SUCCESS.or.status==CHEB_REDUCED) then
     if(rep==CHEB_DIRECT) then
        model%retained_min=minval(u); model%retained_max=maxval(u)
     else
        model%retained_min=lower; model%retained_max=upper
     endif
  endif
  model%fit_status=status
end subroutine cheb_fit


subroutine cheb_fit_direct(u,values,requested_degree,model,status,rcond)
  !! Unweighted rectangular SVD fit on actual, possibly repeated/unsorted nodes.
  !! The default RCOND is max(number of samples,initial columns)*epsilon.
  !! A degenerate interval has center=u(1), scale=1 and is fitted as a constant.
  !! Rank loss causes a new lower-degree least-squares solve, never truncation.
  !! Input contract (also cheb_evaluate): finite data are supported; quiet NaNs
  !! and infinities are rejected by status. Signaling NaNs are not trap-safe:
  !! they may raise IEEE invalid before a status can be returned. No FP modes
  !! or exception traps are changed by this kernel.
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  real(cheb_dp), intent(in) :: u(:),values(:)
  integer, intent(in) :: requested_degree
  type(cheb_model), intent(out) :: model
  integer, intent(out) :: status
  real(cheb_dp), intent(in), optional :: rcond
  real(cheb_dp), allocatable :: x(:),scaled_values(:),coeff(:)
  real(cheb_dp) :: lower,upper,value_scale,t,y
  integer :: m,degree,i,k,istat,solver_status
  logical :: ok

  model%requested_degree=requested_degree
  model%residual_norm=ieee_value(0.d0,ieee_quiet_nan)
  model%relative_residual=model%residual_norm; model%rms_residual=model%residual_norm
  m=size(u); model%sample_count=m
  status=CHEB_SUCCESS
  fitting: block
     if(requested_degree<0) then
        status=CHEB_INVALID_DEGREE; exit fitting
     endif
     if(size(values)/=m) then
        status=CHEB_INVALID_VALUES; exit fitting
     endif
     if(m==0) then
        status=CHEB_NO_SAMPLES; exit fitting
     endif
     if(.not.all(ieee_is_finite(u))) then
        status=CHEB_INVALID_COORDINATES; exit fitting
     endif
     if(.not.all(ieee_is_finite(values))) then
        status=CHEB_INVALID_VALUES; exit fitting
     endif
     degree=min(requested_degree,m-1)
     model%sample_limited=degree<requested_degree
     model%svd_tolerance=dble(max(m,degree+1))*epsilon(1.d0)
     if(present(rcond)) model%svd_tolerance=rcond
     if(.not.ieee_is_finite(model%svd_tolerance)) then
        status=CHEB_INVALID_TOLERANCE; exit fitting
     endif
     if(model%svd_tolerance<0.d0.or.model%svd_tolerance>=1.d0) then
        status=CHEB_INVALID_TOLERANCE; exit fitting
     endif
     lower=minval(u); upper=maxval(u)
     ! Halve before adding/subtracting, avoiding overflow for large endpoints.
     model%center=0.5d0*lower+0.5d0*upper
     model%scale=0.5d0*upper-0.5d0*lower
     model%degenerate_interval=lower==upper
     if(model%degenerate_interval) then
        model%center=lower; model%scale=1.d0
     endif
     if(.not.ieee_is_finite(model%center).or..not.ieee_is_finite(model%scale)) then
        status=CHEB_INVALID_NORMALIZATION; exit fitting
     endif
     if(model%scale<=0.d0) then
        status=CHEB_INVALID_NORMALIZATION; exit fitting
     endif
     allocate(x(m),scaled_values(m),stat=istat)
     if(istat/=0) then
        status=CHEB_ALLOCATION_FAILURE; exit fitting
     endif
     do i=1,m
        call cheb_add(u(i),-model%center,t,ok)
        if(.not.ok) exit
        call cheb_divide(t,model%scale,x(i),ok)
        if(.not.ok) exit
     enddo
     if(.not.ok) then
        status=CHEB_INVALID_NORMALIZATION; exit fitting
     endif
     ! Scale the RHS, not the basis columns, to keep LAPACK and diagnostics safe.
     value_scale=maxval(abs(values))
     if(value_scale==0.d0) value_scale=1.d0
     scaled_values=values/value_scale
     do
        model%fit_attempts=model%fit_attempts+1
        call cheb_svd(x,scaled_values,degree,model%svd_tolerance,coeff,model%rank, &
             model%singular_value_min,model%singular_value_max,model%lapack_info,solver_status)
        if(solver_status/=CHEB_SUCCESS) then
           status=solver_status; exit fitting
        endif
        if(model%fit_attempts==1) model%initial_rank=model%rank
        if(model%rank==degree+1) exit
        model%rank_reduced=.true.
        if(degree==0.or.model%rank<1) then
           status=CHEB_RANK_FAILURE; exit fitting
        endif
        degree=min(degree-1,model%rank-1)
     enddo
     allocate(model%coefficients(0:degree),model%derivative_coefficients(0:max(0,degree-1)),stat=istat)
     if(istat/=0) then
        status=CHEB_ALLOCATION_FAILURE; exit fitting
     endif
     do k=0,degree
        call cheb_multiply(coeff(k),value_scale,model%coefficients(k),ok)
        if(.not.ok) exit
     enddo
     if(.not.ok) then
        status=CHEB_NONFINITE_COEFFICIENTS; exit fitting
     endif
     ! Differentiate the full-weight-c(0) convention. Form d(0) already halved
     ! to avoid an unnecessary overflow in 2*c(1), including degree one.
     model%derivative_coefficients=0.d0
     do k=degree-1,1,-1
        call cheb_multiply(2.d0*dble(k+1),model%coefficients(k+1),t,ok)
        if(.not.ok) exit
        if(k+2<=degree-1) then
           call cheb_add(t,model%derivative_coefficients(k+2),y,ok)
           if(.not.ok) exit
           t=y
        endif
        model%derivative_coefficients(k)=t
     enddo
     if(degree>=1.and.ok) then
        model%derivative_coefficients(0)=model%coefficients(1)
        if(degree>=3) then
           call cheb_add(model%coefficients(1),0.5d0*model%derivative_coefficients(2),y,ok)
           if(ok) model%derivative_coefficients(0)=y
        endif
     endif
     if(.not.ok) then
        status=CHEB_NONFINITE_COEFFICIENTS; exit fitting
     endif
     model%effective_degree=degree
     if(degree<requested_degree) status=CHEB_REDUCED
     ! Reporting failure never changes coefficients, rank, degree, or fit status.
     call cheb_residual_diagnostics(x,values,model)
  end block fitting
  model%fit_status=status
  if(status/=CHEB_SUCCESS.and.status/=CHEB_REDUCED) then
     if(allocated(model%coefficients)) deallocate(model%coefficients)
     if(allocated(model%derivative_coefficients)) deallocate(model%derivative_coefficients)
     model%effective_degree=-1
  endif
end subroutine cheb_fit_direct


subroutine cheb_residual_diagnostics(x,values,model)
  !! Diagnose the stored coefficients against the original data. Scale before
  !! evaluation/norm accumulation and reconstruct each diagnostic separately:
  !! an out-of-range L2 norm need not prevent a finite RMS or relative residual.
  use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan
  implicit none
  real(cheb_dp), intent(in) :: x(:),values(:)
  type(cheb_model), intent(inout) :: model
  real(cheb_dp), allocatable :: coeff(:),residual(:),scaled_values(:)
  real(cheb_dp) :: data_scale,value_scale,peak,unit_norm,y_norm,nan
  complex(cheb_dp) :: predicted
  integer :: i,m,istat
  logical :: ok
  nan=ieee_value(0.d0,ieee_quiet_nan)
  model%residual_norm=nan; model%relative_residual=nan; model%rms_residual=nan
  model%residual_norm_status=CHEB_DIAG_FAILED
  model%relative_residual_status=CHEB_DIAG_FAILED; model%rms_residual_status=CHEB_DIAG_FAILED
  m=size(x)
  allocate(coeff(0:model%effective_degree),residual(m),scaled_values(m),stat=istat)
  if(istat/=0) return
  value_scale=maxval(abs(values))
  data_scale=max(value_scale,maxval(abs(model%coefficients)))
  if(data_scale==0.d0) data_scale=1.d0
  coeff=model%coefficients/data_scale; scaled_values=values/data_scale
  ! Do not report an apparent zero caused by an unrepresentable rescaling.
  if(any(coeff==0.d0.and.model%coefficients/=0.d0)) return
  if(any(scaled_values==0.d0.and.values/=0.d0)) return
  do i=1,m
     call cheb_clenshaw(coeff,cmplx(x(i),0.d0,cheb_dp),predicted,ok)
     if(.not.ok) return
     call cheb_add(real(predicted,cheb_dp),-scaled_values(i),residual(i),ok)
     if(.not.ok) return
  enddo
  peak=maxval(abs(residual))
  if(peak==0.d0) then
     model%residual_norm=0.d0; model%relative_residual=0.d0; model%rms_residual=0.d0
     model%residual_norm_status=CHEB_DIAG_VALID
     model%relative_residual_status=CHEB_DIAG_VALID; model%rms_residual_status=CHEB_DIAG_VALID
     return
  endif
  unit_norm=norm2(residual/peak)
  call cheb_diagnostic_product([peak,unit_norm,data_scale],[1.d0], &
       model%residual_norm,model%residual_norm_status)
  call cheb_diagnostic_product([peak,unit_norm/sqrt(dble(m)),data_scale],[1.d0], &
       model%rms_residual,model%rms_residual_status)
  if(value_scale>0.d0) then
     y_norm=norm2(values/value_scale)
     call cheb_diagnostic_product([peak,unit_norm/y_norm,data_scale],[value_scale], &
          model%relative_residual,model%relative_residual_status)
  endif
end subroutine cheb_residual_diagnostics


pure subroutine cheb_diagnostic_product(numerator,denominator,value,status)
  !! Positive products/ratios in fraction/exponent form avoid unnecessary
  !! intermediate overflow/underflow. RANGE means overflow or underflow to zero
  !! at working precision; FAILED means computation/inputs were unavailable.
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  real(cheb_dp), intent(in) :: numerator(:),denominator(:)
  real(cheb_dp), intent(out) :: value
  integer, intent(out) :: status
  real(cheb_dp) :: mantissa,t
  integer :: power,i
  value=ieee_value(0.d0,ieee_quiet_nan); status=CHEB_DIAG_FAILED
  if(.not.all(ieee_is_finite(numerator)).or..not.all(ieee_is_finite(denominator))) return
  if(any(numerator<0.d0).or.any(denominator<=0.d0)) return
  if(any(numerator==0.d0)) then
     value=0.d0; status=CHEB_DIAG_VALID; return
  endif
  mantissa=1.d0; power=0
  do i=1,size(numerator)
     t=mantissa*fraction(numerator(i))
     power=power+exponent(numerator(i))+exponent(t)
     mantissa=fraction(t)
  enddo
  do i=1,size(denominator)
     t=mantissa/fraction(denominator(i))
     power=power-exponent(denominator(i))+exponent(t)
     mantissa=fraction(t)
  enddo
  status=CHEB_DIAG_RANGE
  if(power>maxexponent(value)) return
  t=scale(mantissa,power)
  if(.not.ieee_is_finite(t)) return
  if(t==0.d0) return
  value=t; status=CHEB_DIAG_VALID
end subroutine cheb_diagnostic_product


subroutine cheb_svd(x,y,degree,tolerance,coeff,rank,smin,smax,lapack_info,status)
  !! Private LAPACK adapter; default INTEGER matches the existing LP64 calls.
  !! DGELSD integer workspace: 3*n*NLVL+11*n, NLVL from ILAENV(9).
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  real(cheb_dp), intent(in) :: x(:),y(:),tolerance
  integer, intent(in) :: degree
  real(cheb_dp), allocatable, intent(out) :: coeff(:)
  integer, intent(out) :: rank,lapack_info,status
  real(cheb_dp), intent(out) :: smin,smax
  real(cheb_dp), allocatable :: a(:,:),b(:,:),singular(:),work(:)
  integer, allocatable :: iwork(:)
  real(cheb_dp) :: query(1),product,twice
  integer :: m,n,smlsiz,nlvl,liwork,lwork,k,i,istat
  logical :: ok
  m=size(x); n=degree+1; rank=0; lapack_info=0; smin=0.d0; smax=0.d0
  status=CHEB_LAPACK_FAILURE
  smlsiz=ilaenv(9,'DGELSD',' ',0,0,0,0)
  if(smlsiz<1.or.smlsiz>=huge(0)-1) return
  nlvl=max(0,int(log(dble(n)/dble(smlsiz+1))/log(2.d0))+1)
  status=CHEB_ALLOCATION_FAILURE
  if(n>huge(0)/(3*nlvl+11)) return
  liwork=max(1,n*(3*nlvl+11))
  allocate(a(m,n),b(m,1),singular(n),iwork(liwork),stat=istat)
  if(istat/=0) return
  a(:,1)=1.d0
  if(n>=2) a(:,2)=x
  do k=2,degree
     do i=1,m
        call cheb_multiply(x(i),a(i,k),product,ok)
        if(.not.ok) exit
        call cheb_add(product,product,twice,ok)
        if(.not.ok) exit
        call cheb_add(twice,-a(i,k-1),a(i,k+1),ok)
        if(.not.ok) exit
     enddo
     if(.not.ok) then
        status=CHEB_RANK_FAILURE; return
     endif
  enddo
  if(.not.all(ieee_is_finite(a))) then
     status=CHEB_INVALID_NORMALIZATION; return
  endif
  b(:,1)=y
  call dgelsd(m,n,1,a,m,b,m,singular,tolerance,rank,query,-1,iwork,lapack_info)
  status=CHEB_LAPACK_FAILURE
  if(lapack_info/=0) return
  if(.not.ieee_is_finite(query(1))) return
  if(query(1)<1.d0.or.query(1)>dble(huge(0)-1)) return
  lwork=ceiling(query(1))
  allocate(work(lwork),coeff(0:degree),stat=istat)
  status=CHEB_ALLOCATION_FAILURE
  if(istat/=0) return
  call dgelsd(m,n,1,a,m,b,m,singular,tolerance,rank,work,lwork,iwork,lapack_info)
  status=CHEB_LAPACK_FAILURE
  if(lapack_info/=0) return
  if(rank<0.or.rank>n) then
     status=CHEB_RANK_FAILURE; return
  endif
  if(.not.all(ieee_is_finite(singular))) return
  coeff=b(1:n,1)
  status=CHEB_NONFINITE_COEFFICIENTS
  if(.not.all(ieee_is_finite(coeff))) return
  smax=singular(1); smin=singular(n)
  status=CHEB_SUCCESS
end subroutine cheb_svd


pure subroutine cheb_evaluate_real(model,u,value,status,derivative,value_status,derivative_status)
  implicit none
  type(cheb_model), intent(in) :: model
  real(cheb_dp), intent(in) :: u
  real(cheb_dp), intent(out) :: value
  integer, intent(out) :: status
  real(cheb_dp), intent(out), optional :: derivative
  integer, intent(out), optional :: value_status,derivative_status
  complex(cheb_dp) :: z,dz
  if(present(derivative)) then
     call cheb_evaluate_complex(model,cmplx(u,0.d0,cheb_dp),z,status,dz,value_status,derivative_status)
     derivative=real(dz,cheb_dp)
  else
     call cheb_evaluate_complex(model,cmplx(u,0.d0,cheb_dp),z,status, &
          value_status=value_status,derivative_status=derivative_status)
  endif
  value=real(z,cheb_dp)
end subroutine cheb_evaluate_real


pure subroutine cheb_evaluate_complex(model,u,value,status,derivative,value_status,derivative_status)
  !! Default evaluation reconstructs the physical value for LOG10 models.
  !! cheb_evaluate_polynomial always returns the raw series P and dP/du.
  !! UNDERFLOW is a nonfatal evaluation condition: affected nonzero components
  !! below the subnormal range become zero, while other components survive.
  !! Optional per-output statuses distinguish value and derivative underflow.
  !! Hard failures publish only quiet NaNs; the fitted model remains valid.
  use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan
  implicit none
  type(cheb_model), intent(in) :: model
  complex(cheb_dp), intent(in) :: u
  complex(cheb_dp), intent(out) :: value
  integer, intent(out) :: status
  complex(cheb_dp), intent(out), optional :: derivative
  integer, intent(out), optional :: value_status,derivative_status
  complex(cheb_dp) :: p,dpdu,v,d
  real(cheb_dp) :: nan
  integer :: vs,ds
  nan=ieee_value(0.d0,ieee_quiet_nan)
  value=cmplx(nan,nan,cheb_dp)
  if(present(derivative)) derivative=value
  if(present(derivative)) then
     call cheb_polynomial_complex(model,u,p,status,dpdu)
  else
     call cheb_polynomial_complex(model,u,p,status)
  endif
  vs=status; ds=status
  if(status==CHEB_SUCCESS) then
     if(model%representation==CHEB_DIRECT) then
        value=p
        if(present(derivative)) derivative=dpdu
     else
        if(present(derivative)) then
           call cheb_reconstruct(p,v,vs,dpdu,d,ds)
        else
           call cheb_reconstruct(p,v,vs)
           ds=CHEB_SUCCESS
        endif
        status=vs
        if(vs==CHEB_SUCCESS.or.vs==CHEB_LOG_UNDERFLOW) then
           if(ds==CHEB_SUCCESS.or.ds==CHEB_LOG_UNDERFLOW) then
              value=v
              if(present(derivative)) derivative=d
              if(ds==CHEB_LOG_UNDERFLOW) status=CHEB_LOG_UNDERFLOW
           else
              status=CHEB_LOG_DERIVATIVE_FAILURE
           endif
        endif
     endif
  endif
  if(present(value_status)) value_status=vs
  if(present(derivative_status)) then
     derivative_status=CHEB_NOT_REQUESTED
     if(present(derivative)) derivative_status=ds
  endif
end subroutine cheb_evaluate_complex


pure subroutine cheb_reconstruct(p,value,status,dpdu,derivative,derivative_status)
  !! Exp is applied only to a small binary remainder; SCALE restores range.
  !! Derivative components combine ln(f), ln(ln(10)), and ln(|P_u|) before
  !! exponentiation, so an underflowing f can have a representable derivative.
  !! Trig routines perform argument reduction. Reject phases above 1/sqrt(eps)
  !! radians: multiplication by ln(10) then loses O(sqrt(eps)) absolute phase
  !! accuracy. This numerical reliability limit is unrelated to physical support.
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  complex(cheb_dp), intent(in) :: p
  complex(cheb_dp), intent(out) :: value
  integer, intent(out) :: status
  complex(cheb_dp), intent(in), optional :: dpdu
  complex(cheb_dp), intent(out), optional :: derivative
  integer, intent(out), optional :: derivative_status
  real(cheb_dp), parameter :: ln10=log(10.d0)
  real(cheb_dp) :: nan,t,phase,magnitude,re,im,s,logd,dr,di,tmp
  complex(cheb_dp) :: direction,rotated,scaled
  integer :: sr,si,ds
  logical :: ok
  nan=ieee_value(0.d0,ieee_quiet_nan); value=cmplx(nan,nan,cheb_dp)
  if(present(derivative)) derivative=value
  if(present(derivative_status)) derivative_status=CHEB_LOG_DERIVATIVE_FAILURE
  status=CHEB_PHASE_FAILURE
  ! Check before multiplication: an extreme but finite phase must not trap.
  if(abs(aimag(p))>1.d0/(sqrt(epsilon(1.d0))*ln10)) return
  phase=aimag(p)*ln10
  direction=cmplx(cos(phase),sin(phase),cheb_dp)
  call cheb_multiply(real(p,cheb_dp),ln10,t,ok)
  if(.not.ok) then
     if(real(p,cheb_dp)>0.d0) then
        status=CHEB_LOG_OVERFLOW
     else
        ! No finite P_u can compensate for an exponent of this magnitude.
        value=cmplx(0.d0,0.d0,cheb_dp); status=CHEB_LOG_UNDERFLOW
        if(present(derivative)) derivative=value
        if(present(derivative_status)) then
           derivative_status=CHEB_LOG_UNDERFLOW
           if(present(dpdu)) then
              if(dpdu==cmplx(0.d0,0.d0,cheb_dp)) derivative_status=CHEB_SUCCESS
           endif
        endif
     endif
     return
  endif
  call cheb_exp_component(t,1.d0,magnitude,status)
  if(status/=CHEB_SUCCESS.and.status/=CHEB_LOG_UNDERFLOW) return
  re=magnitude*real(direction,cheb_dp); im=magnitude*aimag(direction)
  if(magnitude/=0.d0) then
     if((re==0.d0.and.real(direction,cheb_dp)/=0.d0).or.(im==0.d0.and.aimag(direction)/=0.d0)) &
          status=CHEB_LOG_UNDERFLOW
  endif
  if(.not.ieee_is_finite(re).or..not.ieee_is_finite(im)) then
     status=CHEB_LOG_VALUE_FAILURE; return
  endif
  value=cmplx(re,im,cheb_dp)
  if(.not.present(derivative)) return
  ds=CHEB_SUCCESS; dr=0.d0; di=0.d0
  s=max(abs(real(dpdu,cheb_dp)),abs(aimag(dpdu)))
  if(s/=0.d0) then
     ! Each rotated, s-normalized component is bounded by two. Below this
     ! threshold even the largest finite P_u component cannot rescue a result.
     if(t<log(nearest(0.d0,1.d0))-log(huge(1.d0))-log(2.d0*ln10)) then
        ds=CHEB_LOG_UNDERFLOW
     else
        scaled=cmplx(real(dpdu,cheb_dp)/s,aimag(dpdu)/s,cheb_dp)
        call cheb_complex_multiply(direction,scaled,rotated,ok)
        if(.not.ok) then
           derivative_status=CHEB_LOG_DERIVATIVE_FAILURE; return
        endif
        call cheb_add(t,log(ln10),tmp,ok)
        if(.not.ok) then
           derivative_status=CHEB_LOG_DERIVATIVE_FAILURE; return
        endif
        call cheb_add(tmp,log(s),logd,ok)
        if(.not.ok) then
           derivative_status=CHEB_LOG_DERIVATIVE_FAILURE; return
        endif
        call cheb_exp_component(logd,real(rotated,cheb_dp),dr,sr)
        call cheb_exp_component(logd,aimag(rotated),di,si)
        if(sr==CHEB_LOG_OVERFLOW.or.si==CHEB_LOG_OVERFLOW) then
           ds=CHEB_LOG_OVERFLOW
        elseif(sr/=CHEB_SUCCESS.and.sr/=CHEB_LOG_UNDERFLOW) then
           ds=CHEB_LOG_DERIVATIVE_FAILURE
        elseif(si/=CHEB_SUCCESS.and.si/=CHEB_LOG_UNDERFLOW) then
           ds=CHEB_LOG_DERIVATIVE_FAILURE
        elseif(sr==CHEB_LOG_UNDERFLOW.or.si==CHEB_LOG_UNDERFLOW) then
           ds=CHEB_LOG_UNDERFLOW
        endif
     endif
  endif
  if(ds==CHEB_SUCCESS.or.ds==CHEB_LOG_UNDERFLOW) then
     if(ieee_is_finite(dr).and.ieee_is_finite(di)) then
        derivative=cmplx(dr,di,cheb_dp)
     else
        ds=CHEB_LOG_DERIVATIVE_FAILURE
     endif
  endif
  if(present(derivative_status)) derivative_status=ds
end subroutine cheb_reconstruct


pure subroutine cheb_exp_component(logamp,factor,value,status)
  !! Reconstruct factor*exp(logamp) without an out-of-range exp intermediate.
  !! A nonzero magnitude below the least subnormal is explicitly zero/UNDERFLOW.
  !! Near the upper limit, rounded logarithms can conservatively reject a value;
  !! no boundary value is clamped and no overflowing SCALE is ever attempted.
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  real(cheb_dp), intent(in) :: logamp,factor
  real(cheb_dp), intent(out) :: value
  integer, intent(out) :: status
  real(cheb_dp), parameter :: ln2=log(2.d0)
  real(cheb_dp) :: t,mantissa
  integer :: power
  logical :: ok
  value=ieee_value(0.d0,ieee_quiet_nan); status=CHEB_LOG_VALUE_FAILURE
  if(.not.ieee_is_finite(logamp).or..not.ieee_is_finite(factor)) return
  if(factor==0.d0) then
     value=0.d0; status=CHEB_SUCCESS; return
  endif
  call cheb_add(logamp,log(abs(factor)),t,ok)
  if(.not.ok) return
  status=CHEB_LOG_OVERFLOW
  if(t>log(huge(1.d0))) return
  if(t<log(nearest(0.d0,1.d0))) then
     value=sign(0.d0,factor); status=CHEB_LOG_UNDERFLOW; return
  endif
  power=floor(t/ln2)
  mantissa=0.5d0*exp(t-dble(power)*ln2); power=power+1
  ! Renormalize: remainder rounding can place the mantissa at exactly one.
  power=power+exponent(mantissa); mantissa=fraction(mantissa)
  if(power>maxexponent(1.d0)) return
  value=sign(scale(mantissa,power),factor)
  status=CHEB_SUCCESS
  if(value==0.d0) status=CHEB_LOG_UNDERFLOW
  if(.not.ieee_is_finite(value)) then
     value=ieee_value(0.d0,ieee_quiet_nan); status=CHEB_LOG_VALUE_FAILURE
  endif
end subroutine cheb_exp_component


pure subroutine cheb_polynomial_real(model,u,value,status,derivative)
  implicit none
  type(cheb_model), intent(in) :: model
  real(cheb_dp), intent(in) :: u
  real(cheb_dp), intent(out) :: value
  integer, intent(out) :: status
  real(cheb_dp), intent(out), optional :: derivative
  complex(cheb_dp) :: z,dz
  if(present(derivative)) then
     call cheb_polynomial_complex(model,cmplx(u,0.d0,cheb_dp),z,status,dz)
     derivative=real(dz,cheb_dp)
  else
     call cheb_polynomial_complex(model,cmplx(u,0.d0,cheb_dp),z,status)
  endif
  value=real(z,cheb_dp)
end subroutine cheb_polynomial_real


pure subroutine cheb_polynomial_complex(model,u,value,status,derivative)
  !! No support/unit-disk cutoff. Failure outputs are quiet NaNs, not zeros.
  !! Checked arithmetic works even when the caller enables overflow traps.
  !! Quiet NaNs/infinities are rejected; signaling NaNs may raise IEEE invalid
  !! before a status is returned. The caller's FP environment is never changed.
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  type(cheb_model), intent(in) :: model
  complex(cheb_dp), intent(in) :: u
  complex(cheb_dp), intent(out) :: value
  integer, intent(out) :: status
  complex(cheb_dp), intent(out), optional :: derivative
  complex(cheb_dp) :: x,p,dpdx
  real(cheb_dp) :: t,xr,xi,dr,di,nan
  logical :: ok
  nan=ieee_value(0.d0,ieee_quiet_nan)
  value=cmplx(nan,nan,cheb_dp)
  if(present(derivative)) derivative=value
  status=CHEB_INVALID_REPRESENTATION
  if(model%representation/=CHEB_DIRECT.and.model%representation/=CHEB_LOG10) return
  status=CHEB_INVALID_MODEL
  if(.not.cheb_model_valid(model)) return
  status=CHEB_INVALID_COORDINATES
  if(.not.ieee_is_finite(real(u,cheb_dp)).or..not.ieee_is_finite(aimag(u))) return
  status=CHEB_EVALUATION_FAILURE
  if(model%effective_degree==0) then
     p=cmplx(model%coefficients(0),0.d0,cheb_dp); dr=0.d0; di=0.d0
  else
     call cheb_add(real(u,cheb_dp),-model%center,t,ok)
     if(.not.ok) return
     call cheb_divide(t,model%scale,xr,ok)
     if(.not.ok) return
     call cheb_divide(aimag(u),model%scale,xi,ok)
     if(.not.ok) return
     x=cmplx(xr,xi,cheb_dp)
     call cheb_clenshaw(model%coefficients,x,p,ok)
     if(.not.ok) return
     if(present(derivative)) then
        call cheb_clenshaw(model%derivative_coefficients,x,dpdx,ok)
        if(.not.ok) return
        call cheb_divide(real(dpdx,cheb_dp),model%scale,dr,ok)
        if(.not.ok) return
        call cheb_divide(aimag(dpdx),model%scale,di,ok)
        if(.not.ok) return
     endif
  endif
  ! Publish outputs atomically, only after all requested results are finite.
  if(.not.ieee_is_finite(real(p,cheb_dp)).or..not.ieee_is_finite(aimag(p))) return
  if(present(derivative)) then
     if(.not.ieee_is_finite(dr).or..not.ieee_is_finite(di)) return
     derivative=cmplx(dr,di,cheb_dp)
  endif
  value=p; status=CHEB_SUCCESS
end subroutine cheb_polynomial_complex


pure logical function cheb_model_valid(model) result(valid)
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  type(cheb_model), intent(in) :: model
  integer :: degree
  valid=.false.; degree=model%effective_degree
  if(model%fit_status/=CHEB_SUCCESS.and.model%fit_status/=CHEB_REDUCED) return
  if(degree<0.or.degree>model%requested_degree) return
  if(degree==huge(0)) return
  if(model%rank/=degree+1.or.model%sample_count<degree+1) return
  if(.not.allocated(model%coefficients).or..not.allocated(model%derivative_coefficients)) return
  if(lbound(model%coefficients,1)/=0.or.ubound(model%coefficients,1)/=degree) return
  if(lbound(model%derivative_coefficients,1)/=0.or. &
       ubound(model%derivative_coefficients,1)/=max(0,degree-1)) return
  if(.not.ieee_is_finite(model%center).or..not.ieee_is_finite(model%scale)) return
  if(model%scale<=0.d0) return
  if(.not.all(ieee_is_finite(model%coefficients))) return
  if(.not.all(ieee_is_finite(model%derivative_coefficients))) return
  valid=.true.
end function cheb_model_valid


pure subroutine cheb_clenshaw(coeff,x,value,ok)
  implicit none
  real(cheb_dp), intent(in) :: coeff(0:)
  complex(cheb_dp), intent(in) :: x
  complex(cheb_dp), intent(out) :: value
  logical, intent(out) :: ok
  complex(cheb_dp) :: b1,b2,b,t,twice
  integer :: k
  b1=cmplx(0.d0,0.d0,cheb_dp); b2=b1; ok=.true.
  do k=ubound(coeff,1),1,-1
     call cheb_complex_multiply(x,b1,t,ok)
     if(.not.ok) return
     call cheb_complex_add(t,t,twice,ok)
     if(.not.ok) return
     call cheb_complex_add(cmplx(coeff(k),0.d0,cheb_dp),twice,t,ok)
     if(.not.ok) return
     call cheb_complex_add(t,-b2,b,ok)
     if(.not.ok) return
     b2=b1; b1=b
  enddo
  call cheb_complex_multiply(x,b1,t,ok)
  if(.not.ok) return
  call cheb_complex_add(cmplx(coeff(0),0.d0,cheb_dp),t,b,ok)
  if(.not.ok) return
  call cheb_complex_add(b,-b2,value,ok)
end subroutine cheb_clenshaw


pure subroutine cheb_add(a,b,value,ok)
  !! For like signs, correct the rounded H-small bound with its exact residual.
  !! small<=H/2 makes H-bound exact (Sterbenz); the residual's sign identifies
  !! an upward-rounded bound. Unlike a blanket ULP margin, exact boundary sums
  !! such as (H-spacing(H))+spacing(H) remain accepted. Opposite signs are safe.
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  real(cheb_dp), intent(in) :: a,b
  real(cheb_dp), intent(out) :: value
  logical, intent(out) :: ok
  real(cheb_dp) :: large,small,bound,error
  ok=.false.; value=0.d0
  if(.not.ieee_is_finite(a).or..not.ieee_is_finite(b)) return
  if((a>0.d0.and.b>0.d0).or.(a<0.d0.and.b<0.d0)) then
     large=max(abs(a),abs(b)); small=min(abs(a),abs(b))
     if(small>0.5d0*huge(a)) return
     bound=huge(a)-small
     error=(huge(a)-bound)-small
     if(large>bound) return
     if(large==bound.and.error<0.d0) return
  endif
  value=a+b; ok=ieee_is_finite(value)
end subroutine cheb_add

pure subroutine cheb_multiply(a,b,value,ok)
  !! Test the rounded product in a unit-size binade before multiplying large
  !! operands. Power-of-radix scaling preserves rounding at the overflow end;
  !! no potentially upward-rounded H/abs(b) threshold is used.
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  real(cheb_dp), intent(in) :: a,b
  real(cheb_dp), intent(out) :: value
  logical, intent(out) :: ok
  real(cheb_dp) :: product,bound
  integer :: power
  ok=.false.; value=0.d0
  if(.not.ieee_is_finite(a).or..not.ieee_is_finite(b)) return
  if(abs(a)>1.d0.and.abs(b)>1.d0) then
     power=exponent(a)+exponent(b)
     if(power>maxexponent(a)+1) return
     if(power>=maxexponent(a)) then
        product=abs(fraction(a))*abs(fraction(b))
        bound=scale(fraction(huge(a)),maxexponent(a)-power)
        if(product>bound) return
     endif
  endif
  value=a*b; ok=ieee_is_finite(value)
end subroutine cheb_multiply

pure subroutine cheb_divide(a,b,value,ok)
  !! b is a previously validated, positive coordinate scale.
  !! A fraction/exponent comparison also protects division without forming a
  !! rounding-sensitive H*b bound. Zero/unit/non-increasing cases remain exact.
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  real(cheb_dp), intent(in) :: a,b
  real(cheb_dp), intent(out) :: value
  logical, intent(out) :: ok
  real(cheb_dp) :: quotient,bound
  integer :: power
  ok=.false.; value=0.d0
  if(.not.ieee_is_finite(a).or..not.ieee_is_finite(b)) return
  if(b<=0.d0) return
  if(b<1.d0.and.a/=0.d0) then
     power=exponent(a)-exponent(b)
     if(power>maxexponent(a)) return
     if(power>=maxexponent(a)-1) then
        quotient=abs(fraction(a))/fraction(b)
        bound=scale(fraction(huge(a)),maxexponent(a)-power)
        if(quotient>bound) return
     endif
  endif
  value=a/b; ok=ieee_is_finite(value)
end subroutine cheb_divide

pure subroutine cheb_complex_add(a,b,value,ok)
  implicit none
  complex(cheb_dp), intent(in) :: a,b
  complex(cheb_dp), intent(out) :: value
  logical, intent(out) :: ok
  real(cheb_dp) :: re,im
  call cheb_add(real(a,cheb_dp),real(b,cheb_dp),re,ok)
  if(.not.ok) return
  call cheb_add(aimag(a),aimag(b),im,ok)
  if(ok) value=cmplx(re,im,cheb_dp)
end subroutine cheb_complex_add

pure subroutine cheb_complex_multiply(a,b,value,ok)
  implicit none
  complex(cheb_dp), intent(in) :: a,b
  complex(cheb_dp), intent(out) :: value
  logical, intent(out) :: ok
  real(cheb_dp) :: ar,ai,br,bi,p,q,re,im
  ar=real(a,cheb_dp); ai=aimag(a); br=real(b,cheb_dp); bi=aimag(b)
  call cheb_multiply(ar,br,p,ok)
  if(.not.ok) return
  call cheb_multiply(ai,bi,q,ok)
  if(.not.ok) return
  call cheb_add(p,-q,re,ok)
  if(.not.ok) return
  call cheb_multiply(ar,bi,p,ok)
  if(.not.ok) return
  call cheb_multiply(ai,br,q,ok)
  if(.not.ok) return
  call cheb_add(p,q,im,ok)
  if(ok) value=cmplx(re,im,cheb_dp)
end subroutine cheb_complex_multiply



double complex function eval_fit(is,iperp,ppar_valC,dfpar,status)
  !! This function evaluates the fit to f0 at and the complex parallel
  !! momentum ppar_valC. It requires the fit parameters that will be determined
  !! by the subroutine [[determine_param_fit(subroutine)]].

        use alps_var, only : fit_type, pp, param_fit, n_fits, gamma_rel,nspec,relativistic
        use alps_var, only : ACmethod, poly_fit_coeffs, poly_order
		use alps_io, only : alps_error,alps_continuation_error
	use alps_distribution_analyt, only : distribution_analyt
	implicit none


	integer, intent(in) :: is
	!! Index of species for which [[eval_fit(function)]] is executed.

	integer, intent(in) :: iperp
	!! Index of perpendicular momentum at which [[eval_fit(function)]] is executed.

	double complex, intent(in) :: ppar_valC
	!! Complex parallel momentum at which [[eval_fit(function)]] is executed.

	double precision :: pperp_val
	!! Perpendicular momentum corrsponding to index iperp.

	integer :: par_ind
	!! Parameter index for the fit parameters.

	integer :: ifit
	!! Index for the fits that add up to the total fit for a given species.

	integer :: n_params
	!! Total number of fit parameters for a given species.

	integer :: is_rel
	!! Index for relativistic species (if any).

	integer :: sproc_rel
	!! MPI process number for relativistic species.

	integer :: is_run
	!! Running variable over species to determine relativistic evaluation.

	double precision,allocatable, dimension(:) :: params
	!! Array of fit parameters.

	double complex, intent(out), optional :: dfpar
	!! Analytic derivative with respect to parallel momentum from polynomial representation. 

integer, intent(out), optional :: status
integer :: cheb_status
! For relativistic method 2, iperp is the existing Gamma row and ppar_valC
! is already u=pparbar. No coordinate transformation or derivative stencil
! belongs in this adapter. Calls without status fail explicitly on any
! non-SUCCESS result, including underflow; no silent zero policy is introduced.
if(present(status)) status=CHEB_SUCCESS
if(ACmethod(is)==2.and.relativistic(is)) then
   sproc_rel=count(relativistic(1:is))
   call eval_rel_cheb(sproc_rel,iperp,ppar_valC,eval_fit,cheb_status,dfpar)
   if(present(status)) then
      status=cheb_status
   elseif(cheb_status/=CHEB_SUCCESS) then
      call alps_continuation_error(is,iperp,cheb_status,ppar_valC)
   endif
   return
endif

if (present(dfpar)) then
   if (ACmethod(is) /= 2) then
	call alps_error(11)
   endif
endif

 select case (ACmethod(is))
 case (0)
    ! Use the pre-coded distribution from distribution/distribution_analyt.f90
    if (relativistic(is)) then

			! Determine your sproc_rel:
			 is_rel=0
			 do is_run=1,nspec
				if (relativistic(is_run)) then
					is_rel=is_rel+1
					if (is_run.EQ.is) sproc_rel=is_rel
				endif
			 enddo

		eval_fit=distribution_analyt(is,gamma_rel(sproc_rel,iperp,1),ppar_valC,1,0.d0,0.d0)

	else !non-relativistic:

		eval_fit=distribution_analyt(is,pp(is,iperp,1,1),ppar_valC,1,0.d0,0.d0)
	endif


    return


 case (1)
    ! Use the 'n_fits' functions described with 'fit_type'
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

case (2)
   ! use the orthogonal polynomials of kind described in &poly_is
   
   if (present(dfpar)) then
	  dfpar = (0.d0,0.d0)      
      eval_fit = fit_function_poly(is,iperp,ppar_valC,poly_order(is),&
           poly_fit_coeffs(is,iperp,0:poly_order(is)),dfpar=dfpar)
   else
      eval_fit = fit_function_poly(is,iperp,ppar_valC,poly_order(is),&
           poly_fit_coeffs(is,iperp,0:poly_order(is)))
   endif
   
end select
end function eval_fit



double complex function fit_function(is,n_params,params,pperp_val,ppar_val)
	!! This function evaluates the fit to f0 at real pperp_val and complex ppar_val,
	!! provided that the one-dimensional fit-parameter array params is fed into the
	!! function. This is only used during the fitting. For the evaluation in ALPS,
	!! use [[eval_fit(function)]].
	use alps_var, only : fit_type, n_fits, ms, vA, perp_correction
	implicit none

	integer, intent(in) :: is
	!! Index of species for which [[eval_fit(function)]] is executed.

	integer, intent(in) :: n_params
	!! Total number of fit parameters for a given species.

	double precision, intent(in) :: params(n_params)
	!! Array of fit parameters.

	double precision, intent(in) :: pperp_val
	!! Perpendicular momentum.

	double complex, intent(in) :: ppar_val
	!! Complex parallel momentum.

	integer :: ifit
	!! Index for the fits that add up to the total fit for a given species.

	integer :: par_ind
	!! Parameter index to loop over fits.

	double complex :: sqrtpart
	!! Square-root part of Juettner distribution.

	double complex :: kappapart
	!! Kappa part of the kappa distribution.


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
end function fit_function

double complex function fit_function_poly(is,iperp,ppar_val,n_poly,fit_coeffs,dfpar)
  !! This function evaluates the orthogonal polynomical fit to f0
  !! at real pperp_val and complex ppar_val,
  !! with the one-dimensional polynomical coefficient array fit_coeffs is fed into the
  !! function.
  !! For the evaluation in ALPS, use [[eval_fit(function)]].
  use alps_var, only : pp, poly_kind, npar, logfit
  use alps_var, only : poly_log_max
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
 !! Index of perpendicular momentum at which [[fit_function_poly(function)]] is executed.

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

 fit_function_poly=cmplx(0.d0,0.d0,kind(1.d0))
 if (present(dfpar)) dfpar = (0.d0,0.d0)
 poly_deriv = (0.d0,0.d0)
 
 select case (poly_kind(is))
 case (1)
    !Chebyshev polynomials range from [-1,1]
    norm_1=5.d-1*(pp(is,iperp,npar,2)+pp(is,iperp,0,2))
    norm_2=5.d-1*(pp(is,iperp,npar,2)-pp(is,iperp,0,2))
    ppar_val_tmp=(ppar_val-norm_1)/norm_2

    if (abs(ppar_val_tmp).gt.1.d0) then !!kgk: should this be .gt. or .ge. ?
       fit_function_poly=0.d0
    else
	!Set n=0
    n=0
    poly_basis(n)=cmplx(1.d0,0.d0,kind(1.d0))
	if (present(dfpar)) &
		poly_basis_deriv(n) = (0.d0,0.d0)
    fit_function_poly=fit_function_poly+fit_coeffs(n)*poly_basis(n)
    
	!Set n=1
	n=1
    poly_basis(n)=ppar_val_tmp
	if (present(dfpar)) &
		poly_basis_deriv(n) = (1.d0,0.d0)
    fit_function_poly=fit_function_poly+fit_coeffs(n)*poly_basis(n)

	!Start Summing Derivatives
	if (present(dfpar)) &
		poly_deriv = fit_coeffs(n)

	!Start recursion
    do n=2,n_poly
       poly_basis(n) = cmplx(2.d0,0.d0,kind(1.d0)) * &
			ppar_val_tmp * poly_basis(n-1) - poly_basis(n-2)

       fit_function_poly=fit_function_poly+&
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

       !if (abs(fit_function_poly).gt.20.) then
       !KGK: how does this impact the 'noise'
       !if ((real(fit_function_poly).gt.1.).or.(aimag(fit_function_poly).gt.1.)&
       !     .or.(real(fit_function_poly).lt.-18.).or.(aimag(fit_function_poly).lt.-18.)) then

       !if ((real(fit_function_poly).gt.10.).or.(aimag(fit_function_poly).gt.10.)&
       !     .or.(real(fit_function_poly).lt.-22.).or.(aimag(fit_function_poly).lt.-22.)) then

       !KGK: Need to investigate the appropriate threshold here for when this functional
       !representation ceases to work.
       !if ((real(fit_function_poly).gt.17.).or.(aimag(fit_function_poly).gt.17.)&
         !if ((real(fit_function_poly).lt.-poly_log_max(is)).or.& !??
              !(aimag(fit_function_poly).lt.-poly_log_max(is))) then

       !if ((real(fit_function_poly).lt.-poly_log_max(is)).or.& !??
       !     (aimag(fit_function_poly).lt.-poly_log_max(is)).or.&
       !     (real(fit_function_poly).gt. poly_log_max(is)).or.& !??
       !     (aimag(fit_function_poly).gt.poly_log_max(is))) then

       !KGK: 2026-10-12
       if ((real(fit_function_poly).gt. poly_log_max(is))) then
       !if ((abs(fit_function_poly).lt.-100.d0).or.(abs(fit_function_poly).gt.100.d0)) then
       !if ((abs(fit_function_poly).lt.-50.d0).or.(abs(fit_function_poly).gt.50.d0)) then
       !if (abs(fit_function_poly).gt.poly_log_max(is)) then
          fit_function_poly=cmplx(0.d0,0.d0,kind(1.d0))
       else
          !write(*,*)iproc,is,iperp,ppar_val,ppar_val_tmp,fit_function_poly
          fit_function_poly=10.d0**(fit_function_poly)
          !write(*,*)iproc,is,iperp,ppar_val,ppar_val_tmp,fit_function_poly
       endif
    endif
	!Calculate the Actual Derivative
	if (present(dfpar)) then
   		if (logfit(is)) then
      		dfpar = log(10.d0)*fit_function_poly*poly_deriv/norm_2
   		else
      		dfpar = poly_deriv/norm_2
   		endif
	endif
 endif
end select

 return
end function fit_function_poly

subroutine determine_param_fit
	!! This is the fitting routine for the hybrid analytic continuation. It determines
	!! the full field [[alps_var(module):param_fit(variable)]].
	use alps_var, only : writeOut, fit_type, param_fit, n_fits, nspec, f0, nperp, npar, logfit, runname
        use alps_var, only : relativistic,npparbar,f0_rel,ngamma, perp_correction, gamma_rel, usebM
        use alps_var, only : ACmethod, poly_fit_coeffs 
	implicit none

	integer :: ifit
	!! Index for the fits that add up to the total fit for a given species.

	integer :: n_params
	!! Total number of fit parameters for a given species.

	integer :: par_ind
	!! Parameter index for the fit parameters.

	integer :: iperp
	!! Index of perpendicular momentum.

	integer :: is
	!! Index of species.

	integer :: is_run
	!! Running variable over species to determine relativistic evaluation.

	integer :: is_rel
	!! Index for relativistic species (if any).

	integer :: sproc_rel
	!! Local process index for relativistic calculation.

	integer :: nJT
	!! First dimension of matrix JT.

	integer :: ipparbar
	!! Index of parallel momentum (relativistic).

	integer :: ipparbar_lower
	!! Lower index of parallel momentum (relativistic).

	integer :: ipparbar_upper
	!! Upper index of parallel momentum (relativistic).

	integer :: upperlimit
	!! Upper limit of iperp space (relativistic and non-relativistic).

	integer :: unit_spec
	!! Unit to write fit parameters to file.

	logical :: found_lower
	!! Check whether lower boundary was found.

	logical :: found_upper
	!! Check whether upper boundary was found.

	double precision :: quality
	!! Quality of the individual fit result.

	double precision :: qualitytotal
	!! Quality of the total fit result.

	logical, allocatable, dimension (:) :: param_mask
	!! Bit mask for required fit parameters.

	double precision,allocatable,dimension (:) :: g
	!! Array of function to be fitted.

	double precision,allocatable,dimension(:) :: params
	!! Array of fit parameters.

	character (10) :: specwrite
	!! File name to write fit parameters to file.

	if (writeOut) then
		write (*,'(a)') 'Determine fit parameters for hybrid analytic continuation...'
	endif

	qualitytotal=0.d0

	do is=1,nspec

		if (usebM(is)) then
     write (*,'(a,i2)') ' Bi-Maxwellian/cold-plasma calculation: no fits necessary for species ',is
     param_fit(is,:,:,:)=0.d0
     fit_type(is,:)=1
     perp_correction(is,:)=1.d0

		 elseif (ACmethod(is).EQ.0) then

			 	write (*,'(a,i2)') ' Using analytical function from distribution/distribution_analyt.f90: no fits necessary for species ',is

				param_fit(is,:,:,:)=0.d0
				fit_type(is,:)=1
    perp_correction(is,:)=1.d0

    		 elseif (ACmethod(is).EQ.2) then

                                ! Relativistic rows were built once by initialize_rel_cheb.
                                if(relativistic(is)) cycle

			 	write (*,'(a,i2)') ' Using polynomial representation for species ',is

				param_fit(is,:,:,:)=0.d0
				fit_type(is,:)=1
                                perp_correction(is,:)=1.d0

                                call set_polynomial_basis(is)

                                call determine_GLLS(is)

                                unit_spec=2500+is
                                write(specwrite,'(i0)') is
                                open(unit = unit_spec,file = 'distribution/'&
                                     //trim(runname)//'.poly_parameters.'//trim(specwrite)&
                                     //'.out', status = 'replace')
                                do iperp=0,nperp
                                   write(unit_spec,*)&
                                        iperp,poly_fit_coeffs(is,iperp,:)
                                enddo

                                close(unit_spec)

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
				ipparbar_lower = 0
				ipparbar_upper = npparbar
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

endif !end if (bM or ACmethod) selection 
enddo	! End loop over is


if(writeOut) then
   call output_fit(qualitytotal)
   write(*,'(a)') '-=-=-=-=-=-=-=-=-'
endif


end subroutine determine_param_fit


subroutine determine_GLLS(is)
  !! This subroutine evaluates the General Linear Least Squares fit
  !! to the distribution function for component 'is' using the selected
  !! polynomial basis functions.
  use alps_var, only : polynomials, nperp, npar, f0, poly_fit_coeffs, poly_order
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
     call least_squares_fit(polynomials(is,:,:),f0_fit,poly_fit_coeffs(is,iperp,:),poly_order(is))
  enddo
  deallocate(f0_fit)
  
end subroutine determine_GLLS

subroutine least_squares_fit(AA,BB,coeffs,npoly)
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
  
end subroutine least_squares_fit

subroutine set_polynomial_basis(is)
  !! This subroutine evaluates the General Linear Least Squares fit
  !! to the distribution function for component 'is' using the selected
  !! polynomial basis functions.
  use alps_var, only : polynomials, poly_kind, poly_order
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
  
end subroutine set_polynomial_basis

subroutine output_fit(qualitytotal)
	!! This subroutine outputs the fit parameters for iperp=0 to stdout to monitor the fit.
	use alps_io, only : isnancheck, alps_error
	use alps_var, only : fit_type, param_fit, n_fits, nspec, nperp, npar, pp, f0, pi, vA, runname
        use alps_var, only : relativistic,gamma_rel,pparbar_rel,ngamma,npparbar,f0_rel, ms, usebM
        use alps_var, only : ACmethod
        use alps_var, only : ns, qs, ms, bMpdrifts
        use alps_var, only : density_int, current_int, poly_order
	implicit none

	double precision, intent(in) :: qualitytotal
	!! Quality of the total fit result.

	integer :: is
	!! Index of species.

	integer :: ifit
	!! Index for the fits that add up to the total fit for a given species.

	integer :: n_params
	!! Total number of fit parameters for a given species.

	integer :: iparam
	!! Variable running over fit parameters.

	integer :: unit_spec
	!! Unit to write fit results to file.

	integer :: ipar
	!! Index of parallel momentum.

	integer :: iperp
	!! Index of perpendicular momentum.

	integer :: is_rel
	!! Index for relativistic species (if any).

	integer :: sproc_rel
	!! Local process index for relativistic calculation.

	integer :: igamma
	!! Index of gamma (relativistic).

	integer :: ipparbar
	!! Index of pparbar (relativistic).

	integer :: is_run
	!! Index for relativistic species (if any).

	double precision :: integrate
 !! Integration of fit result.

 	double precision, dimension(0:nspec) :: charge
  !! Charge Density of species.
  !! Zeroth index is sum over all species

   	double precision, dimension(0:nspec) :: current
  !! Current Density of species.
  !! Zeroth index is sum over all species.

	double precision :: dpperp
	!! Step size in pperp for integration of fit result.

	double precision :: dppar
	!! Step size in ppar for integration of fit result.

	double precision :: dgamma_rel
	!! Step size in gamma for integration of relativistic fit result.

	double precision :: dpparbar
	!! Step size in pparbar for integration of relativistic fit result.

	double complex :: ppar_comp
	!! Complex parallel momentum to evaluate fit function at.

	character (10) :: specwrite
	!! File name to write fit results to file.

	write (*,'(a)') ' Results of the fit for hybrid analytic continuation at iperp = 1:'

 do is=1,nspec
    select case (ACmethod(is))
    case (1)
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
    case (2)
       
    end select
 enddo


	write (*,'(a,es14.4e3)') ' Sum of all least-squares: ', qualitytotal
	write (*,'(a,es14.4e3)') ' Standard error of the estimate: ', sqrt(qualitytotal/(1.d0*(nspec*nperp*npar)))


	if (isnancheck(qualitytotal)) call alps_error(7)

 charge=0.d0
 current=0.d0

	write (*,'(a)') ' Writing fit result to files'
	do is=1,nspec

    if (usebM(is)) then
       charge(is)=ns(is)*qs(is)
       current(is)=ns(is)*qs(is)*&
            bMpdrifts(is)/ms(is)
       write(*,'(a)')&
            '-=-=-=-='
       write(*,'(a,i3,a)')&
            'Bi-Maxwellian/cold-plasma Species ',is,':'
       write(*,'(a, 2es14.4e3)') &
            ' Charge density:           ', charge(is)
       write(*,'(a, 2es14.4e3)') &
            ' Parallel current density: ', current(is)

    else
       select case (ACmethod(is))
	
	  case (0,1)
	
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

		dgamma_rel = gamma_rel(sproc_rel,2,2)-gamma_rel(sproc_rel,1,2)
		dpparbar = pparbar_rel(sproc_rel,2,2)-pparbar_rel(sproc_rel,2,1)

		do igamma=0,ngamma
		 do ipparbar=0,npparbar
			ppar_comp=cmplx(pparbar_rel(sproc_rel,igamma,ipparbar),0.d0,kind=kind(1.d0))
			if (f0_rel(sproc_rel,igamma,ipparbar).EQ.-1.d0) then
				write (unit_spec,*) gamma_rel(sproc_rel,igamma,ipparbar),pparbar_rel(sproc_rel,igamma,ipparbar),&
					"-1.0",	abs(-1.d0-f0_rel(sproc_rel,igamma,ipparbar))
			else

				write (unit_spec,*) gamma_rel(sproc_rel,igamma,ipparbar),pparbar_rel(sproc_rel,igamma,ipparbar),&
					real(eval_fit(is,igamma,ppar_comp)),&
						abs(real(eval_fit(is,igamma,ppar_comp))-f0_rel(sproc_rel,igamma,ipparbar))/f0_rel(sproc_rel,igamma,ipparbar)

						  integrate = integrate + &
                     gamma_rel(is_rel,igamma,ipparbar) * real(eval_fit(is,igamma,ppar_comp)) * &
                     2.d0 * pi * dgamma_rel * dpparbar * (ms(is) / vA)**3

        charge(is) = charge(is) + &
             qs(is)* ns(is)* gamma_rel(is_rel,igamma,ipparbar) * real(eval_fit(is,igamma,ppar_comp)) * &
             2.d0 * pi * dgamma_rel * dpparbar * (ms(is) / vA)**3

        current(is) = current(is) + &
             (real(ppar_comp)/ms(is))*qs(is)* ns(is)* real(eval_fit(is,igamma,ppar_comp)) * &
             2.d0 * pi * dgamma_rel * dpparbar * (ms(is) / vA)**3

			endif
		 enddo
		enddo

	  else 		! non-relativistic:

	    dpperp = pp(is,2,2,1)-pp(is,1,2,1)
	    dppar = pp(is,2,2,2)-pp(is,2,1,2)

		do iperp=0,nperp
			do ipar=0,npar
				ppar_comp=cmplx(pp(is,iperp,ipar,2),0.d0,kind=kind(1.d0))
				write (unit_spec,*) pp(is,iperp,ipar,1), pp(is,iperp,ipar,2), real(eval_fit(is,iperp,ppar_comp)), &
         abs(real(eval_fit(is,iperp,ppar_comp))-f0(is,iperp,ipar))/f0(is,iperp,ipar)

    integrate = integrate + pp(is,iperp,ipar,1) * real(eval_fit(is,iperp,ppar_comp)) * &
             2.d0 * pi * dpperp * dppar
    charge(is) = charge(is) + qs(is)* ns(is)* pp(is,iperp,ipar,1) * real(eval_fit(is,iperp,ppar_comp)) * &
         2.d0 * pi * dpperp * dppar
    current(is) = current(is) + (qs(is)* ns(is)/ms(is))*&
         pp(is,iperp,ipar,2) *pp(is,iperp,ipar,1) * real(eval_fit(is,iperp,ppar_comp)) * &
         2.d0 * pi * dpperp * dppar

			enddo
enddo

	  endif
	  close(unit_spec)

   write(*,'(a)')&
        '-=-=-=-='
   write(*,'(a,i3,a)')&
        'Species ',is,':'
   write(*,'(a, 2es14.4e3)') &
        ' Integration of fit/analytical function:              ', integrate
   write(*,'(a, 2es14.4e3)') &
        ' Charge density of fit/analytical function:           ', charge(is)
   write(*,'(a, 2es14.4e3)') &
        ' Parallel current density of fit/analytical function: ', current(is)

case (2)

   integrate = 0.d0
   dpperp = pp(is,2,2,1)-pp(is,1,2,1)
   dppar = pp(is,2,2,2)-pp(is,2,1,2)

   
     unit_spec=2000+is
     write(specwrite,'(i0)') is
     open(unit = unit_spec,file = 'distribution/'//trim(runname)//'.poly_fit_result.'//trim(specwrite)//'.out', status = 'replace')
     write (*,'(a,i0)') ' Writing GLLS fit result to files: ispec',is
     do iperp=0,nperp
        do ipar=0,npar
           ppar_comp=cmplx(pp(is,iperp,ipar,2),0.d0,kind=kind(1.d0))
           if (f0(is,iperp,ipar).gt.0.d0) then
              write (unit_spec,*) pp(is,iperp,ipar,1), pp(is,iperp,ipar,2), &
                   real(eval_fit(is,iperp,ppar_comp)), &
                   abs(real(eval_fit(is,iperp,ppar_comp))-f0(is,iperp,ipar))/f0(is,iperp,ipar)
           else
              write (unit_spec,*) pp(is,iperp,ipar,1), pp(is,iperp,ipar,2), &
                   real(eval_fit(is,iperp,ppar_comp)), &
                   0.d0
           endif

    integrate = integrate + pp(is,iperp,ipar,1) * real(eval_fit(is,iperp,ppar_comp)) * &
             2.d0 * pi * dpperp * dppar
    charge(is) = charge(is) + qs(is)* ns(is)* pp(is,iperp,ipar,1) * real(eval_fit(is,iperp,ppar_comp)) * &
         2.d0 * pi * dpperp * dppar
    current(is) = current(is) + (qs(is)* ns(is)/ms(is))*&
         pp(is,iperp,ipar,2) *pp(is,iperp,ipar,1) * real(eval_fit(is,iperp,ppar_comp)) * &
         2.d0 * pi * dpperp * dppar
           
        enddo
     enddo
    
     close(unit_spec)

        write(*,'(a)')&
        '-=-=-=-='
   write(*,'(a,i3,a)')&
        'Species ',is,':'
   write(*,'(a, 2es14.4)') &
        ' Integration of polynomial representation:              ', integrate
   write(*,'(a, i0, es14.4)') &
        ' Relative difference in density for Order:', poly_order(is), &
        2.d0*(abs(density_int(is))-abs(integrate))/(abs(density_int(is))+abs(integrate))
   write(*,'(a, 2es14.4)') &
        ' Charge density of polynomial representation:           ', charge(is)
   if (relativistic(is)) then
      write(*,'(a)')'Relativistic parallel current density not yet implemented!'
   else
      write(*,'(a, es14.4)') &
           ' Parallel current density of polynomial representation: ', current(is)
      write(*,'(a, i0, es14.4)') &
           ' Relative difference in current density for Order:', poly_order(is), &
           2.d0*(abs(current(is))-abs(current_int(is)))/(abs(current_int(is))+abs(current(is)))
   endif


   
end select
   
endif
enddo
write(*,'(a)')         '-=-=-=-='
write(*,'(a, es14.4e3)') ' Total charge density of fit/analytical function:           ', sum(charge(1:nspec))
write(*,'(a, es14.4e3)') ' Total parallel current density of fit/analytical function: ', sum(current(1:nspec))
write(*,'(a)')         '-=-=-=-=-=-=-=-='


end subroutine output_fit




subroutine determine_JT(is,n_params,nJT,JT,params,iperp,upper_limit,ipparbar_lower)
	!! This subroutine calculates the transposed Jacobian matrix of the fit function with respect to the
	!! fit parameter array.
	use alps_var, only : fit_type, n_fits, pp, ms, vA, perp_correction
	use alps_var, only : gamma_rel,pparbar_rel,nspec,relativistic
	implicit none

	integer, intent(in) :: is
	!! Index of species for which [[determine_JT]] is executed.

	integer, intent(in) :: n_params
	!! Total number of fit parameters for a given species.

	integer, intent(in) :: nJT
	!! First dimension of matrix JT.

	integer, intent(in) :: upper_limit
	!! Upper limit of iperp space (relativistic and non-relativistic).

	double precision, intent(out) :: JT(nJT,0:upper_limit)
	!! Transposed Jacobian matrix of the fit function.

	double precision, intent(in) :: params(n_params)
	!! Array of fit parameters.

	integer, intent(in) :: iperp
	!! Index of perpendicular momentum at which JT is evaluated.

	integer, intent(in) :: ipparbar_lower
	!! Lower index of parallel momentum (relativistic).

	integer :: ifit
	!! Index for the fits that add up to the total fit for a given species.

	integer :: par_ind
	!! Parameter index to loop over fits.

	integer :: ipar
	!! Index of parallel momentum.

	integer :: is_rel
	!! Index for relativistic species (if any).

	integer :: sproc_rel
	!! Local process index for relativistic calculation.

	integer :: is_run
	!! Index for relativistic species (if any).

	integer :: JT_ind
	!! Index running over [[JT]]

	double precision :: ppar_val
	!! Parallel momentum.

	double precision :: pperp_val
	!! Perpendicular momentum.

	double precision :: sqrtpart
	!! Square-root part of Juettner distribution.

	double precision :: expterm
	!! Exponential part of the Maxwellian distribution.

	double precision :: kappapart
	!! Kappa part of the kappa distribution.

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





subroutine LM_nonlinear_fit(is,g,n_params,nJT,params,param_mask,iperp,npar,ipparbar_lower,quality)
	!! This subroutine processes the nonlinear Levenberg-Marquart algorithm and returns
	!! the one-dimensional array params at a given iperp.
	!! The variable quality is the sum of the squares of all residuals.
	use alps_var, only : lambda_initial_fit, pp, lambdafac_fit, logfit
	use alps_var, only : epsilon_fit, maxsteps_fit
	use alps_var, only : gamma_rel,pparbar_rel,relativistic,nspec
	use alps_io, only : alps_error
	use mpi
	implicit none

	integer, intent(in) :: is
	!! Index of species for which [[LM_nonlinear_fit]] is executed.

	double precision, intent(in) :: g(0:npar)
	!! Array of function to be fitted.

	integer, intent(in) :: n_params
	!! Total number of fit parameters for a given species.

	integer, intent(in) :: nJT
	!! First dimension of matrix JT (see [[determine_JT(subroutine)]]).

	double precision, intent(inout) :: params(n_params)
	!! Array of fit parameters.

	logical, intent(in) :: param_mask(n_params)
	!! Bit mask for required fit parameters.

	integer, intent(in) :: iperp
	!! Index of perpendicular momentum.

	integer, intent(in) :: npar
	!! Number of steps in parallel momentum

	integer, intent(in) :: ipparbar_lower
	!! Lower index of parallel momentum (relativistic).

	double precision, intent(out) :: quality
	!! Quality of the individual fit result.

	logical :: converged
	!! Check whether fit has converged.

	integer :: ipar
	!! Index running over parallel momentum.

	integer :: k
	!! Index running over entries of [[JT]].

	integer :: counter
	!! Count of fit iterations.

	integer :: is_rel
	!! Index for relativistic species (if any).

	integer :: is_run
	!! Index for relativistic species (if any).

	integer :: sproc_rel
	!! Local process index for relativistic calculation.

	integer :: l
	!! Index running over the delta's in the L-M fit.

	integer :: ipiv(nJT)
	!! Integer array required for matrix inversion.

	integer :: info
	!! Integer required for matrix inversion.

	double precision :: LSQ
	!! Least squares for L-M fit.

	double precision :: LSQnew
	!! Next iteration of least squares for L-M fit.

	double precision :: lambda_fit
	!! Lambda in L-M fit.

	double precision :: pperp_val
	!! Perpendicular momentum.

	double precision :: ppar_val
	!! Parallel momentum.

	double precision :: residuals(0:npar)
	!! Array of residuals.

	double precision :: deltaparam_fit(nJT)
	!! delta's in the L-M fit.

	double precision :: JTJ(nJT,nJT)
	!! Matrix product of JT and J.

	double precision :: JT(nJT,0:npar)
	!! Transposed Jacobian matrix of the fit function.

	double precision :: diagmat(nJT,nJT)
	!! Diagonal matrix of JT.

	double precision :: Amat(nJT,nJT)
	!! Matrix to be inverted.

	double precision :: work_array(nJT)
	!! Work array for matrix inversion.

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

end subroutine LM_nonlinear_fit

  integer function factorial(n)
    implicit none
    integer, intent(in) :: n
    integer :: i, Ans

    Ans=1
    do i=1,n
       Ans=Ans*i
    enddo
    Factorial=Ans
  end function factorial

end module alps_analyt
