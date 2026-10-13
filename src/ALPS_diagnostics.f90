! Optional observations of the existing MAP calculation. No solver decisions live here.
module alps_diagnostics
  use, intrinsic :: ieee_arithmetic
  use alps_diagnostic_math
  use alps_var, only: nspec, nperp, pp, ns, qs, poly_order, ACmethod, logfit, &
       usebM, relativistic, proc0, iproc, nproc, runname, kperp, kpar, vA, kperp_norm, &
       npar, ms, pi, poly_kind, poly_fit_coeffs, poly_log_max
  use mpi
  implicit none
  private
  public :: diagnose_analytic_continuation, diagnostic_prefix, diagnostic_species
  public :: diagnostic_all_harmonics, diagnostic_harmonic, diagnostic_verbosity
  public :: diagnostic_iperp_first, diagnostic_iperp_last, diagnostic_iperp_stride
  public :: diagnostic_eval_first, diagnostic_eval_last, diagnostic_eval_stride
  public :: diag_active, diag_begin, diag_finish, diag_close, diag_configure
  public :: diag_context, diag_integral, diag_landau, diag_full, diag_ee
  public :: diag_point, diag_summand, diag_polynomial, diag_capture_matrix
  public :: diag_logabs, diag_abs, diag_finite
  public :: diag_compare_residue, diagnostic_mode

  logical :: diagnose_analytic_continuation=.false., diag_active=.false.
  logical :: diagnostic_all_harmonics=.true.
  character(500) :: diagnostic_prefix=''
  integer :: diagnostic_species=0, diagnostic_harmonic=0, diagnostic_verbosity=1
  integer :: diagnostic_iperp_first=0, diagnostic_iperp_last=-1, diagnostic_iperp_stride=1
  integer :: diagnostic_eval_first=1, diagnostic_eval_last=-1, diagnostic_eval_stride=1
  integer :: diagnostic_mode=0
  integer, parameter :: width=96, nfiles=10
  integer, parameter :: ci(7)=[1,2,3,1,1,2,3], cj(7)=[1,2,3,2,3,3,3]

  type :: table_buffer
     double precision, allocatable :: rows(:,:)
     integer :: count=0, columns=0, unit=0
  end type
  type :: harmonic_summary
     double complex :: integral=(0.d0,0.d0), residue=(0.d0,0.d0)
     double precision :: maximum=0.d0
     integer :: maximum_index=-1, limiter_count=0, discarded=0, failed=0
  end type
  type :: frequency_summary
     double precision :: extrema(5)
     double precision :: location(9,5)
     logical :: extrema_seen(5)=.false.
     integer :: limiter_count=0, invalid_count=0, calls=0
     double precision :: first(9)
  end type
  type :: residue_point
     double complex :: parallel, perpendicular, summand
     logical :: visited=.false.
  end type
  type(residue_point), allocatable, save :: comparison_points(:)
  logical, save :: comparison_ready=.false.
  double complex, save :: alternative_sums(4), reconstructed_residue, pre_discard
  integer, save :: comparison_cutoffs=0, comparison_missing=0
  type(table_buffer), save :: tables(nfiles)
  type(harmonic_summary), save :: harmonic
  type(frequency_summary), save :: frequency
  integer, save :: eval_id=0, map_ir=0, map_ii=0, species=0, harmonic_n=0, tensor_mode=0
  integer, save :: residue_index=-1
  logical, save :: opened=.false., ee_seen=.false., first_written=.false.
  double complex, save :: omega, resonance=(0.d0,0.d0)
  double complex, allocatable, save :: reconstructed(:,:,:), actual_chi(:,:,:)
  double precision, allocatable, save :: absolute_sum(:,:,:)
  double complex, save :: actual_wave(3,3)

contains

  logical function diag_finite(z)
    double complex, intent(in) :: z
    diag_finite=ieee_is_finite(real(z)).and.ieee_is_finite(aimag(z))
  end function

  double precision function nan_value()
    nan_value=ieee_value(0.d0,ieee_quiet_nan)
  end function

  double precision function diag_logabs(z) result(value)
    double complex, intent(in) :: z
    double precision :: a,b
    if (ieee_is_nan(real(z)).or.ieee_is_nan(aimag(z))) then
       value=nan_value()
    elseif (.not.diag_finite(z)) then
       value=ieee_value(0.d0,ieee_positive_inf)
    else
       a=max(abs(real(z)),abs(aimag(z)))
       b=min(abs(real(z)),abs(aimag(z)))
       if (a==0.d0) then
          value=ieee_value(0.d0,ieee_negative_inf)
       else
          value=log10(a)+0.5d0*log10(1.d0+(b/a)**2)
       endif
    endif
  end function

  double precision function diag_abs(z) result(value)
    double complex, intent(in) :: z
    double precision :: a,b,factor
    if (.not.diag_finite(z)) then
       value=diag_logabs(z)
       return
    endif
    a=max(abs(real(z)),abs(aimag(z)))
    b=min(abs(real(z)),abs(aimag(z)))
    value=0.d0
    if (a==0.d0) return
    factor=sqrt(1.d0+(b/a)**2)
    if (a>huge(a)/factor) then
       value=ieee_value(0.d0,ieee_positive_inf)
    else
       value=a*factor
    endif
  end function

  double precision function safe_sum(a,b) result(value)
    double precision, intent(in) :: a,b
    ! All callers supply nonnegative magnitudes. An unrepresentable sum is
    ! reported as infinity rather than raising a diagnostic-only overflow.
    if (ieee_is_nan(a).or.ieee_is_nan(b)) then
       value=nan_value()
    elseif (a>huge(a)-b) then
       value=ieee_value(0.d0,ieee_positive_inf)
    else
       value=a+b
    endif
  end function

  double precision function signed_sum(a,b) result(value)
    double precision, intent(in) :: a,b
    value=nan_value()
    if (ieee_is_nan(a).or.ieee_is_nan(b)) return
    if (.not.ieee_is_finite(a).or..not.ieee_is_finite(b)) then
       if (.not.ieee_is_finite(a).and..not.ieee_is_finite(b)) then
          if (sign(1.d0,a)/=sign(1.d0,b)) return
       endif
       value=a
       if (ieee_is_finite(a)) value=b
    elseif (b>0.d0) then
       if (a>huge(a)-b) then
          value=ieee_value(0.d0,ieee_positive_inf)
       else
          value=a+b
       endif
    elseif (b<0.d0) then
       if (a<-huge(a)-b) then
          value=ieee_value(0.d0,ieee_negative_inf)
       else
          value=a+b
       endif
    else
       value=a
    endif
  end function

  double complex function complex_sum(a,b) result(value)
    double complex, intent(in) :: a,b
    value=cmplx(signed_sum(real(a),real(b)),signed_sum(aimag(a),aimag(b)),kind(1.d0))
  end function

  double complex function normalized(value) result(result_value)
    double complex, intent(in) :: value
    double complex :: normalization
    normalization=omega*omega*vA*vA
    if (.not.kperp_norm) normalization=normalization*kperp*kperp
    result_value=cmplx(nan_value(),nan_value(),kind(1.d0))
    if (.not.diag_finite(value).or..not.diag_finite(normalization)) return
    if (diag_abs(normalization)==0.d0) return
    if (diag_logabs(value)-diag_logabs(normalization)>log10(huge(1.d0))) return
    result_value=value/normalization
  end function

  subroutine diag_configure
    integer :: ierr
    ! Called by every rank after the normal input broadcasts.
    call MPI_BCAST(diagnose_analytic_continuation,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (.not.diagnose_analytic_continuation) return
    call MPI_BCAST(diagnostic_prefix,len(diagnostic_prefix),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(diagnostic_all_harmonics,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(diagnostic_species,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(diagnostic_harmonic,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(diagnostic_verbosity,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(diagnostic_mode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(diagnostic_iperp_first,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(diagnostic_iperp_last,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(diagnostic_iperp_stride,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(diagnostic_eval_first,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(diagnostic_eval_last,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(diagnostic_eval_stride,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if (diagnostic_iperp_stride<1.or.diagnostic_eval_stride<1.or. &
         diagnostic_verbosity<1.or.diagnostic_verbosity>2.or.diagnostic_species<0.or. &
         diagnostic_species>nspec.or.diagnostic_iperp_first<0.or.diagnostic_eval_first<1.or. &
         diagnostic_mode<0.or.diagnostic_mode>6) then
       if (proc0) write(*,*) 'Invalid analytic-continuation diagnostic configuration'
       call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    endif
  end subroutine

  logical function selected(s,n)
    integer, intent(in) :: s,n
    selected=(diagnostic_species==0.or.s==diagnostic_species).and. &
         (diagnostic_all_harmonics.or.n==diagnostic_harmonic)
  end function

  logical function detailed(s,n,p)
    integer, intent(in) :: s,n,p
    detailed=.false.
    if (diagnostic_verbosity/=2.or..not.selected(s,n)) return
    if (diagnostic_mode/=0.and.tensor_mode/=diagnostic_mode) return
    if (eval_id<diagnostic_eval_first.or.p<diagnostic_iperp_first) return
    if (diagnostic_eval_last>=0.and.eval_id>diagnostic_eval_last) return
    if (diagnostic_iperp_last>=0.and.p>diagnostic_iperp_last) return
    detailed=mod(eval_id-diagnostic_eval_first,diagnostic_eval_stride)==0.and. &
         mod(p-diagnostic_iperp_first,diagnostic_iperp_stride)==0
  end function

  subroutine append_row(file,row)
    integer, intent(in) :: file
    double precision, intent(in) :: row(:)
    double precision, allocatable :: larger(:,:)
    integer :: capacity
    if (.not.allocated(tables(file)%rows)) allocate(tables(file)%rows(width,64))
    capacity=size(tables(file)%rows,2)
    if (tables(file)%count==capacity) then
       allocate(larger(width,2*capacity))
       larger(:,1:capacity)=tables(file)%rows
       call move_alloc(larger,tables(file)%rows)
    endif
    tables(file)%count=tables(file)%count+1
    tables(file)%rows(:,tables(file)%count)=0.d0
    tables(file)%rows(1:size(row),tables(file)%count)=row
  end subroutine

  subroutine open_table(file,suffix,columns,header)
    integer, intent(in) :: file,columns
    character(*), intent(in) :: suffix,header
    integer :: ios,ierr
    tables(file)%columns=columns
    if (.not.proc0) return
    open(newunit=tables(file)%unit,file=trim(diagnostic_prefix)//suffix, &
         status='replace',action='write',iostat=ios)
    if (ios/=0) then
       write(*,*) 'Cannot open diagnostic file: ',trim(diagnostic_prefix)//suffix
       call MPI_ABORT(MPI_COMM_WORLD,ios,ierr)
    endif
    write(tables(file)%unit,'(a)') '# ALPS analytic-continuation diagnostics v1'
    write(tables(file)%unit,'(a)') '# omega: Omega_ref; k: 1/d_ref; p: m_ref*vA_ref; P: log10(f)'
    write(tables(file)%unit,'(a)') '# '//header
  end subroutine

  subroutine diag_begin(om,ir,ii)
    double complex, intent(in) :: om
    integer, intent(in) :: ir,ii
    integer :: f
    if (.not.diagnose_analytic_continuation) return
    if (.not.opened) then
       if (proc0.and.len_trim(diagnostic_prefix)==0) diagnostic_prefix='solution/'//trim(runname)//'_analytic'
       call open_table(1,'_determinant.dat',17, &
            'eval ir ii wr wi kperp kpar order det_re det_im logabsD maxD minD finite_det finite_matrix limiters invalid_f')
       call open_table(2,'_matrix.dat',12, &
            'eval kind species i j re im conventional_chi_re conventional_chi_im reconstruction_error cancellation coverage')
       call open_table(3,'_harmonics.dat',25, &
            'eval species n mode order integral_re integral_im residue_re residue_im total_re total_im '// &
            'chi_integral_re chi_integral_im chi_residue_re chi_residue_im chi_total_re chi_total_im '// &
            'max_weighted_integrand iperp_max integrated_residue_abs limiters discarded failed coverage rank')
       call open_table(4,'_glls_extrema.dat',13, &
            'eval metric maximum species n mode iperp ppar_re ppar_im role rank call coverage')
       call open_table(5,'_limiter_events.dat',14, &
            'eval rank count calls invalid_f first_species first_n first_mode first_iperp '// &
            'first_ppar_re first_ppar_im first_role first_call first_scan_crossing')
       if (diagnostic_verbosity==2) then
          call open_table(6,'_residue_integrand.dat',20, &
               'eval species n mode iperp pperp pres_re pres_im h summand_re summand_im '// &
               'dfpar_re dfpar_im dfperp_re dfperp_im logabs_dfpar finite_df finite_summand rank multiplier')
          call open_table(7,'_glls_evaluations.dat',26, &
               'eval species n mode iperp pperp ppar_re ppar_im pres_re pres_im role '// &
               'P_re P_im logabs_f logabs_Pprime logabs_df_est limiter representable underflow '// &
               'finite_f domain_suppressed derivative_available rank call z_re z_im')
          call open_table(8,'_derivative_comparison.dat',32, &
               'eval species n mode iperp pperp pres_re pres_im observed factor h '// &
               'analytic_re analytic_im fd_re fd_im abs_error mixed_error magnitude_ratio phase_delta '// &
               'roundoff_scale convergence_ratio center_flags plus_flags minus_flags contaminated '// &
               'logabs_f logabs_analytic analytic_available resolved actual_dfpar_re actual_dfpar_im rank')
          call open_table(9,'_bracket_comparison.dat',62, &
               'eval species n mode iperp pperp pres_re pres_im observed h center_flags '// &
               'dfpar_re dfpar_im dfperp_re dfperp_im analytic_re analytic_im '// &
               'original_fd_re original_fd_im simplified_fd_re simplified_fd_im '// &
               'original_analytic_re original_analytic_im simplified_analytic_re simplified_analytic_im '// &
               'difference_re difference_im abs_difference normalized_difference '// &
               'parallel_re parallel_im negative_perp_re negative_perp_im omega_perp_re omega_perp_im '// &
               'perp_sum_re perp_sum_im resonance_residual_re resonance_residual_im '// &
               'intended_perp_re intended_perp_im residual_to_parallel tensor_re tensor_im '// &
               'actual_summand_re actual_summand_im weighted_original_fd_re weighted_original_fd_im '// &
               'weighted_simplified_fd_re weighted_simplified_fd_im weighted_original_analytic_re '// &
               'weighted_original_analytic_im weighted_simplified_analytic_re weighted_simplified_analytic_im rank finite '// &
               'abs_parallel abs_negative_perp abs_omega_perp abs_perp_sum abs_resonance_residual abs_intended_perp')
          call open_table(10,'_residue_alternatives.dat',21, &
               'eval species n mode variant value_re value_im actual_re actual_im abs_difference '// &
               'normalized_difference pre_discard_re pre_discard_im discarded cutoff_points '// &
               'missing_analytic_points finite multiplier species_factor rank order')
       endif
       allocate(reconstructed(nspec,3,3),actual_chi(nspec,3,3),absolute_sum(nspec,3,3))
       opened=.true.
    endif
    do f=1,nfiles
       tables(f)%count=0
    enddo
    eval_id=eval_id+1
    map_ir=ir; map_ii=ii; omega=om
    frequency%extrema=ieee_value(0.d0,ieee_negative_inf); frequency%location=0.d0
    frequency%extrema_seen=.false.
    frequency%limiter_count=0; frequency%invalid_count=0; frequency%calls=0; frequency%first=0.d0
    reconstructed=(0.d0,0.d0); absolute_sum=0.d0; ee_seen=.false.
    diag_active=.true.
  end subroutine

  subroutine diag_context(s,n,mode)
    integer, intent(in) :: s,n,mode
    if (.not.diag_active) return
    species=s; harmonic_n=n; tensor_mode=mode
    harmonic=harmonic_summary()
    comparison_ready=.false.
    if (allocated(comparison_points)) comparison_points%visited=.false.
    residue_index=-1; resonance=(0.d0,0.d0)
  end subroutine

  subroutine diag_point(p,pres)
    integer, intent(in) :: p
    double complex, intent(in) :: pres
    if (.not.diag_active) return
    residue_index=p; resonance=pres
  end subroutine

  subroutine diag_integral(value)
    double complex, intent(in) :: value
    if (diag_active) harmonic%integral=value
  end subroutine

  subroutine diag_landau(value,discarded)
    double complex, intent(in) :: value
    logical, intent(in), optional :: discarded
    if (.not.diag_active) return
    harmonic%residue=value
    if (present(discarded)) then
       if (discarded) harmonic%discarded=1
    endif
  end subroutine

  subroutine diag_summand(p,h,term,dfpar,dfperp,pres)
    integer, intent(in) :: p
    double precision, intent(in) :: h
    double complex, intent(in) :: term,dfpar,dfperp,pres
    double precision :: magnitude,multiplier
    if (.not.diag_active) return
    if (comparison_selected()) then
       if (.not.allocated(comparison_points)) allocate(comparison_points(0:nperp))
       comparison_points(p)%parallel=dfpar
       comparison_points(p)%perpendicular=dfperp
       comparison_points(p)%summand=term
       comparison_points(p)%visited=.true.
    endif
    magnitude=diag_abs(term)
    if (.not.diag_finite(term).or..not.diag_finite(dfpar).or..not.diag_finite(dfperp)) harmonic%failed=1
    if (magnitude>harmonic%maximum.or.harmonic%maximum_index==-1) then
       harmonic%maximum=magnitude; harmonic%maximum_index=p
    endif
    multiplier=1.d0
    if (aimag(omega)<0.d0) multiplier=2.d0
    if (detailed(species,harmonic_n,p)) call append_row(6,[dble(eval_id),dble(species), &
         dble(harmonic_n),dble(tensor_mode),dble(p),pp(species,p,1,1),real(pres),aimag(pres),h, &
         real(term),aimag(term),real(dfpar),aimag(dfpar),real(dfperp),aimag(dfperp),diag_logabs(dfpar), &
         dble(merge(1,0,diag_finite(dfpar).and.diag_finite(dfperp))),dble(merge(1,0,diag_finite(term))), &
         dble(iproc),multiplier])
  end subroutine

  logical function comparison_selected()
    comparison_selected=.false.
    if (.not.diag_active.or.diagnostic_verbosity/=2) return
    if (.not.selected(species,harmonic_n)) return
    if (diagnostic_mode/=0.and.tensor_mode/=diagnostic_mode) return
    if (eval_id<diagnostic_eval_first) return
    if (diagnostic_eval_last>=0.and.eval_id>diagnostic_eval_last) return
    if (mod(eval_id-diagnostic_eval_first,diagnostic_eval_stride)/=0) return
    if (.not.allocated(poly_fit_coeffs)) return
    if (ACmethod(species)/=2.or.poly_kind(species)/=1.or..not.logfit(species)) return
    if (usebM(species).or.relativistic(species)) return
    comparison_selected=.true.
  end function

  subroutine comparison_sample(p,arg,sample)
    integer, intent(in) :: p
    double complex, intent(in) :: arg
    type(fitted_sample), intent(out) :: sample
    double precision :: norm1,norm2
    norm1=0.5d0*(pp(species,p,npar,2)+pp(species,p,0,2))
    norm2=0.5d0*(pp(species,p,npar,2)-pp(species,p,0,2))
    call sample_log_chebyshev(poly_fit_coeffs(species,p,0:poly_order(species)), &
         (arg-norm1)/norm2,norm2,poly_log_max(species),sample)
  end subroutine

  subroutine diag_compare_residue(accumulated,discarded,tensor)
    ! This is a shadow experiment. Diagnostic samples never enter eval_fit,
    ! never change its event counts, and never supply the production result.
    double complex, intent(in) :: accumulated
    logical, intent(in) :: discarded
    interface
       double complex function tensor(n,p,pres,mode)
         integer, intent(in) :: n,p,mode
         double complex, intent(in) :: pres
       end function
    end interface
    double precision, parameter :: steps(6)=[2.d0,1.d0,0.5d0,0.25d0,0.125d0,0.0625d0]
    type(fitted_sample) :: center,plus,minus,upper,lower
    double complex :: pres,dfpar,dfperp,fd,b(4),terms(6),weighted(4),t,weight,difference
    double complex :: visited_sum,actual_term,normalization
    double precision :: dpar,dperp,h,metrics(5),previous,ratio,scale,denominator,relative,residual_ratio
    logical :: ok,observed,resolved,contaminated,previous_clean
    integer :: p,position,j,k,cutoff
    if (.not.comparison_selected()) return
    if (.not.allocated(comparison_points)) allocate(comparison_points(0:nperp))
    alternative_sums=(0.d0,0.d0); visited_sum=(0.d0,0.d0)
    comparison_cutoffs=0; comparison_missing=0
    dpar=abs(pp(species,2,2,2)-pp(species,2,1,2))
    dperp=pp(species,2,2,1)-pp(species,1,2,1)
    pres=(ms(species)*omega-dble(harmonic_n)*qs(species))/kpar
    ! Match the solver's traversal, including its existing nperp-1 weight.
    do position=1,nperp+1
       p=position
       if (position==nperp) p=0
       if (position==nperp+1) p=nperp
       h=1.d0
       if (p==0.or.p==nperp-1.or.p==nperp) h=0.5d0
       observed=comparison_points(p)%visited
       call comparison_sample(p,pres,center)
       call comparison_sample(p,pres+dpar,plus)
       call comparison_sample(p,pres-dpar,minus)
       dfpar=derivative_difference(plus%value,minus%value,dpar)
       call comparison_sample(min(p+1,nperp),pres,upper)
       call comparison_sample(max(p-1,0),pres,lower)
       if (p==0.or.p==nperp) then
          dfperp=derivative_difference(upper%value,lower%value,0.5d0*dperp)
       else
          dfperp=derivative_difference(upper%value,lower%value,dperp)
       endif
       actual_term=cmplx(nan_value(),nan_value(),kind(1.d0))
       if (observed) then
          dfpar=comparison_points(p)%parallel; dfperp=comparison_points(p)%perpendicular
          actual_term=comparison_points(p)%summand
       endif
       cutoff=ior(center%flags,ior(plus%flags,minus%flags))
       if (cutoff/=0) comparison_cutoffs=comparison_cutoffs+1
       if (.not.center%analytic_available) comparison_missing=comparison_missing+1
       previous=nan_value(); previous_clean=.false.
       do j=1,6
          call comparison_sample(p,pres+steps(j)*dpar,plus)
          call comparison_sample(p,pres-steps(j)*dpar,minus)
          fd=derivative_difference(plus%value,minus%value,steps(j)*dpar)
          if (j==2.and.observed) fd=dfpar
          scale=max(diag_abs(center%value),diag_abs(plus%value),diag_abs(minus%value))
          call difference_metrics(center%analytic,fd,scale,dpar,steps(j)*dpar,metrics,resolved)
          contaminated=ior(center%flags,ior(plus%flags,minus%flags))/=0
          ratio=nan_value()
          if (previous_clean.and..not.contaminated.and.center%analytic_available) then
             if (metrics(1)>metrics(5).and.previous>metrics(5)) ratio=safe_ratio(previous,metrics(1))
          endif
          previous=metrics(1); previous_clean=.not.contaminated
          if (detailed(species,harmonic_n,p)) call append_row(8,[dble(eval_id),dble(species), &
               dble(harmonic_n),dble(tensor_mode),dble(p),pp(species,p,1,1),real(pres),aimag(pres), &
               dble(merge(1,0,observed)),steps(j),steps(j)*dpar,real(center%analytic),aimag(center%analytic), &
               real(fd),aimag(fd),metrics,ratio,dble(center%flags),dble(plus%flags),dble(minus%flags), &
               dble(merge(1,0,contaminated)),center%log_function,center%log_derivative, &
               dble(merge(1,0,center%analytic_available)),dble(merge(1,0,resolved)),real(dfpar),aimag(dfpar),dble(iproc)])
       enddo
       call compare_brackets(omega,kpar,ms(species),qs(species),harmonic_n,pp(species,p,1,1), &
            pres,dfpar,dfperp,center%analytic,b,terms,ok)
       t=cmplx(nan_value(),nan_value(),kind(1.d0))
       ! An early discard can leave this tensor unevaluated in production.
       ! Do not introduce a diagnostic-only square of an enormous resonance.
       if (diag_abs(pres)<huge(1.d0)**0.125d0) t=tensor(harmonic_n,p,pres,tensor_mode)
       ! Match the production left-to-right weighting, including its minus sign.
       weight=guarded_product(cmplx(h,0.d0,kind(1.d0)),t)
       weight=guarded_product(weight,cmplx(qs(species)/abs(kpar),0.d0,kind(1.d0)))
       do k=1,4
          weighted(k)=-guarded_product(weight,b(k))
          alternative_sums(k)=complex_sum(alternative_sums(k),weighted(k))
       enddo
       if (observed) visited_sum=complex_sum(visited_sum,weighted(1))
       difference=complex_sum(b(1),-b(2))
       denominator=max(diag_abs(b(2)),epsilon(1.d0)* &
            maxval([diag_abs(terms(1)),diag_abs(terms(2)),diag_abs(terms(3))]),tiny(1.d0))
       relative=safe_ratio(diag_abs(difference),denominator)
       residual_ratio=safe_ratio(diag_abs(terms(5)),max(diag_abs(terms(1)),tiny(1.d0)))
       if (detailed(species,harmonic_n,p)) call append_row(9,[dble(eval_id),dble(species), &
            dble(harmonic_n),dble(tensor_mode),dble(p),pp(species,p,1,1),real(pres),aimag(pres), &
            dble(merge(1,0,observed)),h,dble(center%flags),real(dfpar),aimag(dfpar), &
            real(dfperp),aimag(dfperp),real(center%analytic),aimag(center%analytic), &
            real(b(1)),aimag(b(1)),real(b(2)),aimag(b(2)),real(b(3)),aimag(b(3)),real(b(4)),aimag(b(4)), &
            real(difference),aimag(difference),diag_abs(difference),relative, &
            real(terms(1)),aimag(terms(1)),real(terms(2)),aimag(terms(2)),real(terms(3)),aimag(terms(3)), &
            real(terms(4)),aimag(terms(4)),real(terms(5)),aimag(terms(5)),real(terms(6)),aimag(terms(6)), &
            residual_ratio,real(t),aimag(t),real(actual_term),aimag(actual_term), &
            real(weighted(1)),aimag(weighted(1)),real(weighted(2)),aimag(weighted(2)), &
            real(weighted(3)),aimag(weighted(3)),real(weighted(4)),aimag(weighted(4)), &
            dble(iproc),dble(merge(1,0,ok)),[(diag_abs(terms(k)),k=1,6)]])
    enddo
    normalization=cmplx(0.d0,dperp*pi*2.d0*pi,kind(1.d0))
    do k=1,4
       alternative_sums(k)=guarded_product(alternative_sums(k),normalization)
    enddo
    pre_discard=guarded_product(accumulated,normalization)
    reconstructed_residue=guarded_product(visited_sum,normalization)
    if (discarded) reconstructed_residue=(0.d0,0.d0)
    comparison_ready=.true.
  end subroutine

  subroutine write_residue_alternatives(actual,multiplier,factor)
    double complex, intent(in) :: actual
    double precision, intent(in) :: multiplier,factor
    double complex :: value,difference,pre
    double precision :: denominator
    integer :: variant
    pre=guarded_product(pre_discard,cmplx(multiplier*factor,0.d0,kind(1.d0)))
    do variant=0,5
       if (variant==0) then
          value=actual
       elseif (variant==1) then
          value=guarded_product(reconstructed_residue,cmplx(multiplier*factor,0.d0,kind(1.d0)))
       else
          value=guarded_product(alternative_sums(variant-1),cmplx(multiplier*factor,0.d0,kind(1.d0)))
       endif
       difference=complex_sum(value,-actual)
       denominator=max(diag_abs(actual),epsilon(1.d0)*diag_abs(value),tiny(1.d0))
       call append_row(10,[dble(eval_id),dble(species),dble(harmonic_n),dble(tensor_mode),dble(variant), &
            real(value),aimag(value),real(actual),aimag(actual),diag_abs(difference), &
            safe_ratio(diag_abs(difference),denominator),real(pre),aimag(pre),dble(harmonic%discarded), &
            dble(comparison_cutoffs),dble(comparison_missing),dble(merge(1,0,diag_finite(value))), &
            multiplier,factor,dble(iproc),dble(poly_order(species))])
    enddo
  end subroutine

  subroutine diag_full(value,found_res)
    double complex, intent(in) :: value
    logical, intent(in) :: found_res
    double precision :: factor,multiplier
    double complex :: integral,residue,total
    integer :: i,j,coverage
    if (.not.diag_active) return
    factor=ns(species)*qs(species)
    multiplier=0.d0
    if (found_res.and.aimag(omega)<=0.d0) multiplier=1.d0
    if (found_res.and.aimag(omega)<0.d0) multiplier=2.d0
    coverage=merge(0,1,relativistic(species))
    integral=harmonic%integral*factor
    residue=harmonic%residue*multiplier*factor
    if (comparison_ready) call write_residue_alternatives(residue,multiplier,factor)
    total=value*factor
    i=ci(tensor_mode); j=cj(tensor_mode)
    reconstructed(species,i,j)=complex_sum(reconstructed(species,i,j),total)
    absolute_sum(species,i,j)=safe_sum(absolute_sum(species,i,j),diag_abs(total))
    if (.not.diag_finite(total)) harmonic%failed=1
    if (coverage==0) then
       integral=cmplx(nan_value(),nan_value(),kind(1.d0)); residue=integral
    endif
    if (selected(species,harmonic_n)) call harmonic_row(integral,residue,total,coverage)
  end subroutine

  subroutine harmonic_row(integral,residue,total,coverage)
    double complex, intent(in) :: integral,residue,total
    integer, intent(in) :: coverage
    double complex :: conventional_integral,conventional_residue,conventional_total
    conventional_integral=normalized(integral)
    conventional_residue=normalized(residue)
    conventional_total=normalized(total)
    call append_row(3,[dble(eval_id),dble(species),dble(harmonic_n),dble(tensor_mode), &
         dble(poly_order(species)),real(integral),aimag(integral),real(residue),aimag(residue), &
         real(total),aimag(total),real(conventional_integral),aimag(conventional_integral), &
         real(conventional_residue),aimag(conventional_residue), &
         real(conventional_total),aimag(conventional_total),harmonic%maximum, &
         dble(harmonic%maximum_index),diag_abs(residue),dble(harmonic%limiter_count), &
         dble(harmonic%discarded),dble(harmonic%failed),dble(coverage),dble(iproc)])
  end subroutine

  subroutine diag_ee(s,value)
    integer, intent(in) :: s
    double complex, intent(in) :: value
    double complex :: total
    if (.not.diag_active.or.ee_seen) return
    ee_seen=.true.
    total=value*ns(s)*qs(s)
    if (.not.kperp_norm) total=kperp*kperp*total
    reconstructed(s,3,3)=complex_sum(reconstructed(s,3,3),total)
    absolute_sum(s,3,3)=safe_sum(absolute_sum(s,3,3),diag_abs(total))
    species=s; harmonic_n=0; tensor_mode=7; harmonic=harmonic_summary()
    if (selected(s,0)) call harmonic_row(total,cmplx(0.d0,0.d0,kind(1.d0)),total,1)
  end subroutine

  subroutine diag_polynomial(s,p,arg,z,polynomial,coeffs,norm2,value,outside)
    integer, intent(in) :: s,p
    double complex, intent(in) :: arg,z,polynomial,value
    double precision, intent(in) :: coeffs(0:),norm2
    logical, intent(in) :: outside
    double complex :: t0,t1,t2,d0,d1,d2,derivative
    double precision :: logs(5),location(9),lp,ld,cap,representable,underflow,role
    logical :: limited,available
    integer :: n,m
    if (.not.diag_active) return
    frequency%calls=frequency%calls+1
    role=0.d0
    if (p>residue_index) role=3.d0
    if (p<residue_index) role=4.d0
    if (p==residue_index.and.real(arg)>real(resonance)) role=1.d0
    if (p==residue_index.and.real(arg)<real(resonance)) role=2.d0
    limited=.false.
    ! The caller supplies the original P and the actual returned value.
    if (.not.outside.and.logfit(s)) limited=real(polynomial)>diagnostic_limit(s)
    if (limited) then
       harmonic%limiter_count=harmonic%limiter_count+1
       frequency%limiter_count=frequency%limiter_count+1
    endif
    if (.not.diag_finite(value)) frequency%invalid_count=frequency%invalid_count+1
    location=[dble(s),dble(harmonic_n),dble(tensor_mode),dble(p),real(arg),aimag(arg),role, &
         dble(iproc),dble(frequency%calls)]
    if (limited.and.frequency%limiter_count==1) frequency%first=location
    lp=nan_value(); ld=nan_value(); representable=-1.d0; underflow=-1.d0
    available=.not.outside.and.diag_finite(polynomial).and.norm2/=0.d0
    cap=sqrt(huge(1.d0))/100.d0
    derivative=(0.d0,0.d0); t0=(1.d0,0.d0); t1=z; d0=(0.d0,0.d0); d1=(1.d0,0.d0)
    if (available) then
       if (size(coeffs)>1) derivative=cmplx(coeffs(1),0.d0,kind(1.d0))
       do n=2,size(coeffs)-1
          if (max(diag_abs(t0),diag_abs(t1),diag_abs(d0),diag_abs(d1), &
               diag_abs(derivative),abs(coeffs(n)),diag_abs(z))>cap) then
             available=.false.; exit
          endif
          t2=2.d0*z*t1-t0
          d2=2.d0*t1+2.d0*z*d1-d0
          derivative=derivative+coeffs(n)*d2
          t0=t1; t1=t2; d0=d1; d1=d2
       enddo
       if (available) then
          ! Logarithms avoid division overflow for a very small momentum range.
          lp=diag_logabs(derivative)-log10(abs(norm2))
          if (logfit(s)) ld=real(polynomial)+log10(log(10.d0))+lp
       endif
    endif
    if (.not.outside.and.logfit(s)) then
       representable=dble(merge(1,0,ieee_is_finite(real(polynomial)).and. &
            real(polynomial)<=log10(huge(1.d0))))
       underflow=dble(merge(1,0,real(polynomial)<log10(tiny(1.d0))))
       logs=[real(polynomial),abs(aimag(polynomial)),real(polynomial),ld,lp]
       do m=1,5
          if (.not.ieee_is_nan(logs(m))) then
             if (.not.frequency%extrema_seen(m).or.logs(m)>frequency%extrema(m)) then
                frequency%extrema(m)=logs(m); frequency%location(:,m)=location
                frequency%extrema_seen(m)=.true.
             endif
          endif
       enddo
    endif
    if (detailed(s,harmonic_n,p)) call append_row(7,[dble(eval_id),dble(s),dble(harmonic_n), &
         dble(tensor_mode),dble(p),pp(s,p,1,1),real(arg),aimag(arg),real(resonance),aimag(resonance),role, &
         merge(real(polynomial),nan_value(),.not.outside),merge(aimag(polynomial),nan_value(),.not.outside), &
         merge(real(polynomial),nan_value(),.not.outside.and.logfit(s)),lp,ld, &
         dble(merge(1,0,limited)),representable,underflow,dble(merge(1,0,diag_finite(value))), &
         dble(merge(1,0,outside)),dble(merge(1,0,available)),dble(iproc),dble(frequency%calls),real(z),aimag(z)])
  end subroutine

  double precision function diagnostic_limit(s)
    use alps_var, only: poly_log_max
    integer, intent(in) :: s
    diagnostic_limit=poly_log_max(s)
  end function

  subroutine diag_capture_matrix(chi,wave)
    double complex, intent(in) :: chi(:,:,:),wave(3,3)
    if (.not.diag_active) return
    actual_chi=chi; actual_wave=wave
  end subroutine

  subroutine collect_table(f)
    integer, intent(in) :: f
    integer :: counts(nproc),offsets(nproc),i,ierr,total
    double precision, allocatable :: received(:,:)
    if (.not.allocated(tables(f)%rows)) allocate(tables(f)%rows(width,1))
    call MPI_GATHER(width*tables(f)%count,1,MPI_INTEGER,counts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    total=0; offsets=0
    if (proc0) then
       do i=1,nproc
          offsets(i)=total; total=total+counts(i)
       enddo
    endif
    allocate(received(width,max(1,total/width)))
    call MPI_GATHERV(tables(f)%rows,width*tables(f)%count,MPI_DOUBLE_PRECISION,received, &
         counts,offsets,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (proc0) then
       do i=1,total/width
          if (f==5) then
             if (received(3,i)>0.d0.and..not.first_written) then
                received(14,i)=1.d0; first_written=.true.
             endif
          endif
          write(tables(f)%unit,'(*(es25.16e3,1x))') received(1:tables(f)%columns,i)
       enddo
       flush(tables(f)%unit)
    endif
  end subroutine

  subroutine diag_finish(determinant)
    double complex, intent(in) :: determinant
    double complex :: rebuilt(nspec,3,3),total,conventional,rebuild_total
    double complex, allocatable :: gathered_rebuilt(:,:,:,:)
    double precision, allocatable :: gathered_sums(:,:,:,:)
    double precision :: sums(nspec,3,3),magnitudes(3,3),floor,cancellation,error,total_absolute
    double precision :: global_extrema(5),global_location(9,5)
    integer :: counts(2),global_counts(2),i,j,s,f,m,r,ierr,order,coverage,owner,winning_rank
    integer :: seen(5),global_seen(5)
    if (.not.diag_active) return
    ! Gather diagnostic copies: unlike the production reduction, their order
    ! differs from the physical accumulation. Guard against overflow here.
    r=1
    if (proc0) r=nproc
    allocate(gathered_rebuilt(nspec,3,3,r),gathered_sums(nspec,3,3,r))
    call MPI_GATHER(reconstructed,size(rebuilt),MPI_DOUBLE_COMPLEX,gathered_rebuilt, &
         size(rebuilt),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(absolute_sum,size(sums),MPI_DOUBLE_PRECISION,gathered_sums, &
         size(sums),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    rebuilt=(0.d0,0.d0); sums=0.d0
    if (proc0) then
       do r=1,nproc
          do j=1,3
             do i=1,3
                do s=1,nspec
                   rebuilt(s,i,j)=complex_sum(rebuilt(s,i,j),gathered_rebuilt(s,i,j,r))
                   sums(s,i,j)=safe_sum(sums(s,i,j),gathered_sums(s,i,j,r))
                enddo
             enddo
          enddo
       enddo
    endif
    counts=[frequency%limiter_count,frequency%invalid_count]
    call MPI_REDUCE(counts,global_counts,2,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(frequency%extrema,global_extrema,5,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    seen=merge(1,0,frequency%extrema_seen)
    call MPI_ALLREDUCE(seen,global_seen,5,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    do m=1,5
       owner=nproc
       if (frequency%extrema_seen(m).and.frequency%extrema(m)==global_extrema(m)) owner=iproc
       if (global_seen(m)==0) owner=0
       call MPI_ALLREDUCE(owner,winning_rank,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
       global_location(:,m)=frequency%location(:,m)
       call MPI_BCAST(global_location(:,m),9,MPI_DOUBLE_PRECISION,winning_rank,MPI_COMM_WORLD,ierr)
    enddo
    call append_row(5,[dble(eval_id),dble(iproc),dble(frequency%limiter_count),dble(frequency%calls), &
         dble(frequency%invalid_count),frequency%first(1:7),frequency%first(9),0.d0])
    if (proc0) then
       order=-1
       do s=1,nspec
          if (ACmethod(s)/=2.or.usebM(s)) cycle
          if (order==-1) then
             order=poly_order(s)
          elseif (order/=poly_order(s)) then
             order=-2; exit
          endif
       enddo
       do i=1,3
          do j=1,3
             magnitudes(i,j)=diag_abs(actual_wave(i,j))
          enddo
       enddo
       if (any(ieee_is_nan(magnitudes))) magnitudes=nan_value()
       call append_row(1,[dble(eval_id),dble(map_ir),dble(map_ii),real(omega),aimag(omega),kperp,kpar, &
            dble(order),real(determinant),aimag(determinant),diag_logabs(determinant), &
            maxval(magnitudes),minval(magnitudes),dble(merge(1,0,diag_finite(determinant))), &
            dble(merge(1,0,all(ieee_is_finite(real(actual_wave))).and. &
            all(ieee_is_finite(aimag(actual_wave))))),dble(global_counts(1)),dble(global_counts(2))])
       do m=1,5
          coverage=global_seen(m)
          if (coverage==0) global_extrema(m)=nan_value()
          call append_row(4,[dble(eval_id),dble(m),global_extrema(m),global_location(:,m),dble(coverage)])
       enddo
       do i=1,3
          do j=1,3
             call matrix_row(0,0,i,j,actual_wave(i,j),cmplx(nan_value(),nan_value(),kind(1.d0)), &
                  nan_value(),nan_value(),1)
          enddo
       enddo
       do m=1,6
          i=ci(m); j=cj(m)
          total=(0.d0,0.d0)
          rebuild_total=(0.d0,0.d0); total_absolute=0.d0
          do s=1,nspec
             total=complex_sum(total,actual_chi(s,i,j))
             rebuild_total=complex_sum(rebuild_total,rebuilt(s,i,j))
             total_absolute=safe_sum(total_absolute,sums(s,i,j))
             coverage=merge(1,0,.not.usebM(s).and..not.relativistic(s))
             error=nan_value(); cancellation=nan_value()
             if (coverage==1) then
                error=diag_abs(complex_sum(rebuilt(s,i,j),-actual_chi(s,i,j)))
                floor=max(diag_abs(actual_chi(s,i,j)),epsilon(1.d0)*sums(s,i,j),tiny(1.d0))
                cancellation=sums(s,i,j)/floor
             endif
             conventional=normalized(actual_chi(s,i,j))
             call matrix_row(1,s,i,j,actual_chi(s,i,j),conventional,error,cancellation,coverage)
          enddo
          coverage=merge(1,0,.not.any(usebM).and..not.any(relativistic))
          cancellation=nan_value(); error=nan_value()
          if (coverage==1) then
             floor=max(diag_abs(total),epsilon(1.d0)*total_absolute,tiny(1.d0))
             cancellation=total_absolute/floor
             error=diag_abs(complex_sum(rebuild_total,-total))
          endif
          call matrix_row(2,0,i,j,total,normalized(total),error,cancellation,coverage)
       enddo
    endif
    do f=1,nfiles
       if (f>5.and.diagnostic_verbosity/=2) cycle
       call collect_table(f)
    enddo
    diag_active=.false.
  end subroutine

  subroutine matrix_row(kind,s,i,j,value,conventional,error,cancellation,coverage)
    integer, intent(in) :: kind,s,i,j,coverage
    double complex, intent(in) :: value,conventional
    double precision, intent(in) :: error,cancellation
    call append_row(2,[dble(eval_id),dble(kind),dble(s),dble(i),dble(j),real(value),aimag(value), &
         real(conventional),aimag(conventional),error,cancellation,dble(coverage)])
  end subroutine

  subroutine diag_close
    integer :: f
    if (.not.opened) return
    if (proc0) then
       do f=1,nfiles
          if (f>5.and.diagnostic_verbosity/=2) cycle
          close(tables(f)%unit)
       enddo
    endif
    opened=.false.
    deallocate(reconstructed,actual_chi,absolute_sum)
  end subroutine
end module alps_diagnostics
