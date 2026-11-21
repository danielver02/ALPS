! MPI/ALPS adapter for the independent scalar winding kernel.
module alps_nyquist_driver
  use ALPS_nyquist
  use alps_nyquist_growth
  use alps_var, only: proc0, kperp,kpar,vA,nspec,ns,qs,ms,usebM,relativistic, &
       bMbetas,bMalphas,bMpdrifts,runname,foldername,unit_error,nproc,n_scan,scan,scan_option
  use alps_fns, only: disp,scan_wavevector,update_scan_wavevector
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  private
  public :: read_nyquist, run_nyquist, run_nyquist_scan
  public :: nyq_scan_coordinates
  type :: nyq_point_result
    type(nyq_result) :: initial
    type(nyq_growth_result) :: growth
    integer :: logical_requests=0,disp_calls=0,cache_hits=0,min_calls=0,max_calls=0
    real(nyq_dp) :: wall_seconds=0,determinant_seconds=0
  end type
  character(16) :: nyq_geometry='rectangle'
  real(nyq_dp) :: nyq_omega_min=0,nyq_omega_max=0,nyq_gamma_min=1.d-4,nyq_gamma_max=0
  type(nyq_options) :: nyq_sampling
  logical :: nyq_save_samples=.false.,nyq_normalize_rectangle=.false.
  logical :: nyq_find_max_growth=.false.
  real(nyq_dp) :: nyq_growth_tolerance=1.d-4
  integer :: nyq_growth_max_iterations=32
  character(16) :: nyq_growth_spacing='linear'
  type(nyq_path) :: path
  integer, parameter :: command_stop=0,command_evaluate=1
  namelist /nyquist/ nyq_geometry,nyq_omega_min,nyq_omega_max,nyq_gamma_min,nyq_gamma_max, &
       nyq_sampling,nyq_save_samples,nyq_normalize_rectangle,nyq_find_max_growth, &
       nyq_growth_tolerance,nyq_growth_max_iterations,nyq_growth_spacing
contains
  subroutine read_nyquist
    ! Rank 0 only, after ordinary species input. No per-frequency initialization.
    integer :: u,ios,code
    character(512) :: message
    ! Explicit real-frequency/upper bounds are mandatory; zero is a valid supplied real bound.
    nyq_omega_min=ieee_value(0._nyq_dp,ieee_quiet_nan)
    nyq_omega_max=ieee_value(0._nyq_dp,ieee_quiet_nan)
    nyq_gamma_max=ieee_value(0._nyq_dp,ieee_quiet_nan)
    open(newunit=u,file=trim(foldername)//trim(runname)//'.in',status='old',action='read',iostat=ios)
    if (ios/=0) call input_error('Cannot open input for &nyquist.')
    read(u,nml=nyquist,iostat=ios,iomsg=message)
    close(u)
    if (ios/=0) call input_error('Cannot read required &nyquist: '//trim(message))
    if (.not.ieee_is_finite(nyq_growth_tolerance)) call input_error('nyq_growth_tolerance must be finite.')
    if (nyq_growth_tolerance<=0) call input_error('nyq_growth_tolerance must be strictly positive.')
    if (nyq_growth_max_iterations<1) call input_error('nyq_growth_max_iterations must be positive.')
    if (nyq_find_max_growth.and.trim(nyq_geometry)/='rectangle') &
        call input_error('Maximum-growth estimation requires nyq_geometry=rectangle.')
    if (nyq_find_max_growth.and.nyq_growth_spacing/='linear'.and.nyq_growth_spacing/='log') &
        call input_error('nyq_growth_spacing must be linear or log.')
    if (.not.ieee_is_finite(nyq_gamma_min)) call input_error('nyq_gamma_min must be finite.')
    if (nyq_gamma_min<=0) call input_error('nyq_gamma_min must be strictly positive.')
    select case(trim(nyq_geometry))
    case('rectangle')
      call nyq_rectangle(nyq_omega_min,nyq_omega_max,nyq_gamma_min,nyq_gamma_max,path,code)
    case('horizontal')
      if (.not.ieee_is_finite(nyq_gamma_max)) nyq_gamma_max=0 ! unused; preserve horizontal default output
      call nyq_horizontal(nyq_omega_min,nyq_omega_max,nyq_gamma_min,path,code)
      if (nyq_normalize_rectangle) call input_error('nyq_normalize_rectangle applies only to rectangles.')
    case default
      call input_error('nyq_geometry must be rectangle or horizontal.')
    end select
    if (code/=nyq_ok) call input_error('Nyquist bounds must be finite and correctly ordered.')
    if ((scan_option==4.or.scan_option==5).and..not.path%closed) &
        call input_error('Nyquist line/plane scans require rectangle geometry.')
    if (.not.all(ieee_is_finite([kperp,kpar,vA]))) call input_error('Wavevector and vA must be finite.')
    if (kpar==0.or.vA<=0) call input_error('A2 requires kpar /= 0 and vA > 0.')
    ! Physical interpretation is validated only for the warm analytic backend in A2.
    if (.not.all(usebM).or.any(relativistic)) call input_error('A2 requires use_bM=T, relat=F for all species.')
    if (.not.all(ieee_is_finite(ns)).or..not.all(ieee_is_finite(qs)).or. &
        .not.all(ieee_is_finite(ms)).or..not.all(ieee_is_finite(bMbetas)).or. &
        .not.all(ieee_is_finite(bMalphas)).or..not.all(ieee_is_finite(bMpdrifts))) &
        call input_error('Species parameters must be finite.')
    if (any(ns<=0).or.any(ms<=0).or.any(qs==0).or.any(bMbetas<=0).or.any(bMalphas<=0)) &
        call input_error('A2 requires positive densities, masses, beta, anisotropy and nonzero charges.')
  end subroutine

  subroutine input_error(message)
    character(*), intent(in) :: message
    integer :: ierr
    write(*,'(a)') 'ERROR: '//message
    write(unit_error,'(a)') 'ERROR: '//message
    call mpi_abort(MPI_COMM_WORLD,31,ierr)
    error stop 'Invalid Nyquist input'
  end subroutine

  subroutine run_nyquist(operation_status,point,point_basename)
    integer, intent(out) :: operation_status
    type(nyq_point_result), intent(out), optional :: point
    character(*), intent(in), optional :: point_basename
    type(nyq_result) :: result
    type(nyq_growth_result) :: growth
    type(nyq_raw_cache) :: cache
    integer :: command, ierr, local_calls, min_calls, max_calls, u, ios, sample_ios, j, close_ios
    real(nyq_dp) :: started,elapsed,wall_seconds,sample_seconds,bottom_phase,other_phase,sample_started,driver_seconds
    real(nyq_dp) :: determinant_seconds,miss_started
    complex(nyq_dp) :: frequency,determinant
    logical :: normalized,have_closure,cache_active
    integer :: logical_requests,cache_hits,history_ios
    character(:), allocatable :: basename
    local_calls=0; operation_status=0; sample_ios=0; sample_seconds=0
    have_closure=.false.; bottom_phase=0; other_phase=0
    logical_requests=0; cache_hits=0; determinant_seconds=0; history_ios=0
    cache=nyq_raw_cache() ! Never carry values into a later operation/wavevector.
    cache_active=nyq_find_max_growth
    normalized=trim(nyq_geometry)=='horizontal'.or.nyq_normalize_rectangle
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    started=MPI_Wtime()
    if (proc0) then
      if (nyq_find_max_growth) then
        call nyq_growth_search(collective_evaluator,nyq_omega_min,nyq_omega_max,nyq_gamma_min,nyq_gamma_max, &
             nyq_sampling,nyq_growth_tolerance,nyq_growth_max_iterations,result,growth, &
             nyq_save_samples.or.nyq_normalize_rectangle,nyq_growth_spacing)
      else
        call nyq_evaluate(collective_evaluator,path,nyq_sampling,result, &
             nyq_save_samples.or.nyq_normalize_rectangle)
      endif
      ! Always release workers, including invalid values/options and exhausted budgets.
      command=command_stop
      call mpi_bcast(command,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    else
      do
        call mpi_bcast(command,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (command==command_stop) exit
        call mpi_bcast(frequency,1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
        determinant=disp(frequency)
        local_calls=local_calls+1
      end do
    endif
    elapsed=MPI_Wtime()-started
    call mpi_allreduce(local_calls,min_calls,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(local_calls,max_calls,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(elapsed,wall_seconds,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (proc0) then
      basename='solution/'//trim(runname)//'.nyquist'
      if(present(point_basename)) basename=point_basename
      if (allocated(result%z).and.path%closed.and.result%converged) then
        do j=1,size(result%z)
          if (result%z(j)==path%vertices(2)) then
            bottom_phase=result%cumulative_phase(j)
            other_phase=result%phase-bottom_phase
            have_closure=.true.
            exit
          endif
        end do
      endif
      if (nyq_save_samples) then
        sample_started=MPI_Wtime()
        call nyq_write_ascii(basename//'.samples',result,sample_ios)
        sample_seconds=MPI_Wtime()-sample_started
      endif
      if (nyq_find_max_growth.and..not.present(point)) call nyq_growth_history(basename//'.growth',growth,history_ios)
      driver_seconds=MPI_Wtime()-started
      if (min_calls/=max_calls.or.max_calls+cache_hits/=logical_requests.or.sample_ios/=0.or.history_ios/=0) &
          operation_status=1
      if (.not.result%converged) operation_status=1
      if (nyq_find_max_growth) then
        if (growth%logical_requests/=logical_requests) operation_status=1
        if (growth%status/=growth_ok.and.growth%status/=growth_no_roots) operation_status=1
      endif
      if(present(point)) then
        point%initial=result; point%growth=growth
        point%logical_requests=logical_requests; point%disp_calls=local_calls; point%cache_hits=cache_hits
        point%min_calls=min_calls; point%max_calls=max_calls
        point%wall_seconds=wall_seconds; point%determinant_seconds=determinant_seconds
      else
      call summary(6,ios)
      if (ios/=0) operation_status=1
      open(newunit=u,file=basename//'.summary',status='replace',action='write',iostat=ios)
      if (ios==0) then
        call summary(u,ios)
        close(u,iostat=close_ios)
        if (ios==0) ios=close_ios
      endif
      if (ios/=0) then
        write(*,'(a,i0)') 'ERROR: cannot write Nyquist summary; iostat=',ios
        operation_status=1
      endif
      endif
    endif
    call mpi_bcast(operation_status,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  contains
    subroutine collective_evaluator(z,value,status)
      complex(nyq_dp), intent(in) :: z
      complex(nyq_dp), intent(out) :: value
      integer, intent(out) :: status
      integer :: request
      logical :: found
      logical_requests=logical_requests+1
      status=0; found=.false.
      if (cache_active) call nyq_cache_get(cache,z,value,found)
      if (found) then
        cache_hits=cache_hits+1
        if (normalized) call normalize(z,value,status)
        return ! Workers stay in the request loop; no communication on cache hits.
      endif
      miss_started=MPI_Wtime()
      request=command_evaluate
      call mpi_bcast(request,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      frequency=z
      call mpi_bcast(frequency,1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
      value=disp(z)
      local_calls=local_calls+1
      determinant_seconds=determinant_seconds+MPI_Wtime()-miss_started
      ! Pass invalid raw values unchanged to A1: do not conceal determinant failures.
      if (.not.all(ieee_is_finite([real(value,nyq_dp),aimag(value)]))) return
      if (max(abs(real(value,nyq_dp)),abs(aimag(value)))==0) return
      if (cache_active) call nyq_cache_put(cache,z,value,status)
      if (status/=0) return
      if (normalized) call normalize(z,value,status)
    end subroutine

    subroutine summary(out,output_status)
      integer, intent(in) :: out
      integer, intent(out) :: output_status
      integer :: s
      output_status=0
      write(unit=out,iostat=output_status,err=990,fmt='(a)') &
          'NYQUIST single-wavevector result (frequencies in reference gyrofrequency units)'
      write(unit=out,iostat=output_status,err=990,fmt='(a,a)') 'equilibrium_input=',trim(foldername)//trim(runname)//'.in'
      write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'kperp=',kperp
      write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'kpar=',kpar
      write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'vA_over_c=',vA
      write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'mpi_ranks=',nproc
      do s=1,nspec
        write(unit=out,iostat=output_status,err=990,fmt='(a,i0,a,6(es24.16,1x))') 'species_',s,'_n_q_m_beta_alpha_drift=', &
             ns(s),qs(s),ms(s),bMbetas(s),bMalphas(s),bMpdrifts(s)
      enddo
      write(unit=out,iostat=output_status,err=990,fmt='(a,a)') 'geometry=',trim(nyq_geometry)
      write(unit=out,iostat=output_status,err=990,fmt='(a,4(es24.16,1x))') 'omega_min_max_gamma_min_max=', &
           nyq_omega_min,nyq_omega_max,nyq_gamma_min,nyq_gamma_max
      write(unit=out,iostat=output_status,err=990,fmt='(a,l1)') 'normalized_determinant=',normalized
      write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'status=',result%status
      select case(result%status)
      case(nyq_ok)
        write(unit=out,iostat=output_status,err=990,fmt='(a)') 'status_description=sampling converged on prescribed finite path'
      case(nyq_bad_input)
        write(unit=out,iostat=output_status,err=990,fmt='(a)') 'status_description=invalid path or sampling controls'
      case(nyq_bad_value)
        write(unit=out,iostat=output_status,err=990,fmt='(a)') 'status_description=zero or nonfinite determinant value'
      case(nyq_evaluator_failed)
        write(unit=out,iostat=output_status,err=990,fmt='(a)') &
            'status_description=evaluator failure; code 101: normalization range; code 102: cache allocation'
      case(nyq_budget)
        write(unit=out,iostat=output_status,err=990,fmt='(a)') 'status_description=total evaluation budget exhausted'
      case(nyq_depth)
        write(unit=out,iostat=output_status,err=990,fmt='(a)') &
            'status_description=refinement depth or coordinate resolution exhausted'
      case(nyq_no_agreement)
        write(unit=out,iostat=output_status,err=990,fmt='(a)') &
            'status_description=global sampling passes or integer consistency did not agree'
      end select
      write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'evaluator_status=',result%evaluator_status
      write(unit=out,iostat=output_status,err=990,fmt='(a,l1)') 'sampling_converged=',result%converged
      write(unit=out,iostat=output_status,err=990,fmt='(a,l1)') 'integer_accepted=',result%integer_accepted
      write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'total_phase=',result%phase
      write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'phase_over_2pi=',result%raw_winding
      if (result%converged) then
        write(unit=out,iostat=output_status,err=990,fmt='(a)') 'phase_scope=complete prescribed finite path'
      else
        write(unit=out,iostat=output_status,err=990,fmt='(a)') &
            'phase_scope=unaccepted latest path or partial prefix; endpoints may not reach requested bounds'
      endif
      if (result%integer_accepted) then
        write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'accepted_winding=',result%winding
        write(unit=out,iostat=output_status,err=990,fmt='(a)') &
            'interpretation=roots within stated finite rectangle (warm analytic response, no upper poles)'
      else
        write(unit=out,iostat=output_status,err=990,fmt='(a)') 'accepted_winding=NONE'
      endif
      if (.not.path%closed) then
        write(unit=out,iostat=output_status,err=990,fmt='(a)') &
            'asymptotic_closure=UNVERIFIED; finite-interval phase is not an instability count'
      endif
      write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'disp_calls=',local_calls
      write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'logical_requests=',logical_requests
      write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'initial_logical_requests=',result%evaluations
      write(unit=out,iostat=output_status,err=990,fmt='(a,l1)') 'cache_enabled=',cache_active
      write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'cache_hits=',cache_hits
      write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'cache_misses=',local_calls
      write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') &
          'cache_hit_fraction=',real(cache_hits,nyq_dp)/real(max(1,logical_requests),nyq_dp)
      write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'determinant_seconds=',determinant_seconds
      write(unit=out,iostat=output_status,err=990,fmt='(a,2(i0,1x))') 'rank_calls_min_max=',min_calls,max_calls
      write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'wall_seconds=',wall_seconds
      write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'driver_seconds_before_summary=',driver_seconds
      write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'sample_output_seconds=',sample_seconds
      write(unit=out,iostat=output_status,err=990,fmt='(a,l1)') 'sample_output_enabled=',nyq_save_samples
      write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'sample_output_iostat=',sample_ios
      write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'passes=',result%passes
      write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'refinements=',result%refinements
      write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'deepest_level=',result%deepest_level
      write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'max_increment=',result%max_increment
      write(unit=out,iostat=output_status,err=990,fmt='(a,2(es24.16,1x))') 'endpoint_first=',result%endpoints(1)
      write(unit=out,iostat=output_status,err=990,fmt='(a,2(es24.16,1x))') 'endpoint_last=',result%endpoints(2)
      if (allocated(result%pass_phases)) &
          write(unit=out,iostat=output_status,err=990,fmt='(a,*(es24.16,1x))') 'pass_phases=',result%pass_phases(:result%passes)
      if (have_closure) then
        write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'bottom_phase=',bottom_phase
        write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'other_three_phase=',other_phase
      endif
      write(unit=out,iostat=output_status,err=990,fmt='(a,l1)') 'growth_requested=',nyq_find_max_growth
      if (nyq_find_max_growth) then
        write(unit=out,iostat=output_status,err=990,fmt='(a)') 'status_scope=initial contour; growth_status is overall search'
        write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'growth_status=',growth%status
        write(unit=out,iostat=output_status,err=990,fmt='(a,a)') 'growth_spacing=',trim(nyq_growth_spacing)
        write(unit=out,iostat=output_status,err=990,fmt='(a,2(es24.16,1x))') &
            'growth_original_gamma_bounds=',nyq_gamma_min,nyq_gamma_max
        write(unit=out,iostat=output_status,err=990,fmt='(a,a)') &
            'growth_status_description=',trim(nyq_growth_status_text(growth%status))
        if (result%integer_accepted) then
          write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'growth_initial_winding=',result%winding
        else
          write(unit=out,iostat=output_status,err=990,fmt='(a)') 'growth_initial_winding=NONE'
        endif
        write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'growth_tolerance=',nyq_growth_tolerance
        if (nyq_growth_spacing=='log') then
          write(unit=out,iostat=output_status,err=990,fmt='(a)') 'growth_tolerance_type=relative'
        else
          write(unit=out,iostat=output_status,err=990,fmt='(a)') 'growth_tolerance_type=absolute'
        endif
        write(unit=out,iostat=output_status,err=990,fmt='(a,a)') 'growth_history_file=',basename//'.growth'
        write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'growth_failure_A1_status=',growth%failure_status
        write(unit=out,iostat=output_status,err=990,fmt='(a,l1)') 'growth_converged=',growth%converged
        write(unit=out,iostat=output_status,err=990,fmt='(a,l1)') 'growth_no_roots=',growth%no_roots
        write(unit=out,iostat=output_status,err=990,fmt='(a,l1)') 'growth_bracket_valid=',growth%bracket_valid
        write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'growth_iterations=',growth%iterations
        write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'growth_attempts=',growth%attempts
        write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'growth_total_passes=',growth%passes
        write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'growth_total_refinements=',growth%refinements
        write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'growth_max_depth=',growth%depth
        write(unit=out,iostat=output_status,err=990,fmt='(a,i0)') 'growth_history_iostat=',history_ios
        if (growth%no_roots) write(unit=out,iostat=output_status,err=990,fmt='(a)') &
            'growth_interpretation=no roots above selected lower boundary within stated frequency domain'
        if (growth%bracket_valid) then
          write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'growth_lower=',growth%lower
          write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'growth_upper=',growth%upper
          write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') 'growth_width=',growth%upper-growth%lower
          write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') &
              'growth_relative_width=',nyq_growth_relative_width(growth%lower,growth%upper)
        endif
        if (growth%converged) then
          write(unit=out,iostat=output_status,err=990,fmt='(a,es24.16)') &
              'growth_midpoint_estimate=',0.5_nyq_dp*growth%lower+0.5_nyq_dp*growth%upper
          write(unit=out,iostat=output_status,err=990,fmt='(a)') &
              'growth_interpretation=maximum growth bracket within input real-frequency and upper-growth bounds; not a root'
        else
          write(unit=out,iostat=output_status,err=990,fmt='(a)') 'growth_midpoint_estimate=NONE'
        endif
      endif
      write(unit=out,nml=nyquist,iostat=output_status,err=990)
990   continue
    end subroutine
  end subroutine

  subroutine nyq_scan_coordinates(plane,coordinates,segments,indices)
    ! Rank-0 coordinate preflight only. No determinant calls or MPI scheduling.
    logical, intent(in) :: plane
    real(nyq_dp), allocatable, intent(out) :: coordinates(:,:)
    integer, allocatable, intent(out) :: segments(:),indices(:,:)
    real(nyq_dp) :: original(2),base(2),theta,magnitude,theta_inner,magnitude_inner,delta
    integer :: count,s,q,r,j,nouter,ninner,ios
    original=[kperp,kpar]
    if(plane) then
      if(n_scan/=2) call input_error('scan_option=5 requires n_scan=2.')
      if(any(scan%type_s<1).or.any(scan%type_s>4).or.scan(1)%type_s==scan(2)%type_s) &
          call input_error('Nyquist plane requires two different scan types in 1..4.')
      if(any(scan%n_out==huge(count))) call input_error('Too many Nyquist grid points.')
      nouter=scan(1)%n_out+1; ninner=scan(2)%n_out+1
      if(nouter>huge(count)/ninner) call input_error('Too many Nyquist grid points.')
      count=nouter*ninner
    else
      count=1
      do s=1,n_scan
        if(scan(s)%n_out>huge(count)-count) call input_error('Too many Nyquist scan points.')
        count=count+scan(s)%n_out
      enddo
    endif
    allocate(coordinates(2,count),segments(count),indices(2,count),stat=ios)
    if(ios/=0) call input_error('Cannot allocate Nyquist wavevector path.')
    indices=0
    if(plane) then
      theta=atan(original(1)/original(2)); magnitude=sqrt(original(1)**2+original(2)**2)
      j=0
      do q=0,scan(1)%n_out
        kperp=original(1); kpar=original(2)
        call scan_wavevector(1,q*scan(1)%n_res,original(1),original(2),theta,magnitude,delta)
        call validate_coordinate
        base=[kperp,kpar]
        theta_inner=atan(base(1)/base(2)); magnitude_inner=sqrt(base(1)**2+base(2)**2)
        do r=0,scan(2)%n_out
          kperp=base(1); kpar=base(2)
          call scan_wavevector(2,r*scan(2)%n_res,base(1),base(2),theta_inner,magnitude_inner,delta)
          call validate_coordinate
          j=j+1; coordinates(:,j)=[kperp,kpar]
          segments(j)=q+1; indices(:,j)=[q+1,r+1]
        enddo
      enddo
    else
      coordinates(:,1)=original; segments(1)=1; j=1
      do s=1,n_scan
        base=[kperp,kpar]
        theta=atan(base(1)/base(2)); magnitude=sqrt(base(1)**2+base(2)**2)
        do q=1,scan(s)%n_out
          call scan_wavevector(s,q*scan(s)%n_res,base(1),base(2),theta,magnitude,delta)
          call validate_coordinate
          j=j+1; coordinates(:,j)=[kperp,kpar]; segments(j)=s
        enddo
      enddo
    endif
    kperp=original(1); kpar=original(2)
  contains
    subroutine validate_coordinate
      if(.not.all(ieee_is_finite([kperp,kpar]))) call input_error('Nonfinite scan wavevector.')
      if(kperp<0.or.kpar==0) call input_error('Nyquist scan requires kperp >= 0 and kpar /= 0.')
    end subroutine
  end subroutine nyq_scan_coordinates

  subroutine run_nyquist_scan(operation_status,shared_initialization_seconds)
    integer, intent(out) :: operation_status
    real(nyq_dp), intent(in) :: shared_initialization_seconds
    type(nyq_point_result) :: point
    real(nyq_dp), allocatable :: coordinates(:,:)
    integer, allocatable :: segments(:),indices(:,:)
    real(nyq_dp) :: pair(2),update_started,update_seconds
    real(nyq_dp) :: started,nan,lo,hi,mid,width,relative_width,initialization_seconds
    integer :: count,j,u,v,w,ios,close_ios,ierr,point_status,gstatus
    integer :: positive,zero,unresolved,growth_success,growth_failed,total_calls,rebuilds,passes,refinements,depth,iterations
    integer :: negative,growth_search_failed
    character(32) :: winding,tag
    character(:), allocatable :: basename
    character(:), allocatable :: index_header
    character(3) :: row_advance
    logical :: rebuild,scan_open,data_open,plane
    operation_status=0; count=0
    positive=0; zero=0; unresolved=0; growth_success=0; growth_failed=0; total_calls=0; rebuilds=0
    negative=0; growth_search_failed=0; scan_open=.false.; data_open=.false.
    plane=scan_option==5; index_header=''; row_advance='yes'
    if(plane) then
      index_header=' outer_index inner_index'; row_advance='no'
    endif
    nan=ieee_value(0._nyq_dp,ieee_quiet_nan)
    started=MPI_Wtime()
    if(proc0) then
      call nyq_scan_coordinates(plane,coordinates,segments,indices)
      count=size(segments)
      basename='solution/'//trim(runname)//'.nyquist'
      open(newunit=u,file=basename//'.scan',status='replace',action='write',iostat=ios)
      scan_open=ios==0
      if(ios==0) then
        open(newunit=v,file=basename//'.dat',status='replace',action='write',iostat=ios)
        data_open=ios==0
      endif
      if(ios==0) then
        call scan_headers(u,v,ios)
      endif
      if(ios/=0) then
        if(scan_open) close(u,iostat=close_ios)
        if(data_open) close(v,iostat=close_ios)
        write(*,'(a)') 'ERROR: cannot initialize Nyquist scan output files.'
        operation_status=1
      endif
    endif
    call mpi_bcast(operation_status,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if(operation_status/=0) return
    call mpi_bcast(count,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    initialization_seconds=MPI_Wtime()-started
    do j=1,count
      if(proc0) pair=coordinates(:,j)
      call mpi_bcast(pair,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      rebuild=pair(1)/=kperp
      update_started=MPI_Wtime()
      kperp=pair(1); kpar=pair(2)
      call update_scan_wavevector(rebuild)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      update_seconds=MPI_Wtime()-update_started
      if(rebuild) rebuilds=rebuilds+1
      ! A fresh local cache and counters are constructed on every invocation.
      ! All workers return to this outer loop at the point-end command, not MPI_Finalize.
      write(tag,'(a,i6.6)') '.point_',j
      if(proc0) then
        call run_nyquist(point_status,point,basename//trim(tag))
      else
        call run_nyquist(point_status,point)
      endif
      if(proc0) then
        winding='NaN'; lo=nan; hi=nan; mid=nan; width=nan; relative_width=nan
        gstatus=-1; iterations=0
        passes=point%initial%passes; refinements=point%initial%refinements; depth=point%initial%deepest_level
        if(point%initial%integer_accepted.and.point%initial%converged) then
          write(winding,'(i0)') point%initial%winding
          if(point%initial%winding>0) positive=positive+1
          if(point%initial%winding==0) zero=zero+1
          if(point%initial%winding<0) then
            point_status=1; negative=negative+1
          endif
        else
          unresolved=unresolved+1
        endif
        if(nyq_find_max_growth) then
          gstatus=point%growth%status; iterations=point%growth%iterations
          passes=point%growth%passes; refinements=point%growth%refinements; depth=point%growth%depth
          if(point%growth%bracket_valid) then
            lo=point%growth%lower; hi=point%growth%upper
            width=hi-lo; relative_width=nyq_growth_relative_width(lo,hi)
          endif
          if(point%growth%converged) then
            mid=0.5_nyq_dp*lo+0.5_nyq_dp*hi; growth_success=growth_success+1
          else if(point%growth%status/=growth_no_roots) then
            growth_failed=growth_failed+1
            if(point%initial%integer_accepted.and.point%initial%winding>0) growth_search_failed=growth_search_failed+1
          endif
        endif
        total_calls=total_calls+point%disp_calls
        if(point_status/=0) operation_status=1
        write(u,'(2(i0,1x),2(es25.16e3,1x),a,1x,2(i0,1x),5(es25.16e3,1x),3(i0,1x),'// &
            '3(es25.16e3,1x),8(i0,1x))',iostat=ios,advance=row_advance) j,segments(j),kperp,kpar,trim(winding), &
            point%initial%status,gstatus,lo,hi,mid,width,relative_width,point%logical_requests, &
            point%disp_calls,point%cache_hits,point%wall_seconds,update_seconds,point%determinant_seconds, &
            passes,refinements,depth,iterations,point%min_calls,point%max_calls,point_status,point%initial%evaluator_status
        if(plane.and.ios==0) write(u,'(2(1x,i0))',iostat=ios) indices(:,j)
        if(ios==0) write(v,'(2(es25.16e3,1x),a,1x,es25.16e3)',iostat=ios) kperp,kpar,trim(winding),mid
        if(ios/=0) then
          write(*,'(a)') 'ERROR: cannot write Nyquist scan row; stopping scan.'
          operation_status=2
        endif
      endif
      call mpi_bcast(operation_status,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(operation_status==2) exit ! Output failure is global; numerical point failures continue.
    enddo
    if(proc0) then
      close(u,iostat=close_ios)
      if(close_ios/=0) operation_status=1
      close(v,iostat=close_ios)
      if(close_ios/=0) operation_status=1
      open(newunit=w,file=basename//'.summary',status='replace',action='write',iostat=ios)
      if(ios==0) then
        call scan_totals(w,ios)
        close(w,iostat=close_ios)
        if(close_ios/=0) ios=close_ios
      endif
      if(ios/=0) operation_status=1
      call scan_totals(6,ios)
    endif
    call mpi_bcast(operation_status,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  contains
    subroutine scan_headers(out,data,code)
      integer, intent(in) :: out,data
      integer, intent(out) :: code
      integer :: species
      if(plane) then
        write(out,'(a)',iostat=code) '# Nyquist plane scan; wavevectors k*d_ref, frequencies omega/Omega_ref; bounded domain only'
      else
        write(out,'(a)',iostat=code) '# Nyquist line scan; wavevectors k*d_ref, frequencies omega/Omega_ref; bounded domain only'
      endif
      if(code/=0) return
      write(out,'(a,a)',iostat=code) '# equilibrium_input=',trim(foldername)//trim(runname)//'.in'
      if(code/=0) return
      write(out,'(a,4(es25.16e3,1x))',iostat=code) '# frequency_bounds=', &
          nyq_omega_min,nyq_omega_max,nyq_gamma_min,nyq_gamma_max
      if(code/=0) return
      write(out,'(a,l1,a,a,a,es25.16e3)',iostat=code) '# growth_requested=',nyq_find_max_growth, &
          ' spacing=',trim(nyq_growth_spacing),' tolerance=',nyq_growth_tolerance
      if(code/=0) return
      if(nyq_growth_spacing=='log') then
        write(out,'(a)',iostat=code) '# growth_tolerance_type=relative'
      else
        write(out,'(a)',iostat=code) '# growth_tolerance_type=absolute'
      endif
      if(code/=0) return
      write(out,'(a,4(i0,1x),a,i0)',iostat=code) '# sampling intervals depth evaluations passes=', &
          nyq_sampling%initial_intervals,nyq_sampling%max_depth,nyq_sampling%max_evaluations,nyq_sampling%max_passes, &
          ' growth_max_iterations=',nyq_growth_max_iterations
      if(code/=0) return
      write(out,'(a,3(es25.16e3,1x))',iostat=code) '# sampling phase curvature agreement=', &
          nyq_sampling%phase_limit,nyq_sampling%curvature_limit,nyq_sampling%phase_tolerance
      if(code/=0) return
      write(out,'(a,l1,a,es25.16e3)',iostat=code) '# normalized_determinant=',nyq_normalize_rectangle,' vA/c=',vA
      if(code/=0) return
      do species=1,nspec
        write(out,'(a,i0,a,6(es25.16e3,1x))',iostat=code) '# species_',species,' n q m beta alpha drift=', &
            ns(species),qs(species),ms(species),bMbetas(species),bMalphas(species),bMpdrifts(species)
        if(code/=0) return
      enddo
      write(out,'(a)',iostat=code) '# point segment kperp kpar initial_winding A1_status growth_status '// &
          'gamma_lower gamma_upper gamma_midpoint gamma_absolute_width gamma_relative_width '// &
          'logical_requests disp_calls cache_hits wall_seconds update_seconds determinant_seconds '// &
          'passes refinements max_depth growth_iterations rank_calls_min rank_calls_max point_status evaluator_status'//index_header
      if(code/=0) return
      if(plane) then
        write(out,'(a,2(i0,1x),a,2(i0,1x))',iostat=code) '# grid_dimensions=',scan%n_out+1, &
            'scan_types=',scan%type_s
        if(code/=0) return
        write(out,'(a)',iostat=code) '# indices one-based; segment=outer_index; outer-major, inner varies fastest'
        if(code/=0) return
      endif
      write(data,'(a)',iostat=code) '# kperp kpar Wn gamma_max'
      if(code/=0) return
      write(data,'(a)',iostat=code) '# k*d_ref, gamma/Omega_ref; Wn=original rectangle winding; NaN=unavailable'
      if(code/=0) return
      write(data,'(a,4(es25.16e3,1x))',iostat=code) '# original_frequency_bounds=', &
          nyq_omega_min,nyq_omega_max,nyq_gamma_min,nyq_gamma_max
      if(code/=0) return
      if(plane) write(data,'(a,2(i0,1x),a)',iostat=code) '# grid_dimensions=',scan%n_out+1, &
          'outer-major, scan 2 varies fastest; no blank rows'
    end subroutine
    subroutine scan_totals(out,code)
      integer, intent(in) :: out
      integer, intent(out) :: code
      code=0
      if(plane) then
        write(out,'(a)',iostat=code,err=991) 'NYQUIST wavevector-plane scan (bounded frequency domain)'
        write(out,'(a,2(i0,1x))',iostat=code,err=991) 'grid_dimensions=',scan%n_out+1
        write(out,'(a,2(i0,1x))',iostat=code,err=991) 'scan_types=',scan%type_s
        write(out,'(a)',iostat=code,err=991) 'grid_order=outer-major; scan 2 varies fastest; indices one-based'
      else
        write(out,'(a)',iostat=code,err=991) 'NYQUIST wavevector-line scan (bounded frequency domain)'
      endif
      write(out,'(a,i0)',iostat=code,err=991) 'scan_status=',operation_status
      write(out,'(a,i0)',iostat=code,err=991) 'requested_points=',count
      write(out,'(a,i0)',iostat=code,err=991) 'evaluated_points=',min(j,count)
      write(out,'(a,i0)',iostat=code,err=991) 'positive_initial_windings=',positive
      write(out,'(a,i0)',iostat=code,err=991) 'zero_initial_windings=',zero
      write(out,'(a,i0)',iostat=code,err=991) 'negative_initial_windings=',negative
      write(out,'(a,i0)',iostat=code,err=991) 'unresolved_initial_contours=',unresolved
      write(out,'(a,i0)',iostat=code,err=991) 'successful_growth_searches=',growth_success
      write(out,'(a,i0)',iostat=code,err=991) 'failed_growth_calculations=',growth_failed
      write(out,'(a,i0)',iostat=code,err=991) 'failed_growth_searches=',growth_search_failed
      write(out,'(a,i0)',iostat=code,err=991) 'disp_calls=',total_calls
      write(out,'(a,i0)',iostat=code,err=991) 'harmonic_rebuilds=',rebuilds
      write(out,'(a,es25.16e3)',iostat=code,err=991) 'scan_preparation_seconds=',initialization_seconds
      write(out,'(a,es25.16e3)',iostat=code,err=991) 'shared_initialization_seconds=',shared_initialization_seconds
      write(out,'(a,es25.16e3)',iostat=code,err=991) 'wall_seconds=',MPI_Wtime()-started
991   continue
    end subroutine
  end subroutine run_nyquist_scan

  subroutine normalize(z,value,status)
    ! Log-polar division by (z*vA)**6 avoids forming a tiny/huge sixth power.
    ! Preserve magnitude as well as phase; report unrepresentable results rather than clipping.
    complex(nyq_dp), intent(in) :: z
    complex(nyq_dp), intent(inout) :: value
    integer, intent(out) :: status
    real(nyq_dp) :: sf,sz,lf,lz,amplitude_log
    complex(nyq_dp) :: uf,uz
    status=0
    sf=max(abs(real(value,nyq_dp)),abs(aimag(value)))
    sz=max(abs(real(z,nyq_dp)),abs(aimag(z)))
    uf=value/sf; uz=z/sz
    lf=log(sf)+log(abs(uf)); lz=log(sz)+log(abs(uz))+log(vA)
    amplitude_log=lf-6*lz
    if (amplitude_log>log(huge(1._nyq_dp)).or. &
        amplitude_log<log(tiny(1._nyq_dp))+log(epsilon(1._nyq_dp))) then
      status=101
      return
    endif
    uf=uf/abs(uf); uz=uz/abs(uz)
    value=exp(amplitude_log)*uf*conjg(uz)**6
  end subroutine
end module
