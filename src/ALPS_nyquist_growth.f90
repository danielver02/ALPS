! Closed-strip growth search and exact-value storage, independent of ALPS/MPI.
! A1 remains responsible for every contour's full numerical convergence checks.
module alps_nyquist_growth
  use ALPS_nyquist
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: growth_ok=0,growth_no_roots=1,growth_initial_failed=2, &
       growth_strip_failed=3,growth_iteration_limit=4,growth_inconsistent=5,growth_resolution=6,growth_resource=7
  type, public :: nyq_growth_step
    real(nyq_dp) :: trial=0,strip_lower=0,strip_upper=0,lower=0,upper=0,phase=0
    logical :: accepted=.false.
    integer :: status=nyq_bad_input,winding=0,requests=0,passes=0,refinements=0,depth=0
    integer :: above=-1
    character(16) :: anchor='current_upper',spacing='linear'
  end type
  type, public :: nyq_growth_result
    logical :: converged=.false.,bracket_valid=.false.,no_roots=.false.
    integer :: status=growth_initial_failed,failure_status=nyq_ok,iterations=0,attempts=0
    integer :: logical_requests=0,passes=0,refinements=0,depth=0
    real(nyq_dp) :: lower=0,upper=0
    real(nyq_dp) :: original_lower=0,original_upper=0
    character(16) :: spacing='linear'
    type(nyq_growth_step), allocatable :: history(:)
  end type
  type, public :: nyq_raw_cache
    ! Cache lifetime is one operation at one wavevector; values always remain raw.
    complex(nyq_dp), allocatable :: z(:),raw(:)
    logical, allocatable :: occupied(:)
    integer :: entries=0
  end type
  public :: nyq_growth_search,nyq_growth_history,nyq_cache_get,nyq_cache_put,nyq_growth_status_text
  public :: nyq_log_trial,nyq_growth_relative_width
contains
  ! Valid positive, ordered endpoints are required by logarithmic searches.
  pure real(nyq_dp) function nyq_growth_relative_width(lower,upper) result(width)
    real(nyq_dp), intent(in) :: lower,upper
    width=(upper-lower)/(0.5_nyq_dp*lower+0.5_nyq_dp*upper)
  end function

  logical function growth_tolerance_met(lower,upper,tolerance,mode) result(met)
    real(nyq_dp), intent(in) :: lower,upper,tolerance
    character(*), intent(in) :: mode
    if (mode=='log') then
      met=nyq_growth_relative_width(lower,upper)<tolerance
    else
      met=upper-lower<tolerance
    endif
  end function

  subroutine nyq_log_trial(lower,upper,mid,status)
    real(nyq_dp), intent(in) :: lower,upper
    real(nyq_dp), intent(out) :: mid
    integer, intent(out) :: status
    mid=0; status=nyq_bad_input
    if (.not.all(ieee_is_finite([lower,upper]))) return
    if (lower<=0.or.upper<=lower) return
    mid=exp(0.5_nyq_dp*log(lower)+0.5_nyq_dp*log(upper))
    status=growth_resolution
    if (.not.ieee_is_finite(mid)) return
    if (mid<=lower.or.mid>=upper) return
    status=growth_ok
  end subroutine

  function nyq_growth_status_text(code) result(message)
    integer, intent(in) :: code
    character(80) :: message
    select case(code)
    case(growth_ok)
      message='converged bounded maximum-growth bracket'
    case(growth_no_roots)
      message='valid initial zero winding; no growth search or estimate'
    case(growth_initial_failed)
      message='initial contour unresolved; no root count or growth bracket established'
    case(growth_strip_failed)
      message='strip unresolved; last valid bracket retained'
    case(growth_iteration_limit)
      message='iteration limit reached; last valid bracket retained'
    case(growth_inconsistent)
      message='negative or inconsistent count; growth search unresolved'
    case(growth_resolution)
      message='growth-coordinate resolution exhausted; last valid bracket retained'
    case default
      message='growth-search storage allocation failed'
    end select
  end function

  subroutine nyq_growth_search(evaluate,xmin,xmax,gmin,gmax,sampling,tolerance,max_iterations,initial,answer, &
       save_initial,spacing)
    procedure(nyq_evaluator) :: evaluate
    real(nyq_dp), intent(in) :: xmin,xmax,gmin,gmax,tolerance
    type(nyq_options), intent(in) :: sampling
    integer, intent(in) :: max_iterations
    type(nyq_result), intent(out) :: initial
    type(nyq_growth_result), intent(out) :: answer
    logical, intent(in) :: save_initial
    character(*), intent(in), optional :: spacing
    type(nyq_path) :: contour
    type(nyq_result) :: strip
    real(nyq_dp) :: mid,strip_lower,strip_upper,logmid,logbottom,logtop
    integer :: code,j,lower_count,above
    character(16) :: mode,anchor
    answer%failure_status=nyq_bad_input
    answer%original_lower=gmin; answer%original_upper=gmax
    if (.not.ieee_is_finite(tolerance)) return
    if (tolerance<=0.or.max_iterations<1) return
    mode='linear'; logbottom=0; logtop=0
    if (present(spacing)) then
      if (spacing/='linear'.and.spacing/='log') return
      mode=spacing
    endif
    answer%spacing=mode
    if (mode=='log') then
      if (.not.all(ieee_is_finite([gmin,gmax]))) return
      if (gmin<=0.or.gmax<=0) return
      logbottom=log(gmin); logtop=log(gmax)
    endif
    call nyq_rectangle(xmin,xmax,gmin,gmax,contour,code)
    if (code/=nyq_ok) return
    call nyq_evaluate(evaluate,contour,sampling,initial,save_initial)
    answer%failure_status=nyq_ok
    answer%logical_requests=initial%evaluations; answer%passes=initial%passes
    answer%refinements=initial%refinements; answer%depth=initial%deepest_level
    if (.not.initial%integer_accepted.or..not.initial%converged) then
      answer%failure_status=initial%status
      return
    endif
    if (initial%winding<0) then
      answer%status=growth_inconsistent
      return
    endif
    if (initial%winding==0) then
      answer%status=growth_no_roots; answer%no_roots=.true.
      return
    endif
    answer%lower=gmin; answer%upper=gmax; answer%bracket_valid=.true.
    lower_count=initial%winding
    allocate(answer%history(max_iterations),stat=code)
    if (code/=0) then
      answer%status=growth_resource
      return
    endif
    do j=1,max_iterations
      if (growth_tolerance_met(answer%lower,answer%upper,tolerance,mode)) exit
      mid=0.5_nyq_dp*answer%lower+0.5_nyq_dp*answer%upper
      if (mode=='log') then
        call nyq_log_trial(answer%lower,answer%upper,mid,code)
        if (code/=growth_ok) then
          answer%status=growth_resolution
          return
        endif
      endif
      if (mid<=answer%lower.or.mid>=answer%upper) then
        answer%status=growth_resolution
        return
      endif
      strip_lower=mid; anchor='current_upper'
      if (mode=='linear') then
      strip_upper=answer%upper
      else
        logmid=log(mid)
        ! Exact distance ties choose the original upper boundary deterministically.
        if (abs(logmid-logbottom)<abs(logtop-logmid)) then
          strip_lower=gmin; strip_upper=mid; anchor='original_lower'
        else
          strip_upper=gmax; anchor='original_upper'
        endif
      endif
      call nyq_rectangle(xmin,xmax,strip_lower,strip_upper,contour,code)
      if (code/=nyq_ok) then
        answer%status=growth_resolution
        return
      endif
      call nyq_evaluate(evaluate,contour,sampling,strip)
      answer%attempts=j
      answer%logical_requests=answer%logical_requests+strip%evaluations
      answer%passes=answer%passes+strip%passes
      answer%refinements=answer%refinements+strip%refinements
      answer%depth=max(answer%depth,strip%deepest_level)
      associate(step=>answer%history(j))
        step%trial=mid; step%strip_lower=strip_lower; step%strip_upper=strip_upper
        step%anchor=anchor; step%spacing=mode
        step%status=strip%status; step%accepted=strip%integer_accepted.and.strip%converged
        step%winding=strip%winding; step%phase=strip%phase
        step%requests=strip%evaluations; step%passes=strip%passes
        step%refinements=strip%refinements; step%depth=strip%deepest_level
        step%lower=answer%lower; step%upper=answer%upper
        if (.not.step%accepted) then
          answer%status=growth_strip_failed; answer%failure_status=strip%status
          return ! The last valid bracket remains unchanged.
        endif
        if (strip%winding<0.or.strip%winding>initial%winding) then
          answer%status=growth_inconsistent
          return
        endif
        above=strip%winding
        if (anchor=='original_lower') above=initial%winding-strip%winding
        if (above<0.or.above>lower_count) then
          answer%status=growth_inconsistent
          return
        endif
        step%above=above
        if (above>0) then
          answer%lower=mid; lower_count=above
        else
          answer%upper=mid
        endif
        answer%iterations=j
        step%lower=answer%lower; step%upper=answer%upper
      end associate
    enddo
    if (growth_tolerance_met(answer%lower,answer%upper,tolerance,mode)) then
      answer%status=growth_ok; answer%converged=.true.
    else
      answer%status=growth_iteration_limit
    endif
  end subroutine

  subroutine nyq_growth_history(filename,answer,iostat)
    character(*), intent(in) :: filename
    type(nyq_growth_result), intent(in) :: answer
    integer, intent(out) :: iostat
    integer :: u,j,close_code
    character(16) :: count,above
    open(newunit=u,file=filename,status='replace',action='write',iostat=iostat)
    if (iostat/=0) return
    write(u,'(a)',iostat=iostat) '# attempt trial strip_low strip_high winding A1_status new_low new_high '// &
        'logical_requests phase passes refinements depth anchor N_above spacing'
    if(iostat==0) write(u,'(a,2(1x,es25.16e3),1x,a,a)',iostat=iostat) &
        '# original_gamma_bounds=',answer%original_lower,answer%original_upper,'spacing=',trim(answer%spacing)
    do j=1,answer%attempts
      if (iostat/=0) exit
      associate(s=>answer%history(j))
        count='NONE'
        above='NONE'
        if (s%accepted) write(count,'(i0)') s%winding
        if (s%above>=0) write(above,'(i0)') s%above
        write(u,'(i0,1x,3(es25.16e3,1x),a,1x,i0,1x,2(es25.16e3,1x),i0,1x,es25.16e3,3(1x,i0),3(1x,a))', &
            iostat=iostat) &
            j,s%trial,s%strip_lower,s%strip_upper,trim(count),s%status,s%lower,s%upper,s%requests, &
            s%phase,s%passes,s%refinements,s%depth,trim(s%anchor),trim(above),trim(s%spacing)
      end associate
    enddo
    close(u,iostat=close_code)
    if (iostat==0) iostat=close_code
  end subroutine

  integer function cache_slot(z,capacity) result(slot)
    complex(nyq_dp), intent(in) :: z
    integer, intent(in) :: capacity
    integer(int64) :: x,y,h
    real(nyq_dp) :: re,im
    re=real(z,nyq_dp); im=aimag(z)
    ! Numerical equality includes signed zeros, so their hashes must agree too.
    if (re==0) re=0._nyq_dp
    if (im==0) im=0._nyq_dp
    x=transfer(re,x); y=transfer(im,y)
    h=ieor(x,ishftc(y,29))
    h=ieor(h,ishft(h,-17)); h=ieor(h,ishft(h,-31))
    slot=int(iand(h,int(capacity-1,int64)))+1
  end function

  subroutine nyq_cache_get(cache,z,raw,found)
    type(nyq_raw_cache), intent(in) :: cache
    complex(nyq_dp), intent(in) :: z
    complex(nyq_dp), intent(out) :: raw
    logical, intent(out) :: found
    integer :: slot,n
    found=.false.; raw=cmplx(0,0,nyq_dp)
    if (.not.allocated(cache%z)) return
    n=size(cache%z); slot=cache_slot(z,n)
    do while(cache%occupied(slot))
      if (cache%z(slot)==z) then
        found=.true.; raw=cache%raw(slot)
        return
      endif
      slot=mod(slot,n)+1
    enddo
  end subroutine

  subroutine nyq_cache_put(cache,z,raw,status)
    type(nyq_raw_cache), intent(inout) :: cache
    complex(nyq_dp), intent(in) :: z,raw
    integer, intent(out) :: status
    type(nyq_raw_cache) :: grown
    integer :: n,j
    status=0; n=1024
    if (allocated(cache%z)) then
      n=size(cache%z)
      if (cache%entries<n/2) then
        call insert(cache,z,raw)
        return
      endif
      if (n>huge(n)-n) then
        status=102
        return
      endif
      n=2*n
    endif
    allocate(grown%z(n),grown%raw(n),grown%occupied(n),stat=status)
    if (status/=0) then
      status=102
      return
    endif
    grown%occupied=.false.
    if (allocated(cache%z)) then
      do j=1,size(cache%z)
        if (cache%occupied(j)) call insert(grown,cache%z(j),cache%raw(j))
      enddo
    endif
    call insert(grown,z,raw)
    call move_alloc(grown%z,cache%z); call move_alloc(grown%raw,cache%raw)
    call move_alloc(grown%occupied,cache%occupied); cache%entries=grown%entries
  end subroutine

  subroutine insert(cache,z,raw)
    type(nyq_raw_cache), intent(inout) :: cache
    complex(nyq_dp), intent(in) :: z,raw
    integer :: slot,n
    n=size(cache%z); slot=cache_slot(z,n)
    do while(cache%occupied(slot))
      if (cache%z(slot)==z) exit
      slot=mod(slot,n)+1
    enddo
    if (.not.cache%occupied(slot)) cache%entries=cache%entries+1
    cache%z(slot)=z; cache%raw(slot)=raw; cache%occupied(slot)=.true.
  end subroutine
end module
