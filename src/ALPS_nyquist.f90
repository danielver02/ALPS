! General scalar phase winding. No ALPS state, normalization, or MPI dependencies.
module ALPS_nyquist
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: nyq_dp = kind(1.d0)
  real(nyq_dp), parameter :: pi = acos(-1._nyq_dp)
  integer, parameter, public :: nyq_ok=0, nyq_bad_input=1, nyq_bad_value=2, &
       nyq_evaluator_failed=3, nyq_budget=4, nyq_depth=5, nyq_no_agreement=6

  type, public :: nyq_path
    ! Ordered polyline vertices. Closed paths implicitly connect last to first.
    complex(nyq_dp), allocatable :: vertices(:)
    logical :: closed=.false.
  end type

  type, public :: nyq_options
    integer :: initial_intervals=8 ! per geometric segment
    integer :: max_depth=24, max_evaluations=200000, max_passes=6
    real(nyq_dp) :: phase_limit=0.7_nyq_dp, curvature_limit=0.15_nyq_dp
    real(nyq_dp) :: log_magnitude_limit=1.5_nyq_dp, local_ratio=0.05_nyq_dp
    real(nyq_dp) :: phase_tolerance=1.e-8_nyq_dp, integer_tolerance=1.e-8_nyq_dp
  end type

  type, public :: nyq_result
    logical :: converged=.false., integer_accepted=.false.
    integer :: status=nyq_bad_input, evaluator_status=0, winding=0
    integer :: evaluations=0, passes=0, refinements=0, deepest_level=0
    real(nyq_dp) :: phase=0, raw_winding=0, max_increment=0
    complex(nyq_dp) :: endpoints(2)=(0._nyq_dp,0._nyq_dp)
    real(nyq_dp), allocatable :: pass_phases(:)
    complex(nyq_dp), allocatable :: z(:), values(:)
    real(nyq_dp), allocatable :: cumulative_phase(:)
  end type

  abstract interface
    subroutine nyq_evaluator(z, value, status)
      import nyq_dp
      complex(nyq_dp), intent(in) :: z
      complex(nyq_dp), intent(out) :: value
      integer, intent(out) :: status ! zero means success
    end subroutine
  end interface
  public :: nyq_evaluator, nyq_rectangle, nyq_horizontal, nyq_accumulate, nyq_evaluate, nyq_write_ascii
contains
  logical function finite_complex(z)
    complex(nyq_dp), intent(in) :: z
    finite_complex=ieee_is_finite(real(z,nyq_dp)).and.ieee_is_finite(aimag(z))
  end function

  logical function usable(z)
    complex(nyq_dp), intent(in) :: z
    usable=finite_complex(z)
    if (usable) usable=max(abs(real(z,nyq_dp)),abs(aimag(z)))>0
  end function

  function unit(z) result(u)
    complex(nyq_dp), intent(in) :: z
    complex(nyq_dp) :: u
    ! Component scaling also works when abs(z) itself would overflow.
    u=z/max(abs(real(z,nyq_dp)),abs(aimag(z)))
    u=u/abs(u)
  end function

  real(nyq_dp) function increment(a,b) result(d)
    complex(nyq_dp), intent(in) :: a,b
    complex(nyq_dp) :: q
    q=unit(b)*conjg(unit(a))
    d=atan2(aimag(q),real(q,nyq_dp))
  end function

  real(nyq_dp) function log_magnitude(z) result(v)
    complex(nyq_dp), intent(in) :: z
    real(nyq_dp) :: s
    s=max(abs(real(z,nyq_dp)),abs(aimag(z)))
    v=log(s)+log(abs(z/s))
  end function

  subroutine nyq_accumulate(values, phase, status, cumulative, max_increment)
    complex(nyq_dp), intent(in) :: values(:)
    real(nyq_dp), intent(out) :: phase
    integer, intent(out) :: status
    real(nyq_dp), intent(out), optional :: cumulative(:), max_increment
    real(nyq_dp) :: d, compensation, corrected, total, largest
    integer :: j
    phase=0; status=nyq_bad_input; largest=0; compensation=0
    if (present(max_increment)) max_increment=0
    if (present(cumulative)) then
      cumulative=0
      if (size(cumulative)/=size(values)) return
    end if
    if (size(values)<2) return
    status=nyq_bad_value
    do j=1,size(values)
      if (.not.usable(values(j))) return
    end do
    do j=2,size(values)
      d=increment(values(j-1),values(j))
      largest=max(largest,abs(d))
      corrected=d-compensation
      total=phase+corrected
      compensation=(total-phase)-corrected
      phase=total
      if (present(cumulative)) cumulative(j)=phase
    end do
    if (present(max_increment)) max_increment=largest
    status=nyq_ok
  end subroutine

  subroutine nyq_rectangle(xmin,xmax,ymin,ymax,path,status)
    real(nyq_dp), intent(in) :: xmin,xmax,ymin,ymax
    type(nyq_path), intent(out) :: path
    integer, intent(out) :: status
    status=nyq_bad_input
    if (.not.all(ieee_is_finite([xmin,xmax,ymin,ymax]))) return
    if (xmin>=xmax.or.ymin>=ymax) return
    path%vertices=[cmplx(xmin,ymin,nyq_dp),cmplx(xmax,ymin,nyq_dp), &
                   cmplx(xmax,ymax,nyq_dp),cmplx(xmin,ymax,nyq_dp)]
    path%closed=.true.; status=nyq_ok
  end subroutine

  subroutine nyq_horizontal(xmin,xmax,y,path,status)
    real(nyq_dp), intent(in) :: xmin,xmax,y
    type(nyq_path), intent(out) :: path
    integer, intent(out) :: status
    status=nyq_bad_input
    if (.not.all(ieee_is_finite([xmin,xmax,y]))) return
    if (xmin>=xmax) return
    path%vertices=[cmplx(xmin,y,nyq_dp),cmplx(xmax,y,nyq_dp)]
    status=nyq_ok
  end subroutine

  subroutine nyq_evaluate(evaluate,path,options,result,save_samples)
    procedure(nyq_evaluator) :: evaluate ! deliberately neither pure nor concurrent
    type(nyq_path), intent(in) :: path
    type(nyq_options), intent(in) :: options
    type(nyq_result), intent(out) :: result
    logical, intent(in), optional :: save_samples
    complex(nyq_dp), allocatable :: zs(:),fs(:)
    real(nyq_dp), allocatable :: phases(:)
    complex(nyq_dp) :: za,zb,fa,fb,first_f,left,right
    real(nyq_dp) :: t,limit,old_phase,older_phase
    integer :: n,pass,edge,j,nedges,nout,stat,nvertices
    logical :: keep,seeded

    keep=.false.
    if (present(save_samples)) keep=save_samples
    if (.not.allocated(path%vertices)) return
    nvertices=size(path%vertices)
    if (nvertices<2) return
    do j=1,nvertices
      if (.not.finite_complex(path%vertices(j))) return
    end do
    if (options%initial_intervals<1.or.options%max_depth<0.or.options%max_evaluations<3.or. &
        options%max_passes<3) return
    if (.not.all(ieee_is_finite([options%phase_limit,options%curvature_limit, &
        options%log_magnitude_limit,options%local_ratio,options%phase_tolerance,options%integer_tolerance]))) return
    if (options%phase_limit<=0.or.options%phase_limit>=pi.or.options%curvature_limit<=0.or. &
        options%log_magnitude_limit<=0.or.options%local_ratio<=0.or.options%local_ratio>=1.or. &
        options%phase_tolerance<=0.or.options%integer_tolerance<=0) return
    nedges=nvertices-1
    if (path%closed) nedges=nvertices
    do edge=1,nedges
      if (path%vertices(edge)==path%vertices(mod(edge,nvertices)+1)) return
    end do
    allocate(zs(options%max_evaluations+1),fs(options%max_evaluations+1), &
             phases(options%max_evaluations+1),result%pass_phases(options%max_passes))
    result%pass_phases=0
    n=options%initial_intervals; old_phase=0; older_phase=0; nout=0
    result%status=nyq_no_agreement
    do pass=1,options%max_passes
      result%passes=pass; stat=nyq_ok; nout=0
      result%phase=0; result%raw_winding=0; result%max_increment=0
      result%endpoints=cmplx(0,0,nyq_dp)
      limit=options%phase_limit
      if (pass>1) limit=0.7_nyq_dp*limit
      ! Odd passes >=3 use a non-dyadic, nonuniform grid, independent of midpoint subdivision.
      seeded=pass>=3.and.mod(pass,2)==1
      za=path%vertices(1)
      call fetch(za,fa)
      if (stat/=nyq_ok) exit
      first_f=fa
      call append(za,fa)
      do edge=1,nedges
        left=path%vertices(edge); right=path%vertices(mod(edge,nvertices)+1)
        do j=1,n
          t=real(j,nyq_dp)/real(n,nyq_dp)
          if (seeded.and.j<n) t=t+0.2718281828459_nyq_dp*sin(pi*t)/real(n,nyq_dp)
          ! Convex interpolation avoids overflowing right-left for large bounds.
          zb=(1-t)*left+t*right
          if (path%closed.and.edge==nedges.and.j==n) then
            fb=first_f
          else
            call fetch(zb,fb)
          end if
          if (stat/=nyq_ok) exit
          call refine(za,fa,zb,fb,0)
          if (stat/=nyq_ok) exit
          za=zb; fa=fb
        end do
        if (stat/=nyq_ok) exit
      end do
      if (stat/=nyq_ok) exit
      call nyq_accumulate(fs(:nout),result%phase,stat,phases(:nout),result%max_increment)
      if (stat/=nyq_ok) exit
      result%raw_winding=result%phase/(2*pi)
      result%endpoints=[fs(1),fs(nout)]
      result%pass_phases(pass)=result%phase
      if (pass>=3) then
        if (abs(result%phase-old_phase)<=options%phase_tolerance.and. &
            abs(result%phase-older_phase)<=options%phase_tolerance) then
          if (path%closed) then
            ! Avoid integer overflow on conversion, even for exotic user grids.
            if (abs(result%raw_winding)>=real(huge(0)-1,nyq_dp)) exit
            j=nint(result%raw_winding)
            if (abs(result%raw_winding-real(j,nyq_dp))>options%integer_tolerance) exit
            result%winding=j; result%integer_accepted=.true.
          end if
          result%converged=.true.; result%status=nyq_ok
          exit
        end if
      end if
      older_phase=old_phase; old_phase=result%phase
      if (n>options%max_evaluations/2) then
        stat=nyq_budget
        exit
      end if
      n=n*2
    end do
    if (stat/=nyq_ok) result%status=stat
    ! On failure phase/samples describe only the completed prefix of the last pass.
    if (.not.result%converged.and.nout>=2) then
      call nyq_accumulate(fs(:nout),result%phase,j,phases(:nout),result%max_increment)
      result%raw_winding=result%phase/(2*pi)
      result%endpoints=[fs(1),fs(nout)]
    end if
    if (keep.and.nout>=2) then
      result%z=zs(:nout); result%values=fs(:nout); result%cumulative_phase=phases(:nout)
    end if
  contains
    subroutine fetch(z,f)
      complex(nyq_dp), intent(in) :: z
      complex(nyq_dp), intent(out) :: f
      integer :: code
      f=cmplx(0,0,nyq_dp)
      if (result%evaluations>=options%max_evaluations) then
        stat=nyq_budget
        return
      end if
      result%evaluations=result%evaluations+1
      call evaluate(z,f,code)
      if (code/=0) then
        result%evaluator_status=code; stat=nyq_evaluator_failed
      else if (.not.usable(f)) then
        stat=nyq_bad_value
      end if
    end subroutine

    subroutine append(z,f)
      complex(nyq_dp), intent(in) :: z,f
      nout=nout+1
      zs(nout)=z; fs(nout)=f
    end subroutine

    recursive subroutine refine(a,va,b,vb,depth)
      complex(nyq_dp), intent(in) :: a,va,b,vb
      integer, intent(in) :: depth
      complex(nyq_dp) :: mid,vm,sa,sb,sm
      real(nyq_dp) :: d1,d2,direct,scale,curve,logs(3)
      logical :: split
      result%deepest_level=max(result%deepest_level,depth)
      mid=0.5_nyq_dp*a+0.5_nyq_dp*b
      if (mid==a.or.mid==b) then
        stat=nyq_depth
        return
      end if
      call fetch(mid,vm)
      if (stat/=nyq_ok) return
      d1=increment(va,vm); d2=increment(vm,vb); direct=increment(va,vb)
      scale=max(abs(real(va,nyq_dp)),abs(aimag(va)),abs(real(vm,nyq_dp)), &
                abs(aimag(vm)),abs(real(vb,nyq_dp)),abs(aimag(vb)))
      sa=va/scale; sb=vb/scale; sm=vm/scale
      curve=abs(sm-(0.5_nyq_dp*sa+0.5_nyq_dp*sb))/max(abs(sa),abs(sb),abs(sm))
      logs=[log_magnitude(va),log_magnitude(vm),log_magnitude(vb)]
      split=max(abs(d1),abs(d2))>limit.or.abs(d1+d2-direct)>options%phase_tolerance.or. &
            curve>options%curvature_limit.or. &
            max(abs(logs(2)-logs(1)),abs(logs(3)-logs(2)))>options%log_magnitude_limit.or. &
            minval(logs)-maxval(logs)<log(options%local_ratio)
      if (split) then
        if (depth>=options%max_depth) then
          stat=nyq_depth
          return
        end if
        result%refinements=result%refinements+1
        call refine(a,va,mid,vm,depth+1)
        if (stat/=nyq_ok) return
        call refine(mid,vm,b,vb,depth+1)
      else
        call append(mid,vm)
        call append(b,vb)
      end if
    end subroutine
  end subroutine

  subroutine nyq_write_ascii(filename,result,iostat)
    character(*), intent(in) :: filename
    type(nyq_result), intent(in) :: result
    integer, intent(out) :: iostat
    integer :: u,j,close_status
    iostat=nyq_bad_input
    if (.not.allocated(result%z)) return
    open(newunit=u,file=filename,status='replace',action='write',iostat=iostat)
    if (iostat/=0) return
    write(u,'(a)',iostat=iostat) '# Re(z) Im(z) Re(F) Im(F) cumulative_phase'
    do j=1,size(result%z)
      if (iostat/=0) exit
      write(u,'(5(es25.16e3,1x))',iostat=iostat) real(result%z(j),nyq_dp),aimag(result%z(j)), &
           real(result%values(j),nyq_dp),aimag(result%values(j)),result%cumulative_phase(j)
    end do
    close(u,iostat=close_status)
    if (iostat==0) iostat=close_status
  end subroutine
end module ALPS_nyquist
