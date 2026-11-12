! Analytic row-layer tests; no original-grid interpolation or continuation wrapper.
program test_rel_chebyshev_rows
  use alps_analyt
  use relativistic_references, only : juttner,anisotropic,drifting
  use, intrinsic :: iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter :: dp=kind(1.d0),ng=6,np=16
  real(dp) :: gamma(3,0:ng,0:np),u(3,0:ng,0:np),f(3,0:ng,0:np),original(3,0:ng,0:np)
  real(dp) :: badg(1,0:9,0:4),badu(1,0:9,0:4),badf(1,0:9,0:4)
  real(dp) :: repeatg(1,0:0,0:4),repeatu(1,0:0,0:4),repeatf(1,0:0,0:4)
  real(dp) :: mass(3)=[1.d0,2.d0,.6d0],c(3)=[.7d0,1.3d0,2.1d0],va=.6d0
  real(dp) :: abs_error(3),rel_error(3),deriv_error(3),complex_error(3),residual(3)
  real(dp) :: amplitude,gval,x,nan,inf,v,d,old_coefficient
  integer :: species(3)=[1,3,5],orders(5),degrees(3),status,s,g,j,k,n,ns,last,calls,initial_calls
  integer :: passed=0,failed=0,unit,ios,lines
  integer(int64) :: bytes
  type(cheb_model) :: model,copy
  type(rel_cheb_row_info) :: info
  complex(dp) :: z,value,derivative,want,dwant,fg
  logical :: modes_before(3),modes_after(3)
  character(512) :: line

  call ieee_get_halting_mode(ieee_invalid,modes_before(1))
  call ieee_get_halting_mode(ieee_divide_by_zero,modes_before(2))
  call ieee_get_halting_mode(ieee_overflow,modes_before(3))
  nan=ieee_value(0.d0,ieee_quiet_nan); inf=ieee_value(0.d0,ieee_positive_inf)
  call get_rel_cheb_row(1,0,model,info,status)
  call check(status==CHEB_INVALID_MODEL.and.info%state==REL_CHEB_NOT_BUILT,'initial collection unavailable')
  orders=[0,17,2,9,1]; degrees=orders(species) ! Existing poly_order(original species) convention.
  do s=1,3
     do g=0,ng
        do j=0,np
           gamma(s,g,j)=1.d0+.03d0*g
           u(s,g,j)=(va/mass(s))*(-1.3d0+2.6d0*j/np)
           call reference(s,gamma(s,g,j),cmplx(u(s,g,j),0.d0,dp),want,fg,dwant)
           f(s,g,j)=real(want,dp)
           if(u(s,g,j)**2>gamma(s,g,j)**2-1.d0) f(s,g,j)=-1.d0
        enddo
     enddo
  enddo
  original=f
  call build_rel_cheb_rows(gamma,u,f,species,degrees,status)
  call check(status==CHEB_SUCCESS,'analytic collection construction')
  call check(all(f==original),'builder preserves distribution and masks')
  call rel_cheb_collection_info(ns,last,calls,bytes)
  call check(ns==3.and.last==ng.and.calls==3*(ng+1).and.bytes>0,'collection bounds and fit count')
  initial_calls=calls
  abs_error=0.d0; rel_error=0.d0; deriv_error=0.d0; complex_error=0.d0; residual=0.d0
  do s=1,3
     amplitude=maxval(f(s,:,:))
     do g=0,ng
        call get_rel_cheb_row(s,g,model,info,status)
        call check(status==CHEB_SUCCESS.or.status==CHEB_REDUCED,'available analytic row')
        call check(info%species==species(s).and.info%packed_species==s.and.info%gamma_index==g,'original/packed mapping')
        call check(info%geometric_samples==count(f(s,g,:)/=-1.d0).and. &
             info%masked_samples==count(f(s,g,:)==-1.d0),'geometric counts')
        call check(model%representation==CHEB_LOG10.and.model%requested_degree==degrees(s),'representation and requested degree')
        call check(model%rank==model%effective_degree+1.and.model%retained_sample_count==info%geometric_samples, &
             'accepted rank and retained samples')
        call check(model%lower_tail_zeros==0.and.model%upper_tail_zeros==0,'no excluded analytic tail zeros')
        residual(s)=max(residual(s),model%residual_norm)
        do j=0,np
           if(f(s,g,j)==-1.d0) cycle
           call cheb_evaluate(model,u(s,g,j),v,status,d)
           call check(status==CHEB_SUCCESS,'retained-node evaluation')
           abs_error(s)=max(abs_error(s),abs(v-f(s,g,j)))
           rel_error(s)=max(rel_error(s),abs(v-f(s,g,j))/max(f(s,g,j),1.d-14*amplitude))
        enddo
        if(g==0) then
           call check(model%effective_degree==0,'cone apex singleton constant')
           call cheb_evaluate(model,(.4d0,.2d0),value,status,derivative)
           call check(status==CHEB_SUCCESS.and.derivative==(0.d0,0.d0),'singleton derivative exactly zero')
           cycle ! A singleton cannot identify nonconstant analytic derivatives.
        endif
        do k=1,3
           select case(k)
           case(1); z=cmplx(.08d0,-.06d0,dp)
           case(2); z=cmplx(.8d0,.3d0,dp) ! Beyond sampled intervals; no support cutoff.
           case(3); z=cmplx(.95d0*sqrt(gamma(s,g,0)**2-1.d0),.04d0,dp)
           end select
           call reference(s,gamma(s,g,0),z,want,fg,dwant)
           call cheb_evaluate(model,z,value,status,derivative)
           call check(status==CHEB_SUCCESS,'analytic complex extrapolation status')
           complex_error(s)=max(complex_error(s),abs(value-want)/max(abs(want),1.d-14*amplitude))
           deriv_error(s)=max(deriv_error(s),abs(derivative-dwant)/amplitude)
           call close(value,want,3.d-12,'analytic complex value')
           call close(derivative,dwant,3.d-12,'analytic physical-u derivative without mass/vA factor')
        enddo
     enddo
     call check(rel_error(s)<3.d-13.and.deriv_error(s)<3.d-12,'analytic real/derivative error budget')
     write(*,'(a,i0,a,5es13.4)') 'ROW analytic species=',species(s), &
          ' abs,rel,log_L2,du_scaled,complex_rel=',abs_error(s),rel_error(s),residual(s),deriv_error(s),complex_error(s)
  enddo
  ! Neighbor fixture: current Gamma=1.03, u=.23 is outside Gamma=1's interval.
  call get_rel_cheb_row(1,0,model,info,status)
  call cheb_evaluate(model,(.23d0,-.02d0),value,status)
  call reference(1,1.d0,(.23d0,-.02d0),want,fg,dwant)
  call close(value,want,3.d-14,'neighbor-row continuation is not masked')
  call get_rel_cheb_row(2,3,model,info,status)
  old_coefficient=model%coefficients(0); model%coefficients(0)=99.d0
  call get_rel_cheb_row(2,3,copy,info,status)
  call check(copy%coefficients(0)==old_coefficient,'accessor is an independent deep copy')
  do j=1,20
     call cheb_evaluate(copy,cmplx(.01d0*j,.02d0,dp),value,status,derivative)
  enddo
  call rel_cheb_collection_info(ns,last,calls,bytes)
  call check(calls==initial_calls,'retrieval and repeated evaluations never refit')

  ! Explicit selection keeps the original packed indexing; unselected slots do not fit.
  call build_rel_cheb_rows(gamma,u,f,species,degrees,status,selected=[.false.,.true.,.false.])
  call rel_cheb_collection_info(ns,last,calls,bytes)
  call check(status==CHEB_SUCCESS.and.calls==ng+1,'selected slots only')
  call get_rel_cheb_row(1,2,model,info,status)
  call check(status==CHEB_INVALID_MODEL.and.info%state==REL_CHEB_NOT_SELECTED.and. &
       .not.allocated(model%coefficients),'unselected row has no invented model')
  call get_rel_cheb_row(2,2,model,info,status)
  call check(status==CHEB_SUCCESS.and.info%species==3,'selected packed species retained')

  ! Expected unavailable states and actual sample failures in one collection.
  badg=1.5d0
  do g=0,9
     badu(1,g,:)=[-.4d0,-.2d0,0.d0,.2d0,.4d0]
  enddo
  badf=1.d0
  badf(1,0,:)=-1.d0; badf(1,1,:)=0.d0
  badf(1,2,:)=[-1.d0,0.d0,2.d0,0.d0,-1.d0]
  badf(1,3,:)=[1.d0,0.d0,1.d0,-1.d0,-1.d0]
  badf(1,4,:)=[1.d0,-2.d0,1.d0,-1.d0,-1.d0]
  badf(1,5,1)=nan; badf(1,6,1)=inf; badu(1,7,1)=nan; badg(1,8,1)=1.6d0
  badf(1,9,:)=[1.d0,10.d0,100.d0,-1.d0,-1.d0]
  open(newunit=unit,file='row-failures.log',status='new',action='readwrite')
  call build_rel_cheb_rows(badg,badu,badf,[9],[4],status,diagnostic_unit=unit)
  call check(status==CHEB_INTERIOR_ZERO,'collection reports first unexpected failure without abort')
  do g=0,9
     call get_rel_cheb_row(1,g,model,info,status)
     select case(g)
     case(0)
        call check(info%state==REL_CHEB_NO_GEOMETRY.and.status==CHEB_NO_SAMPLES,'no geometry classification')
     case(1)
        call check(info%state==REL_CHEB_ZERO.and.status==CHEB_NO_POSITIVE_SAMPLES,'all-zero classification')
        call check(info%geometric_samples==5.and.info%positive_samples==0,'all-zero counts')
     case(2)
        call check(info%state==REL_CHEB_REDUCED.and.status==CHEB_REDUCED.and.model%effective_degree==0,'zero-tail singleton reduction')
        call check(model%lower_tail_zeros==1.and.model%upper_tail_zeros==1.and. &
             model%original_sample_count==3.and.model%retained_sample_count==1,'mask exclusion before zero filtering')
     case(3:8)
        call check(info%state==REL_CHEB_INVALID_SAMPLES,'invalid values/coordinates classified')
        if(g==3) call check(status==CHEB_INTERIOR_ZERO,'interior zero not discarded')
        if(g==4) call check(status==CHEB_NEGATIVE_SAMPLES,'negative beyond sentinel not discarded')
        if(g==5.or.g==6) call check(status==CHEB_INVALID_VALUES,'nonfinite value reaches validation')
        if(g==7.or.g==8) call check(status==CHEB_INVALID_COORDINATES,'invalid coordinate reaches validation')
     case(9)
        call check(info%state==REL_CHEB_REDUCED.and.model%effective_degree==2,'later valid row survives earlier errors')
     end select
     if(g==0.or.g==1.or.(g>=3.and.g<=8)) then
        call check(.not.allocated(model%coefficients),'unavailable row has no coefficients')
        call cheb_evaluate(model,0.d0,v,status,d)
        call check(status==CHEB_INVALID_MODEL.and.ieee_is_nan(v).and.ieee_is_nan(d),'unavailable evaluation fails explicitly')
     endif
  enddo
  rewind(unit); lines=0
  do
     read(unit,'(a)',iostat=ios) line
     if(ios/=0) exit
     lines=lines+1
     call check(index(line,'species=9')>0.and.index(line,'packed_species=1')>0.and.index(line,'igamma=')>0.and. &
          index(line,'kernel_status=')>0.and.index(line,'unmasked=')>0.and.index(line,'positive=')>0,'failure diagnostic identifies row/status/counts')
  enddo
  close(unit)
  call check(lines==6,'every unexpected sample failure diagnosed')
  repeatg=1.5d0; repeatu(1,0,:)=[0.d0,.2d0,0.d0,.2d0,0.d0]; repeatf=10.d0**repeatu
  call build_rel_cheb_rows(repeatg,repeatu,repeatf,[3],[4],status)
  call get_rel_cheb_row(1,0,model,info,status)
  call check(status==CHEB_REDUCED.and.model%rank_reduced.and.model%effective_degree==1.and.model%rank==2, &
       'rank-limited row is refitted')
  call check(model%fit_attempts==2.and..not.model%sample_limited,'rank versus sample reduction reason')
  call cheb_evaluate(model,.1d0,v,status,d)
  call check(abs(d-log(10.d0)*10.d0**.1d0)<1.d-13,'rank-reduced physical derivative')
  repeatu(1,0,:)=[1.d0,nearest(1.d0,1.d0),1.d0,nearest(1.d0,1.d0),1.d0]
  repeatf(1,0,:)=[1.d0,10.d0,1.d0,10.d0,1.d0]
  call build_rel_cheb_rows(repeatg,repeatu,repeatf,[3],[1],status)
  call get_rel_cheb_row(1,0,model,info,status)
  call check(status==CHEB_SUCCESS.and..not.model%degenerate_interval,'narrow nondegenerate physical interval')
  call cheb_evaluate_polynomial(model,1.d0,v,status,d)
  call check(abs(d*spacing(1.d0)-1.d0)<1.d-13,'narrow-coordinate derivative scaling')
  repeatu=0.d0; repeatf=10.d0
  call build_rel_cheb_rows(repeatg,repeatu,repeatf,[3],[5],status)
  call get_rel_cheb_row(1,0,model,info,status)
  call check(status==CHEB_REDUCED.and.model%degenerate_interval.and.model%effective_degree==0,'degenerate row constant')
  repeatf=nearest(0.d0,1.d0)
  call build_rel_cheb_rows(repeatg,repeatu,repeatf,[3],[0],status)
  call get_rel_cheb_row(1,0,model,info,status)
  call check(status==CHEB_SUCCESS.and.info%positive_samples==5.and.model%retained_sample_count==5,'positive subnormals retained')
  repeatu(1,0,:)=[.4d0,-.4d0,0.d0,.2d0,-.2d0]; repeatf=1.d-250*10.d0**repeatu
  call build_rel_cheb_rows(repeatg,repeatu,repeatf,[3],[1],status)
  call get_rel_cheb_row(1,0,model,info,status)
  call cheb_evaluate(model,.1d0,v,status,d)
  call check(status==CHEB_SUCCESS.and.abs(v/(1.d-250*10.d0**.1d0)-1.d0)<1.d-12,'low-amplitude relative reconstruction')
  call check(abs(d/(log(10.d0)*1.d-250*10.d0**.1d0)-1.d0)<1.d-12,'low-amplitude relative derivative')

  ! Reinitialization replaces shapes and clears stale models, including on invalid input.
  call get_rel_cheb_row(1,1,model,info,status)
  call check(status==CHEB_INVALID_MODEL,'old Gamma bounds not retained')
  call get_rel_cheb_row(0,0,model,info,status)
  call check(status==CHEB_INVALID_MODEL,'invalid packed index')
  call get_rel_cheb_row(1,-1,model,info,status)
  call check(status==CHEB_INVALID_MODEL,'invalid negative row index')
  call build_rel_cheb_rows(gamma,u,f,[1,3],degrees,status)
  call check(status==CHEB_INVALID_VALUES,'invalid metadata shape')
  call rel_cheb_collection_info(ns,last,calls,bytes)
  call check(ns==0.and.calls==0.and.bytes==0,'invalid rebuild releases old collection')
  call build_rel_cheb_rows(gamma,u,f,[1,1,5],degrees,status)
  call check(status==CHEB_INVALID_VALUES,'duplicate species mapping rejected')
  call build_rel_cheb_rows(repeatg,repeatu,repeatf,[3],[-1],status)
  call check(status==CHEB_INVALID_DEGREE,'invalid requested degree reported')
  call get_rel_cheb_row(1,0,model,info,status)
  call check(info%state==REL_CHEB_FIT_FAILURE.and..not.allocated(model%coefficients),'configuration failure classification')
  call release_rel_cheb_rows
  call release_rel_cheb_rows
  call get_rel_cheb_row(1,0,model,info,status)
  call check(status==CHEB_INVALID_MODEL.and.info%state==REL_CHEB_NOT_BUILT,'idempotent release')

  call benchmark(100,200)
  call benchmark(500,500)
  call release_rel_cheb_rows
  call ieee_get_halting_mode(ieee_invalid,modes_after(1))
  call ieee_get_halting_mode(ieee_divide_by_zero,modes_after(2))
  call ieee_get_halting_mode(ieee_overflow,modes_after(3))
  call check(all(modes_before.eqv.modes_after),'row builder preserves IEEE halting modes')
  write(*,'(a,i0,a,i0)') 'Relativistic Chebyshev rows: passed=',passed,', failed=',failed
  if(failed/=0) error stop 1
contains
  subroutine benchmark(ngrid,nparallel)
    integer, intent(in) :: ngrid,nparallel
    real(dp), allocatable :: bg(:,:,:),bu(:,:,:),bf(:,:,:)
    real(dp) :: start,finish
    integer :: i,j
    allocate(bg(1,0:ngrid,0:nparallel),bu(1,0:ngrid,0:nparallel),bf(1,0:ngrid,0:nparallel))
    do i=0,ngrid
       do j=0,nparallel
          bg(1,i,j)=1.d0+dble(i)/ngrid; bu(1,i,j)=-1.d0+2.d0*j/nparallel
          bf(1,i,j)=exp(-2.d0*bg(1,i,j)-.15d0*bu(1,i,j)**2)
          if(bu(1,i,j)**2>bg(1,i,j)**2-1.d0) bf(1,i,j)=-1.d0
       enddo
    enddo
    call cpu_time(start)
    call build_rel_cheb_rows(bg,bu,bf,[1],[30],status)
    call cpu_time(finish)
    call check(status==CHEB_SUCCESS,'benchmark construction without unexpected failures')
    call rel_cheb_collection_info(ns,last,calls,bytes)
    call check(ns==1.and.last==ngrid.and.calls==ngrid+1,'one fit call per benchmark row')
    call get_rel_cheb_row(1,ngrid/2,model,info,status)
    call cheb_evaluate(model,(.15d0,.04d0),value,status,derivative)
    call anisotropic(bg(1,ngrid/2,0),(.15d0,.04d0),2.d0,.15d0,1.d0,want,fg,dwant)
    call close(value,want,1.d-11,'benchmark complex reconstruction')
    call close(derivative,dwant,1.d-10,'benchmark derivative')
    write(*,'(a,2i5,a,f10.5,a,i0,a,i0)') 'ROW benchmark intervals=',ngrid,nparallel, &
         ' build_CPU_seconds=',finish-start,' coefficient_bytes=',bytes,' approximate_stored_bytes=', &
         bytes+int(ngrid+1,int64)*((storage_size(model)+storage_size(info))/8)
  end subroutine
  subroutine reference(s,g,z,f,fg,fu)
    integer, intent(in) :: s
    real(dp), intent(in) :: g
    complex(dp), intent(in) :: z
    complex(dp), intent(out) :: f,fg,fu
    select case(s)
    case(1); call juttner(g,z,.8d0,c(s),f,fg,fu)
    case(2); call anisotropic(g,z,2.1d0,.35d0,c(s),f,fg,fu)
    case(3); call drifting(g,z,.7d0,.3d0,c(s),f,fg,fu)
    end select
  end subroutine
  subroutine close(actual,expected,tolerance,label)
    complex(dp), intent(in) :: actual,expected
    real(dp), intent(in) :: tolerance
    character(*), intent(in) :: label
    if(.not.ieee_is_finite(real(actual,dp)).or..not.ieee_is_finite(aimag(actual))) then
       call check(.false.,label)
    else
       call check(abs(actual-expected)<=tolerance*max(1.d0,abs(expected)),label)
    endif
  end subroutine
  subroutine check(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if(condition) then
       passed=passed+1
    else
       failed=failed+1; write(*,'(a)') 'FAIL: '//label
    endif
  end subroutine
end program test_rel_chebyshev_rows
