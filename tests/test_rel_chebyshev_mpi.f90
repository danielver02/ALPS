program test_rel_chebyshev_mpi
  use mpi
  use alps_var, only : proc0,iproc,nproc,ierror,nspec,relativistic,ACmethod,usebM,poly_order, &
       gamma_rel,pparbar_rel,f0_rel,df0_rel,ms,vA,pi,nperp,npar,ngamma,npparbar,nspec_rel, &
       pp,f0,writeOut,profile_startup,rel_cache,arrayName,unit_error
  use alps_analyt
  use alps_com, only : pass_rel_cheb_models,rel_cheb_distribution_info
  use alps_fns_rel, only : derivative_f0_rel
  use, intrinsic :: iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter :: dp=kind(1.d0)
  type(cheb_model) :: model
  type(rel_cheb_row_info) :: info
  type(cheb_model), allocatable :: serial(:,:)
  integer :: status,code,ns,ng,calls,builds,transfers,collectives,old_calls,old_builds,old_collectives, &
       old_transfers,s,g,j,k,nmax,row,at,n,valid,failed=0,passed=0,total_failures,phase
  integer(int64) :: bytes,volume,evaluations,old_evaluations
  integer(int64), allocatable :: signature(:),minimum(:),maximum(:)
  integer(int64) :: output_bits(10),output_min(10),output_max(10)
  real(dp) :: start,construction,distribution,dt,max_distribution,evaluation_time,max_evaluation
  real(dp) :: x,gamma,h,du,r,a
  complex(dp) :: z,value,derivative,expected,expected_d,actual,plus,minus,forward,centered
  complex(dp) :: want_forward,want_centered,fu_fd,want_fd
  logical :: available
  character(32) :: mode

  call mpi_init(ierror)
  call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror)
  call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)
  proc0=iproc==0
  call get_command_argument(1,mode)
  nspec=5; nspec_rel=4
  allocate(relativistic(5),ACmethod(5),usebM(5),poly_order(5),ms(5))
  relativistic=[.true.,.false.,.true.,.true.,.true.]
  ACmethod=[2,2,1,2,2]; usebM=[.false.,.false.,.false.,.false.,.true.]
  poly_order=[1,3,2,6,6]; ms=[1.2d0,.3d0,1.d0,2.d0,1.d0]
  vA=.6d0; pi=4.d0*atan(1.d0); writeOut=.false.; profile_startup=.false.; rel_cache='off'

  do phase=1,3
     ngamma=4; npparbar=8
     if(phase==2) then
        ngamma=12; npparbar=24
     elseif(phase==3) then
        ngamma=100; npparbar=200; poly_order(1)=30; poly_order(4)=30
     endif
     if(proc0) then
        call prepare_tables(phase)
        ! Separate serial oracle, outside the measured initialization. This
        ! does not fit on workers or supply the distributed collection.
        call build_rel_cheb_rows(gamma_rel,pparbar_rel,f0_rel,[1,3,4,5],poly_order([1,3,4,5]), &
             status,selected=[.true.,.false.,.true.,.false.])
        call check(status==CHEB_SUCCESS,'serial reference construction')
        allocate(serial(4,0:ngamma))
        do g=0,ngamma
           do s=1,4
              call get_rel_cheb_row(s,g,serial(s,g),info,status)
           enddo
        enddo
     endif
     start=MPI_Wtime()
     call initialize_rel_cheb(status) ! Internally root-gated; same production hook.
     construction=MPI_Wtime()-start
     call check(status==CHEB_SUCCESS,'conditional initialization')
     call rel_cheb_collection_info(ns,ng,calls,bytes,builds,evaluations)
     if(proc0) then
        call check(builds==1.and.ns==4.and.ng==ngamma,'one root builder invocation')
        call check(calls==2*(ngamma+1)-merge(1,0,phase==1),'selected rows only')
        do g=0,ngamma
           do s=1,4
              call get_rel_cheb_row(s,g,model,info,status)
              call check(model%fit_status==serial(s,g)%fit_status,'initializer versus serial row status')
              if(.not.allocated(model%coefficients)) cycle
              call check(all(model%coefficients==serial(s,g)%coefficients).and. &
                   all(model%derivative_coefficients==serial(s,g)%derivative_coefficients), &
                   'initializer versus serial coefficients')
           enddo
        enddo
        deallocate(serial)
     else
        call check(ns==0.and.calls==0.and.builds==0,'worker initialization never fits')
     endif
     call mpi_barrier(MPI_COMM_WORLD,ierror)
     start=MPI_Wtime()
     call pass_rel_cheb_models(status)
     dt=MPI_Wtime()-start
     call mpi_reduce(dt,max_distribution,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierror)
     call check(status==CHEB_SUCCESS,'model distribution')
     call rel_cheb_collection_info(ns,ng,calls,bytes,builds,evaluations)
     call rel_cheb_distribution_info(transfers,collectives,volume)
     call check(ns==4.and.ng==ngamma.and.transfers==1.and.collectives==6,'single fixed-collective transfer')
     if(.not.proc0) call check(calls==0.and.builds==0,'worker collection imported without fitting')

     ! Independent bitwise cross-rank comparison through test-only MIN/MAX
     ! reductions, not through the production packer and not per-row broadcasts.
     nmax=maxval(poly_order([1,3,4,5])); n=4*(ngamma+1)*(11+2*(nmax+1))
     allocate(signature(n),minimum(n),maximum(n)); signature=0_int64; at=0; valid=0
     do g=0,ngamma
        do s=1,4
           call get_rel_cheb_row(s,g,model,info,status)
           signature(at+1:at+8)=int([info%species,info%state,info%kernel_status,model%requested_degree, &
                model%effective_degree,model%rank,model%sample_count,model%representation],int64)
           signature(at+9:at+11)=transfer([model%center,model%scale,info%gamma],[0_int64,0_int64,0_int64])
           if(allocated(model%coefficients)) then
              valid=valid+1
              do j=0,model%effective_degree
                 signature(at+12+j)=transfer(model%coefficients(j),0_int64)
              enddo
              do j=0,max(0,model%effective_degree-1)
                 signature(at+12+nmax+1+j)=transfer(model%derivative_coefficients(j),0_int64)
              enddo
              if(.not.proc0) call check(model%residual_norm_status==CHEB_DIAG_UNAVAILABLE.and. &
                   ieee_is_nan(model%residual_norm).and.info%input_samples==-1,'omitted diagnostics explicitly unavailable')
           endif
           at=at+11+2*(nmax+1)
        enddo
     enddo
     call mpi_allreduce(signature,minimum,n,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,ierror)
     call mpi_allreduce(signature,maximum,n,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,ierror)
     call check(all(minimum==maximum),'all required metadata/mappings/coefficients bit-identical')
     deallocate(signature,minimum,maximum)
     ! Direct cross-rank comparison of outputs and statuses as well as storage.
     ! The second row is unavailable in phase 1, exercising NaN/status transfer.
     call eval_rel_cheb(1,ngamma/2,(.15d0,-.02d0),value,status,derivative)
     output_bits(1)=int(status,int64)
     output_bits(2:5)=transfer([value,derivative],[0_int64,0_int64,0_int64,0_int64])
     call eval_rel_cheb(3,3,(.15d0,-.02d0),value,status,derivative)
     output_bits(6)=int(status,int64)
     output_bits(7:10)=transfer([value,derivative],[0_int64,0_int64,0_int64,0_int64])
     call mpi_allreduce(output_bits,output_min,10,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,ierror)
     call mpi_allreduce(output_bits,output_max,10,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,ierror)
     call check(all(output_min==output_max),'cross-rank values/derivatives/statuses bit-identical')
     if(trim(mode)=='fatal') then
        if(proc0) open(newunit=unit_error,file='fatal.log',status='new')
        value=eval_fit(4,3,(.1d0,.02d0)) ! No status: expected communicator abort.
        error stop 'unavailable continuation silently accepted'
     endif

     ! Removing the transformed tables proves neither local evaluator nor the
     ! continuation adapter reconstructs/accesses them during the hot path.
     if(proc0) deallocate(gamma_rel,pparbar_rel,f0_rel,df0_rel)
     old_calls=calls; old_builds=builds; old_collectives=collectives; old_transfers=transfers
     do g=0,ngamma
        do s=1,4
           call get_rel_cheb_row(s,g,model,info,status)
           available=info%state==REL_CHEB_FITTED.or.info%state==REL_CHEB_REDUCED
           do k=1,4
              select case(k)
              case(1); z=(0.d0,0.d0)
              case(2); z=(.2d0,-.01d0)
              case(3); z=(.8d0,.2d0)
              case(4); z=(.65d0,.02d0) ! Neighbor-row extrapolation.
              end select
              if(phase==3) z=cmplx(model%center+.3d0*model%scale,.02d0*model%scale,dp)
              call eval_rel_cheb(s,g,z,value,status,derivative)
              if(available) then
                 call cheb_evaluate(model,z,expected,code,expected_d)
                 call check(status==code.and.value==expected.and.derivative==expected_d,'local stored evaluation equivalence')
                 actual=eval_fit(info%species,g,z,dfpar=expected_d,status=code)
                 call check(code==status.and.actual==value.and.expected_d==derivative,'production continuation adapter equivalence')
                 ! Legacy signature succeeds for valid rows as well.
                 actual=eval_fit(info%species,g,z)
                 call check(actual==value,'source-compatible continuation call')
              else
                 call check(status==info%kernel_status.and.ieee_is_nan(real(value,dp)).and. &
                      ieee_is_nan(real(derivative,dp)),'unavailable local row status and NaNs')
                 if(s==1.or.s==3) then
                    actual=eval_fit(info%species,g,z,status=code)
                    call check(code==status.and.ieee_is_nan(real(actual,dp)),'unavailable continuation propagates status')
                 endif
              endif
           enddo
        enddo
     enddo

     ! Preserve the existing fixed-u Gamma and parallel finite-difference forms.
     z=(.25d0,-.03d0); du=2.d0/npparbar
     call get_rel_cheb_row(1,1,model,info,status); h=info%gamma
     call get_rel_cheb_row(1,2,model,info,status); h=info%gamma-h
     actual=eval_fit(1,1,z,status=status)
     plus=eval_fit(1,2,z,status=code)
     forward=(plus-actual)/h
     call eval_rel_cheb(1,2,z,expected,status)
     call eval_rel_cheb(1,1,z,value,status)
     want_forward=(expected-value)/h
     plus=eval_fit(1,3,z,status=status); minus=eval_fit(1,1,z,status=code)
     centered=(plus-minus)/(2.d0*h)
     call eval_rel_cheb(1,3,z,expected,status); call eval_rel_cheb(1,1,z,value,status)
     want_centered=(expected-value)/(2.d0*h)
     plus=eval_fit(1,1,z+du,status=status); minus=eval_fit(1,1,z-du,status=code)
     fu_fd=(plus-minus)/(2.d0*du)
     call eval_rel_cheb(1,1,z+du,expected,status); call eval_rel_cheb(1,1,z-du,value,status)
     want_fd=(expected-value)/(2.d0*du)
     call check(forward==want_forward.and.centered==want_centered.and.fu_fd==want_fd,'unchanged finite-difference stencil values')
     call rel_cheb_collection_info(ns,ng,calls,bytes,builds,old_evaluations)
     start=MPI_Wtime()
     do k=1,10000
        value=eval_fit(1,ngamma/2,(.15d0,-.02d0),status=status)
     enddo
     evaluation_time=MPI_Wtime()-start
     call mpi_reduce(evaluation_time,max_evaluation,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierror)
     call rel_cheb_collection_info(ns,ng,calls,bytes,builds,evaluations)
     call rel_cheb_distribution_info(transfers,collectives,volume)
     call check(calls==old_calls.and.builds==old_builds.and.transfers==old_transfers.and. &
          collectives==old_collectives.and.evaluations==old_evaluations+10000_int64,'repeated evaluations do not fit or communicate')
     if(proc0) then
        write(*,'(a,3i5,a,3es12.4)') 'MPI timing ranks/phase/rows=',nproc,phase,valid, &
             ' construction/transfer/10000_evaluations_seconds=',construction,max_distribution,max_evaluation
        write(*,'(a,3i5,a,3i12)') 'MPI root builders/fits/collectives=',builds,calls,collectives, &
             ' transfer_bytes/coefficient_bytes/local_evaluations=',volume,bytes,evaluations
     endif
  enddo

  ! Failure is broadcast before any payload and leaves no rank with a model.
  if(proc0) then
     ngamma=4; npparbar=8
     call prepare_tables(1)
     f0_rel(3,1,4)=-2.d0
  endif
  call initialize_rel_cheb(status)
  if(proc0) call check(status==CHEB_NEGATIVE_SAMPLES,'root unexpected construction failure')
  call pass_rel_cheb_models(status)
  call check(status==CHEB_NEGATIVE_SAMPLES,'root failure agreed by every rank')
  call rel_cheb_collection_info(ns,ng,calls,bytes,builds,evaluations)
  call check(ns==0.and.bytes==0,'failed initialization clears every rank')
  call eval_rel_cheb(1,0,(0.d0,0.d0),value,status)
  call check(status==CHEB_INVALID_MODEL.and.ieee_is_nan(real(value,dp)),'failed initialization cannot evaluate')
  ! No selection must require neither transformed data nor polynomial setup.
  ACmethod=1
  call initialize_rel_cheb(status)
  call pass_rel_cheb_models(status)
  call rel_cheb_collection_info(ns,ng,calls,bytes,builds,evaluations)
  call rel_cheb_distribution_info(transfers,collectives,volume)
  call check(status==CHEB_SUCCESS.and.ns==0.and.builds==0.and.transfers==0.and.collectives==1,'unselected path stays empty')
  call mpi_allreduce(failed,total_failures,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)
  write(*,'(a,i0,a,i0,a,i0)') 'MPI checks rank=',iproc,' passed=',passed,' failed=',failed
  if(total_failures/=0) call mpi_abort(MPI_COMM_WORLD,99,ierror)
  if(proc0) write(*,'(a)') 'PASS: MPI construction, distribution, continuation and lifecycle'
  call mpi_finalize(ierror)
contains
  subroutine prepare_tables(kind)
    integer, intent(in) :: kind
    integer :: s,i,j,original
    real(dp) :: u,g,r,b
    if(allocated(gamma_rel)) deallocate(gamma_rel,pparbar_rel,f0_rel,df0_rel)
    allocate(gamma_rel(4,0:ngamma,0:npparbar),pparbar_rel(4,0:ngamma,0:npparbar), &
         f0_rel(4,0:ngamma,0:npparbar),df0_rel(4,0:ngamma,0:npparbar,2))
    gamma_rel=1.d0; pparbar_rel=0.d0; f0_rel=-1.d0; df0_rel=0.d0
    if(kind==2) then
       nperp=4; npar=8; arrayName='mpi-transformed'
       allocate(pp(5,0:nperp,0:npar,2),f0(5,0:nperp,0:npar)); pp=0.d0; f0=1.d0
       do s=1,3,2
          original=1
          if(s==3) original=4
          do i=0,nperp
             do j=0,npar
                r=1.2d0*i/nperp; u=-1.d0+2.d0*j/npar; g=sqrt(1.d0+r*r+u*u)
                pp(original,i,j,1)=r*ms(original)/vA; pp(original,i,j,2)=u*ms(original)/vA
                if(s==1) then
                   f0(original,i,j)=exp(-2.d0*g+.3d0*u)
                else
                   f0(original,i,j)=exp(-2.d0*g-.15d0*u*u)
                endif
             enddo
          enddo
          call derivative_f0_rel(original,s)
       enddo
       deallocate(pp,f0)
       return
    endif
    do i=0,ngamma
       do j=0,npparbar
          gamma=1.02d0+.2d0*i
          if(kind==3) gamma=1.d0+dble(i)/ngamma
          u=-1.d0+2.d0*j/npparbar
          gamma_rel(:,i,j)=gamma; pparbar_rel(:,i,j)=u
          f0_rel(1,i,j)=exp(-2.d0*gamma+.3d0*u)
          f0_rel(3,i,j)=exp(-2.d0*gamma-.15d0*u*u)
          if(u*u>gamma*gamma-1.d0) then
             f0_rel(1,i,j)=-1.d0; f0_rel(3,i,j)=-1.d0
          endif
       enddo
    enddo
    if(kind==1) then
       f0_rel(1,ngamma,:)=-1.d0
       do j=0,npparbar
          if(f0_rel(3,3,j)/=-1.d0) f0_rel(3,3,j)=0.d0
       enddo
    endif
  end subroutine prepare_tables
  subroutine check(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if(condition) then
       passed=passed+1
    else
       failed=failed+1
       write(*,'(a,i0,2a)') 'FAIL rank=',iproc,' ',label
    endif
  end subroutine check
end program test_rel_chebyshev_mpi
