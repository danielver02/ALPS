! Task 0: exercise the real namelist reader without loading any distributions.
program test_input_validation
  use mpi
  use alps_var, only : iproc,nproc,ierror,proc0,unit_error,ACmethod,relativistic,param_fit
  use alps_io, only : alps_error_init,init_param
  implicit none
  character(32) :: mode
  integer :: out,is,method
  logical :: rel
  call mpi_init(ierror)
  call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror)
  call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)
  proc0=iproc==0
  call get_command_argument(2,mode)
  if(proc0) then
     open(newunit=out,file='validation.in',status='new')
     write(out,'(a)') '&system nspec=2,nroots=0,use_map=F,n_scan=0,nperp=4,npar=8,'
     write(out,'(a)') ' kperp=0.001,kpar=0.1,vA=1,arrayName="MUST_NOT_BE_READ",writeOut=F /'
     do is=1,2
        method=1; rel=.true.
        select case(trim(mode))
        case('invalid_first')
           if(is==1) method=2
        case('invalid_second')
           if(is==1) rel=.false.
           if(is==2) method=2
        case('valid_mixed')
           if(is==1) then
              method=2; rel=.false.
           endif
        case('valid_rel4','valid_rel5')
        case default
           error stop 'unknown validation fixture'
        end select
        write(out,'(a,i0)') '&spec_',is
        write(out,'(a,i0,a,l1,a)') ' nn=1,qq=1,mm=1,ff=1,AC_method=',method,',relat=',rel, &
             ',log_fit=T,use_bM=F /'
        write(out,'(a,i0)') '&poly_spec_',is
        write(out,'(a)') ' kind=1,order=2,log_max=18,derivative_type=0 /'
        write(out,'(a,i0,a)') '&ffit_',is,'_1'
        method=4
        if(trim(mode)=='valid_rel5') method=5
        write(out,'(a,i0,a)') ' fit_type_in=',method, &
             ',fit_1=1,fit_2=0.2,fit_3=0,fit_4=0,fit_5=0,perpcorr=2 /'
     enddo
     close(out)
     call alps_error_init
     call init_param
     if(.not.allocated(param_fit)) error stop 'supported input did not finish initialization'
     if(trim(mode)=='valid_mixed') then
        if(ACmethod(1)/=2.or.relativistic(1)) error stop 'nonrelativistic method changed'
        if(ACmethod(2)/=1.or..not.relativistic(2)) error stop 'relativistic method changed'
     endif
     write(*,'(a)') 'PASS: supported species selections preserved'
     close(unit_error)
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierror)
  call mpi_finalize(ierror)
end program
