!=============================================================================!
!=============================================================================!
!!*ALPS                                                                  *!!
!!                                                                           !!
!!Kristopher Klein, Daniel Verscharen                                        !!
!!
!!Space Science Center, University of New Hampshire                          !!
!!                                                                           !!
!!*MAIN PROGRAM                                                             *!!
!!
!! Current Version: 2016-05-13
!=============================================================================!
!=============================================================================!
!

program alps
  use alps_var,    only : nproc, iproc, ierror, proc0,sproc
  use alps_var,    only : use_map, unit_error, scan_option
  use alps_var,    only : writeOut, f0, nperp, npar, nspec,n_scan
  use alps_io,     only : init_param, read_f0, get_unused_unit
  use alps_io,     only : alps_error, alps_error_init, output_time, display_credits
  use alps_fns,    only : derivative_f0, disp, secant
  use alps_fns,    only : refine_guess
  use alps_fns,    only : determine_nmax, split_processes, determine_bessel_array
  use alps_fns,    only : map_search
  use alps_fns,    only : om_scan, om_double_scan
  use alps_com,    only : pass_instructions, pass_distribution
  use alps_analyt, only : determine_param_fit

  implicit none
  include 'mpif.h'      !Include MPI library variables
  !Local
  integer :: ik
  
  !Initialize MPI message passing---------------------------------------------
  call mpi_init (ierror)
  call mpi_comm_size (mpi_comm_world, nproc, ierror)
  call mpi_comm_rank (mpi_comm_world, iproc, ierror)

  !Set logical proc0=true if iproc=0
  proc0= (iproc == 0)

  if (proc0) then
     call alps_error_init
     write(*,'(a)')'Starting ALPS===================================' 
     write(*,'(a)') "Time:"
     call output_time
     write(*,'(a)')'================================================' 
     call display_credits
  endif
  	
  call mpi_barrier(mpi_comm_world,ierror)
  if (proc0) write(*,'(a)') 'All processes are up and running.'	

  !Check to be sure nproc is even, otherwise shutdown
  if (mod(nproc,2) .ne. 0) call alps_error(0)

  !Read parameters------------------------------------------------------------
  if (proc0) call init_param !alps_io
       
  !Split Problem Amongst Processors
  !Pass relevant information from readin
  call pass_instructions !alps_com

  if (proc0) then     

     !Allocate background distribution function f0
     allocate(f0(1:nspec,0:nperp,0:npar)); f0=0.d0

     !Read in f0
     call read_f0 !alps_io

     !Calculate pperp, ppar derivatives of f0
     call derivative_f0 !alps_fns

     call determine_param_fit
     
     !f0 not needed for dispersion calculation
     !Deallocate to save space.
     deallocate(f0)

  endif

  call pass_distribution

  ! Once we know kperp, we can determine nmax and split the processes:
  ! The following three routines will also be called when kperp changes
  call determine_nmax
  call split_processes

  ! All processes determine their Bessel function array:
  if(.NOT.(sproc.EQ.0)) call determine_bessel_array
  
  call mpi_barrier(mpi_comm_world,ierror)
  
  !Either:
  !use_map=.true.  : 
  !    use a scan over (omega,gamma) to local dispersion solutions
  !OR
  !use_map=.false. : 
  !    use user input roots as initial guesses for dispersion solutions
  if (use_map) then
     if (writeOut.and.proc0) &
          write(*,'(a)')'Map Search'           
     call map_search     
  else       
     if (writeOut.and.proc0) &
          write(*,'(a)') 'Starting Secant Method'     
     call refine_guess
  endif

  if (n_scan.gt.0) then !setting n_scan=0 turns off wavevector scanning
     select case (scan_option)
     case (1) !scan along perscribed paths in wavevector space
        do ik = 1, n_scan
           call om_scan(ik)
        enddo
     case(2)
        if (n_scan==2) then
           call om_double_scan
        else
           call alps_error(4) !
        endif
     case default
        call alps_error(3) !scan_option not selected
     end select
  endif


  !Finalize MPI message passing
  call mpi_finalize (ierror)
  if (proc0) then
     write(*,'(a)')'Finishing ALPS===================================' 
     write(*,'(a)') "Time:"
     call output_time
     write(*,'(a)')'================================================' 
     write(unit_error,'(a)')'Run completed without fatal error.'
     close(unit_error)
  endif
end program alps
