!Copyright (c) 2023, Kristopher G. Klein and Daniel Verscharen
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
!I                              Version 1.0                                    I
!I                                                                             I
!I  Kristopher Klein   (kgklein@arizona.edu)                                   I
!I  Daniel Verscharen  (d.verscharen@ucl.ac.uk)                                I
!I                                                                             I
!===============================================================================

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
  use alps_check,  only : check_parameters
  use mpi
  implicit none

  integer :: ik
  !! Index for iterating through wavevector scans with [[om_scan(subroutine)]]
  !! or [[om_double_scan(subroutine)]]

  !Initialize MPI message passing:
  call mpi_init (ierror)
  call mpi_comm_size (mpi_comm_world, nproc, ierror)
  call mpi_comm_rank (mpi_comm_world, iproc, ierror)

  !Set logical proc0=true if iproc=0:
  proc0= (iproc == 0)

  if (proc0) then
     call alps_error_init !alps_io
     write(*,'(a)')'Starting ALPS==================================='
     write(*,'(a)') "Time:"
     call output_time !alps_io
     write(*,'(a)')'================================================'
     call display_credits !alps_io
  endif

  call mpi_barrier(mpi_comm_world,ierror)
  if (proc0) write(*,'(a)') 'All processes are up and running.'

  !Check to be sure nproc is even and greater than 2, otherwise shutdown:
  if ((mod(nproc,2) .ne. 0).or.(nproc .le. 2)) call alps_error(0) !alps_io

  !Read parameters:
  if (proc0) call init_param !alps_io

  !Split Problem Amongst Processors.
  !Pass relevant information from readin:
  call pass_instructions !alps_com

  if (proc0) then

     !Allocate background distribution function f0:
     allocate(f0(1:nspec,0:nperp,0:npar)); f0=0.d0

     !Read in f0:
     call read_f0 !alps_io

     !Calculate pperp, ppar derivatives of f0:
     call derivative_f0 !alps_fns
     !Calculate best fit to f0:
     call determine_param_fit !alps_analyt

     !f0 not needed for dispersion calculation.
     !Deallocate to save space:
     deallocate(f0)

     !Check Input Parameters and Distributions:
     call check_parameters !alps_check

  endif

  !Distribute input and derived parameters:
  call pass_distribution ! alps_com

  ! Once we know kperp, we can determine nmax and split the processes.
  ! The following three routines will also be called when kperp changes:
  call determine_nmax ! alps_fns
  call split_processes ! alps_fns
  ! All processes determine their Bessel function array:
  if(.NOT.(sproc.EQ.0)) call determine_bessel_array ! alps_fns

  call mpi_barrier(mpi_comm_world,ierror)

  !Either:
  !use_map=.true.  :
  !    use a scan over (omega,gamma) to local dispersion solutions:
  !OR
  !use_map=.false. :
  !    use user input roots as initial guesses for dispersion solutions:
  if (use_map) then
     if (writeOut.and.proc0) &
          write(*,'(a)')'Map Search'
     call map_search !alps_fns
  else
     if (writeOut.and.proc0) &
          write(*,'(a)') 'Starting Secant Method'
     call refine_guess !alps_fns
  endif

  if (n_scan.gt.0) then !setting n_scan=0 turns off wavevector scanning
     select case (scan_option)
     case (1)
        !scan along perscribed paths in wavevector space:
        do ik = 1, n_scan
           call om_scan(ik) !alps_fns
        enddo
     case(2)
        !scan along a plane in wavevector space:
        if (n_scan==2) then
           call om_double_scan !alps_fns
        else
           call alps_error(4) !alps_io
        endif
     case default
        !scan_option not selected
        call alps_error(3) !alps_io
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
