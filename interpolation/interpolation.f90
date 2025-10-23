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

program interpolate
  !! This program includes the interpolation routine used by ALPS to fill a grid in momentum space.
  !! It is based on the polyharmonic spline algorithm as described in the ALPS code paper.
implicit none

double precision, allocatable, dimension (:) :: grid_coarse
!! Coarse input grid for interpolation.
!! (1:n_coarse)

double precision, allocatable, dimension (:) :: pperp_coarse
!! Coordinates of perpendicular momentum on coarse grid.
!! (1:n_coarse)

double precision, allocatable, dimension (:) :: ppar_coarse
!! Coordinates of parallel momentum on coarse grid.
!! (1:n_coarse)

double precision, allocatable, dimension (:,:) :: grid_fine
!! Fine output grid after interpolation.
!! (0:nperp,0:npar)

double precision, allocatable, dimension (:,:) :: pperp
!! Coordinates of perpendicular momentum on fine output grid.
!! (0:nperp,0:npar)

double precision, allocatable, dimension (:,:) :: ppar
!! Coordinates of parallel momentum on fine output grid.
!! (0:nperp,0:npar)

double precision :: smoothing
!! Smoothing parameter for spline interpolation.

double precision :: pperp_min
!! Minimum perpendicuar momentum.

double precision :: pperp_max
!! Maximum perpendicular momentum.

double precision :: ppar_min
!! Minimum parallel momentum.

double precision :: ppar_max
!! Maximum parallel momentum.

double precision :: threshold
!! Lower treshold for f0-values (coarse grid) to be included.

double precision :: mult_pperp
!! Scaling factor for perpendicular momentum.

double precision :: mult_ppar
!! Scaling factor for parallel momentum.

double precision :: mult_f
!! Scaling factor for f0.

double precision :: r_pperp
!! Read-in variable for perpendicular momentum.

double precision :: r_ppar
!! Read-in variable for parallel momentum.

double precision :: r_f
!! Read-in variable for f0.

double precision :: pperp_max_set
!! Forced maximum perpendicular momentum for output.

double precision :: pperp_min_set
!! Forced minimum perpendicular momentum for output.

double precision :: ppar_max_set
!! Forced maximum parallel momentum for output.

double precision :: ppar_min_set
!! Forced minimum perpendicular momentum for output.

integer :: nperp
!! Number of perpendicular steps on fine output grid.

integer :: npar
!! Number of parallel steps on fine output grid.

integer :: n_coarse
!! Number of entries in coarse grid.

integer :: mode
!! Format of input grid (order pperp/ppar).

integer :: i
!! Index for loops.

integer :: j
!! Index for loops.

integer :: io_error
!! Error flag for i/o.

integer :: status_read
!! Status flag for i/o.

integer :: i_coarse
!! Index to loop over the coarse grid entries.

character(LEN=256) :: filename
!! File name of input file for interpolation.

character(LEN=256) :: output_file
!! File name of output file in ALPS distribution format.

logical :: out_to_file
!! Check whether output should be written to file.

logical :: do_normalize
!! Check whether normalisation should be applied.

!I/O values for namelist:
integer :: unit
!! Unit variable for opening namelist.

integer, parameter :: stdout_unit=6
!! Stdout unit for opening namelist.

integer, save :: input_unit_no
!! Unit index for opening namelist.

integer, save :: error_unit_no=stdout_unit
!! Error unit for opening namelist.

character(500) :: runname
!! String for parameter input file.

character(500) :: foldername
!! String for parameter input folder.


call read_in_params

! Determine length of input file and min/max values:
n_coarse=0
pperp_max=-999999.d0
pperp_min=999999.d0
ppar_max=-999999.d0
ppar_min=999999.d0
open(unit=10,file=filename,status='old',action='read',iostat=io_error)
   do
   	  if (mode.EQ.1) read (10,*,iostat=status_read) r_pperp,r_ppar,r_f
   	  if (mode.EQ.2) read (10,*,iostat=status_read) r_ppar,r_pperp,r_f
	 if (r_f.GT.threshold) then
	     if (status_read.LT.0) exit
    	pperp_max=max(pperp_max,r_pperp)
    	pperp_min=min(pperp_min,r_pperp)
    	ppar_max=max(ppar_max,r_ppar)
    	ppar_min=min(ppar_min,r_ppar)
    	  n_coarse=n_coarse+1
      endif
  if (status_read.LT.0) exit
	enddo
close(10)

if (n_coarse.EQ.0) then
	write (*,*) "No data points in selection."
	stop
endif

pperp_max=pperp_max*mult_pperp
pperp_min=pperp_min*mult_pperp
ppar_max=ppar_max*mult_ppar
ppar_min=ppar_min*mult_ppar

if (pperp_max_set.NE.-9999.d0) pperp_max=pperp_max_set
if (ppar_max_set.NE.-9999.d0) ppar_max=ppar_max_set
if (pperp_min_set.NE.-9999.d0) pperp_min=pperp_min_set
if (ppar_min_set.NE.-9999.d0) ppar_min=ppar_min_set

write (*,'(a,i12)')    "Number of points in the coarse grid: ",n_coarse
write (*,*) " "
write (*,'(a)') "Properties of the fine grid:"
write (*,'(a,i12)')     "Number of grid points in nperp:      ",nperp
write (*,'(a,i12)')     "Number of grid points in npar:       ",npar
write (*,'(a,es12.4)') "Maximum Pperp:                       ", pperp_max
write (*,'(a,es12.4)') "Minimum Pperp:                       ", pperp_min
write (*,'(a,es12.4)') "Maximum Ppar:                        ", ppar_max
write (*,'(a,es12.4)') "Minimum Ppar:                        ", ppar_min
write (*,*) " "


! Allocate grids:
allocate(grid_coarse(n_coarse))
allocate(pperp_coarse(n_coarse))
allocate(ppar_coarse(n_coarse))
allocate(grid_fine(0:nperp,0:npar))
allocate(pperp(0:nperp,0:npar))
allocate(ppar(0:nperp,0:npar))


! Read in file:
i_coarse=1
open(unit=10,file=filename,status='old',action='read',iostat=io_error)
   do
   	  if (mode.EQ.1) read (10,*,iostat=status_read) r_pperp,r_ppar,r_f
   	  if (mode.EQ.2) read (10,*,iostat=status_read) r_ppar,r_pperp,r_f
      if (status_read.LT.0) exit

      if (r_f.GT.threshold) then
		pperp_coarse(i_coarse)=r_pperp
		ppar_coarse(i_coarse)=r_ppar
		grid_coarse(i_coarse)=log(r_f)
       i_coarse=i_coarse+1
      endif

	enddo
close(10)

pperp_coarse=pperp_coarse*mult_pperp
ppar_coarse=ppar_coarse*mult_ppar


! Create the pperp and ppar grids for the fine grid:
do i=0,nperp
do j=0,npar
	 pperp(i,j)=pperp_min+(pperp_max-pperp_min)*i/(1.d0*nperp)
	 ppar(i,j)=ppar_min+(ppar_max-ppar_min)*j/(1.d0*npar)
enddo
enddo



call polyharmonic_spline(grid_coarse,pperp_coarse,ppar_coarse,n_coarse,pperp,ppar,nperp,npar,smoothing,grid_fine)

if (do_normalize) call normalize (pperp,ppar,grid_fine,nperp,npar)


! Output:
if (out_to_file) then
	output_file=trim(runname)//".array"

	write (*,'(a,a)') "Writing output to file ", output_file
	open(unit=20,file=output_file,status='unknown',action='write',iostat=io_error)
	do i=0,nperp
	do j=0,npar
		write (20,*) pperp(i,j),ppar(i,j),mult_f*exp(grid_fine(i,j))
	enddo
	enddo
	close(20)
else
	do i=0,nperp
	do j=0,npar
		write (*,*) pperp(i,j),ppar(i,j),mult_f*exp(grid_fine(i,j))
	enddo
	enddo
endif


contains

	subroutine read_in_params
		!! This subroutine reads in system parameters input file (namelist) as argument after executable:
		!!  `./interpolate input.in`
		implicit none

		nameList /system/ &
				 filename, smoothing, nperp, npar, ppar_min_set, &
				 ppar_max_set, pperp_min_set, pperp_max_set, mode, &
				 threshold, out_to_file, do_normalize, mult_pperp, &
				 mult_ppar, mult_f

		nperp=50
		npar=100
		smoothing=1.d0
		mode=1
		pperp_max_set=-9999.d0
		ppar_max_set=-9999.d0
		pperp_min_set=-9999.d0
		ppar_min_set=-9999.d0
		filename=""
		out_to_file=.FALSE.
		do_normalize=.FALSE.
		threshold=0.d0
		mult_pperp=1.d0
		mult_ppar=1.d0
		mult_f=1.d0


		call get_unused_unit (input_unit_no)
		call get_runname(runname,foldername)
		runname=trim(runname)//".in"
		unit=input_unit_no
		open (unit=unit,file=runname,status='old',action='read')
		read (unit=unit,nml=system)

		unit=input_unit_no
		close (unit)

	end subroutine read_in_params



	  function input_unit_exist (nml,exist)
      !! This function checks whether a unit number exists. It is taken from the AstroGK code.
	    implicit none

	    character(*), intent (in) :: nml
      !! Namelist identifier.

	    logical, intent(out) :: exist
      !! Flags whether input unit exists.

	    integer :: input_unit_exist
      !! Unit number.

      integer :: iostat
      !! Flag for i/o.

	    character(500) :: line
      !! Read-in variable for namelist entries.

	    intrinsic adjustl, trim

	    input_unit_exist = input_unit_no
	    exist = .true.
	    if (input_unit_no > 0) then
	       rewind (unit=input_unit_no)
	       do
	          read (unit=input_unit_no, fmt="(a)", iostat=iostat) line
	          if (iostat /= 0) then
	             rewind (unit=input_unit_no)
	             exit
	          end if
	          if (trim(adjustl(line)) == "&"//nml) then
	             backspace (unit=input_unit_no)
	             return
	          end if
	       end do
	    end if
	    exist = .false.
	  end function input_unit_exist



	  function input_unit (nml)
      !! This function returns a unit number for a namelist. It is taken from the AstroGK code.
      implicit none

      character(*), intent (in) :: nml
      !! Namelist identifier.

	    integer :: input_unit
      !! Unit number for namelist.

      integer :: iostat
      !! Flag for i/o.

	    character(500) :: line
      !! Read-in variable for namelist entries.

	    intrinsic adjustl, trim
	    input_unit = input_unit_no
	    if (input_unit_no > 0) then
	       rewind (unit=input_unit_no)
	       do
	          read (unit=input_unit_no, fmt="(a)", iostat=iostat) line
	          if (iostat /= 0) then
	             rewind (unit=input_unit_no)
	             exit
	          end if
	          if (trim(adjustl(line)) == "&"//nml) then
	             backspace (unit=input_unit_no)
	             return
	          end if
	       end do
	    end if
	    write (unit=error_unit_no, fmt="('Couldn''t find namelist: ',a)") nml
	    write (unit=*, fmt="('Couldn''t find namelist: ',a)") nml
	  end function input_unit



	  subroutine get_unused_unit (unit)
      !! This subroutine returns an available unit number. It is taken from the AstroGK code.

	    implicit none

	    integer, intent (out) :: unit
      !! Unit number.

	    logical :: od
      !! Check whether unit has been opened.

	    unit = 50
	    do
	       inquire (unit=unit, opened=od)
	       if (.not.od) return
	       unit = unit + 1
	    end do
	  end subroutine get_unused_unit



  subroutine get_runname(runname,foldername)
    !! Get runname for output files from input argument.
    implicit none

    integer       :: l
    !!Dummy Length.

    integer       :: pathend
    !!Directory divider.

    character(500) :: arg
    !!Input Argument.

    character(500), intent(out) :: runname
    !!Basename for file I/O.

    character(500), intent(out) :: foldername
    !!Directory in which input file is stored.

    !Get the first argument of the program execution command:
    call getarg(1,arg)
    pathend=0

    !Check if this is the input file and trim .in extension to get runname.
    !Also remove any folder structure from the runname:
    l = len_trim (arg)
    pathend = scan(arg, "/", .true.)
    if (l > 3 .and. arg(l-2:l) == ".in") then
       runname = arg(pathend+1:l-3)
       foldername = arg(1:pathend)
    end if

    end subroutine get_runname

end program





subroutine normalize (pperp,ppar,grid_fine,nperp,npar)
  !! This subroutine normalises the fine interpolation grid.
implicit none

double precision, intent(in) :: pperp(0:nperp,0:npar)
!! Coordinates of perpendicular momentum on fine grid.

double precision, intent(in) :: ppar(0:nperp,0:npar)
!! Coordinates of parallel momentum on fine grid.

double precision, intent(inout) :: grid_fine(0:nperp,0:npar)
!! Fine output grid from interpolation.

integer, intent(in) :: nperp
!! Number of perpendicular steps on fine output grid.

integer, intent(in) :: npar
!! Number of parallel steps on fine output grid.

integer :: i
!! Index to loop over perpendicular momentum.

integer :: j
!! Index to loop over parallel momentum.

double precision :: dpperp
!! Inifinitesimal step in perpendicular momentum.

double precision :: dppar
!! Inifinitesimal step in parallel momentum.

double precision :: integral
!! Integral of the distribution function.

double precision :: M_PI=2.d0*acos(0.d0)
!! Pi

dpperp=pperp(2,1)-pperp(1,1)
dppar=ppar(1,2)-ppar(1,1)
integral=0.d0
! Integrate over the full grid:
do i=0,nperp
do j=0,npar
	integral=integral+exp(grid_fine(i,j))*pperp(i,j)
enddo
enddo
integral=integral*dpperp*dppar

grid_fine=grid_fine-log(integral*2.d0*M_PI)
end subroutine


subroutine polyharmonic_spline(grid_coarse,pperp_coarse,ppar_coarse,n_coarse, &
                               pperp,ppar,nperp,npar,smoothing,grid_fine)
!! This subroutine interpolates the grid with a polyharmonic thin-plate spline.
!! Needs BLAS/LAPACK (dgetrf, dgetrs, dgemv, ddot). We solve A*[w;a] = [y;0].
!! Then we optionally estimate the effective DoF via Hutchinson:
!!   edf ≈ q + trace( K (K + λ I)^(-1) ), with q=3 and λ = smoothing.

  implicit none

  ! ---- Inputs/outputs ----
  integer, intent(in) :: n_coarse, nperp, npar
  double precision, intent(in) :: grid_coarse(n_coarse)
  double precision, intent(in) :: pperp_coarse(n_coarse), ppar_coarse(n_coarse)
  double precision, intent(in) :: pperp(0:nperp,0:npar), ppar(0:nperp,0:npar)
  double precision, intent(in) :: smoothing
  double precision, intent(out) :: grid_fine(0:nperp,0:npar)

  ! ---- Locals (existing) ----
  integer :: i, j, k
  integer :: permutation_index(n_coarse+3)
  double precision :: fullmatrix(n_coarse+3,n_coarse+3)  ! A = [[K+λI, P],[P^T, 0]]
  double precision :: grid_vector(n_coarse+3)
  double precision :: weight_param(n_coarse+3)
  double precision :: r
  integer :: INFO   ! <-- LAPACK INFO should be INTEGER

  ! ---- Build RHS y with 3 appended zeros ----
  grid_vector = 0.d0
  do i=1,n_coarse
    grid_vector(i) = grid_coarse(i)
  end do
  ! grid_vector(n_coarse+1:n_coarse+3) remain 0

  ! ---- Assemble A = [[K+λI, P],[P^T, 0]] ----
  fullmatrix = 0.d0
  do i=1,n_coarse
    do j=1,n_coarse
      r = sqrt( (pperp_coarse(i)-pperp_coarse(j))**2 + (ppar_coarse(i)-ppar_coarse(j))**2 )
      if (r .ge. 1.d0) then
        fullmatrix(i,j) = r*r*log(r)
      else
        fullmatrix(i,j) = r*log(r**r)
      end if
    end do

    ! add smoothing (λ) to K-diagonal
    fullmatrix(i,i) = fullmatrix(i,i) + smoothing

    ! P columns
    fullmatrix(i,n_coarse+1) = 1.d0
    fullmatrix(i,n_coarse+2) = pperp_coarse(i)
    fullmatrix(i,n_coarse+3) = ppar_coarse(i)

    ! P^T rows
    fullmatrix(n_coarse+1,i) = 1.d0
    fullmatrix(n_coarse+2,i) = pperp_coarse(i)
    fullmatrix(n_coarse+3,i) = ppar_coarse(i)
  end do
  ! lower-right 3x3 block is already 0

  ! ---- Solve A * weight_param = [y; 0] ----
  weight_param = grid_vector
  call dgesv(n_coarse+3, 1, fullmatrix, n_coarse+3, permutation_index, &
             weight_param, n_coarse+3, INFO)

  ! ---- (OPTIONAL) Estimate effective DoF: edf ≈ q + trace(K (K+λI)^(-1)) ----
  ! This uses Hutchinson’s estimator with Rademacher probes. Cheap for N ~ few hundred.

  block
    integer :: N, mtrials
    integer, allocatable :: ipivK(:)
    double precision, allocatable :: K(:,:), Klam(:,:), z(:), x(:), t(:), urand(:)
    double precision :: acc, edf_eff, dotval, lambda_local
    external ddot
    double precision ddot

    N = n_coarse
    lambda_local = smoothing
    mtrials = 20        ! increase to ~50 for tighter estimate if desired

    if (N >= 3) then
      allocate(K(N,N), Klam(N,N), ipivK(N), z(N), x(N), t(N), urand(N))

      ! Build K with the SAME kernel as above
      do i=1,N
        do j=1,N
          r = sqrt( (pperp_coarse(i)-pperp_coarse(j))**2 + (ppar_coarse(i)-ppar_coarse(j))**2 )
          if (r .ge. 1.d0) then
            K(i,j) = r*r*log(r)
          else
            K(i,j) = r*log(r**r)
          end if
        end do
      end do

      ! Klam = K + λ I
      Klam = K
      do i=1,N
        Klam(i,i) = Klam(i,i) + lambda_local
      end do

      ! LU factorization of (K + λ I)
      call dgetrf(N, N, Klam, N, ipivK, INFO)
      if (INFO .ne. 0) then
        write(*,*) 'dgetrf failed for (K+lambda I), info=', INFO
      else
        acc = 0.d0

        do j=1,mtrials
          ! Rademacher z ∈ {±1}^N
          call random_number(urand)
          do i=1,N
            if (urand(i) .lt. 0.5d0) then
              z(i) = -1.d0
            else
              z(i) =  1.d0
            end if
          end do

          ! Solve (K + λ I) x = z
          x = z
          call dgetrs('N', N, 1, Klam, N, ipivK, x, N, INFO)
          if (INFO .ne. 0) then
            write(*,*) 'dgetrs failed, info=', INFO
            exit
          end if

          ! t = K * x
          t      = matmul(K, x)
          dotval = dot_product(z, t)
          !call dgemv('N', N, N, 1.0d0, K, N, x, 1, 0.0d0, t, 1)

          ! accumulate z^T t
          !dotval = ddot(N, z, 1, t, 1)
          acc = acc + dotval
        end do

        if (mtrials > 0) then
          acc = acc / dble(mtrials)
          edf_eff = 3.d0 + acc    ! q = 3 (affine tail: 1, pperp, ppar)
          write(*,'(a,1pe12.5)') 'Estimated effective DoF (edf) = ', edf_eff

          if (edf_eff .lt. 3.d0 .or. edf_eff .gt. dble(N) + 1.d-8) then
            write(*,'(a,i0,a,1pe12.5)') 'Warning: edf outside [3,', N, ']: ', edf_eff
          end if
        end if
      end if

      deallocate(K, Klam, ipivK, z, x, t, urand)
    end if
  end block

  ! ---- Evaluate on fine grid as before ----
  grid_fine = 0.d0
  do i=0,nperp
    do j=0,npar
      do k=1,n_coarse
        r = sqrt( (pperp(i,j)-pperp_coarse(k))**2 + (ppar(i,j)-ppar_coarse(k))**2 )
        if (r .ge. 1.d0) then
          grid_fine(i,j) = grid_fine(i,j) + weight_param(k) * r*r*log(r)
        else
          grid_fine(i,j) = grid_fine(i,j) + weight_param(k) * r*log(r**r)
        end if
      end do

      grid_fine(i,j) = grid_fine(i,j) + weight_param(n_coarse+1) &
                     + weight_param(n_coarse+2)*pperp(i,j) &
                     + weight_param(n_coarse+3)*ppar(i,j)
    end do
  end do

end subroutine

