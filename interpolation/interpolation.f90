program interpolate
implicit none

double precision, allocatable, dimension (:) :: grid_coarse,pperp_coarse,ppar_coarse
double precision, allocatable, dimension (:,:) :: grid_fine,pperp,ppar
double precision :: smoothing,pperp_min,pperp_max,ppar_min,ppar_max,threshold
double precision :: mult_pperp,mult_ppar,mult_f
double precision :: r_pperp,r_ppar,r_f,pperp_max_set,pperp_min_set,ppar_max_set,ppar_min_set

integer :: nperp,npar,n_coarse,mode
integer :: i,j,io_error,status_read,i_coarse
character(LEN=256) :: filename,output_file
logical :: out_to_file,do_normalize

!I/O values for namelist
  integer :: unit
  integer, parameter :: stdout_unit=6
  integer, save :: input_unit_no, error_unit_no=stdout_unit


!string for parameter input file
character(50) :: runname

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
		!Read in system parameters
		!input file is argument after executable:
		!$ ./interpolate input.in
		implicit none
		!Dummy values for reading in species parameters

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
		call get_runname(runname)
		runname=trim(runname)//".in"
		unit=input_unit_no
		open (unit=unit,file=runname,status='old',action='read')
		read (unit=unit,nml=system)

		unit=input_unit_no
		close (unit)

	end subroutine read_in_params


	!-=-=-=-=-=-
	!The following routines:
	!    input_unit_exist
	!    get_unused_unit
	!    input_unit
	!were all adopted from the Astrophysical Gyrokinetic Code (AGK)
	!as a means of allowing arbitrary namelist group name input.
	!A bit of hassle, but worth the effort.
	!-=-=-=-=-=-


	  function input_unit_exist (nml,exist)
	    implicit none
	    character(*), intent (in) :: nml
	    logical, intent(out) :: exist
	    integer :: input_unit_exist, iostat
	    character(500) :: line
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
	    implicit none
	    character(*), intent (in) :: nml
	    integer :: input_unit, iostat
	    character(500) :: line
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
	    implicit none
	    integer, intent (out) :: unit
	    logical :: od
	    unit = 50
	    do
	       inquire (unit=unit, opened=od)
	       if (.not.od) return
	       unit = unit + 1
	    end do
	  end subroutine get_unused_unit
	!-=-=-=-=-=-

	!---------------------------------------------------------------
	! Get runname for output files from input argument
	  subroutine get_runname(runname)
	    implicit none
	    integer       :: l
	    character(50) :: arg
	    character(50), intent(out) :: runname

	    !Get the first argument of the program execution command
	    call getarg(1,arg)

	    !Check if this is the input file and trim .in extension to get runname
	    l = len_trim (arg)
	    if (l > 3 .and. arg(l-2:l) == ".in") then
	       runname = arg(1:l-3)
	    end if
	  end subroutine get_runname
	!------------------------------------------------------------------------------


end program





subroutine normalize (pperp,ppar,grid_fine,nperp,npar)
implicit none

integer :: nperp,npar,i,j
double precision :: pperp(0:nperp,0:npar),ppar(0:nperp,0:npar),grid_fine(0:nperp,0:npar)
double precision :: dpperp,dppar,integral
double precision :: M_PI=2.d0*acos(0.d0)

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



! Polyharmonic Spline:

subroutine polyharmonic_spline(grid_coarse,pperp_coarse,ppar_coarse,n_coarse,pperp,ppar,nperp,npar,smoothing,grid_fine)
!
! This soubroutine interpolates the grid with a polyharmonic thin-plate spline.
!
! Input:
! grid_coarse is a vector of length n_coarse that includes the values of f at each point i
! pperp_coarse is a vector of length n_coarse that includes the values of pperp at each point i
! ppar_coarse is a vector of length n_coarse that includes the values of ppar at each point i
! n_coarse is the total number of points in the coarse grid (nperp_coarse * npar_coarse)
!
! pperp is the value of pperp in the fine grid. It is a field of rank (nperp, npar)
! ppar is the value of ppar in the fine grid. It is a field of rank (nperp, npar)
! nperp is the number of perpendicular data points in the fine grid
! npar is the number of parallel data points in the fine grid
!
!
! Output:
! grid_fine is the interpolated grid. It is a field of rank (nperp, npar)
!
! This subroutine needs the LUPACK and BLAS libraries to evoke the dgesv subroutine
!

! This is the Thin Plate Spline:
! We use these resources:
! http://cseweb.ucsd.edu/~sjb/eccv_tps.pdf
! http://www.univie.ac.at/nuhag-php/bibtex/open_files/po94_M%20J%20D%20Powell%2003%2093.pdf
! http://vision.ucsd.edu/sites/default/files/fulltext(4).pdf
implicit none

integer :: i,j,k,permutation_index(n_coarse+3)
integer :: nperp,npar,n_coarse
double precision :: pperp_coarse(n_coarse),ppar_coarse(n_coarse)
double precision :: grid_coarse(n_coarse),grid_fine(0:nperp,0:npar)

double precision :: fullmatrix(n_coarse+3,n_coarse+3)
double precision :: grid_vector(n_coarse+3),weight_param(n_coarse+3)
double precision :: pperp(0:nperp,0:npar),ppar(0:nperp,0:npar)
double precision :: r,smoothing
double precision :: INFO

grid_vector=0.d0
do i=1,n_coarse
	grid_vector(i)=grid_coarse(i)
enddo
! grid_vector has three additional entries. The last three entries are all zero.

fullmatrix=0.d0
do i=1,n_coarse
	do j=1,n_coarse

		! Do the K-matrix part first:
		r=sqrt((pperp_coarse(i)-pperp_coarse(j))**2+(ppar_coarse(i)-ppar_coarse(j))**2)
		if(r.GE.1.d0) then
			fullmatrix(i,j)=r*r*log(r)
		else
			fullmatrix(i,j)=r*log(r**r)
		endif

	enddo

	fullmatrix(i,i)=fullmatrix(i,i)+smoothing

	! Now the P-matrix parts:
	fullmatrix(i,n_coarse+1)=1.d0
	fullmatrix(i,n_coarse+2)=pperp_coarse(i)
	fullmatrix(i,n_coarse+3)=ppar_coarse(i)

	! and the transposed P-matrix:
	fullmatrix(n_coarse+1,i)=1.d0
	fullmatrix(n_coarse+2,i)=pperp_coarse(i)
	fullmatrix(n_coarse+3,i)=ppar_coarse(i)
enddo

weight_param=grid_vector
call dgesv(n_coarse+3,1,fullmatrix,n_coarse+3,permutation_index,weight_param,n_coarse+3,INFO)


grid_fine=0.d0
do i=0,nperp
do j=0,npar

 do k=1,n_coarse
	r=sqrt((pperp(i,j)-pperp_coarse(k))**2+(ppar(i,j)-ppar_coarse(k))**2)
	if (r.GE.1.d0) then
		grid_fine(i,j)=grid_fine(i,j)+weight_param(k)*r*r*log(r)
	else
		grid_fine(i,j)=grid_fine(i,j)+weight_param(k)*r*log(r**r)
	endif
 enddo


 grid_fine(i,j)=grid_fine(i,j)+weight_param(n_coarse+1)+weight_param(n_coarse+2)*pperp(i,j)+&
 weight_param(n_coarse+3)*ppar(i,j)

enddo
enddo

end subroutine




!
