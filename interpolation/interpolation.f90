program interpolate
implicit none

double precision, allocatable, dimension (:) :: grid_coarse,pperp_coarse,ppar_coarse
double precision, allocatable, dimension (:,:) :: grid_fine,pperp,ppar
double precision :: smoothing,pperp_min,pperp_max,ppar_min,ppar_max,threshold
double precision :: mult_pperp,mult_ppar,mult_f
double precision :: r_pperp,r_ppar,r_f,pperp_max_set,pperp_min_set,ppar_max_set,ppar_min_set

integer :: nperp,npar,n_coarse,narg,run_arg,mode
integer :: i,j,io_error,status_read,i_coarse
character(LEN=256) :: filename,argument,output_file
logical :: skip_next,out_to_file,do_normalize

nperp=-1		
npar=-1			
smoothing=-1.d0 
mode=1			
pperp_max_set=-9999.d0
ppar_max_set=-9999.d0
pperp_min_set=-9999.d0
ppar_min_set=-9999.d0
filename=" " 
out_to_file=.FALSE.
do_normalize=.FALSE.
threshold=0.d0
mult_pperp=1.d0
mult_ppar=1.d0
mult_f=1.d0

narg=command_argument_count() 
if (narg.EQ.0) then
	write (*,*) "Specify argument. Use --help or -h for help."
	stop
endif

skip_next=.FALSE.

 do run_arg=1,narg
  if (skip_next.EQV..FALSE.) then
  call get_command_argument(run_arg,argument)
   
	select case(adjustl(argument))
	case("--help","-h")
		call write_help
    	stop
    case("-f","-F")
		call get_command_argument(run_arg+1,argument)
		filename=adjustl(argument)
		skip_next=.TRUE.
   case("-o","-O")
		call get_command_argument(run_arg+1,argument)
		output_file=adjustl(argument)
		skip_next=.TRUE.
		out_to_file=.TRUE.		
	case("-s","-S")
		call get_command_argument(run_arg+1,argument)
		read (argument,*) smoothing
		skip_next=.TRUE.
	case("-nperp")
		call get_command_argument(run_arg+1,argument)
		read (argument,*) nperp
		skip_next=.TRUE.
	case("-npar")
		call get_command_argument(run_arg+1,argument)
		read (argument,*) npar
		skip_next=.TRUE.	
	case("-ppar_min")
		call get_command_argument(run_arg+1,argument)
		read (argument,*) ppar_min_set
		skip_next=.TRUE.	
	case("-ppar_max")
		call get_command_argument(run_arg+1,argument)
		read (argument,*) ppar_max_set
		skip_next=.TRUE.	
	case("-pperp_min")
		call get_command_argument(run_arg+1,argument)
		read (argument,*) pperp_min_set
		skip_next=.TRUE.	
	case("-pperp_max")
		call get_command_argument(run_arg+1,argument)
		read (argument,*) pperp_max_set
		skip_next=.TRUE.	
	case("-m")
		call get_command_argument(run_arg+1,argument)
		read (argument,*) mode
		skip_next=.TRUE.	
	case("-thr")
		call get_command_argument(run_arg+1,argument)
		read (argument,*) threshold
		skip_next=.TRUE.	
	case("-mult_pperp")
		call get_command_argument(run_arg+1,argument)
		read (argument,*) mult_pperp
		skip_next=.TRUE.	
	case("-mult_ppar")
		call get_command_argument(run_arg+1,argument)
		read (argument,*) mult_ppar
		skip_next=.TRUE.	
	case("-mult_f")
		call get_command_argument(run_arg+1,argument)
		read (argument,*) mult_f
		skip_next=.TRUE.	
	case("-n","-N")
		do_normalize=.TRUE.
	case default
		write(*,*) "Option ",adjustl(argument),"unknown."
	end select
 else
 	skip_next=.FALSE.
 endif	
end do
 
 
 
 
 if (len_trim(filename).EQ.0) then
 	write (*,*) "No filename selected."
 	stop
 endif

 if (smoothing.LT.0.d0) then
 	write (*,*) "# Using standard smoothing parameter: 0.1"
 	smoothing=0.1d0
 endif
 
  if (nperp.LT.0) then
 	write (*,*) "# Using standard nperp: 50"
 	nperp=50
 endif
 
if (npar.LT.0) then
 	write (*,*) "# Using standard nperp: 50"
 	npar=50
 endif

if (pperp_max_set.EQ.-9999.d0) then
	write (*,*) "# Using maximum pperp as upper bound."
endif

if (pperp_min_set.EQ.-9999.d0) then
	write (*,*) "# Using minimum pperp as lower bound."
endif

if (ppar_max_set.EQ.-9999.d0) then
	write (*,*) "# Using maximum ppar as upper bound."
endif

if (ppar_min_set.EQ.-9999.d0) then
	write (*,*) "# Using minimum ppar as lower bound."
endif

write (*,*) "# Mode:",mode
write (*,*) "# Threshold for f:",threshold



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
	write (*,*) "# No data points in selection."
	stop
endif

pperp_max=pperp_max*mult_pperp
pperp_min=pperp_min*mult_pperp
ppar_max=ppar_max*mult_ppar
ppar_min=ppar_min*mult_ppar


write (*,*) "# Number of points on coarse grid: ",n_coarse
write (*,*) "# Maximum Pperp:", pperp_max
write (*,*) "# Minimum Pperp:", pperp_min
write (*,*) "# Maximum Ppar:", ppar_max
write (*,*) "# Minimum Ppar:", ppar_min

if (pperp_max_set.NE.-9999.d0) pperp_max=pperp_max_set
if (ppar_max_set.NE.-9999.d0) ppar_max=ppar_max_set
if (pperp_min_set.NE.-9999.d0) pperp_min=pperp_min_set
if (ppar_min_set.NE.-9999.d0) ppar_min=ppar_min_set


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
		grid_coarse(i_coarse)=r_f
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

grid_fine=grid_fine*mult_f

! Output:
if (out_to_file) then
	open(unit=20,file=output_file,status='unknown',action='write',iostat=io_error) 
	do i=0,nperp
	do j=0,npar
		write (20,*) pperp(i,j),ppar(i,j),grid_fine(i,j)
	enddo
	enddo
	close(20)
else
	do i=0,nperp
	do j=0,npar
		write (*,*) pperp(i,j),ppar(i,j),grid_fine(i,j)
	enddo
	enddo
endif


end program





! Help output
subroutine write_help
implicit none

write(*,*) "This is the interpolation auxiliary program for ALPS."
write (*,*) " "
write (*,*) " " 
write (*,*) "Command line options:"
write (*,*) " "
write (*,*) "   -f, -F <file>         Specify file that contains the coarse grid."
write (*,*) "   -o, -O <file>         Output file for the fine grid."
write (*,*) "   -s, -S <arg>          Set smoothing parameter to <arg>. If not selected, smoothing = 0.1"
write (*,*) " " 
write (*,*) "   -nperp <arg>          Set the number of steps in nperp for the fine grid. Standard: 50"
write (*,*) "   -npar <arg>           Set the number of steps in npar for the fine grid. Standard: 50"
write (*,*) " "
write (*,*) "   -pperp_min <arg>      Set the minium value of pperp for the fine grid. If not selected, "
write (*,*) "                         the minimum value of the fine grid is used."
write (*,*) "   -pperp_max <arg>      Set the maximum value of pperp for the fine grid. If not selected, "
write (*,*) "                         the maximum value of the fine grid is used."
write (*,*) "   -ppar_min <arg>       Set the minium value of ppar for the fine grid. If not selected, "
write (*,*) "                         the minimum value of the fine grid is used."
write (*,*) "   -ppar_max <arg>       Set the maximum value of ppar for the fine grid. If not selected, "
write (*,*) "                         the maximum value of the fine grid is used."
write (*,*) " " 
write (*,*) "   -thr <arg>            Set a threshold on f0. Every data point with f0 below this"
write (*,*) "                         threshold is ignored (Standard: 0)."
write (*,*) " " 
write (*,*) "   -m <arg>              Set the mode for the input file (Standard: 1):"
write (*,*) "                            1:    column 1: pperp, column 2: ppar,  column 3: f0"
write (*,*) "                            2:    column 1: ppar,  column 2: pperp, column 3: f0"
write (*,*) " " 
write (*,*) "   -n, -N                Normalize the output to unity in cylindrical coordinates."
write (*,*) " "
write (*,*) "   -mult_pperp <arg>     Multiply every value of pperp in the grid with a constant factor <arg>."
write (*,*) "   -mult_ppar <arg>      Multiply every value of ppar in the grid with a constant factor <arg>."
write (*,*) "   -mult_f <arg>         Multiply every value of f0 with a constant factor <arg> (after normalization)."
write (*,*) " "
write (*,*) "   --help, -h            Displays this help screen."

end subroutine




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
	integral=integral+grid_fine(i,j)*pperp(i,j)
enddo
enddo
integral=integral*dpperp*dppar

grid_fine=grid_fine/(integral*2.d0*M_PI)

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
! This subroutine needs the subroutines LUDCMP and LUBKSB -- see below
!

! This is the Thin Plate Spline:
! We use these resources:
! http://cseweb.ucsd.edu/~sjb/eccv_tps.pdf
! http://www.univie.ac.at/nuhag-php/bibtex/open_files/po94_M%20J%20D%20Powell%2003%2093.pdf
! http://vision.ucsd.edu/sites/default/files/fulltext(4).pdf
implicit none

integer :: i,j,k,permutation_index(n_coarse+3),odd_even,code
integer :: nperp,npar,n_coarse
double precision :: pperp_coarse(n_coarse),ppar_coarse(n_coarse)
double precision :: grid_coarse(n_coarse),grid_fine(0:nperp,0:npar)

double precision :: fullmatrix(n_coarse+3,n_coarse+3)
double precision :: grid_vector(n_coarse+3),weight_param(n_coarse+3)
double precision :: pperp(0:nperp,0:npar),ppar(0:nperp,0:npar)
double precision :: r,smoothing


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
call LUDCMP(fullmatrix,n_coarse+3,permutation_index,odd_even,code)
call LUBKSB(fullmatrix,n_coarse+3,permutation_index,weight_param)



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






!-=-=-=-=-=-=
! The following two functions are from pre-written sources:
!-=-=-=-=-=-=
!*******************************************************
!*    LU decomposition routines 				       *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************


!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
subroutine LUDCMP(a,n,indx,d,code)
implicit none
 
integer, parameter :: nmax=10000
double precision, parameter :: tiny=1.5d-16
double precision :: amax, dum, sum, a(n,n), vv(nmax)
integer :: n, indx(n), code, d, i, j, k, imax
 

 d=1; code=0


 do i=1,n
   amax=0.d0
   do j=1,n
     if (dabs(a(i,j)).GT.amax) amax=DABS(a(i,j))
   enddo ! j loop
   if(amax.LT.tiny) then
     code = 1
     return
   endif
   vv(i) = 1.d0 / amax
 enddo ! i loop

 do j=1,n
   do i=1,j-1
     sum = a(i,j)
     do k=1,i-1
       sum = sum - a(i,k)*A(k,j) 
     enddo ! k loop
     a(i,j) = sum
   enddo ! i loop
   amax = 0.d0
   do i=j,n
     sum = a(i,j)
     do k=1,j-1
       sum = sum - a(i,k)*a(k,j) 
     enddo ! k loop
     a(i,j) = sum
     dum = vv(i)*dabs(sum)
     if(dum.GE.amax) then
       imax = i
       amax = dum
     endif
   enddo ! i loop  
   
   if(j.NE.imax) then
     do k=1,n
       dum = a(imax,k)
       a(imax,k) = a(j,k)
       a(j,k) = dum
     enddo ! k loop
     d = -d
     vv(imax) = vv(j)
   endif

   indx(j) = imax
   if(dabs(a(j,j)) < tiny) a(j,j) = tiny

   if(j.NE.n) then
     dum = 1.d0 / a(j,j)
     do i=j+1,n
       a(i,j) = a(i,j)*dum
     enddo ! i loop
   endif 
 enddo ! j loop

 return
 end subroutine LUDCMP


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
subroutine LUBKSB(a,n,indx,b)
implicit none

integer :: n, indx(n), ii, i, ll, j 
double precision :: sum, a(n,n),b(n)


ii = 0

 do i=1,n
   ll = indx(i)
   sum = b(ll)
   b(ll) = b(i)
   if(ii.NE.0) then
     do j=ii,i-1
       sum = sum - a(i,j)*b(j)
     enddo ! j loop
   else if(sum.NE.0.d0) then
     ii = i
   endif
   b(i) = sum
 enddo ! i loop

 do i=n,1,-1
   sum = b(i)
   if(i < n) then
     do j=i+1,n
       sum = sum - a(i,j)*b(j)
     enddo ! j loop
   endif
   b(i) = sum / a(i,i)
 enddo ! i loop

 return
 end subroutine LUBKSB

