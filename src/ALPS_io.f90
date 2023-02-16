!=============================================================================!
!=============================================================================!
!!*ALPS                                                                  *!!
!!                                                                           !!
!!Kristopher Klein, Daniel Verscharen                                        !!
!!
!!Space Science Center, University of New Hampshire                          !!
!!                                                                           !!
!!*INPUT/OUTPUT FUNCTIONS                                                   *!!
!=============================================================================!
!=============================================================================!
!

module alps_io
  implicit none

  !Numbers for identifying input/output files
  integer :: unit
  integer, parameter :: stdout_unit=6
  integer, save :: input_unit_no, error_unit_no=stdout_unit

  private :: get_runname, get_indexed_namelist_unit
  private :: input_unit_exist, input_unit
  private :: map_read, solution_read, spec_read, scan_read, bM_read
  public  :: init_param, read_f0, get_unused_unit, alps_error
  public :: output_time, display_credits, isnancheck

contains

!-=-=-=-=-
!Read in System parameters from *.in file
!Only processor 0 calls this routine
!-=-=-=-=-
  subroutine init_param
    use alps_var, only : runname, kperp, kpar, option, nroots, D_prec, D_gap
    use alps_var, only : kperp_last, kpar_last, kperp_0, kpar_0
    use alps_var, only : use_map, writeOut, wroots, nspec, numroots
    use alps_var, only : nperp, npar, arrayName, fit_check, param_fit, fit_type, perp_correction
    use alps_var, only : ns, qs, ms, vA, Bessel_zero, numiter, D_threshold,positions_principal
    use alps_var, only : determine_minima, n_resonance_interval, ngamma, npparbar, Tlim
    use alps_var, only : scan_option, n_scan, scan, use_secant, relativistic, logfit, usebM
    use alps_var, only : maxsteps_fit, n_fits, lambda_initial_fit, lambdafac_fit, epsilon_fit
    use alps_var, only : bMnmaxs, bMBessel_zeros, bMbetas, bMalphas, bMpdrifts
    implicit none

    integer :: ik, is, ifit, ip !Solution, species, fit, additional index

    !Namelist read in from input file.
    nameList /system/ &
         kperp, kpar, nspec, option, nroots, use_map, writeOut,&
         nperp, npar, ngamma, npparbar, vA, arrayName, Bessel_zero, numiter, D_threshold, &
         D_prec, D_gap, positions_principal, Tlim, &
         maxsteps_fit, lambda_initial_fit, lambdafac_fit, epsilon_fit, fit_check, &
         determine_minima, n_resonance_interval, scan_option, n_scan,use_secant

    !Get a unassigned unit number for input/output
    call get_unused_unit (input_unit_no)

    !Read in system parameters
    !runname called earlier in alps_error_init
    unit=input_unit_no
    open (unit=unit,file=trim(runname)//".in",status='old',action='read')
    read (unit=unit,nml=system)

    if (writeOut) &
         write(*,'(2a)') &
         'Reading from Input File: ', trim(runname)

    !save initial kperp,kpar values
    kperp_last=kperp; kpar_last=kpar
    kperp_0=kperp; kpar_0=kpar

    !Allocate solution space for nroots dispersion solutions
    allocate(wroots(1:numroots));wroots=cmplx(0.d0,0.d0,kind(1.d0))

    !Read in map scan parameters
    if (use_map) then
       if (writeOut) &
            write(*,'(a)') 'MAP ROUTINE'
       ik = 1
       !Either read in parameters for complex frequency map
       !Some FORTRAN compliers require namelists to be the first
       !declared item, thus requiring a separate subroutine
       unit=input_unit_no
       call get_indexed_namelist_unit (unit, "maps", ik)
       call map_read
       close(unit)
    else
    !OR
    !read in guesses for solutions
       if (writeOut) &
            write(*,'(a)') 'GUESS ROUTINE: '
       do ik = 1, nroots
          unit=input_unit_no
          call get_indexed_namelist_unit (unit, "guess", ik)
          !Some FORTRAN compliers require namelists to be the first
          !declared item, thus requiring a separate subroutine)
          call solution_read(ik)
          write(*,'(a,i3,a,2es14.4)')'Intial Guess ',ik,' : ', wroots(ik)
          close(unit)
       enddo
    endif

    !read in species density, charge, and mass from input file
    if (writeOut) &
         write(*,'(a)') 'SPECIES PARAMETERS: '
    allocate(ns(1:nspec)); ns = 0.d0
    allocate(qs(1:nspec)); qs = 0.d0
    allocate(ms(1:nspec)); ms = 0.d0
    allocate(n_fits(1:nspec)); n_fits = 1
    allocate(relativistic(1:nspec)); relativistic=.FALSE.
    allocate(logfit(1:nspec)); logfit=.TRUE.

    allocate(usebM(1:nspec)); usebM=.TRUE.
    allocate(bMnmaxs(1:nspec)); bMnmaxs=500
    allocate(bMBessel_zeros(1:nspec)); bMBessel_zeros=1.d-50
    allocate(bMbetas(1:nspec)); bMbetas=1.d0
    allocate(bMalphas(1:nspec)); bMalphas=1.d0
    allocate(bMpdrifts(1:nspec)); bMpdrifts=0.d0


    do is = 1, nspec
       !READ IN SPECIES PARAMETERS
       unit=input_unit_no
       call get_indexed_namelist_unit (unit, "spec", is)
          !Some FORTRAN compliers require namelists to be the first
          !declared item, thus requiring a separate subroutine)
       call spec_read(is)
       write(*,'(a,i3,a)')'Species ',is,' : '
       write(*,'(a,es11.4,a,es11.4,a,es11.4)')&
            ' ns/nREF = ',ns(is),' | qs/qREF = ',qs(is),' | ms/mREF = ',ms(is)
       write(*,'(a,i4)')&
            ' Number of fitted functions = ',n_fits(is)
       write(*,'(a,l1)')&
            ' Relativistic effects = ',relativistic(is)

       close(unit)
    enddo

    !I assume that n_fits for any species won't be more than 10 smaller
    !than the n_fits(ref). KGK -------> If we do it here, we don't have to assume that. DV 2016-05-18
    !Quite right. A much cleaner solution. KGK- 2016-05-18
    allocate(param_fit(1:nspec,0:max(nperp,ngamma),5,maxval(n_fits)))
    allocate(fit_type(1:nspec,maxval(n_fits)))
    allocate(perp_correction(1:nspec,maxval(n_fits)))

    write(*,'(a)')'-=-=-=-=-=-=-=-=-'

    do is = 1, nspec
       !READ IN SPECIES FIT PARAMETERS

       if (usebM(is)) then
             write (*,'(a,i2,a)') 'Species ',is,' uses bi-Maxwellian calculation...skipping fits. Parameters:'

             call get_indexed_namelist_unit (unit, "bM_spec", is)
             call bM_read(is)
             write(*,'(a,i4,a,es11.4)')&
                  '  nmax = ',bMnmaxs(is),',        Bessel_zero = ',bMBessel_zeros(is)
             write(*,'(a,es11.4,a,es11.4,a,es11.4)')&
                   '  beta = ',bMbetas(is),', alpha = ',bMalphas(is),', drift momentum = ',bMpdrifts(is)
      else

       do ifit=1,n_fits(is)
          !call get_indexed_namelist_unit (unit, "ffit", ifit)
          call get_indexed_double_namelist_unit (unit, "ffit", is, ifit)
          call fit_read(is,ifit)
          select case(fit_type(is,ifit))
          case(1)
             write(*,'(a,i2,a,i2,a)')&
                  'Species ',is,', function ',ifit,', Maxwellian fit: '
          case(2)
             write(*,'(a,i2,a,i2,a)')&
                  'Species ',is,', function ',ifit,', kappa fit: '
          case(3)
             write(*,'(a,i2,a,i2,a)')&
                  'Species ',is,', function ',ifit,', Juettner fit (pperp and ppar): '
		      case(4)
             write(*,'(a,i2,a,i2,a)')&
                  'Species ',is,', function ',ifit,', Juettner fit (gamma-dependent only): '
          case(5)
             write(*,'(a,i2,a,i2,a)')&
                  'Species ',is,', function ',ifit,', Juettner fit (gamma and pparbar): '
          case(6)
             write(*,'(a,i2,a,i2,a)')&
                  'Species ',is,', function ',ifit,', bi-Moyal fit: '
          case default
             write(*,'(a)')&
                  'Function fit undefined'
             stop
          end select

          do ip = 1, 5
             write(*,'(a,i2,a,es14.4)')&
                  ' Initial fit parameter ',ip,' = ',param_fit(is,0,ip,ifit)
          enddo
           write(*,'(a,es14.4)')&
                  ' Perpendicular correction:  ',perp_correction(is,ifit)
          close(unit)
       enddo
     endif
    enddo

    !Read in selection for scan paramter
    if (n_scan.gt.0) then
       write(*,'(a)')'-=-=-=-=-=-=-=-=-'
       allocate(scan(n_scan))
       do ip = 1, n_scan
          unit=input_unit_no
          call get_indexed_namelist_unit (unit, "scan_input", ip)
          call scan_read(ip)
          close(unit)
       enddo
       write(*,'(a)')'-=-=-=-=-=-=-=-=-'
    endif

    close(unit)

  end subroutine init_param
!-=-=-=-=-
!-=-=-=-=-

!-=-=-=-=-
!Subroutine for reading in map parameters
!-=-=-=-=-
subroutine map_read!KGK
  use alps_var, only : loggridw, loggridg, omi, omf, gami, gamf, ni, nr
  implicit none

  nameList /maps/ &
       loggridw, loggridg, omi, omf, gami, gamf, ni, nr
  read (unit=unit,nml=maps)

end subroutine map_read
!-=-=-=-=-
!-=-=-=-=-

!-=-=-=-=-
!Subroutine for reading in initial guesses for dispersion solutions
!-=-=-=-=-
subroutine solution_read(ik)
  use alps_var, only : wroots
  implicit none
  !Passed
  integer :: ik !solution index
  !Local
  double precision    :: g_om,g_gam !Guesses for (real, imaginary) solution

  nameList /guess/ &
       g_om,g_gam
  read (unit=unit,nml=guess)
  wroots(ik)=cmplx(g_om,g_gam,kind(1.d0))
end subroutine solution_read
!-=-=-=-=-
!-=-=-=-=-

!-=-=-=-=-
!Subroutine for reading in species parameters
!-=-=-=-=-
subroutine spec_read(is)
  use alps_var, only : ns, qs, ms, n_fits, relativistic, logfit, usebM
  implicit none
  !Passed
  integer :: is !species index
  !Local
  double precision :: nn,qq,mm    ! Read in values for ns, qs, ms
  integer :: ff !read in value for number of fitted functions for species
  logical :: relat=.false.
  logical :: log_fit=.true.
  logical :: use_bM=.false.

  nameList /spec/ &
       nn,qq,mm,ff,relat,log_fit,use_bM
  read (unit=unit,nml=spec)
  ns(is) = nn; qs(is) = qq; ms(is) = mm; n_fits(is)=ff; relativistic(is)=relat; &
   logfit(is)=log_fit; usebM(is)=use_bM
end subroutine spec_read
!-=-=-=-=-

!-=-=-=-=-
!Subroutine for reading in bi-Maxwellian parameters
!-=-=-=-=-
subroutine bM_read(is)
  use alps_var, only : bMnmaxs,bMBessel_zeros,bMbetas,bMalphas,bMpdrifts
  implicit none
  !Passed
  integer :: is !species index
  !Local

  integer :: bM_nmaxs
  double precision :: bM_Bessel_zeros, bM_betas, bM_alphas, bM_pdrifts

  nameList /bM_spec/ &
       bM_nmaxs,bM_Bessel_zeros,bM_betas,bM_alphas,bM_pdrifts

  read (unit=unit,nml=bM_spec)
    bMnmaxs(is)=bM_nmaxs; bMBessel_zeros(is)=bM_Bessel_zeros; bMbetas(is)=bM_betas; &
    bMalphas(is)=bM_alphas; bMpdrifts(is)=bM_pdrifts
end subroutine bM_read

!-=-=-=-=-

!-=-=-=-=-
!Subroutine for reading in wavevector scan parameters
!-=-=-=-=-
subroutine scan_read(is)
  use alps_var, only : scan,kperp_last,kpar_last
  implicit none
  !Passed
  integer :: is !scan index
  !Local
  integer :: scan_type,ns,nres
  double precision    :: swi,swf
  logical :: swlog,heating,eigen
  double precision :: theta_0
  double precision :: k_0

  nameList /scan_input/ &
       scan_type,swi,swf,swlog,ns,nres,&
       heating,eigen
  read (unit=unit,nml=scan_input)

  scan(is)%range_i =swi
  scan(is)%range_f =swf
  scan(is)%log_scan=swlog
  scan(is)%type_s  =scan_type
  scan(is)%n_out   =ns
  scan(is)%n_res   =nres
  scan(is)%eigen_s =eigen
  scan(is)%heat_s  =heating

  !calculate step size
  select case (scan_type)
  case(0)
     !scan from k_0 to k_1
     write(*,'(a,i0,a,es12.3,a,es12.3,a,es12.3,a,es12.3,a)')&
          'Scan ',is,': (kpar,kperp) from (',&
          kperp_last,',',kpar_last,') to (',swi,',',swf,')'
     if (swlog) then
        scan(is)%diff =(log10(swi)-log10(kperp_last))/(1.d0*ns*nres)
        scan(is)%diff2=(log10(swf)-log10(kpar_last))/(1.d0*ns*nres)
     else
        scan(is)%diff =(swi-kperp_last)/(1.d0*ns*nres)
        scan(is)%diff2=(swf-kpar_last )/(1.d0*ns*nres)
     endif
     kperp_last=swi;kpar_last=swf
  case(1)
     !scan from theta_0 to theta_1
     theta_0=atan(kperp_last/kpar_last); k_0=sqrt(kperp_last**2+kpar_last**2)
     write(*,'(a,i0,a,es12.3,a,es12.3)')&
          'Scan ',is,': theta from ',&
          theta_0*180.d0/(4.d0*atan(1.d0)),' to ',swf
     if (swlog) then
        scan(is)%diff=(log10(swf*4.d0*atan(1.d0)/180.)-log10(theta_0))/(1.d0*ns*nres)
     else
        scan(is)%diff=((swf*4.d0*atan(1.d0)/180.)-theta_0)/(1.d0*ns*nres)
     endif
     kpar_last=k_0*cos(swf*4.d0*atan(1.d0)/180.d0)
     kperp_last=k_0*sin(swf*4.d0*atan(1.d0)/180.d0)
  case(2)
     !scan from |k_0| to |k_1| at constant theta
     theta_0=atan(kperp_last/kpar_last); k_0=sqrt(kperp_last**2+kpar_last**2)
     write(*,'(a,i0,a,es12.3,a,es12.3,a,es12.3,a,es12.3,a)')&
          'Scan ',is,': |k| from ',&
          k_0,' to ',swf,' at theta=',theta_0*180.d0/(4.d0*atan(1.d0))
     if (swlog) then
        scan(is)%diff =(log10(swf*sin(theta_0))-log10(kperp_last))/(1.d0*ns*nres)
        scan(is)%diff2=(log10(swf*cos(theta_0))-log10(kpar_last))/(1.d0*ns*nres)
     else
        scan(is)%diff =(swf*sin(theta_0)-kperp_last)/(1.d0*ns*nres)
        scan(is)%diff2=(swf*cos(theta_0)-kpar_last )/(1.d0*ns*nres)
     endif
     kpar_last=swf*cos(theta_0)
     kperp_last=swf*sin(theta_0)
  case(3)
     !scan of kperp; kpar constant
     write(*,'(a,i0,a,es12.3,a,es12.3)')&
          'Scan ',is,': kperp from ',kperp_last,' to ',swf
     if (swlog) then
        scan(is)%diff=(log10(swf)-log10(kperp_last))/(1.d0*ns*nres)
     else
        scan(is)%diff=(swf-kperp_last)/(1.d0*ns*nres)
     endif
     kperp_last=swf
  case(4)
     !scan of kpar; kperp constant
     write(*,'(a,i0,a,es12.3,a,es12.3)')&
          'Scan ',is,': kpar from ',kpar_last,' to ',swf
     if (swlog) then
        scan(is)%diff=(log10(swf)-log10(kpar_last))/(1.d0*ns*nres)
     else
        scan(is)%diff=(swf-kpar_last)/(1.d0*ns*nres)
     endif
     kpar_last=swf
  end select

end subroutine scan_read
!-=-=-=-=-
!-=-=-=-=-

!-=-=-=-=-
!Subroutine for reading in species parameters
!-=-=-=-=-
subroutine fit_read(is,ifit)
  use alps_var, only : fit_type, param_fit, perp_correction
  implicit none
  !Passed
  integer :: is, ifit !species,  index
  !Local
  double precision :: fit_1, fit_2, fit_3, fit_4, fit_5    ! Read in values for fit functions
  double precision :: perpcorr
  integer :: fit_type_in !read in value for type of analytic function

  nameList /ffit/ &
       fit_type_in, fit_1, fit_2, fit_3, fit_4, fit_5, perpcorr

  fit_1=0.d0;fit_2=0.d0;fit_3=0.d0;fit_4=0.d0;fit_5=0.d0;perpcorr=0.d0
  read (unit=unit,nml=ffit)
  fit_type(is,ifit)=fit_type_in
  param_fit(is,0,1,ifit)=fit_1
  param_fit(is,0,2,ifit)=fit_2
  param_fit(is,0,3,ifit)=fit_3
  param_fit(is,0,4,ifit)=fit_4
  param_fit(is,0,5,ifit)=fit_5
  perp_correction(is,ifit)=perpcorr

end subroutine fit_read
!-=-=-=-=-
!-=-=-=-=-

!-=-=-=-=-
! Get runname for output files from input argument
!-=-=-=-=-
  subroutine get_runname(runname)
    implicit none
    integer       :: l
    character(500) :: arg
    character(500), intent(out) :: runname

    !Get the first argument of the program execution command
    call getarg(1,arg)

    !Check if this is the input file and trim .in extension to get runname
    l = len_trim (arg)
    if (l > 3 .and. arg(l-2:l) == ".in") then
       runname = arg(1:l-3)
    end if
  end subroutine get_runname
!-=-=-=-=-
!-=-=-=-=-

!-=-=-=-=-
! Subroutine for reading in background distribution function
!-=-=-=-=-
subroutine read_f0
  use alps_var, only : nperp, npar, arrayName, f0, pp, nspec
  use alps_var, only : writeOut, usebM
  implicit none
  !Local
  integer :: iperp,ipar !parallel and perpendicular indices
  integer :: is         !species index
  character (100) :: readname


     if (writeOut) &
          write(*,'(2a)')&
          'Attempting to read f0 array from file: ',trim(arrayName)

     call get_unused_unit (input_unit_no)
     unit=input_unit_no
     !The f0 arrays are stored in the distribution folder
     !arrayName is read in from *.in input file
     !each species is has a unique file for f0 and pp
     do is = 1, nspec

       if (usebM(is)) then
          write (*,'(a,i2)') ' Bi-Maxwellian calcuation: not reading f0 array for species ',is
          f0(is,:,:)=0.d0
          pp(is,:,:,:)=0.d0

       else

        write(readname, '(3a,i0,a)') &
             "distribution/",trim(arrayName),'.',is,".array"
        open (unit=unit, file=trim(readname), status='old', action='read')

        do iperp=0,nperp
           do ipar=0,npar
              !read in: pperp, ppar, f0
              read(unit,*) pp(is,iperp,ipar,1),pp(is,iperp,ipar,2),f0(is,iperp,ipar)
           enddo
        enddo
        close(unit)
      endif
     enddo


end subroutine read_f0
!-=-=-=-=-
!-=-=-=-=-

!-=-=-=-=-=-
!The following routines:
!    get_indexed_namelist_unit
!    input_unit_exist
!    get_unused_unit
!    input_unit
!were all adopted from the Astrophysical Gyrokinetic Code (AGK)
!as a means of allowing arbitrary namelist group name input.
!A bit of hassle, but worth the effort.
!-=-=-=-=-=-
  subroutine get_indexed_namelist_unit (unit, nml, index_in)
    use alps_var, only : runname
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: nml
    integer, intent (in) :: index_in
    character(500) :: line
    integer :: iunit, iostat, in_file
    integer :: ind_slash
    logical :: exist

    call get_unused_unit (unit)
    ind_slash=index(runname,"/",.True.)
    if (ind_slash.EQ.0) then !No slash in name
        !Original behaviour
        open (unit=unit, file='.'//trim(runname)//'.scratch')
    else
        !General behaviour
        open (unit=unit, file=trim(runname(1:ind_slash))//"."//trim(runname(ind_slash+1:))//".scratch")
    endif

    write (line, *) index_in
    line = nml//"_"//trim(adjustl(line))
    in_file = input_unit_exist(trim(line), exist)

    if (exist) then
       iunit = input_unit(trim(line))
    else
       call alps_error(1)
       !return
    end if

    read (unit=iunit, fmt="(a)") line
    write (unit=unit, fmt="('&',a)") nml

    do
       read (unit=iunit, fmt="(a)", iostat=iostat) line
       if (iostat /= 0 .or. trim(adjustl(line)) == "/") exit
       write (unit=unit, fmt="(a)") trim(line)
    end do
    write (unit=unit, fmt="('/')")
    rewind (unit=unit)
  end subroutine get_indexed_namelist_unit

!KGK- a work around to allow fit parameter readins for
!an arbitrary number of species for an arbitrary number of fitted functions
  subroutine get_indexed_double_namelist_unit (unit, nml, spec_in, index_in)
    use alps_var, only : runname
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: nml
    integer, intent (in) :: index_in
    integer, intent (in) :: spec_in
    character(500) :: line,lines
    integer :: iunit, iostat, in_file
    integer :: ind_slash
    logical :: exist

    call get_unused_unit (unit)
    ind_slash=index(runname,"/",.True.)
    if (ind_slash.EQ.0) then !No slash in name
        !Original behaviour
        open (unit=unit, file='.'//trim(runname)//'.scratch')
    else
        !General behaviour
        open (unit=unit, file=trim(runname(1:ind_slash))//"."//trim(runname(ind_slash+1:))//".scratch")
    endif

    write (line, *) index_in
    write (lines, *) spec_in
    line = nml//"_"//trim(adjustl(lines))//"_"//trim(adjustl(line))
    in_file = input_unit_exist(trim(line), exist)

    if (exist) then
       iunit = input_unit(trim(line))
    else
       call alps_error(1)
    end if

    read (unit=iunit, fmt="(a)") line
    write (unit=unit, fmt="('&',a)") nml

    do
       read (unit=iunit, fmt="(a)", iostat=iostat) line
       if (iostat /= 0 .or. trim(adjustl(line)) == "/") exit
       write (unit=unit, fmt="(a)") trim(line)
    end do
    write (unit=unit, fmt="('/')")
    rewind (unit=unit)
  end subroutine get_indexed_double_namelist_unit

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
!Open a file for the error log
  subroutine alps_error_init
    use alps_var, only : unit_error, runname
    implicit none
    call get_unused_unit(unit_error)
    !Get the run name, which comes from the name
    !of the input file appended after the executable:
    !mpirun -np 8 ./alps.e sample.in
    !yields a run name of 'sample'
    call get_runname(runname)
    open (unit=unit_error,file=trim(runname)//".log",status='replace')
  end subroutine alps_error_init

!Error catching subroutine
  subroutine alps_error(error_id)
    use alps_var, only : ierror,unit_error,nproc,scan_option
    use mpi

    implicit none
    integer :: error_id

!    if (proc0) then
       select case(error_id)
       case(0) !seen by all processors
          write(*,'(a,i6)')'ERROR: Number of processors must be even and greater than 2: nproc= ',nproc
          write(unit_error,'(a,i6)')'ERROR: Number of processors must be even and greater than 2: nproc= ',nproc
       case(1) !seen by proc0
          write(*,'(2a)') "get_indexed_namelist: required input namelist not found "
          write(unit_error,'(2a)') "get_indexed_namelist: required input namelist not found "
       case(2) !seen by proc0
          write (*,'(a)') "ERROR: More fit parameters than data points."
          write(unit_error,'(a)') "ERROR: More fit parameters than data points."
       case(3) !seen by all processors
          write(*,'(a,i6)')&
               'ERROR: scan_option not set to allowable value:',scan_option
          write(unit_error,'(a,i6)')&
               'ERROR: scan_option not set to allowable value:',scan_option
       case(4) !seen by all processors
          write(*,'(a)')&
               'ERROR: n_scan .ne.2 for scan_option=2'
          write(unit_error,'(a)')&
               'ERROR: n_scan .ne.2 for scan_option=2'
       case(5) !seen by all processors
          write(*,'(a)')&
               'ERROR: scan(1)%type_s==scan(2)%type_s for double k scan'
          write(unit_error,'(a)')&
               'ERROR: scan(1)%type_s==scan(2)%type_s for double k scan'
       case(6) !seen by all processors
          write(*,'(a)')&
               'ERROR: scan(*)%type_s=0 not allowed for double k scan'
          write(unit_error,'(a)')&
               'ERROR: scan(*)%type_s=0 not allowed for double k scan'
      case(7) !seen by all processors
          write(*,'(a)')&
               'ERROR: Fit for analytical continuation failed. Adjustment of perpcorr may resolve this problem.'
          write(unit_error,'(a)')&
               'ERROR: Fit for analytical continuation failed. Adjustment of perpcorr may resolve this problem.'
      case(8) !seen by all processors
          write(*,'(a)')&
               'ERROR: Resonance integration covers entire subluminal cone.'
          write(*,'(a)')&
               '       Adjustment of positions_principal, npar, or a non-relativistic run may resolve this problem.'
          write(unit_error,'(a)')&
               'ERROR: Resonance integration covers entire subluminal cone.'
          write(unit_error,'(a)')&
               '       Adjustment of positions_principal, npar, or a non-relativistic run may resolve this problem.'
	  case(9) !seen by all processors
          write(*,'(a)')&
               'ERROR: All roots diverged. Adjustment of initial guesses or step width may resolve this problem.'
          write(unit_error,'(a)')&
               'ERROR: All roots diverged. Adjustment of initial guesses or step width may resolve this problem.'
       case default
          write(*,'(a)')'ERROR: Unspecified...'
          write(unit_error,'(a)')'ERROR: Unspecified...'
       end select
       write(*,'(a)') 'Finishing ALPS=================================='
       write(*,'(a)') "Time:"
       call output_time
       write(*,'(a)')'================================================'
 !   endif

    close(unit_error)
    call mpi_abort(MPI_COMM_WORLD,error_id, ierror)
  stop

  end subroutine alps_error


!checks if double precision number input is NaN
logical function isnancheck(input)
implicit none
double precision :: input

isnancheck=.FALSE.
if (abs(input).GE.huge(1.d0)) isnancheck=.TRUE.
!if ((input+1.d0).EQ.input) isnancheck=.TRUE.
if (input.NE.input) isnancheck=.TRUE.

end function isnancheck


!-=-=-=-=-=-=
!Output the date and time in a given format
!-=-=-=-=-=-=

subroutine output_time
	implicit none
	character(8)  :: date
	character(10) :: time
	character(5)  :: zone
	integer,dimension(8) :: value
	call date_and_time(date,time,zone,value)
	write (*,'(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') value(1),"-",value(2),"-",value(3)," -- ",value(5),":",value(6),":",value(7)

end subroutine


!-=-=-=-=-=-=
!write the opening credits
!-=-=-=-=-=-=
subroutine display_credits
implicit none

 write (*,*) "==========================================================="
 write (*,*) "I                                                         I"
 write (*,*) "I                       A  L  P  S                        I"
 write (*,*) "I              Arbitrary Linear Plasma Solver             I"
 write (*,*) "I                                                         I"
 write (*,*) "I                       Version 1.0                       I"
 write (*,*) "I                                                         I"
 write (*,*) "I  Kristopher Klein   (kris.klein@gmail.com)              I"
 write (*,*) "I  Daniel Verscharen  (d.verscharen@ucl.ac.uk)            I"
 write (*,*) "I                                                         I"
 write (*,*) "==========================================================="


end subroutine

end module alps_io
