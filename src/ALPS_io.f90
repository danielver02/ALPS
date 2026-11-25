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
!I  Daniel Verscharen  (daniel.verscharen@rub.de)                                I
!I                                                                             I
!===============================================================================

module alps_io
  !!Controls input and output functions to and from main program.
  implicit none

  integer :: unit
  !!Index for file I/O.

  integer, parameter :: stdout_unit=6
  !! Standard unit for I/O.

  integer, save :: input_unit_no
  !! Saved input unit for use with multiple read in calls.

  integer, save :: error_unit_no=stdout_unit
  !! Error output unit.

  private :: get_runname, get_indexed_namelist_unit
  private :: input_unit_exist, input_unit
  private :: map_read, solution_read, spec_read, scan_read, bM_read
  private :: poly_read, relativistic_params_read, fitting_params_read
  public  :: init_param, read_f0, get_unused_unit, alps_error
  public :: output_time, display_credits, isnancheck
  public :: alps_continuation_error
  public :: canonicalize_nyquist_grid,decimal_coordinate_radius,fit_rounded_axis

contains


  subroutine alps_continuation_error(species,row,status,u,nonrelativistic)
    !! Report a failed Chebyshev continuation evaluation and terminate all MPI ranks.
    !! Called by eval_fit when its caller omits the optional status output and the
    !! evaluator reports a non-success status, such as an unavailable row, invalid
    !! model, or overflow/underflow while reconstructing the distribution or derivative.
    !! The diagnostic identifies the species, fitted row, kernel status code, and
    !! complex parallel coordinate u at which the evaluation failed. By default,
    !! row is a relativistic Gamma-grid index (igamma); nonrelativistic=.true.
    !! instead identifies a perpendicular-momentum row (iperp).
    !! The calling rank writes to stdout. MPI_Abort stops the entire communicator
    !! so other ranks cannot remain blocked in collectives, and a failed evaluation
    !! cannot silently enter the dispersion calculation as a zero contribution.
    use alps_var, only : proc0,unit_error,ierror
    use mpi
    implicit none
    integer, intent(in) :: species,row,status
    double complex, intent(in) :: u
    logical, intent(in), optional :: nonrelativistic
    character(80) :: label
    character(8) :: row_label
    label='ERROR: relativistic Chebyshev continuation:';row_label=' igamma='
    if(present(nonrelativistic)) then
       if(nonrelativistic) then
          label='ERROR: nonrelativistic Chebyshev continuation:';row_label=' iperp='
       endif
    endif
    write(*,'(a,3(a,i0),a,2es24.15)') trim(label), &
         ' species=',species,trim(row_label),row,' status=',status,' u=',u
    if(proc0) write(unit_error,'(a,3(a,i0),a,2es24.15)') trim(label), &
         ' species=',species,trim(row_label),row,' status=',status,' u=',u
    call mpi_abort(MPI_COMM_WORLD,max(1,status),ierror)
    error stop 'relativistic Chebyshev continuation failed'
  end subroutine alps_continuation_error



  subroutine init_param
    !!Read in system parameters from *.in file.
    !!Only processor 0 calls this routine:
    use alps_var, only : runname, foldername, kperp, kpar, nroots, D_prec, D_gap
    use alps_var, only : profile_startup,unit_error
    use alps_var, only : kperp_last, kpar_last, kperp_0, kpar_0
    use alps_var, only : use_map, writeOut, wroots, nspec, numroots
    use alps_var, only : nperp, npar, arrayName, param_fit, fit_type, perp_correction
    use alps_var, only : ns, qs, ms, vA, Bessel_zero, numiter
    use alps_var, only : D_threshold, D_tol, positions_principal
    use alps_var, only : determine_minima, n_resonance_interval, ngamma, Tlim
    use alps_var, only : scan_option, n_scan, scan, relativistic, logfit, usebM
    use alps_var, only : n_fits
    use alps_var, only : bMnmaxs, bMBessel_zeros, bMbetas, bMalphas, bMpdrifts
    use alps_var, only : ACmethod, poly_kind, poly_order
    use alps_var, only : secant_method, poly_derivative_type
    implicit none

    integer :: ik
    !!Solution index for iterating [[(solution_read(subroutine)]].

    integer :: read_status
    !! I/O status for the checked system namelist read.
    character(1024) :: read_message
    !! Processor-provided diagnostic accompanying a failed system namelist read.

    integer :: is
    !!Species index for iterating [[spec_read(subroutine)]],
    !![[read_f0(subroutine)]], [[bM_read(subroutine)]], and
    !![[fit_read(subroutine)]]

    integer :: ifit
    !!Fit index for iterating [[fit_read(subroutine)]].

    integer :: ip
    !!Index for iterating [[scan_read(subroutine)]] and
    !!internal looping in [[init_param(subroutine)]].

    !Namelist read in from input file:
    nameList /system/ &
         kperp, kpar, nspec, nroots, use_map, writeOut, profile_startup,&
         nperp, npar, vA, arrayName, Bessel_zero, &
         secant_method, numiter, D_threshold, &
         D_prec, D_gap, D_tol, positions_principal, Tlim, &
         determine_minima, n_resonance_interval, scan_option, n_scan

    !Get a unassigned unit number for input/output:
    call get_unused_unit (input_unit_no)

    !Read in system parameters.
    !runname called earlier in alps_error_init:
    unit=input_unit_no
    open (unit=unit,file=trim(foldername)//trim(runname)//".in",status='old',action='read')
    read (unit=unit,nml=system,iostat=read_status,iomsg=read_message)
    if(read_status/=0) then
       write(*,'(a)') 'Cannot read &system: '//trim(read_message)
       write(unit_error,'(a)') 'Cannot read &system: '//trim(read_message)
       call alps_error(16)
    endif

    if (scan_option==3) then
       ! No root/map/scan instructions are consumed by a single-wavevector operation.
       nroots=0; use_map=.false.; n_scan=0
    else if (scan_option==4.or.scan_option==5) then
       nroots=0; use_map=.false.
       if(scan_option==5.and.n_scan/=2) then
          write(*,'(a)') 'ERROR: scan_option=5 requires n_scan=2.'
          call alps_error(20)
       endif
       if(n_scan<1) then
          write(*,'(a)') 'ERROR: scan_option=4 requires n_scan >= 1.'
          call alps_error(20)
       endif
    endif

    call fitting_params_read

    if (writeOut) &
         write(*,'(2a)') &
         'Reading from Input File: ', trim(runname)

    !save initial kperp,kpar values:
    kperp_last=kperp; kpar_last=kpar
    kperp_0=kperp; kpar_0=kpar

    !Allocate solution space for nroots dispersion solutions;
    allocate(wroots(1:numroots));wroots=cmplx(0.d0,0.d0,kind(1.d0))

    !Read in map scan parameters:
    if (scan_option==3.or.scan_option==4.or.scan_option==5) then
       ! Nyquist parameters are read by the separate driver after species input.
    else if (use_map) then
       if (writeOut) &
            write(*,'(a)') 'MAP ROUTINE'
       ik = 1
       !Read in parameters for complex frequency map:
       unit=input_unit_no
       call get_indexed_namelist_unit (unit, "maps", ik)
       call map_read
       close(unit)
    else
    !OR
    !read in guesses for solutions:
       if (writeOut) &
            write(*,'(a)') 'GUESS ROUTINE: '
       do ik = 1, nroots
          unit=input_unit_no
          call get_indexed_namelist_unit (unit, "guess", ik)
          call solution_read(ik)
          write(*,'(a,i3,a,2es14.4e3)')&
               'Intial Guess ',ik,' : ', wroots(ik)
          close(unit)
       enddo
    endif

    !Read in species density, charge, and mass from input file:
    if (writeOut) &
         write(*,'(a)') 'SPECIES PARAMETERS: '
    allocate(ns(1:nspec)); ns = 0.d0
    allocate(qs(1:nspec)); qs = 0.d0
    allocate(ms(1:nspec)); ms = 0.d0

    !Fitting Parameters:
    allocate(n_fits(1:nspec)); n_fits = 1
    allocate(relativistic(1:nspec)); relativistic=.FALSE.
    allocate(logfit(1:nspec)); logfit=.TRUE.

    !Bi-Maxwellian/Cold-plasma Parameters:
    allocate(usebM(1:nspec)); usebM=.TRUE.
    allocate(bMnmaxs(1:nspec)); bMnmaxs=500
    allocate(bMBessel_zeros(1:nspec)); bMBessel_zeros=1.d-50
    allocate(bMbetas(1:nspec)); bMbetas=1.d0
    allocate(bMalphas(1:nspec)); bMalphas=1.d0
    allocate(bMpdrifts(1:nspec)); bMpdrifts=0.d0

    !Basis Function Parameters:
    allocate(ACmethod(1:nspec)); ACmethod = 1
    allocate(poly_kind(1:nspec)); poly_kind = 0
    allocate(poly_order(1:nspec)); poly_order = 0
    ! Methods 0/1 require finite differences; poly_read sets the method-2 default.
    allocate(poly_derivative_type(1:nspec)); poly_derivative_type = 0

    !READ IN SPECIES PARAMETERS:
    do is = 1, nspec
       unit=input_unit_no
       call get_indexed_namelist_unit (unit, "spec", is)
       call spec_read(is)
       if (usebM(is).and.relativistic(is)) &
            call rel_input_error(is,'use_bM=T with relat=T is unsupported; analytic bimaxwellian species require relat=F.')
       write(*,'(a,i3,a)')'Species ',is,' : '
       write(*,'(a,es14.4e3,a,es14.4e3,a,es14.4e3)')&
            ' ns/nREF = ',ns(is),' | qs/qREF = ',qs(is),' | ms/mREF = ',ms(is)
       if (usebM(is)) then
          write(*,'(a)') ' Using bi-Maxwellian/cold-plasma representation'
       else
          select case (ACmethod(is))
          case (0)
             write(*,'(a)') ' Using function defined in distribution/distribution_analyt.f90'
          case (1)
             write(*,'(a,i4)')&
                  ' Number of fitted functions = ',n_fits(is)
          case (2)
             write(*,'(a)')&
                  ' Using a polynomial basis representation'
          end select
          write(*,'(a,l1)')&
               ' Relativistic effects = ',relativistic(is)
       endif
       close(unit)
    enddo

    ! Species flags must be known, and ngamma finalized before sizing fits.
    if(any(relativistic)) call relativistic_params_read

    allocate(param_fit(1:nspec,0:max(nperp,ngamma),5,maxval(n_fits)))
    allocate(fit_type(1:nspec,maxval(n_fits)))
    allocate(perp_correction(1:nspec,maxval(n_fits)))

    write(*,'(a)')'-=-=-=-=-=-=-=-=-'

    !READ IN SPECIES FIT PARAMETERS
    do is = 1, nspec
       if (usebM(is)) then
          write (*,'(a,i2,a)') 'Species ',is,&
               ' uses bi-Maxwellian/cold-plasma calculation... skipping fits. Parameters:'
          call get_indexed_namelist_unit (unit, "bM_spec", is)
          call bM_read(is)

          write(*,'(a,es14.4e3,a,es14.4e3,a,es14.4e3)')&
               '  beta = ',bMbetas(is),', alpha = ',bMalphas(is),', drift momentum = ',bMpdrifts(is)

          if (bMbetas(is).EQ.0.d0) then
              write (*,'(a)') '  Cold-plasma calculation.'
          else
            write(*,'(a,i4,a,es14.4e3)')&
                 '  nmax = ',bMnmaxs(is),',        Bessel_zero = ',bMBessel_zeros(is)
          endif

          close(unit)
       else
          !Read in initial guesses for LM fits.
          select case(ACmethod(is))
          case (1)
             do ifit=1,n_fits(is)
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
          case (2)
             call get_indexed_namelist_unit (unit, "poly_spec", is)
             call poly_read(is)
             if (relativistic(is)) then
                if (.not.logfit(is)) then
                   call rel_input_error(is, &
                        'relativistic method 2 requires log_fit=T; direct fitting is unsupported.')
                endif

                if (poly_order(is)<=0) then
                   call rel_input_error(is, &
                        'relativistic method 2 requires a positive requested polynomial degree.')
                endif
             endif
             select case (poly_kind(is))
             case (1)
                write(*,'(a,i0,a,i0)')&
                     'Chebyshev representation of order ',poly_order(is),' for species ',is
             case default
                call alps_error(10)
             end select
             close(unit)
          end select
       endif
    enddo


    !Read in selection for scan paramter:
    if (n_scan.gt.0) then
       write(*,'(a)')'-=-=-=-=-=-=-=-=-'
       allocate(scan(n_scan))
       do ip = 1, n_scan
          unit=input_unit_no
          call get_indexed_namelist_unit (unit, "scan_input", ip)
          call scan_read(ip)
          close(unit)
       enddo
       if(scan_option==5) then
          if(scan(1)%type_s==scan(2)%type_s) then
             write(*,'(a)') 'ERROR: Nyquist plane requires two different scan types.'
             call alps_error(20)
          endif
       endif
       write(*,'(a)')'-=-=-=-=-=-=-=-=-'
    endif
    close(unit)
  end subroutine init_param



  subroutine fitting_params_read
    !! Read optional global LM settings on rank 0.
    !! No group means module defaults. Missing and malformed groups are distinct; duplicates are errors.
    !! A quote/comment-aware header scan permits arbitrary record length, placement, case and inline data.
    use alps_var, only : runname,foldername,unit_error
    use alps_var, only : maxsteps_fit,lambda_initial_fit,lambdafac_fit,epsilon_fit
    use iso_fortran_env, only : iostat_end,iostat_eor
    implicit none
    integer :: scan_unit
    !! Private read-only input unit used for header scanning and the subsequent namelist read.
    integer :: ios
    !! I/O status for opening, scanning, rewinding, or reading; includes end-of-record and end-of-file.
    integer :: close_status
    !! I/O status returned when closing the private input unit.
    integer :: groups
    !! Number of matching fitting_params headers outside comments and quoted strings.
    !! Zero retains defaults; more than one is rejected as a duplicate active group.
    integer :: token_length
    !! Full length of the header name being scanned, including characters beyond the token buffer.
    integer :: code
    !! ASCII character code used to recognize header-name characters.
    character :: ch
    !! Current input character read without advancing to the next record.
    character :: quote
    !! Active single- or double-quote delimiter; a blank indicates scanning outside a quoted string.
    character(32) :: token
    !! Lowercase header-name buffer following an ampersand; token_length prevents truncated-name matches.
    character(1024) :: message
    !! Processor-provided explanation of an input opening, scanning, rewind, or namelist-read failure.
    logical :: comment
    !! True after an unquoted exclamation mark, until the end of the current record.
    logical :: group_open
    !! Whether a previous namelist header has not yet reached its terminator.
    logical :: in_name
    !! True while collecting a namelist header name after an unquoted, uncommented ampersand.
    namelist /fitting_params/ maxsteps_fit,lambda_initial_fit,lambdafac_fit,epsilon_fit

    open(newunit=scan_unit,file=trim(foldername)//trim(runname)//'.in', &
         status='old',action='read',iostat=ios,iomsg=message)
    if(ios/=0) then
       write(*,'(a)') trim(message)
       write(unit_error,'(a)') trim(message)
       call alps_error(19)
    endif
    groups=0; quote=' '; comment=.false.; in_name=.false.; token=''; token_length=0
    group_open=.false.
    do
       read(scan_unit,'(a)',advance='no',iostat=ios,iomsg=message) ch
       if(ios==iostat_eor) then
          call finish_name
          comment=.false.
          cycle
       endif
       if(ios/=0) exit
       if(comment) cycle
       if(quote/=' ') then
          if(ch==quote) quote=' '
          cycle
       endif
       if(in_name) then
          code=iachar(ch)
          if((code>=iachar('a').and.code<=iachar('z')).or. &
               (code>=iachar('A').and.code<=iachar('Z')).or. &
               (code>=iachar('0').and.code<=iachar('9')).or.ch=='_') then
             token_length=token_length+1
             if(code>=iachar('A').and.code<=iachar('Z')) ch=achar(code+32)
             if(token_length<=len(token)) token(token_length:token_length)=ch
             cycle
          endif
          call finish_name
       endif
       select case(ch)
       case('!')
          comment=.true.
       case("'",'"')
          quote=ch
       case('/')
          group_open=.false.
       case('&')
          in_name=.true.; token=''; token_length=0
       end select
    enddo
    call finish_name
    if(ios/=iostat_end) then
       write(*,'(a)') 'Cannot scan fitting parameters: '//trim(message)
       write(unit_error,'(a)') 'Cannot scan fitting parameters: '//trim(message)
       close(scan_unit,iostat=close_status)
       call alps_error(19)
    endif
    if(groups>1) then
       write(*,'(a)') 'Duplicate &fitting_params blocks are not allowed.'
       write(unit_error,'(a)') 'Duplicate &fitting_params blocks are not allowed.'
       close(scan_unit,iostat=close_status)
       call alps_error(19)
    endif
    if(groups==1) then
       rewind(scan_unit,iostat=ios,iomsg=message)
       if(ios==0) read(scan_unit,nml=fitting_params,iostat=ios,iomsg=message)
       if(ios/=0) then
          write(*,'(a)') 'Cannot read &fitting_params: '//trim(message)
          write(unit_error,'(a)') 'Cannot read &fitting_params: '//trim(message)
          close(scan_unit,iostat=close_status)
          call alps_error(19)
       endif
    endif
    close(scan_unit,iostat=close_status)
    if(close_status/=0) call alps_error(19)

  contains
    subroutine finish_name
      !! Count a complete matching header name, then reset the host-associated token-scanning state.
      if(in_name.and.token_length>=len('fitting_params_')) then
         if(index(token,'fitting_params_')==1) then
            write(*,'(a)') 'ERROR: use one global &fitting_params block without a species index.'
            write(unit_error,'(a)') 'ERROR: use one global &fitting_params block without a species index.'
            call alps_error(19)
         endif
      endif
      if(in_name.and.token_length==len('fitting_params')) then
         if(token=='fitting_params') then
            if(group_open) then
               write(*,'(a)') 'ERROR: &fitting_params must appear between complete namelist blocks.'
               write(unit_error,'(a)') 'ERROR: &fitting_params must appear between complete namelist blocks.'
               call alps_error(19)
            endif
            groups=groups+1
         endif
      endif
      if(in_name) group_open=(token/='end')
      in_name=.false.; token=''; token_length=0
    end subroutine
  end subroutine fitting_params_read

  subroutine relativistic_params_read
    !! Read the optional shared relativistic group, only after any relativistic species is identified.
    !! No group means module defaults. Missing and malformed groups are distinct; duplicates are errors.
    !! A quote/comment-aware header scan permits arbitrary record length, placement, case and inline data.
    use alps_var, only : runname,foldername,unit_error,ngamma,npparbar,rel_cache,rel_cache_dir
    use iso_fortran_env, only : iostat_end,iostat_eor
    implicit none
    integer :: scan_unit
    !! Private read-only input unit used for header scanning and the subsequent namelist read.
    integer :: ios
    !! I/O status for opening, scanning, rewinding, or reading; includes end-of-record and end-of-file.
    integer :: close_status
    !! I/O status returned when closing the private input unit.
    integer :: groups
    !! Number of matching relativistic_params headers outside comments and quoted strings.
    !! Zero retains defaults; more than one is rejected as a duplicate active group.
    integer :: token_length
    !! Full length of the header name being scanned, including characters beyond the token buffer.
    integer :: i
    !! Character index used when normalizing whitespace and case in rel_cache.
    integer :: code
    !! ASCII character code used to recognize header-name characters and normalize the cache mode.
    character :: ch
    !! Current input character read without advancing to the next record.
    character :: quote
    !! Active single- or double-quote delimiter; a blank indicates scanning outside a quoted string.
    character(32) :: token
    !! Lowercase header-name buffer following an ampersand; token_length prevents truncated-name matches.
    character(1024) :: message
    !! Processor-provided explanation of an input opening, scanning, rewind, or namelist-read failure.
    logical :: comment
    !! True after an unquoted exclamation mark, until the end of the current record.
    logical :: in_name
    !! True while collecting a namelist header name after an unquoted, uncommented ampersand.
    namelist /relativistic_params/ ngamma,npparbar,rel_cache,rel_cache_dir

    open(newunit=scan_unit,file=trim(foldername)//trim(runname)//'.in', &
         status='old',action='read',iostat=ios,iomsg=message)
    if(ios/=0) then
       write(*,'(a)') trim(message)
       write(unit_error,'(a)') trim(message)
       call alps_error(17)
    endif
    groups=0; quote=' '; comment=.false.; in_name=.false.; token=''; token_length=0
    do
       read(scan_unit,'(a)',advance='no',iostat=ios,iomsg=message) ch
       if(ios==iostat_eor) then
          call finish_name
          comment=.false.
          cycle
       endif
       if(ios/=0) exit
       if(comment) cycle
       if(quote/=' ') then
          if(ch==quote) quote=' '
          cycle
       endif
       if(in_name) then
          code=iachar(ch)
          if((code>=iachar('a').and.code<=iachar('z')).or. &
               (code>=iachar('A').and.code<=iachar('Z')).or. &
               (code>=iachar('0').and.code<=iachar('9')).or.ch=='_') then
             token_length=token_length+1
             if(code>=iachar('A').and.code<=iachar('Z')) ch=achar(code+32)
             if(token_length<=len(token)) token(token_length:token_length)=ch
             cycle
          endif
          call finish_name
       endif
       select case(ch)
       case('!')
          comment=.true.
       case("'",'"')
          quote=ch
       case('&')
          in_name=.true.; token=''; token_length=0
       end select
    enddo
    call finish_name
    if(ios/=iostat_end) then
       write(*,'(a)') 'Cannot scan relativistic parameters: '//trim(message)
       write(unit_error,'(a)') 'Cannot scan relativistic parameters: '//trim(message)
       close(scan_unit,iostat=close_status)
       call alps_error(17)
    endif
    if(groups>1) then
       write(*,'(a)') 'Duplicate &relativistic_params blocks are not allowed.'
       write(unit_error,'(a)') 'Duplicate &relativistic_params blocks are not allowed.'
       close(scan_unit,iostat=close_status)
       call alps_error(17)
    endif
    if(groups==1) then
       rewind(scan_unit,iostat=ios,iomsg=message)
       if(ios==0) read(scan_unit,nml=relativistic_params,iostat=ios,iomsg=message)
       if(ios/=0) then
          write(*,'(a)') 'Cannot read &relativistic_params: '//trim(message)
          write(unit_error,'(a)') 'Cannot read &relativistic_params: '//trim(message)
          close(scan_unit,iostat=close_status)
          call alps_error(17)
       endif
    endif
    close(scan_unit,iostat=close_status)
    if(close_status/=0) call alps_error(17)

    do i=1,len(rel_cache)
       code=iachar(rel_cache(i:i))
       if(code>=9.and.code<=13) rel_cache(i:i)=' '
    enddo
    rel_cache=adjustl(rel_cache)
    do i=1,len_trim(rel_cache)
       code=iachar(rel_cache(i:i))
       if(code>=iachar('A').and.code<=iachar('Z')) rel_cache(i:i)=achar(code+32)
    enddo
    select case(trim(rel_cache))
    case('off','auto','refresh')
    case default
       call alps_error(15)
    end select
    if(ngamma<2.or.npparbar<4) call alps_error(14)
  contains
    subroutine finish_name
      !! Count a complete matching header name, then reset the host-associated token-scanning state.
      if(in_name.and.token_length==len('relativistic_params')) then
         if(token=='relativistic_params') groups=groups+1
      endif
      in_name=.false.; token=''; token_length=0
    end subroutine
  end subroutine relativistic_params_read


  subroutine map_read
    !!Reads in complex frequency map parameters.
    use alps_var, only : loggridw, loggridg
    use alps_var, only : omi, omf, gami, gamf, ni, nr
    implicit none

    nameList /maps/ &
         loggridw, loggridg, omi, omf, gami, gamf, ni, nr
    read (unit=unit,nml=maps)

  end subroutine map_read



  subroutine solution_read(ik)
    !!Reads in initial guesses for dispersion solutions.
    use alps_var, only : wroots
    implicit none

    integer, intent(in) :: ik
    !!Solution index.

    double precision :: g_om
    !!Guess for real solution \(\omega_{\textrm{r}/\Omega_p \).

    double precision :: g_gam
    !!Guess for imaginary solution \(\gamma/\Omega_p \).

    nameList /guess/ &
         g_om,g_gam
    read (unit=unit,nml=guess)
    wroots(ik)=cmplx(g_om,g_gam,kind(1.d0))
  end subroutine solution_read



  subroutine spec_read(is)
    !!Subroutine for reading in species parameters
    use alps_var, only : ns, qs, ms, n_fits
    use alps_var, only : relativistic, logfit, usebM
    use alps_var, only : ACmethod
    implicit none

    integer,intent(in) :: is
    !!Species index.

    double precision :: nn
    !! Read in value for relative density for \(f_j\).

    double precision :: qq
    !! Read in value for charge for \(f_j\).

    double precision :: mm
    !! Read in value for mass for \(f_j\).

    double precision :: AC_method
    !! Read in value for Analytic Continuation Method
    
    integer :: ff
    !!Read in value for number of fitted functions.

    logical :: relat=.false.
    !! Treat species as non-relativistic or relativistic.

    logical :: log_fit=.true.
    !! Use linear or \(\log_{10}\) fitting routine.

    logical :: use_bM=.false.
    !! Use actual numerical integration or bi-Maxwellian/cold-plasma proxy via NHDS.

    nameList /spec/ &
         nn,qq,mm,AC_method,ff,relat,log_fit,use_bM
    read (unit=unit,nml=spec)
    ns(is) = nn; qs(is) = qq; ms(is) = mm
    n_fits(is)=ff; relativistic(is)=relat
    logfit(is)=log_fit; usebM(is)=use_bM
    ACmethod(is)=AC_method
  end subroutine spec_read

  subroutine poly_read(is)
    !!Reads in Polynomial Basis Function Parameters
    use alps_var, only : poly_kind, poly_order
    use alps_var, only : poly_derivative_type, unit_error
    implicit none

    integer,intent(in) :: is
    !!Species index.

    integer :: kind
    !! Selection of Orthogonal Basis Function
    !! 1) Chebyshev polynomials

    integer :: order
    !! Maximum order of Polynomial

    integer :: read_status
    character(256) :: read_message

    integer :: derivative_type
    !! Kind of derivative for analytic continuation evaluation

    nameList /poly_spec/ &
         kind, order, derivative_type

   derivative_type=1
   
    read(unit=unit,nml=poly_spec,iostat=read_status,iomsg=read_message)
    if(read_status/=0) then
       if(index(read_message,'log_max')>0) then
          write(*,'(a)') 'ERROR: log_max is no longer supported in &poly_spec.'
          write(unit_error,'(a)') 'ERROR: log_max is no longer supported in &poly_spec.'
       else
          write(*,'(a,i0,2a)') 'ERROR: invalid &poly_spec for species ',is,': ',trim(read_message)
          write(unit_error,'(a,i0,2a)') 'ERROR: invalid &poly_spec for species ',is,': ',trim(read_message)
       endif
       call alps_error(10)
    endif
    poly_kind(is)=kind
    poly_order(is)=order
    poly_derivative_type(is)=derivative_type

  end subroutine poly_read

  subroutine bM_read(is)
    !!Reads in bi-Maxwellian/cold-plasma parameters.
    use alps_var, only : bMnmaxs,bMBessel_zeros,bMbetas
    use alps_var, only : bMalphas,bMpdrifts
    implicit none

    integer,intent(in) :: is
    !!Species index.

    integer :: bM_nmaxs
    !!Maximum number of resonances to consider.

    double precision :: bM_Bessel_zeros
    !!Precision threshold for \(I_n\).

    double precision :: bM_betas
    !!\(\beta_{\parallel,j}\) of biMaxwellian distribution \(f_j\).
    !! If bM_betas=0.d0, this species is treated with the cold-plasma susceptibilities.

    double precision :: bM_alphas
    !!\(T_{\perp,j}/T_{\parallel,j} \) of biMaxwellian distribution \(f_j\).

    double precision :: bM_pdrifts
    !!Relative drift of biMaxwellian distribution \(f_j\),
    !!in units of \(m_p v_{A,p}\). Also used as the drift of the cold-plasma species
    !! if bM_betas is set to 0.d0.

    nameList /bM_spec/ &
         bM_nmaxs,bM_Bessel_zeros,bM_betas,bM_alphas,bM_pdrifts

    read (unit=unit,nml=bM_spec)
    bMnmaxs(is)=bM_nmaxs; bMBessel_zeros(is)=bM_Bessel_zeros
    bMbetas(is)=bM_betas;
    bMalphas(is)=bM_alphas; bMpdrifts(is)=bM_pdrifts
  end subroutine bM_read



  subroutine scan_read(is)
    !!The most important subroutine.
    !!Reads in wavevector scan parameters.
    !!Defines [[scanner(type)]], which controls the
    !!behavior of the wavevector scan.
  use alps_var, only : scan,kperp_last,kpar_last,scan_option
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none

  integer, intent(in) :: is
  !!Scan index.

  integer :: scan_type
  !!Determine kind of wavevector scan.

  integer :: ns
  !!Number of output scan values.

  integer :: nres
  !!Resolution between output scan values.

  double precision :: swi
  !!Initial Scan Value.

  double precision :: swf
  !!Final Scan Value.

  logical :: swlog
  !!\(\log_{10}\) or linear scan spacing.

  logical :: heating
  !!Activate heating calculation.

  logical :: eigen
  !!Activate eigenfunction calculation.

  double precision :: theta_0
  !!\(\atan(k_\perp/k_\parallel)\)

  double precision :: k_0
  !!\(\sqrt{k_\perp^2+k_\parallel^2}\)

  integer :: scan_ios

  nameList /scan_input/ &
       scan_type,swi,swf,swlog,ns,nres,&
       heating,eigen
  if(scan_option==4.or.scan_option==5) then
     ! Defaults/validation are confined to the new operation; legacy readers retain behavior.
     scan_type=-1; ns=0; nres=1; swi=0; swf=ieee_value(0.d0,ieee_quiet_nan)
     swlog=.false.; heating=.false.; eigen=.false.
     scan(is)%diff2=0
  endif
  read (unit=unit,nml=scan_input,iostat=scan_ios)
  if(scan_ios/=0) then
     write(*,'(a,i0)') 'ERROR: cannot read scan_input_',is
     call alps_error(20)
  endif
  if(scan_option==4.or.scan_option==5) then
     if(scan_type<0.or.scan_type>4.or.ns<1.or.nres<1) then
        write(*,'(a)') 'ERROR: Nyquist scan requires scan_type=0..4 and ns,nres >= 1.'
        call alps_error(20)
     endif
     if(ns>huge(ns)/nres) call alps_error(20)
     if(scan_option==5.and.scan_type==0) then
        write(*,'(a)') 'ERROR: Nyquist plane does not support scan_type=0.'
        call alps_error(20)
     endif
     if(.not.all(ieee_is_finite([swi,swf,kperp_last,kpar_last]))) call alps_error(20)
     if(kperp_last<0.or.kpar_last==0) call alps_error(20)
     if(swlog) then
        select case(scan_type)
        case(0)
           if(min(swi,swf,kperp_last,kpar_last)<=0) call alps_error(20)
        case(1)
           if(swf<=0.or.atan(kperp_last/kpar_last)<=0) call alps_error(20)
        case(2)
           if(swf<=0) call alps_error(20)
        case(3)
           if(min(swf,kperp_last)<=0) call alps_error(20)
        case(4)
           if(min(swf,kpar_last)<=0) call alps_error(20)
        end select
     endif
  endif

  scan(is)%range_i =swi
  scan(is)%range_f =swf
  scan(is)%log_scan=swlog
  scan(is)%type_s  =scan_type
  scan(is)%n_out   =ns
  scan(is)%n_res   =nres
  scan(is)%eigen_s =eigen
  scan(is)%heat_s  =heating

  !Calculate step size:
  select case (scan_type)
  case(0)
     !Scan from k_0 to k_1:
     write(*,'(a,i0,a,es14.4e3,a,es14.4e3,a,es14.4e3,a,es14.4e3,a)')&
          'Scan ',is,': (kpar,kperp) from (',&
          kperp_last,',',kpar_last,') to (',swi,',',swf,')'
     if (swlog) then
        scan(is)%diff =(log10(swi)-log10(kperp_last))/&
             (1.d0*ns*nres)
        scan(is)%diff2=(log10(swf)-log10(kpar_last))/&
             (1.d0*ns*nres)
     else
        scan(is)%diff =(swi-kperp_last)/(1.d0*ns*nres)
        scan(is)%diff2=(swf-kpar_last )/(1.d0*ns*nres)
     endif
     kperp_last=swi;kpar_last=swf
  case(1)
     !Scan from theta_0 to theta_1:
     theta_0=atan(kperp_last/kpar_last)
     k_0=sqrt(kperp_last**2+kpar_last**2)
     write(*,'(a,i0,a,es14.4e3,a,es14.4e3)')&
          'Scan ',is,': theta from ',&
          theta_0*180.d0/(4.d0*atan(1.d0)),' to ',swf
     if (swlog) then
        scan(is)%diff=(log10(swf*4.d0*atan(1.d0)/180.)-&
             log10(theta_0))/(1.d0*ns*nres)
     else
        scan(is)%diff=((swf*4.d0*atan(1.d0)/180.)-theta_0)/&
             (1.d0*ns*nres)
     endif
     kpar_last=k_0*cos(swf*4.d0*atan(1.d0)/180.d0)
     kperp_last=k_0*sin(swf*4.d0*atan(1.d0)/180.d0)
  case(2)
     !Scan from |k_0| to |k_1| at constant theta.
     theta_0=atan(kperp_last/kpar_last)
     k_0=sqrt(kperp_last**2+kpar_last**2)
     write(*,'(a,i0,a,es14.4e3,a,es14.4e3,a,es14.4e3,a,es14.4e3,a)')&
          'Scan ',is,': |k| from ',&
          k_0,' to ',swf,' at theta=',theta_0*180.d0/&
          (4.d0*atan(1.d0))
     if (swlog) then
        scan(is)%diff =(log10(swf)-log10(k_0))/&
             (1.d0*ns*nres)
     else
        scan(is)%diff =(swf-k_0)/(1.d0*ns*nres)
     endif
     kpar_last=k_0*cos(theta_0)
     kperp_last=k_0*sin(theta_0)
     if(scan_option==4.or.scan_option==5) then
        ! The following segment starts at the requested magnitude endpoint.
        ! Leave the historical option-1/2 bookkeeping untouched.
        kpar_last=swf*cos(theta_0)
        kperp_last=swf*sin(theta_0)
     endif
  case(3)
     !Scan of kperp; kpar constant:
     write(*,'(a,i0,a,es14.4e3,a,es14.4e3)')&
          'Scan ',is,': kperp from ',kperp_last,' to ',swf
     if (swlog) then
        scan(is)%diff=(log10(swf)-log10(kperp_last))/(1.d0*ns*nres)
     else
        scan(is)%diff=(swf-kperp_last)/(1.d0*ns*nres)
     endif
     kperp_last=swf
  case(4)
     !Scan of kpar; kperp constant:
     write(*,'(a,i0,a,es14.4e3,a,es14.4e3)')&
          'Scan ',is,': kpar from ',kpar_last,' to ',swf
     if (swlog) then
        scan(is)%diff=(log10(swf)-log10(kpar_last))/(1.d0*ns*nres)
     else
        scan(is)%diff=(swf-kpar_last)/(1.d0*ns*nres)
     endif
     kpar_last=swf
  end select

end subroutine scan_read



subroutine fit_read(is,ifit)
  !!Reads in fit parameters for component is.
  use alps_var, only : fit_type, param_fit, perp_correction
  implicit none

  integer, intent(in) :: is
  !!Species index.

  integer, intent(in) :: ifit
  !Fit index.

  double precision :: fit_1
  !! Read in values for fit component 1.

  double precision :: fit_2
  !! Read in values for fit component 2.

  double precision :: fit_3
  !! Read in values for fit component 3.

  double precision :: fit_4
  !! Read in values for fit component 4.

  double precision :: fit_5
  !! Read in values for fit component 5.

  double precision :: perpcorr
  !!Perpendicular correction to compensate for exponential
  !!dependency of drift to make fit more reliable
  !![\(y\) in Eqn B1 of Verscharen et al 2018].

  integer :: fit_type_in
  !!Read in value for type of analytic function.

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



subroutine get_runname(runname,foldername)
  !! Get runname for output files from input argument.
  implicit none

  integer       :: l
  !!Dummy Length.

  integer :: pathend
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



  subroutine read_f0
    !! Subroutine for reading in background distribution function
    use alps_var, only : nperp, npar, arrayName, f0, pp, nspec
    use alps_var, only : writeOut, usebM,relativistic,ACmethod,unit_error
    implicit none

    integer :: ios
    integer :: ipar
    !Parallel index.

    integer :: iperp
    !!Perpendicular index.

    integer :: is
    !!Species index.

    character (100) :: readname
    !!Base I/O name.

     if (writeOut) &
          write(*,'(2a)')&
          'Attempting to read f0 array from file: ',trim(arrayName)

     call get_unused_unit (input_unit_no)
     unit=input_unit_no
     !The f0 arrays are stored in the distribution folder.
     !arrayName is read in from *.in input file.
     !each species is has a unique file for f0 as a function of pp.
     do is = 1, nspec

       if (usebM(is)) then
          write (*,'(a,i2)') ' Bi-Maxwellian/cold-plasma calculation: not reading f0 array for species ',is
          f0(is,:,:)=0.d0
          pp(is,:,:,:)=0.d0

       else

        write(readname, '(3a,i0,a)') &
             "distribution/",trim(arrayName),'.',is,".array"
        if (relativistic(is).and.ACmethod(is)==2) then
           open(unit=unit,file=trim(readname),status='old',action='read',iostat=ios)
           if(ios/=0) call rel_input_error(is,'cannot open source distribution: '//trim(readname))
        else
           open (unit=unit, file=trim(readname), status='old', action='read')
        endif

        do iperp=0,nperp
           do ipar=0,npar
              if(relativistic(is).and.ACmethod(is)==2) then
                 read(unit,*,iostat=ios) pp(is,iperp,ipar,1),pp(is,iperp,ipar,2),f0(is,iperp,ipar)
                 if(ios/=0) then
                    write(*,'(a,2i8)') ' Missing/malformed source sample at (iperp,ipar)=',iperp,ipar
                    write(unit_error,'(a,2i8)') ' Missing/malformed source sample at (iperp,ipar)=',iperp,ipar
                    call rel_input_error(is,'missing or malformed source sample.')
                 endif
              else
                 read(unit,*) pp(is,iperp,ipar,1),pp(is,iperp,ipar,2),f0(is,iperp,ipar)
              endif
           enddo
        enddo
        close(unit)
        if(relativistic(is).and.ACmethod(is)==2) call validate_rel_source(is)
      endif
     enddo

end subroutine read_f0


subroutine rel_input_error(is,message)
  !! Predictable input failures are distinct from kernel continuation statuses.
  use alps_var, only : unit_error
  implicit none
  integer,intent(in) :: is
  character(*),intent(in) :: message
  write(*,'(a,i0,2a)') 'ERROR: relativistic input species ',is,': ',trim(message)
  write(unit_error,'(a,i0,2a)') 'ERROR: relativistic input species ',is,': ',trim(message)
  call alps_error(18)
end subroutine rel_input_error

subroutine validate_rel_source(is)
  !! Before derivative preparation, TPS logarithms and cache lookup. The cache
  !! identity includes these exact source arrays, coordinates, mass and vA.
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use alps_var, only : pp,f0,ms,vA,ns,qs,nperp,npar,unit_error
  implicit none
  integer,intent(in) :: is
  integer :: i,j
  double precision :: tol
  if(.not.ieee_is_finite(ms(is)).or..not.ieee_is_finite(vA)) &
       call rel_input_error(is,'mass and vA must be finite.')
  if(ms(is)<=0.d0.or.vA<=0.d0) call rel_input_error(is,'mass and vA must be positive.')
  if(.not.ieee_is_finite(ns(is)).or..not.ieee_is_finite(qs(is))) &
       call rel_input_error(is,'density and charge normalization parameters must be finite.')
  if(nperp<2.or.npar<2) call rel_input_error(is,'source grid requires nperp>=2 and npar>=2 for original-grid preparation.')
  do i=0,nperp
     do j=0,npar
        if(.not.ieee_is_finite(f0(is,i,j))) then
           write(*,'(a,2i8)') ' Nonfinite source sample at (iperp,ipar)=',i,j
           write(unit_error,'(a,2i8)') ' Nonfinite source sample at (iperp,ipar)=',i,j
           call rel_input_error(is,'nonfinite source distribution sample.')
        endif
        if(f0(is,i,j)<=0.d0) then
           write(*,'(a,2i8,a,es24.15)') ' Source sample (iperp,ipar)=',i,j,' value=',f0(is,i,j)
           write(unit_error,'(a,2i8,a,es24.15)') ' Source sample (iperp,ipar)=',i,j,' value=',f0(is,i,j)
           call rel_input_error(is,'source distribution must be strictly positive for log-TPS; no zero floor is applied.')
        endif
        if(.not.all(ieee_is_finite(pp(is,i,j,:)))) call rel_input_error(is,'nonfinite source momentum coordinate.')
     enddo
  enddo
  tol=128.d0*epsilon(1.d0)*max(1.d0,maxval(abs(pp(is,:,:,:))))
  if(minval(pp(is,:,:,1))<0.d0) call rel_input_error(is,'perpendicular source momentum must be nonnegative.')
  do i=0,nperp
     if(i>0) then
        if(pp(is,i,0,1)<=pp(is,i-1,0,1)) call rel_input_error(is,'perpendicular grid must be strictly increasing.')
     endif
     do j=0,npar
        if(abs(pp(is,i,j,1)-pp(is,i,0,1))>tol.or.abs(pp(is,i,j,2)-pp(is,0,j,2))>tol) &
             call rel_input_error(is,'source coordinates must form a rectangular tensor-product grid.')
        if(j>0) then
           if(pp(is,i,j,2)<=pp(is,i,j-1,2)) call rel_input_error(is,'parallel grid must be strictly increasing.')
        endif
     enddo
  enddo
end subroutine validate_rel_source





subroutine get_indexed_namelist_unit (unit, nml, index_in)
  !!Determines unused I/O unit.
  use alps_var, only : runname
  implicit none

  integer, intent (out) :: unit
  !!Unit to be defined.

  character (*), intent (in) :: nml
  !!Character string for namelist to be read in.

  integer, intent (in) :: index_in
  !!Index of namelist to be read in.

  character(500) :: line
  !!I/O dummy variable.

  integer :: iunit, iostat, in_file
  !!I/O dummy indices.

  integer :: ind_slash
  !!I/O dummy index.

  logical :: exist
  !!Check if namelist is open.

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






subroutine get_indexed_double_namelist_unit (unit, nml, spec_in, index_in)
  !!A version of [[get_indexed_namelist_unit(subroutine)]], extended
  !!to allow for double indexing in order to read in multiple fits
  !!for a single species.
  use alps_var, only : runname

  implicit none
  integer, intent (out) :: unit
  !!Unit to be defined.

  character (*), intent (in) :: nml
  !!Character string for namelist to be read in.

  integer, intent (in) :: spec_in
  !!First index of namelist to be read in.

  integer, intent (in) :: index_in
  !!Second index of namelist to be read in.

  character(500) :: line,lines
  !!I/O dummy variable.

  integer :: iunit, iostat, in_file
  !!I/O dummy indices.

  integer :: ind_slash
  !!I/O dummy index.

  logical :: exist
  !!Check if namelist is open.

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
  !!Determine if a particular namelist already opened.
  implicit none

  character(*), intent (in) :: nml
  !!Namelist to be opened.

  logical, intent(out) :: exist
  !!Determination if namelist is open.

  integer :: input_unit_exist, iostat
  !!I/O dummy indices.

  character(500) :: line
  !!I/O dummy variable.

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
    !!Assigns input unit for namelist opening.
    implicit none

    character(*), intent (in) :: nml
    !! Namelist string.

    integer :: input_unit, iostat
    !!I/O dummy indices.

    character(500) :: line
    !!I/O dummy variable.

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
    !!Determine unused number for I/O index.
    implicit none

    integer, intent (out) :: unit
    !!Unit to be assigned.

    logical :: od
    unit = 50
    do
       inquire (unit=unit, opened=od)
       if (.not.od) return
       unit = unit + 1
    end do
  end subroutine get_unused_unit

  subroutine alps_error_init
    !!Open a file for the error log.
    use alps_var, only : unit_error, runname, foldername
    implicit none

    call get_unused_unit(unit_error)
    !Get the run name, which comes from the name
    !of the input file appended after the executable:
    !mpirun -np 8 ./alps.e sample.in
    !yields a run name of 'sample'
    call get_runname(runname,foldername)
    open (unit=unit_error,file=trim(foldername)//trim(runname)//".log",status='replace')
  end subroutine alps_error_init



  subroutine alps_error(error_id)
    !!Error catching subroutine.
    use alps_var, only : ierror,unit_error,nproc,scan_option
    use mpi

    implicit none
    integer :: error_id
    !!Index of error message.

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
       case(10) !seen by proc0
          write (*,'(a)') "ERROR: Unspecified Orthogonal Representation."
          write(unit_error,'(a)') "ERROR: Unspecified Orthogonal Representation."
      case(11) !seen by proc0
          write (*,'(a)') &
            "ERROR: Polynomial Analytic derivative unavailable for this AC method "
          write(unit_error,'(a)') &
            "ERROR: Polynomial Analytic derivative unavailable for this AC method"
       case(12)
          write(*,'(a)') 'ERROR: Relativistic spline DGESV received an invalid argument.'
          write(unit_error,'(a)') 'ERROR: Relativistic spline DGESV received an invalid argument.'
       case(13)
          write(*,'(a)') 'ERROR: Relativistic spline matrix is singular; interpolation stopped.'
          write(unit_error,'(a)') 'ERROR: Relativistic spline matrix is singular; interpolation stopped.'
       case(14)
          write(*,'(a)') 'ERROR: Relativistic derivative stencils require ngamma>=2 and npparbar>=4.'
          write(unit_error,'(a)') 'ERROR: Relativistic derivative stencils require ngamma>=2 and npparbar>=4.'
       case(15)
          write(*,'(a)') 'ERROR: rel_cache must be off, auto, or refresh.'
          write(unit_error,'(a)') 'ERROR: rel_cache must be off, auto, or refresh.'
       case(16)
          write(*,'(a)') 'ERROR: Invalid &system. Move ngamma, npparbar, rel_cache and rel_cache_dir to &relativistic_params.'
          write(unit_error,'(a)') 'Move ngamma, npparbar, rel_cache and rel_cache_dir to &relativistic_params.'
          write(*,'(a)') 'Move maxsteps_fit, lambda_initial_fit, lambdafac_fit, and epsilon_fit to &fitting_params.'
          write(unit_error,'(a)') &
               'Move maxsteps_fit, lambda_initial_fit, lambdafac_fit, and epsilon_fit to &fitting_params.'
       case(17)
          write(*,'(a)') 'ERROR: Cannot read one valid optional &relativistic_params block.'
          write(unit_error,'(a)') 'ERROR: Cannot read one valid optional &relativistic_params block.'
       case(18)
          write(*,'(a)') 'ERROR: Unsupported or invalid relativistic input; see the species diagnostic above.'
          write(unit_error,'(a)') 'ERROR: Unsupported or invalid relativistic input; see the species diagnostic above.'
       case(19)
          write(*,'(a)') 'ERROR: Cannot read one valid optional global &fitting_params block.'
          write(unit_error,'(a)') 'ERROR: Cannot read one valid optional global &fitting_params block.'
       case(20)
          write(*,'(a)') 'ERROR: Invalid wavevector path; check n_scan, scan_input_N, counts and finite spacing endpoints.'
          write(unit_error,'(a)') &
               'ERROR: Invalid wavevector path; check n_scan, scan_input_N, counts and finite spacing endpoints.'
       case default
          write(*,'(a)')'ERROR: Unspecified...'
          write(unit_error,'(a)')'ERROR: Unspecified...'
       end select
       write(*,'(a)') 'Finishing ALPS=================================='
       write(*,'(a)') "Time:"
       call output_time
       write(*,'(a)')'================================================'

    close(unit_error)
    call mpi_abort(MPI_COMM_WORLD,error_id, ierror)
  stop

  end subroutine alps_error





  logical function isnancheck(input)
    !!Checks if double precision number input is NaN.
    implicit none

    double precision :: input
    !! Variable to be checked.

    isnancheck=.FALSE.
    if (abs(input).GE.huge(1.d0)) isnancheck=.TRUE.
    if (input.NE.input) isnancheck=.TRUE.

  end function isnancheck




  subroutine output_time
    !!Outputs the date and time in a given format using intrinsic
    !!FORTRAN function.
    implicit none

    character(8)  :: date
    !!Date.

    character(10) :: time
    !!Time.

    character(5)  :: zone
    !!Time Zone.

    integer,dimension(8) :: value
    !!Output Time Values.

    call date_and_time(date,time,zone,value)
    write (*,'(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') value(1),"-",value(2),"-",value(3)," -- ",value(5),":",value(6),":",value(7)

  end subroutine output_time



  subroutine display_credits
    !!Writes the opening credits.
    implicit none

    write (*,*) "==========================================================="
    write (*,*) "I                                                         I"
    write (*,*) "I                       A  L  P  S                        I"
    write (*,*) "I              Arbitrary Linear Plasma Solver             I"
    write (*,*) "I                                                         I"
    write (*,*) "I                       Version 1.0                       I"
    write (*,*) "I                                                         I"
    write (*,*) "I  Kristopher Klein   (kgklein@arizona.edu)               I"
    write (*,*) "I  Daniel Verscharen  (daniel.verscharen@rub.de)            I"
    write (*,*) "I                                                         I"
    write (*,*) "==========================================================="


  end subroutine display_credits

  subroutine decimal_coordinate_radius(token,radius,status)
    ! Decimal serialization quantum: one half unit in the last written place.
    ! Integers have exact integer semantics; unsupported spellings are rejected.
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    implicit none
    character(*),intent(in)::token
    real(8),intent(out)::radius
    integer,intent(out)::status
    character(:),allocatable::s,mant
    integer::i,e,d,exponent_value,ios,n,start
    s=trim(adjustl(token));status=1;radius=0;exponent_value=0
    n=len(s);if(n==0)return
    e=scan(s,'eEdD');mant=s
    if(e>0)then
      if(e==n)return
      read(s(e+1:),*,iostat=ios)exponent_value
      if(ios/=0)return
      mant=s(:e-1)
    endif
    d=index(mant,'.');start=1
    if(len(mant)==0)return
    if(mant(1:1)=='+'.or.mant(1:1)=='-')start=2
    if(len(mant)<start)return
    do i=start,len(mant)
      if(i==d)cycle
      if(mant(i:i)<'0'.or.mant(i:i)>'9')return
    enddo
    if(d>0)then
      if(len(mant)==1.or.(len(mant)==2.and.start==2))return
      exponent_value=exponent_value-(len(mant)-d)
      if(exponent_value>300.or.exponent_value< -300)return
      radius=.5d0*10.d0**exponent_value
    elseif(e>0)then
      if(exponent_value>300.or.exponent_value< -300)return
      radius=.5d0*10.d0**exponent_value
    endif
    if(.not.ieee_is_finite(radius))return
    status=0
  end subroutine

  subroutine fit_rounded_axis(x,lower,upper,nonnegative,canonical,step,residual,unchanged,status,reason)
    ! An affine model must intersect EVERY supplied decimal rounding interval.
    ! Pairwise interval inequalities give a feasible slope interval, then intercept.
    ! No least-squares residual inflation is used to make an infeasible grid pass.
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    implicit none
    real(8),intent(in)::x(0:),lower(0:),upper(0:)
    logical,intent(in)::nonnegative
    real(8),intent(out)::canonical(0:),step,residual
    logical,intent(out)::unchanged
    integer,intent(out)::status
    character(*),intent(out)::reason
    integer::n,i,j
    real(8)::centre,mean,b,a,den,lo,hi,alo,ahi,fp,zlo,zhi
    n=ubound(x,1);status=1;reason='invalid axis';unchanged=.false.;step=0;residual=0
    if(n<1.or..not.all(ieee_is_finite(x)))return
    if(any(x(1:)<=x(:n-1)))then
      reason='axis is not strictly increasing';return
    endif
    if(nonnegative.and.x(0)<0)then
      reason='negative perpendicular momentum';return
    endif
    fp=32*epsilon(1.d0)*max(maxval(abs(x)),tiny(1.d0))
    centre=.5d0*n;mean=sum(x)/(n+1);den=0;b=0
    do i=0,n
      den=den+(i-centre)**2;b=b+(i-centre)*(x(i)-mean)
    enddo
    b=b/den;a=mean-b*centre
    do i=0,n
      canonical(i)=a+i*b
    enddo
    residual=maxval(abs(x-canonical));step=b
    ! Preserve already uniform high-precision input bits and hence legacy science.
    if(residual<=fp.and.all(x>=lower-fp).and.all(x<=upper+fp))then
      canonical=x;unchanged=.true.;status=0;reason='uniform to floating-point roundoff';return
    endif
    if(any(lower>upper))then
      reason='inconsistent rectangular coordinate intervals';return
    endif
    lo=0;hi=huge(1.d0)
    do i=0,n-1
      do j=i+1,n
        lo=max(lo,(lower(j)-upper(i))/(j-i))
        hi=min(hi,(upper(j)-lower(i))/(j-i))
      enddo
    enddo
    ! Nonnegative origin is a feasibility constraint on spacing, not a
    ! post-fit rejection of one otherwise arbitrary least-squares candidate.
    if(nonnegative)then
      if(upper(0)<0)then
        reason='no nonnegative perpendicular origin';return
      endif
      do i=1,n
        hi=min(hi,upper(i)/i)
      enddo
    endif
    if(hi<=0.or.lo>hi)then
      reason='no positive affine spacing fits all rounding intervals';return
    endif
    b=max(lo,min(hi,b));alo=-huge(1.d0);ahi=huge(1.d0)
    do i=0,n
      alo=max(alo,lower(i)-i*b);ahi=min(ahi,upper(i)-i*b)
    enddo
    if(nonnegative)alo=max(0.d0,alo)
    if(alo>ahi)then
      reason='no affine origin fits all rounding intervals';return
    endif
    a=max(alo,min(ahi,a))
    ! Preserve a written zero perpendicular origin when a feasible line permits it.
    if(nonnegative.and.x(0)==0.and.lower(0)<=0.and.upper(0)>=0)then
      zlo=0;zhi=huge(1.d0)
      do i=1,n
        zlo=max(zlo,lower(i)/i);zhi=min(zhi,upper(i)/i)
      enddo
      if(zlo<=zhi.and.zhi>0)then
        a=0;b=max(zlo,min(zhi,b))
      endif
    endif
    do i=0,n
      canonical(i)=a+i*b
    enddo
    step=b;residual=maxval(abs(x-canonical))
    if(any(canonical<lower-fp).or.any(canonical>upper+fp))then
      reason='canonical coordinate exceeds a rounding interval';return
    endif
    status=0;reason='decimal rounding model accepted'
  end subroutine

  subroutine canonicalize_nyquist_grid(species)
    ! Rank zero only, before derivative_f0 and all continuation initialization.
    ! Re-read coordinate TOKENS for their precision, never alter any VDF value.
    use alps_var,only:pp,nperp,npar,arrayName,unit_error
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    implicit none
    integer,intent(in)::species
    integer::u,ios,i,j,d,n,code,k,beg,finish,idx
    real(8),allocatable::x(:),lower(:),upper(:),canon(:),radius_grid(:,:)
    real(8)::r(2),value(2),fp,step,residual,maxradius(2),local_error,max_error,local_bound
    logical::unchanged
    character(1024)::source,line,token(3)
    character(160)::reason
    write(source,'(3a,i0,a)')'distribution/',trim(arrayName),'.',species,'.array'
    do d=1,2
      n=nperp;if(d==2)n=npar
      allocate(x(0:n),lower(0:n),upper(0:n),canon(0:n),radius_grid(0:nperp,0:npar))
      lower=-huge(1.d0);upper=huge(1.d0);maxradius=0
      if(d==1)x=pp(species,:,0,1)
      if(d==2)x=pp(species,0,:,2)
      open(newunit=u,file=trim(source),status='old',action='read',iostat=ios)
      if(ios/=0)call grid_error('cannot reopen VDF coordinate tokens')
      do i=0,nperp
      do j=0,npar
        ! Ordinary array format: one record per node, three whitespace/comma fields.
        do
          read(u,'(a)',iostat=ios)line
          if(ios/=0)call grid_error('missing coordinate record')
          if(len_trim(line)>0)exit
        enddo
        beg=1
        do k=1,3
          do while(beg<=len_trim(line))
            if(index(' ,'//achar(9),line(beg:beg))==0)exit
            beg=beg+1
          enddo
          finish=beg
          do while(finish<=len_trim(line))
            if(index(' ,'//achar(9),line(finish:finish))>0)exit
            finish=finish+1
          enddo
          if(finish==beg)call grid_error('expected three fields in a VDF record')
          token(k)=line(beg:finish-1);beg=finish+1
        enddo
        read(token(d),*,iostat=ios)value(d)
        if(ios/=0)call grid_error('invalid coordinate token')
        if(value(d)/=pp(species,i,j,d))call grid_error('coordinate token/read mismatch')
        call decimal_coordinate_radius(trim(token(d)),r(d),code)
        if(code/=0)call grid_error('unsupported decimal coordinate spelling or exponent')
        fp=32*epsilon(1.d0)*max(abs(value(d)),maxval(abs(x)),tiny(1.d0))
        idx=i;if(d==2)idx=j
        lower(idx)=max(lower(idx),value(d)-r(d)-fp)
        upper(idx)=min(upper(idx),value(d)+r(d)+fp)
        maxradius(d)=max(maxradius(d),r(d));radius_grid(i,j)=r(d)
        if(i>0.and.d==1)then
          if(pp(species,i,j,d)<=pp(species,i-1,j,d))call grid_error('nonmonotonic perpendicular row')
        endif
        if(j>0.and.d==2)then
          if(pp(species,i,j,d)<=pp(species,i,j-1,d))call grid_error('nonmonotonic parallel row')
        endif
      enddo
      enddo
      close(u)
      ! Rectangular intersections are checked even on the high-precision fast path.
      if(any(lower>upper))call grid_error('inconsistent rectangular coordinates')
      call fit_rounded_axis(x,lower,upper,d==1,canon,step,residual,unchanged,code,reason)
      write(*,'(a,i0,a,i0,a,3es24.16,a,l1)')'NYQUIST_GRID species=',species,' axis=',d, &
        ' spacing/max_axis_residual/max_decimal_radius=',step,residual,maxradius(d),' unchanged=',unchanged
      if(code/=0)call grid_error(trim(reason))
      ! Every coordinate and local increment must satisfy its own rounding bound.
      fp=32*epsilon(1.d0)*max(maxval(abs(x)),tiny(1.d0));max_error=0;local_error=0
      do i=0,nperp
      do j=0,npar
        idx=i;if(d==2)idx=j
        max_error=max(max_error,abs(pp(species,i,j,d)-canon(idx)))
        if(abs(pp(species,i,j,d)-canon(idx))>radius_grid(i,j)+2*fp) &
          call grid_error('coordinate residual exceeds its decimal rounding bound')
        if(d==1.and.i>0)then
          local_bound=radius_grid(i,j)+radius_grid(i-1,j)+2*fp
          local_error=max(local_error,abs(pp(species,i,j,d)-pp(species,i-1,j,d)-step))
          if(abs(pp(species,i,j,d)-pp(species,i-1,j,d)-step)>local_bound) &
            call grid_error('perpendicular increment exceeds summed rounding bounds')
        endif
        if(d==2.and.j>0)then
          local_bound=radius_grid(i,j)+radius_grid(i,j-1)+2*fp
          local_error=max(local_error,abs(pp(species,i,j,d)-pp(species,i,j-1,d)-step))
          if(abs(pp(species,i,j,d)-pp(species,i,j-1,d)-step)>local_bound) &
            call grid_error('parallel increment exceeds summed rounding bounds')
        endif
      enddo
      enddo
      write(*,'(a,3es24.16)')'NYQUIST_GRID max_all_coordinate/max_local_step_residual/fp_allowance=', &
        max_error,local_error,fp
      do i=0,nperp
      do j=0,npar
        if(d==1)pp(species,i,j,1)=canon(i)
        if(d==2)pp(species,i,j,2)=canon(j)
      enddo
      enddo
      write(*,'(a,a)')'NYQUIST_GRID ',trim(reason)
      deallocate(x,lower,upper,canon,radius_grid)
    enddo
  contains
    subroutine grid_error(message)
      use mpi
      character(*),intent(in)::message
      integer::ierr
      write(*,'(a,i0,a,i0,2a)')'ERROR: numerical grid species=',species,' axis=',d,': ',trim(message)
      write(unit_error,'(a)')trim(message)
      call mpi_abort(MPI_COMM_WORLD,31,ierr)
      error stop 'Invalid numerical Nyquist grid'
    end subroutine
  end subroutine

end module alps_io
