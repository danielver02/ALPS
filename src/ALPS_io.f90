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
  private :: poly_read
  public  :: init_param, read_f0, get_unused_unit, alps_error
  public :: output_time, display_credits, isnancheck

contains



  subroutine init_param
    !!Read in system parameters from *.in file.
    !!Only processor 0 calls this routine:
    use alps_var, only : runname, foldername, kperp, kpar, nroots, D_prec, D_gap
    use alps_var, only : kperp_last, kpar_last, kperp_0, kpar_0
    use alps_var, only : use_map, writeOut, wroots, nspec, numroots
    use alps_var, only : nperp, npar, arrayName, fit_check, param_fit, fit_type, perp_correction
    use alps_var, only : ns, qs, ms, vA, Bessel_zero, numiter, D_threshold,positions_principal
    use alps_var, only : determine_minima, n_resonance_interval, ngamma, npparbar, Tlim
    use alps_var, only : scan_option, n_scan, scan, relativistic, logfit, usebM
    use alps_var, only : maxsteps_fit, n_fits, lambda_initial_fit, lambdafac_fit, epsilon_fit
    use alps_var, only : bMnmaxs, bMBessel_zeros, bMbetas, bMalphas, bMpdrifts
    use alps_var, only : ACmethod, poly_kind, poly_order, polynomials, poly_fit_coeffs
    use alps_var, only : poly_log_max
    implicit none

    integer :: ik
    !!Solution index for iterating [[(solution_read(subroutine)]].

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
         kperp, kpar, nspec, nroots, use_map, writeOut,&
         nperp, npar, ngamma, npparbar, vA, arrayName, Bessel_zero, numiter, D_threshold, &
         D_prec, D_gap, positions_principal, Tlim, &
         maxsteps_fit, lambda_initial_fit, lambdafac_fit, epsilon_fit, fit_check, &
         determine_minima, n_resonance_interval, scan_option, n_scan

    !Get a unassigned unit number for input/output:
    call get_unused_unit (input_unit_no)

    !Read in system parameters.
    !runname called earlier in alps_error_init:
    unit=input_unit_no
    open (unit=unit,file=trim(foldername)//trim(runname)//".in",status='old',action='read')
    read (unit=unit,nml=system)

    if (writeOut) &
         write(*,'(2a)') &
         'Reading from Input File: ', trim(runname)

    !save initial kperp,kpar values:
    kperp_last=kperp; kpar_last=kpar
    kperp_0=kperp; kpar_0=kpar

    !Allocate solution space for nroots dispersion solutions;
    allocate(wroots(1:numroots));wroots=cmplx(0.d0,0.d0,kind(1.d0))

    !Read in map scan parameters:
    if (use_map) then
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
    allocate(poly_log_max(1:nspec)); poly_log_max = 0

    !READ IN SPECIES PARAMETERS:
    do is = 1, nspec
       unit=input_unit_no
       call get_indexed_namelist_unit (unit, "spec", is)
       call spec_read(is)
       write(*,'(a,i3,a)')'Species ',is,' : '
       write(*,'(a,es14.4e3,a,es14.4e3,a,es14.4e3)')&
            ' ns/nREF = ',ns(is),' | qs/qREF = ',qs(is),' | ms/mREF = ',ms(is)
       select case (ACmethod(is))
       case (0)
          write(*,'(a)') ' Using function defined in distribution/distribution_analyt.f90'
       case (1)
          write(*,'(a,i4)')&
               ' Number of fitted functions = ',n_fits(is)
       case (2)
          write(*,'(a)')&
               ' Using a Polynomial Basis Representation'
       end select
       write(*,'(a,l1)')&
            ' Relativistic effects = ',relativistic(is)
       close(unit)
    enddo

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
             select case (poly_kind(is))
             case (1)
                write(*,'(a,i0,a,i0)')&
                     'Chebyshev Representation of Order ',poly_order(is),' for species ',is
             case default
                call alps_error(10)
             end select
             close(unit)
          end select
       endif
    enddo

    allocate(polynomials(1:nspec,0:npar,0:maxval(poly_order(:)))); polynomials=0.d0
    allocate(poly_fit_coeffs(1:nspec,0:nperp,0:maxval(poly_order(:)))); poly_fit_coeffs=0.d0

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
       write(*,'(a)')'-=-=-=-=-=-=-=-=-'
    endif
    close(unit)
  end subroutine init_param



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
    use alps_var, only : poly_kind, poly_order, poly_log_max
    implicit none

    integer,intent(in) :: is
    !!Species index.

    integer :: kind
    !! Selection of Orthogonal Basis Function
    !! 1) Chebyshev polynomials

    integer :: order
    !! Maximum order of Polynomial

    double precision :: log_max
    !! Maximum Value of Polynomial Evaluation

    nameList /poly_spec/ &
         kind, order, log_max

    read(unit=unit,nml=poly_spec)
    poly_kind(is)=kind
    poly_order(is)=order
    poly_log_max(is)=log_max

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
  use alps_var, only : scan,kperp_last,kpar_last
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
        scan(is)%diff =(log10(swf*sin(theta_0))-&
             log10(kperp_last))/(1.d0*ns*nres)
        scan(is)%diff2=(log10(swf*cos(theta_0))-&
             log10(kpar_last))/(1.d0*ns*nres)
     else
        scan(is)%diff =(swf*sin(theta_0)-kperp_last)/(1.d0*ns*nres)
        scan(is)%diff2=(swf*cos(theta_0)-kpar_last )/(1.d0*ns*nres)
     endif
     kpar_last=swf*cos(theta_0)
     kperp_last=swf*sin(theta_0)
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
    use alps_var, only : writeOut, usebM
    implicit none

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
        open (unit=unit, file=trim(readname), status='old', action='read')

        do iperp=0,nperp
           do ipar=0,npar
              read(unit,*) pp(is,iperp,ipar,1),pp(is,iperp,ipar,2),f0(is,iperp,ipar)
           enddo
        enddo
        close(unit)
      endif
     enddo

end subroutine read_f0





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
    write (*,*) "I  Kristopher Klein   (kgklein@arizone.edu)               I"
    write (*,*) "I  Daniel Verscharen  (d.verscharen@ucl.ac.uk)            I"
    write (*,*) "I                                                         I"
    write (*,*) "==========================================================="


  end subroutine display_credits

end module alps_io
