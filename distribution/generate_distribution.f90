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

program generate_distribution
  !! This program generates distribution files for input into ALPS based on pre-defined functions.
  use alps_distribution_analyt, only : distribution_analyt
  implicit none

  integer :: iperp
  !! Index to loop over perpendicular momentum.

  integer :: ipar
  !! Index to loop over parallel momentum.

  integer :: nperp
  !! Number of steps in perpendicular momentum.

  integer :: npar
  !! Number of steps in parallel momentum.

  integer :: nspec
  !! Total number of plasma species.

  integer :: is
  !! Index to loop over plasma species.

  integer, dimension(:), allocatable :: unit_out
  !! Index for unit for output. (1:nspec)

  character (500) :: writeName
  !! String for part of file name.

  character (500) :: outputName
  !! String for part of file name.

  double precision :: dpperp
  !! Inifinitesimal step in perpendicular momentum.

  double precision :: dppar
  !! Inifinitesimal step in parallel momentum.

  double precision :: ppar_max
  !! Maximum parallel momentum.

  double precision :: pperp_max
  !! Maximum perpendicular momentum.

  double precision :: iperpcorr
  !! Perpendicular correction \(y\).

  double precision :: ifit_1
  !! Ideal first fit parameter \(u_1\).

  double precision :: ifit_2
  !! Ideal second fit parameter \(u_2\).

  double precision :: ifit_3
  !! Ideal third fit parameter \(u_3\).

  double precision :: ifit_4
  !! Ideal fourth fit parameter \(u_4\).

  double precision :: ifit_5
  !! Ideal fifth fit parameter \(u_5\).

  double precision :: ppar
  !! Parallel momentum.

  double precision :: pperp
  !! Perpendicular momentum.

  double precision :: f0
  !! Distribution function.

  double precision, dimension(:), allocatable :: ms
  !! Mass of species. (1:nspec)

  double precision, dimension(:), allocatable :: tau
  !! Parallel temperature ratio of species. (1:nspec)

  double precision, dimension(:), allocatable :: alph
  !! Temperature anisotropy of species. (1:nspec)

  double precision, dimension(:), allocatable :: p_drift
  !! Drift momentum of species. (1:nspec)

  double precision, dimension(:), allocatable :: kappa
  !! Kappa index of species. (1:nspec)

  double precision, dimension(:), allocatable :: maxPperp
  !! Maximum perpendicular momentum of species. (1:nspec)

  double precision, dimension(:), allocatable ::  maxPpar
  !! Maximum parallel momentum of species. (1:nspec)

  logical, dimension (:), allocatable :: autoscale
  !! Flag for autoscaling of grid for species. (1:nspec)

  integer, dimension(:), allocatable :: distribution
  !! Type of distribution for species. (1:nspec)

  double precision :: beta
  !! Plasma beta.

  double precision :: vA
  !! Ratio of Alfven speed to c.

  double precision :: maxP
  !! Maximum momentum.

  double precision :: pi
  !! Pi.

  double precision :: integrate
  !! Integration variable.

  double precision :: norm
  !! Normalisation factor for f0.

  double precision :: a
  !! Helper variable for parts of distributions.

  double precision :: BESSK
  !! Modified Bessel function.

  double complex :: ppar_C
  !! Complex version of parallel momentum.

  integer :: unit
  !! Unit for file i/o.

  integer, parameter :: stdout_unit=6
  !! Flag for file i/o.

  integer, save :: input_unit_no
  !! Saved input unit for use with multiple read in calls.

  integer, save :: error_unit_no=stdout_unit
  !! Error output unit.

  character(500) :: runname
  !! String parameter for input file.


  pi = atan(1.d0)*4.d0

  call read_in_params

  allocate(unit_out(1:nspec)); unit_out = 0

  do is = 1 , nspec
	integrate =0.d0

     unit_out(is) = is
     write(outputName,'(2a,i0)')&
          trim(writeName),'.',is
     open(unit=unit_out(is),file=trim(outputName)//".array",status='replace')


	! Determine the normalization:
	select case(distribution(is))
    case (0) ! The code will use distribution_analyt, no normalisation necessary
      norm=1.d0

      pperp_max = maxPperp(is)
      ppar_max = maxPpar(is)


	  case (1) ! bi-Maxwellian
	  	norm = pi**(-1.5d0) /((ms(is) * beta * tau(is))**(3.d0/2.d0)*alph(is) )

      if (autoscale(is)) then
		      pperp_max = maxP*sqrt(ms(is) * tau(is) * alph(is))
     	    ppar_max  = maxP*sqrt(ms(is) * tau(is))
        else
          pperp_max = maxPperp(is)
          ppar_max = maxPpar(is)
      endif

      ! These are the ideal fitting parameters for this configuration:
      ifit_1=norm
      ifit_2=1.d0/( beta * ms(is) * tau(is))
      ifit_3=p_drift(is)
      ifit_4=0.d0
      iperpcorr=1.d0/( tau(is) *beta * ms(is) * alph(is))

	  case (2) ! bi-kappa
	    a = sqrt((2.d0*kappa(is)-3.d0)/(2.d0*kappa(is)))
	    norm = 1.d0/((ms(is) * beta * tau(is) * pi * kappa(is))**(3.d0/2.d0)*alph(is) )
       ! gamma(x) is the GAMMA-function:
	    norm = norm*gamma(kappa(is)+1.d0)/(gamma(kappa(is)-0.5d0)*a**3)

      if (autoscale(is)) then
		      pperp_max = maxP*sqrt(ms(is) * tau(is) * alph(is)) * sqrt((kappa(is)-1.5d0)/(kappa(1)-1.5d0))
     	    ppar_max  = maxP*sqrt(ms(is) * tau(is)) * sqrt((kappa(is)-1.5d0)/(kappa(1)-1.5d0))
        else
          pperp_max = maxPperp(is)
          ppar_max = maxPpar(is)
      endif

      ! These are the ideal fitting parameters for this configuration:
      ifit_1=norm
      ifit_2=1.d0/(beta * ms(is) * kappa(is) * a * a * tau(is))
      ifit_3=p_drift(is)
      ifit_4=(-1.d0-kappa(is))
      ifit_5=1.d0
      iperpcorr=1.d0/( tau(is) *beta * ms(is) * kappa(is) * a * a * alph(is))

	  case (3) ! Juettner
	    norm = vA/(2.d0*pi*sqrt(alph(is))*ms(is)**2 * beta * tau(is) )
	    norm = norm / BESSK(2,2.d0 * ms(is)/(vA*vA*alph(is)*beta*tau(is)))

      if (autoscale(is)) then
        	pperp_max = sqrt(maxP*maxP*tau(is) + (tau(is)-ms(is)*ms(is))/(vA*vA)) * sqrt(alph(is))
         	ppar_max  = sqrt(maxP*maxP*tau(is) + (tau(is)-ms(is)*ms(is))/(vA*vA))
       else
         pperp_max = maxPperp(is)
         ppar_max = maxPpar(is)
      endif

      ! These are the ideal fitting parameters for this configuration:
      ifit_1=norm
      ifit_2=vA*vA*alph(is)/(ms(is)*ms(is))
      ifit_5=(1.d0-tau(is))/( 2.d0*beta * ms(is))
      ifit_3=p_drift(is)
      ifit_4=0.d0
      ifit_5=0.d0
      iperpcorr=(2.d0*ms(is)/(vA*vA*beta*tau(is)*alph(is)))

    case (4) ! bi-Moyal
      norm = 1.d0

      if (autoscale(is)) then
		      pperp_max = maxP*sqrt(ms(is) * tau(is) * alph(is))
     	    ppar_max  = maxP*sqrt(ms(is) * tau(is))
        else
          pperp_max = maxPperp(is)
          ppar_max = maxPpar(is)
      endif

      ! These are the ideal fitting parameters for this configuration:

      ifit_1=1.d0
      ifit_2=1.d0/( beta * ms(is) * tau(is))
      ifit_3=p_drift(is)
      ifit_4=1.d0
      ifit_5=0.d0
      iperpcorr=1.d0/( tau(is) *beta * ms(is) * alph(is))

	end select

	   dpperp = pperp_max/real(nperp)
     dppar  = 2.d0*ppar_max/real(npar)

     do iperp=0,nperp
        pperp = real(iperp)*dpperp
        do ipar=0,npar
           ppar = real(ipar)*dppar - ppar_max + p_drift(is)

          select case(distribution(is))

          case (0) ! use function from distribution_analyt

            ppar_C=cmplx(ppar,0.d0,kind(1.d0))
            f0 = real(distribution_analyt(is,pperp,ppar_C))

	        case (1) ! bi-Maxwellian
        	   f0 =  exp( -(( (ppar-p_drift(is))**2.d0)/&
                	( beta * ms(is) * tau(is))&
             	   + (pperp**2.d0)/( tau(is) *beta * ms(is) * alph(is)) ) )

    		case (2) ! bi-kappa
    		    f0 = (1.d0+((ppar-p_drift(is))**2.d0/&
                	(beta * ms(is) * kappa(is) * a * a * tau(is)))&
             	   + (pperp**2.d0)/( tau(is) *beta * ms(is) * kappa(is) * a * a * alph(is)))**(-1.d0-kappa(is))

  			case (3) ! Juettner
  			   f0 = exp( -(2.d0*ms(is)/(vA*vA*beta*tau(is)*alph(is))) *&
  			   		sqrt(1.d0+pperp**2*vA*vA/(ms(is)*ms(is))&
  			   		     +(ppar-p_drift(is))**2*vA*vA*alph(is)/(ms(is)*ms(is))))

        case (4) ! bi-Moyal
          f0 = exp(0.5d0*((( (ppar-p_drift(is))**2.d0)/&
               ( beta * ms(is) * tau(is))&
              + (pperp**2.d0)/( tau(is) *beta * ms(is) * alph(is)) ) - &
               exp( (( (ppar-p_drift(is))**2.d0)/&
               ( beta * ms(is) * tau(is))&
              + (pperp**2.d0)/( tau(is) *beta * ms(is) * alph(is)) ) ) ))

    	  end select

    	    f0 = f0 * norm
           integrate = integrate + dpperp*dppar*2.d0*pi*pperp*f0

           if (distribution(is).NE.4) then
             write(unit_out(is),*)  pperp, ppar, f0
           endif

        enddo
     enddo

     ! Numerically re-normalise the Moyal distribution:
     if (distribution(is).EQ.4) then
       norm = 1.d0/integrate
       ifit_1=norm
       integrate=0.d0
       do iperp=0,nperp
          pperp = real(iperp)*dpperp
          do ipar=0,npar
             ppar = real(ipar)*dppar - ppar_max + p_drift(is)

             f0 = norm*exp(0.5d0*((( (ppar-p_drift(is))**2.d0)/&
                  ( beta * ms(is) * tau(is))&
                 + (pperp**2.d0)/( tau(is) *beta * ms(is) * alph(is)) ) - &
                  exp( (( (ppar-p_drift(is))**2.d0)/&
                  ( beta * ms(is) * tau(is))&
                 + (pperp**2.d0)/( tau(is) *beta * ms(is) * alph(is)) ) ) ))

            integrate = integrate + dpperp*dppar*2.d0*pi*pperp*f0

            write(unit_out(is),*)  pperp, ppar, f0

           enddo
        enddo
     endif





     write (*,'(a,i3)')   "Species", is
     if (distribution(is).EQ.0) &
       write (*,*) "Using analytical function from distribution_analyt.f90. Ignoring autoscaling."
     write (*,'(a,es14.4)') " Integration:   ", integrate
     write (*,'(a,es14.4)') " Norm:          ", norm
     write (*,'(a,es14.4)') " pperp_max:     ", pperp_max
     write (*,'(a,es14.4)') " ppar_max:      ", ppar_max
     write (*,*) " "
     select case(distribution(is))
      case (0)
                write (*,*) "Fit type:          0"
      case (1)
                write (*,*) "Fit type:          1"
                write (*,'(a,es14.4)') " ideal fit_1:   ", ifit_1
                write (*,'(a,es14.4)') " ideal fit_2:   ", ifit_2
                write (*,'(a,es14.4)') " ideal fit_3:   ", ifit_3
                write (*,'(a,es14.4)') " ideal perpcorr:", iperpcorr
      case (2)
                write (*,*) "Fit type:          2"
                write (*,'(a,es14.4)') " ideal fit_1:   ", ifit_1
                write (*,'(a,es14.4)') " ideal fit_2:   ", ifit_2
                write (*,'(a,es14.4)') " ideal fit_3:   ", ifit_3
                write (*,'(a,es14.4)') " ideal fit_4:   ", ifit_4
                write (*,'(a,es14.4)') " ideal fit_5:   ", ifit_5
                write (*,'(a,es14.4)') " ideal perpcorr:", iperpcorr
      case (3)
                write (*,*) "Fit types:         3, 4, or 5"
                write (*,*) "If fit type 3 (in pperp and ppar):"
                write (*,'(a,es14.4)') " ideal fit_1:   ", ifit_1
                write (*,'(a,es14.4)') " ideal fit_2:   ", ifit_2
                write (*,'(a,es14.4)') " ideal fit_3:   ", ifit_3
                write (*,*) " "
                write (*,*) "If fit type 4 (gamma-only):"
                write (*,'(a,es14.4)') " ideal fit_1:   ", ifit_1
                write (*,'(a,es14.4)') " ideal perpcorr:", iperpcorr
                write (*,*) " "
                write (*,*) "If fit type 5 (gamma and pparbar):"
                write (*,'(a,es14.4)') " ideal fit_1:   ", ifit_1
                write (*,'(a,es14.4)') " ideal fit_2:   ", ifit_5
                write (*,'(a,es14.4)') " ideal fit_3:   ", ifit_3
                write (*,'(a,es14.4)') " ideal perpcorr:", iperpcorr
      case (4)
                write (*,*) "Fit type:          6"
                write (*,'(a,es14.4)') " ideal fit_1:   ", ifit_1
                write (*,'(a,es14.4)') " ideal fit_2:   ", ifit_2
                write (*,'(a,es14.4)') " ideal fit_3:   ", ifit_3
                write (*,'(a,es14.4)') " ideal fit_4:   ", ifit_4
                write (*,'(a,es14.4)') " ideal perpcorr:", iperpcorr
                write (*,'(a,es14.4)') " WARNING: The Moyal distribution is always normalised to 1."
                write (*,'(a,es14.4)') "          Do not use Integration as a measure for P-range."
     end select
     write (*,'(a,es14.4)') "============================"

  enddo

  do is = 1, nspec
     close(unit_out(is))
  enddo

contains



  subroutine read_in_params
!! This subroutine reads in the parameters for the generation of distribution functions for ALPS.
    implicit none


    nameList /system/ &
         nspec, beta, vA, nperp, npar, maxP, writeName

    call get_unused_unit (input_unit_no)
    call get_runname(runname)
    runname=trim(runname)//".in"
    unit=input_unit_no
    open (unit=unit,file=runname,status='old',action='read')
    read (unit=unit,nml=system)

    allocate(ms(1:nspec));   ms = 0.d0
    allocate(tau(1:nspec));  tau = 0.d0
    allocate(alph(1:nspec)); alph = 0.d0
    allocate(p_drift(1:nspec));   p_drift = 0.d0
    allocate(distribution(1:nspec));   distribution = 0
    allocate(kappa(1:nspec));   kappa = 0.d0
    allocate(autoscale(1:nspec));  autoscale=.TRUE.
    allocate(maxPperp(1:nspec)); maxPperp = 0.d0
    allocate(maxPpar(1:nspec)); maxPpar = 0.d0


    !Read in species parameters:
    do is = 1,nspec
       unit=input_unit_no
       call get_indexed_namelist_unit (unit, "spec", is)
       call spec_read(is)
       close (unit)
    enddo
    unit=input_unit_no
    close (unit)

  end subroutine read_in_params



subroutine spec_read(is)
!!Subroutine for reading in species parameters.
  implicit none

  integer, intent(in) :: is
  !!Species index.

  double precision :: tauS
  !! Parallel temperature ratio of species.

  double precision :: mS_read
  !! Mass of species.

  double precision :: alphS
  !! Temperature anisotropy of species.

  double precision :: pS
  !! Drift momentum of species.

  double precision :: kappaS
  !! Kappa index of species.

  double precision :: maxPperpS
  !! Maximum perpendicular momentum of species.

  double precision :: maxPparS
  !! Maximum parallel momentum of species.

  integer :: distributionS
  !! Type of distribution for species.

  logical :: autoscaleS
  !! Flag for autoscaling of grid for species.

  nameList /spec/ &
       mS_read, tauS, alphS, pS, kappaS, distributionS, autoscaleS, maxPperpS, maxPparS

  read (unit=unit,nml=spec)
  tau(is)=tauS
  mS(is)=mS_read
  alph(is)=alphS
  p_drift(is)=pS
  kappa(is)=kappaS
  distribution(is)=distributionS
  autoscale(is)=autoscaleS
  maxPperp(is)=maxPperpS
  maxPpar(is)=maxPparS

end subroutine spec_read




  subroutine get_indexed_namelist_unit (unit, nml, index_in)
!! Determines unused I/O unit.
    implicit none

    integer, intent (out) :: unit
    !!Unit to be defined.

    character (*), intent (in) :: nml
    !!Character string for namelist to be read in.

    integer, intent (in) :: index_in
    !!Index of namelist to be read in.

    character(500) :: line
    !!I/O dummy variable.

    integer :: iunit
    !!I/O dummy index.

    integer :: iostat
    !!I/O dummy index.

    integer :: in_file
    !!I/O dummy index.

    integer :: ind_slash
    !!I/O dummy index.

    logical :: exist
    !!Check if namelist is open.


    call get_unused_unit (unit)
    ind_slash=index(runname,"/",.True.)
    if (ind_slash.EQ.0) then !No slash in name
        !Original behaviour
        open (unit=unit, file="."//trim(runname)//".scratch")
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
       write(6,*) "get_indexed_namelist: following namelist not found ",trim(line)
       return
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




  function input_unit_exist (nml,exist)
    !! Determine if a particular namelist already opened.
    implicit none

    character(*), intent (in) :: nml
    !!Namelist to be opened.

    logical, intent(out) :: exist
    !!Determination if namelist is open.

    integer :: input_unit_exist
    !!I/O dummy index.

    integer :: iostat
    !!I/O dummy index.

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

    integer :: input_unit
    !!I/O dummy index.

    integer :: iostat
    !!I/O dummy index.

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
  !! Determine unused number for I/O index.
    implicit none

    integer, intent (out) :: unit
    !!Unit to be assigned.

    logical :: od
    !! Check whether unit is opened.


    unit = 50
    do
       inquire (unit=unit, opened=od)
       if (.not.od) return
       unit = unit + 1
    end do
  end subroutine get_unused_unit




  subroutine get_runname(runname)
  !! Get runname for output files from input argument.
    implicit none

    integer :: l
    !! Dummy Length.

    character(500) :: arg
    !! Input Argument.

    character(500), intent(out) :: runname
    !! Basename for file I/O.

    !Get the first argument of the program execution command
    call getarg(1,arg)

    !Check if this is the input file and trim .in extension to get runname
    l = len_trim (arg)
    if (l > 3 .and. arg(l-2:l) == ".in") then
       runname = arg(1:l-3)
    end if
  end subroutine get_runname


end program generate_distribution





      FUNCTION BESSK(N,X)
        !! This function calculates the modified Bessel function of the third kind
        !! of integer order N, for any REAL X. The classical recursion formula is used. Reference:
        !! C.W.Clenshaw, Chebyshev Series for Mathematical Functions, Mathematical Tables, Vol. 5, 1962.
      IMPLICIT NONE

      integer, intent(in) :: N
      !! Order of Bessel function.

      double precision, intent(in) :: X
      !! Argument of the Bessel function.

      INTEGER J
      REAL *8 BESSK,BESSK0,BESSK1,TOX,BK,BKM,BKP

      IF (N.EQ.0) THEN
      BESSK = BESSK0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSK = BESSK1(X)
      RETURN
      ENDIF
      IF (X.EQ.0.D0) THEN
      BESSK = 1.D30
      RETURN
      ENDIF
      TOX = 2.D0/X
      BK  = BESSK1(X)
      BKM = BESSK0(X)
      DO 11 J=1,N-1
      BKP = BKM+DFLOAT(J)*TOX*BK
      BKM = BK
      BK  = BKP
   11 CONTINUE
      BESSK = BK
      RETURN
      END



      FUNCTION BESSK0(X)
!! This function calculates the modified Bessel function of the third kind of order zero for any positive real argument x.
      IMPLICIT NONE

      double precision, intent(in) :: X
      !! Argument of the Bessel function.

      REAL*8 BESSK0,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,    &
      BESSI0
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0, &
      0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1, &
      -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      IF(X.EQ.0.D0) THEN
      BESSK0=1.D30
      RETURN
      ENDIF
      IF(X.LE.2.D0) THEN
      Y=X*X/4.D0
      AX=-LOG(X/2.D0)*BESSI0(X)
      BESSK0=AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      Y=(2.D0/X)
      AX=EXP(-X)/DSQRT(X)
      BESSK0=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END





      FUNCTION BESSK1(X)
!! This function calculates the modified Bessel function of the third kind of order zero for any positive real argument x.
      IMPLICIT NONE

      double precision, intent(in) :: X
      !! Argument of the Bessel function.

      REAL*8 BESSK1,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,BESSI1
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,0.15443144D0,-0.67278579D0,  &
      -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1, &
      0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
      IF(X.EQ.0.D0) THEN
      BESSK1=1.D32
      RETURN
      ENDIF
      IF(X.LE.2.D0) THEN
      Y=X*X/4.D0
      AX=LOG(X/2.D0)*BESSI1(X)
      BESSK1=AX+(1.D0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      Y=(2.D0/X)
      AX=EXP(-X)/DSQRT(X)
      BESSK1=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END





      FUNCTION BESSI0(X)
        !! This function calculates the modified Bessel function of the first kind of order zero for any positive real argument x.
      IMPLICIT NONE

      double precision, intent(in) :: X
      !! Argument of the Bessel function.

      REAL *8 BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067429D0,  &
      0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,  &
      0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
      0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/DSQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI0=AX*BX
      ENDIF
      RETURN
      END





      FUNCTION BESSI1(X)
!! This function calculates the modified Bessel function of the first kind of order one for any positive real argument x.
      IMPLICIT NONE

      double precision, intent(in) :: X
      !! Argument of the Bessel function.

      REAL *8 BESSI1,Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
      0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
      -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,        &
      -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/DSQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI1=AX*BX
      ENDIF
      RETURN
      END
