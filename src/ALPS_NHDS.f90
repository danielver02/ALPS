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

module alps_nhds
!! Module including the NHDS implementation for bi-Maxwellian reference cases.
!! The original NHDS code can be found under github.com/danielver02/NHDS.
  implicit none

  private :: calc_ypsilon, besselI, BESSI, BESSI0, BESSI1, WOFZ, dispfunct
  public :: calc_chi

contains

  ! This file is part of NHDS
  ! Copyright (C) 2020 Daniel Verscharen (d.verscharen@ucl.ac.uk)
  !All rights reserved.
  !
  !Redistribution and use in source and binary forms, with or without
  !modification, are permitted provided that the following conditions are met:
  !
  !1. Redistributions of source code must retain the above copyright notice, this
  !   list of conditions and the following disclaimer.
  !2. Redistributions in binary form must reproduce the above copyright notice,
  !   this list of conditions and the following disclaimer in the documentation
  !   and/or other materials provided with the distribution.
  !
  !THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  !ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  !WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  !DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
  !ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  !(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  !LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  !ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  !(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  !SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  !
  !The views and conclusions contained in the software and documentation are those
  !of the authors and should not be interpreted as representing official policies,
  !either expressed or implied, of the NHDS project.



  subroutine calc_chi(chi,j,kz,kperp,x)
  !! Subroutine that calculates the susceptibility of species j based on NHDS.
  use alps_var, only : bMnmaxs, bMBessel_zeros, bMbetas, bMalphas,bMpdrifts
  use alps_var, only : ms, qs, ns
  implicit none



  double complex, intent(out) :: chi(3,3)
  !! Susceptibility tensor of species j.

  integer, intent(in) :: j
  !! Index for species.

  double precision, intent(in) :: kz
  !! Normalised parallel wavenumber.

  double precision, intent(in) :: kperp
  !! Normalised perpendicular wavenumber.

  double complex, intent(in) :: x
  !! Normalised complex frequency.

  double complex :: Y(3,3)
  !! Y-tensor according to Stix.

  double complex :: Ynew(3,3)
  !! Iteration of Y-tensor according to Stix.

  double precision :: z
  !! Argument of dispersion function.

  double precision :: Omega
  !! Normalised gyro-frequency.

  double precision :: ell
  !! Normalised inertial length.

  double precision :: vtherm
  !! Normalised thermal speed.

  double precision :: Bessel_zero
  !! Limit of Bessel-function calculation.

  double precision :: vdrift
  !! Normalised drift speed of bi-Maxwellian.

  integer :: n
  !! Index of sum over Bessel functions.

  integer :: i
  !! Index running over tensor.

  integer :: k
  !! Index running over tensor.

  integer :: nmaxrun
  !! Running to maximum index of n for Bessel functions.

  integer :: nmax
  !! Maximum index of n for Bessel functions.

  logical :: Bessel_run
  !! Check whether maximum n has been achieved.


  Omega=qs(j)/ms(j)
  ell=sqrt(ms(j)/(ns(j)*qs(j)*qs(j)))
  vtherm=sqrt(bMbetas(j)/(ns(j)*ms(j)))
  nmax=bMnmaxs(j)
  Bessel_zero=bMBessel_zeros(j)
  vdrift=bMpdrifts(j)/ms(j)

  z=0.5d0*(kperp*vtherm/Omega)*(kperp*vtherm/Omega)*bMalphas(j)

  do i=1,3
   do k=1,3
    Y(i,k)=0.d0
   enddo
  enddo


  !determine maximum n for Besselfunction:
  nmaxrun=nmax
  n=0
  Bessel_run=.TRUE.
  do while (Bessel_run)
    if ((n.GE.nmax).OR.(besselI(n,z).LT.Bessel_zero)) then
       nmaxrun=n
       Bessel_run=.FALSE.
    endif
    n=n+1
  enddo


  do n=-nmaxrun,nmaxrun
     call calc_ypsilon(Ynew,j,n,kz,kperp,x)
     do i=1,3
      do k=1,3
        Y(i,k)=Y(i,k)+exp(-z)*Ynew(i,k)
      enddo
     enddo
  enddo

  chi(1,1)=Y(1,1)/(ell*ell)
  chi(1,2)=Y(1,2)/(ell*ell)
  chi(1,3)=Y(1,3)/(ell*ell)
  chi(2,1)=Y(2,1)/(ell*ell)
  chi(2,2)=Y(2,2)/(ell*ell)
  chi(2,3)=Y(2,3)/(ell*ell)
  chi(3,1)=Y(3,1)/(ell*ell)
  chi(3,2)=Y(3,2)/(ell*ell)
  chi(3,3)=2.d0*x*vdrift/(ell*ell*kz*vtherm*vtherm*bMalphas(j))+Y(3,3)/(ell*ell)



  end subroutine







  subroutine calc_ypsilon(Y,j,n,kz,kperp,x)
  !! Calculates the Y-tensor according to Stix for a bi-Maxwelling, using the NHDS calculation.
  use alps_var, only : bMbetas, bMalphas,bMpdrifts
  use alps_var, only : ms, qs, ns
  implicit none

  double complex, parameter ::  uniti=(0.d0,1.d0)
  !! Imaginary unit.

  double complex, intent(out) :: Y(3,3)

  integer, intent(in) :: j
  !! Index for species.

  integer, intent(in) :: n
  !! Index of sum over Bessel functions.

  double precision, intent(in) :: kz
  !! Normalised parallel wavenumber.

  double precision, intent(in) :: kperp
  !! Normalised perpendicular wavenumber.

  double complex, intent(in) :: x
  !! Normalised complex frequency.

  double complex :: zeta
  !! Argument of dispersion function.

  double complex :: An
  !! An function according to Stix.

  double complex :: Bn
  !! Bn function according to Stix.

  double complex :: resfac
  !! Resonance factor.

  double precision :: BInz
  !! Modified Bessel function, evaluated.

  double precision :: z
  !! Argument of Bessel function.

  double precision :: Omega
  !! Normalised gyro-frequency.

  double precision :: ell
  !! Normalised inertial length.

  double precision :: vtherm
  !! Normalised thermal speed.


  double precision :: vdrift
  !! Normalised drift speed of bi-Maxwellian.

  double precision :: dBInzdz
  !! Derivative of Bessel function.

  logical :: kpos
  !! Check whether kpar is positive.

  kpos=.TRUE.
  if (kz.LT.0.d0) kpos=.FALSE.

  Omega=qs(j)/ms(j)
  ell=sqrt(ms(j)/(ns(j)*qs(j)*qs(j)))
  vtherm=sqrt(bMbetas(j)/(ns(j)*ms(j)))
  vdrift=bMpdrifts(j)/ms(j)


  zeta=(x-kz*vdrift-1.d0*n*Omega)/(kz*vtherm)
  resfac=x-kz*vdrift-1.d0*n*Omega
  z=0.5d0*(kperp*vtherm/Omega)*(kperp*vtherm/Omega)*bMalphas(j)

  An=(bMalphas(j)-1.d0)
  An=An+(1.d0/(kz*vtherm)) *( bMalphas(j)*resfac + 1.d0*n*Omega)*dispfunct(zeta,kpos)

  Bn=(bMalphas(j)*(x-1.d0*n*Omega)-(kz*vdrift-1.d0*n*Omega))/kz
  Bn=Bn+((x-1.d0*n*Omega)*(bMalphas(j)*resfac+1.d0*n*Omega)/(kz*kz*vtherm) )*dispfunct(zeta,kpos)

  if (n.GE.0) then
  	BInz=1.d0*besselI(n,z)
  	dBInzdz=besselI(n+1,z)+1.d0*n*BInz/z
  else
  	BInz=1.d0*besselI(-n,z)
  	dBInzdz=besselI(-n-1,z)+1.d0*n*BInz/z
  endif


  ! The tensor in Stix's (10-57)
  Y(1,1)=1.d0*(n*n)*BInz*An/z
  Y(1,2)=-uniti*n*(BInz-dBInzdz)*An
  Y(1,3)=kperp*n*BInz*Bn/(Omega*z)
  Y(2,1)=uniti*n*(BInz-dBInzdz)*An
  Y(2,2)=(1.d0*(n*n)*BInz/z+2.d0*z*BInz-2.d0*z*dBInzdz)*An
  Y(2,3)=uniti*kperp*(BInz-dBInzdz)*Bn/Omega
  Y(3,1)=kperp*BInz*n*Bn/(Omega*z)
  Y(3,2)=-uniti*kperp*(BInz-dBInzdz)*Bn/Omega
  Y(3,3)=2.d0*(x-1.d0*n*Omega)*BInz*Bn/(kz*vtherm*vtherm*bMalphas(j))


  end subroutine





double precision function besselI(n,x)
  !! Calculates the modified Bessel function of argument x and order n.
  implicit none

  double precision, intent(in) :: x
  !! Argument of the modified Bessel function.

  integer, intent(in) :: n
  !! Order of the modified Bessel function.

  if (n.LT.0) then
  	besselI=BESSI(-n,x)
  else
  	besselI=BESSI(n,x)
  endif

  end function


!
! Calculate the dispersion function:
!
double complex function dispfunct(zeta,kpos)
  !! Calculates the dispersion function based on the complex error function.
  implicit none

  double complex, parameter ::  uniti=(0.d0,1.d0)
  !! Imaginary init.

  double complex, intent(in) :: zeta
  !! Argument of dispersion function.

  logical, intent(in) :: kpos
  !! Check whether kpar is positive.

  logical :: flag
  !! Check whether [[WOFZ(subroutine)]] executes successfully.

  double precision :: U
  !! Real part of the output for [[WOFZ(subroutine)]].

  double precision :: V
  !! Imaginary part of the output for [[WOFZ(subroutine)]].

  double precision :: XI
  !! Real part of the argument for [[WOFZ(subroutine)]].

  double precision :: YI
  !! Imaginary part of the argument for [[WOFZ(subroutine)]].

  double precision :: M_PI = 4.d0*atan(1.d0)
  !! Pi.

  XI=1.d0*real(zeta)
  YI=1.d0*aimag(zeta)

  if (kpos) then
     call WOFZ(XI,YI,U,V,flag)
     dispfunct = uniti*sqrt(M_PI)*(U+uniti*V)
  else
     call WOFZ(-XI,-YI,U,V,flag)
     dispfunct = -uniti*sqrt(M_PI)*(U+uniti*V)
  endif
  end function







  SUBROUTINE WOFZ (XI, YI, U, V, FLAG)
    !!  Given a complex number Z = (XI,YI), this subroutine computes
    !!  the value of the Faddeeva function W(Z) = exp(-Z**2)*erfc(-I*Z),
    !!  where erfc is the complex complementary error function and I
    !! is the imaginary unit.

    !Based on
        !G.P.M. Poppe, C.M.J. Wijers, "More Efficient Computation of the Complex Error-Function",
        !ACM Trans. Math. Software 16, 47 (1990).

        !      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
        !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
        !      VOL. 16, NO. 1, PP. 47


  !  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
  !  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
  !  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
  !  OF THE FUNCTION.
  !  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.

  !  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
  !     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
  !                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
  !                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
  !                FLOATING-POINT ARITHMETIC
  !     RMAXEXP  = LN(RMAX) - LN(2)
  !     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
  !                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
  !  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
  !  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
  !
  !
  !  PARAMETER LIST
  !     XI     = REAL      PART OF Z
  !     YI     = IMAGINARY PART OF Z
  !     U      = REAL      PART OF W(Z)
  !     V      = IMAGINARY PART OF W(Z)
  !     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
  !              OCCUR OR NOT; TYPE LOGICAL;
  !              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
  !              MEANING :
  !              FLAG=.FALSE. : NO ERROR CONDITION
  !              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
  !                             BECOMES INACTIVE
  !  XI, YI      ARE THE INPUT-PARAMETERS
  !  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
  !
  !  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
  !
  !  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
  !  PUT TO 0 UPON UNDERFLOW;
  !

        IMPLICIT DOUBLE PRECISION (A-H, O-Z)
        IMPLICIT INTEGER (I-N)

        LOGICAL A, B, FLAG
  PARAMETER (FACTOR=1.12837916709551257388D0,RMAXREAL = 0.5D+154,RMAXEXP  = 708.503061461606D0,&
  &RMAXGONI = 3.53711887601422D+15)

        FLAG = .FALSE.

        XABS = DABS(XI)
        YABS = DABS(YI)
        X    = XABS/6.3
        Y    = YABS/4.4

  !     THE FOLLOWING IF-STATEMENT PROTECTS
  !     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
  !
        IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) GOTO 100

        QRHO = X**2 + Y**2

        XABSQ = XABS**2
        XQUAD = XABSQ - YABS**2
        YQUAD = 2*XABS*YABS

       A     = QRHO.LT.0.085264D0

        IF (A) THEN

  !  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
  !  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
  !  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
  !  ACCURACY

          QRHO  = (1-0.85*Y)*DSQRT(QRHO)
          N     = IDNINT(6 + 72*QRHO)
          J     = 2*N+1
          XSUM  = 1.0/J
          YSUM  = 0.0D0
          DO 10 I=N, 1, -1
            J    = J - 2
            XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I
            YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I
            XSUM = XAUX + 1.0/J
   10     CONTINUE
          U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0
          V1   =  FACTOR*(XSUM*XABS - YSUM*YABS)
          DAUX =  DEXP(-XQUAD)
          U2   =  DAUX*DCOS(YQUAD)
          V2   = -DAUX*DSIN(YQUAD)

          U    = U1*U2 - V1*V2
          V    = U1*V2 + V1*U2

        ELSE

  !  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
  !  CONTINUED FRACTION
  !  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
  !  ACCURACY
  !
  !  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
  !  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
  !  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
  !  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
  !  TO OBTAIN THE REQUIRED ACCURACY
  !  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
  !  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY

          IF (QRHO.GT.1.0) THEN
            H    = 0.0D0
            KAPN = 0
            QRHO = DSQRT(QRHO)
            NU   = IDINT(3 + (1442/(26*QRHO+77)))
          ELSE
            QRHO = (1-Y)*DSQRT(1-QRHO)
            H    = 1.88*QRHO
            H2   = 2*H
            KAPN = IDNINT(7  + 34*QRHO)
            NU   = IDNINT(16 + 26*QRHO)
          ENDIF

          B = (H.GT.0.0)

          IF (B) QLAMBDA = H2**KAPN

          RX = 0.0
          RY = 0.0
          SX = 0.0
          SY = 0.0

          DO 11 N=NU, 0, -1
            NP1 = N + 1
            TX  = YABS + H + NP1*RX
            TY  = XABS - NP1*RY
            C   = 0.5/(TX**2 + TY**2)
            RX  = C*TX
            RY  = C*TY
            IF ((B).AND.(N.LE.KAPN)) THEN
              TX = QLAMBDA + SX
              SX = RX*TX - RY*SY
              SY = RY*TX + RX*SY
              QLAMBDA = QLAMBDA/H2
            ENDIF
   11     CONTINUE

          IF (H.EQ.0.0) THEN
            U = FACTOR*RX
            V = FACTOR*RY
          ELSE
            U = FACTOR*SX
            V = FACTOR*SY
          END IF

          IF (YABS.EQ.0.0) U = DEXP(-XABS**2)

        END IF

  !  EVALUATION OF W(Z) IN THE OTHER QUADRANTS

        IF (YI.LT.0.0) THEN
          IF (A) THEN
            U2    = 2*U2
            V2    = 2*V2
          ELSE
            XQUAD =  -XQUAD

  !         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
  !        AGAINST OVERFLOW

            IF ((YQUAD.GT.RMAXGONI).OR. (XQUAD.GT.RMAXEXP)) GOTO 100

            W1 =  2*DEXP(XQUAD)
            U2  =  W1*DCOS(YQUAD)
            V2  = -W1*DSIN(YQUAD)
          END IF

          U = U2 - U
          V = V2 - V
          IF (XI.GT.0.0) V = -V
        ELSE
          IF (XI.LT.0.0) V = -V
        END IF

        RETURN

    100 FLAG = .TRUE.
        RETURN

        END





  ! ----------------------------------------------------------------------
        double precision FUNCTION BESSI(N,X)
        !! Function to calculate the first kind modified Bessel function of integer order N
        !! for any real X.
        ! -------------------------------------------------------------------- *
        !   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.   *
        !                                                                      *
        !                               F90 Release 1.1 By J-P Moreau, Paris.  *
        !                                                                      *
        !   Version 1.1: corected value of P4 in BESSIO (P4=1.2067492 and not  *
        !                1.2067429) Aug. 2011.                                 *
  
        IMPLICIT DOUBLE PRECISION (A-H, O-Z)
        IMPLICIT INTEGER (I-N)
  !     This subroutine calculates the first kind modified Bessel function
  !     of integer order N, for any REAL X. We use here the classical
  !     recursion formula, when X > N. For X < N, the Miller's algorithm
  !     is used to avoid overflows.
  !     REFERENCE:
  !     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
  !     MATHEMATICAL TABLES, VOL.5, 1962.
  !
        PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
        REAL *8 X,TOX,BIM,BI,BIP
        IF (N.EQ.0) THEN
        BESSI = BESSI0(X)
        RETURN
        ENDIF
        IF (N.EQ.1) THEN
        BESSI = BESSI1(X)
        RETURN
        ENDIF
        IF(X.EQ.0.D0) THEN
        BESSI=0.D0
        RETURN
        ENDIF
        TOX = 2.D0/X
        BIP = 0.D0
        BI  = 1.D0
        BESSI = 0.D0
        M = 2*((N+INT(SQRT(FLOAT(IACC*N)))))
        DO 12 J = M,1,-1
        BIM = BIP+DFLOAT(J)*TOX*BI
        BIP = BI
        BI  = BIM
        IF (ABS(BI).GT.BIGNO) THEN
        BI  = BI*BIGNI
        BIP = BIP*BIGNI
        BESSI = BESSI*BIGNI
        ENDIF
        IF (J.EQ.N) BESSI = BIP
     12 CONTINUE
        BESSI = BESSI*BESSI0(X)/BI
        RETURN
        END




        double precision FUNCTION BESSI0(X)
        !! Auxiliary Bessel functions for N=0, N=1
        REAL *8 X,Y,P1,P2,P3,P4,P5,P6,P7,  &
        Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
        DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067492D0,  &
        0.2659732D0,0.360768D-1,0.45813D-2/
        DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1, &
        0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
        0.2635537D-1,-0.1647633D-1,0.392377D-2/
        IF(ABS(X).LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
        ELSE
        AX=ABS(X)
        Y=3.75D0/AX
        BX=EXP(AX)/SQRT(AX)
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
        BESSI0=AX*BX
        ENDIF
        RETURN
        END



        double precision FUNCTION BESSI1(X)
        !! Modified Bessel function of order 1.

        REAL *8 X,Y,P1,P2,P3,P4,P5,P6,P7,  &
        Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
        DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
        0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
        DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
        -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1, &
        -0.2895312D-1,0.1787654D-1,-0.420059D-2/
        IF(ABS(X).LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
        ELSE
        AX=ABS(X)
        Y=3.75D0/AX
        BX=EXP(AX)/SQRT(AX)
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
        BESSI1=AX*BX
        ENDIF
        RETURN
        END





end module alps_nhds
