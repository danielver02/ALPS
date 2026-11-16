! Strict record/label matching; tolerances come only from the reviewed reference.
program compare_baselines
  use iso_fortran_env, only : dp=>real64
  use ieee_arithmetic, only : ieee_is_finite
  implicit none
  character(1024) :: refpath,actualpath,line,other
  character(128) :: label,actual_label
  real(dp) :: re,im,ar,ai,atol,rtol,unused1,unused2,err,limit,maxratio
  integer :: ref,actual,ios,other_ios,n
  call get_command_argument(1,refpath)
  call get_command_argument(2,actualpath)
  open(newunit=ref,file=trim(refpath),status='old',action='read')
  open(newunit=actual,file=trim(actualpath),status='old',action='read')
  n=0; maxratio=0
  do
     read(ref,'(a)',iostat=ios) line
     read(actual,'(a)',iostat=other_ios) other
     if(ios<0.and.other_ios<0) exit
     if(ios/=0.or.other_ios/=0) error stop 'baseline record count mismatch'
     read(line,*,iostat=ios) label,re,im,atol,rtol
     read(other,*,iostat=other_ios) actual_label,ar,ai,unused1,unused2
     if(ios/=0.or.other_ios/=0) error stop 'malformed baseline record'
     if(label/=actual_label) error stop 'baseline label mismatch'
     if(.not.all(ieee_is_finite([re,im,ar,ai,atol,rtol]))) error stop 'nonfinite baseline record'
     if(atol<=0.or.rtol<0) error stop 'invalid reference tolerance'
     err=abs(cmplx(ar-re,ai-im,dp))
     limit=atol+rtol*abs(cmplx(re,im,dp))
     maxratio=max(maxratio,err/limit)
     if(err>limit) then
        print '(a,3es24.15e3)', 'FAIL '//trim(label)//': actual,error,limit ',ar,err,limit
        error stop 1
     endif
     n=n+1
  enddo
  if(n==0) error stop 'empty baseline'
  print '(a,i0,a,es12.4)', 'PASS: ',n,' numerical records; maximum tolerance fraction ',maxratio
end program
