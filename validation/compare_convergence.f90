! Report empirical integration sensitivity; this is not a pass/fail accuracy oracle.
program compare_convergence
  use iso_fortran_env, only : dp=>real64
  implicit none
  complex(dp) :: chi(3,3,2,3,2),det(3,2)
  real(dp) :: re,im,unused1,unused2,err,norm
  character(1024) :: path,line
  character(128) :: label
  integer :: file,ios,run,s,i,j,w
  logical :: seen(3,3,2,3,2),seen_det(3,2)
  chi=0; det=0; seen=.false.; seen_det=.false.
  do run=1,2
     call get_command_argument(run,path)
     open(newunit=file,file=trim(path),status='old',action='read')
     do
        read(file,'(a)',iostat=ios) line
        if(ios<0) exit
        if(ios/=0) error stop 'reading convergence data'
        read(line,*) label,re,im,unused1,unused2
        if(index(label,'chi')==1) then
           read(label(4:4),*) i
           read(label(5:5),*) j
           read(label(8:8),*) s
           read(label(11:),*) w
           chi(i,j,s,w,run)=cmplx(re,im,dp)
           chi(j,i,s,w,run)=chi(i,j,s,w,run)
           seen(i,j,s,w,run)=.true.; seen(j,i,s,w,run)=.true.
        else if(index(label,'determinant_s0_i')==1) then
           read(label(17:),*) w
           det(w,run)=cmplx(re,im,dp); seen_det(w,run)=.true.
        endif
     enddo
     close(file)
  enddo
  if(.not.all(seen).or..not.all(seen_det)) error stop 'incomplete convergence observations'
  print '(a)', '# omega_index tensor_scaled_Frobenius determinant_absolute_delta determinant_relative_delta'
  do w=1,3
     err=sqrt(sum(abs(chi(:,:,:,w,2)-chi(:,:,:,w,1))**2))
     norm=max(1._dp,sqrt(sum(abs(chi(:,:,:,w,2))**2)))
     print '(i3,3es24.15e3)', w,err/norm,abs(det(w,2)-det(w,1)), &
          abs(det(w,2)-det(w,1))/max(abs(det(w,2)),tiny(1._dp))
  enddo
end program
