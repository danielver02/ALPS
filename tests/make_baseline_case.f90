! Derive small, isolated runs from the established examples; never edit templates.
program make_baseline_case
  use relativistic_references, only : dp
  implicit none
  character(1024) :: repo,path,line,arg
  character(32) :: mode
  integer :: input,out,ios,ng,nu,nres,is,i,j
  real(dp) :: p(0:30,0:60,2),f(0:30,0:60),norm,step1,step2
  call get_command_argument(1,repo)
  call get_command_argument(2,mode)
  ng=100; nu=200; nres=100
  call get_command_argument(3,arg)
  if(len_trim(arg)>0) read(arg,*) ng
  call get_command_argument(4,arg)
  if(len_trim(arg)>0) read(arg,*) nu
  call get_command_argument(5,arg)
  if(len_trim(arg)>0) read(arg,*) nres
  if(mode(1:3)=='rel') then
     path=trim(repo)//'/tests/test_relativistic.in'
  else
     path=trim(repo)//'/tests/test_chebyshev.in'
  endif
  open(newunit=input,file=trim(path),status='old',action='read')
  open(newunit=out,file='case.in',status='new')
  do
     read(input,'(a)',iostat=ios) line
     if(ios<0) exit
     if(ios/=0) error stop 'reading baseline template'
     line=adjustl(line)
     if(index(line,'n_scan=')==1) line='n_scan=0'
     if(index(line,'profile_startup=')==1) line='profile_startup=F'
     if(index(line,'writeOut=')==1) line='writeOut=F'
     if(index(line,'rel_cache=')==1) line="rel_cache='off'"
     if(index(line,'ngamma=')==1) write(line,'(a,i0)') 'ngamma=',ng
     if(index(line,'npparbar=')==1) write(line,'(a,i0)') 'npparbar=',nu
     if(index(line,'n_resonance_interval=')==1) write(line,'(a,i0)') 'n_resonance_interval=',nres
     if(index(line,'numiter')==1) line='numiter=150'
     if(index(line,'D_tol=')==1) line='D_tol=1.d-11'
     if(trim(mode)=='rel5') then
        if(index(line,'fit_type_in=')==1) line='fit_type_in=5'
        if(index(line,'fit_2=')==1) line='fit_2=0.15d0'
        if(index(line,'fit_3=')==1) line='fit_3=0.d0'
     endif
     write(out,'(a)') trim(line)
  enddo
  close(input); close(out)
  if(mode(1:2)=='nr') then
     ! Insert the derivative selector inside each existing poly namelist.
     open(newunit=input,file='case.in',status='old',action='read')
     open(newunit=out,file='derivatives.tmp',status='new')
     do
        read(input,'(a)',iostat=ios) line
        if(ios<0) exit
        if(ios/=0) error stop 'reading derivative fixture'
        write(out,'(a)') trim(line)
        if(index(line,'&poly_spec_')==1) then
           if(trim(mode)=='nr_analytic') then
              write(out,'(a)') 'derivative_type=1'
           else
              write(out,'(a)') 'derivative_type=0'
           endif
        endif
     enddo
     close(input); close(out)
     open(newunit=input,file='derivatives.tmp',status='old',action='read')
     open(newunit=out,file='case.in',status='replace')
     do
        read(input,'(a)',iostat=ios) line
        if(ios<0) exit
        if(ios/=0) error stop 'copying derivative fixture'
        write(out,'(a)') trim(line)
     enddo
     close(input,status='delete'); close(out)
  endif
  if(trim(mode)=='rel5') then
     ! Test-only F=C exp(-2 Gamma-0.15 u^2), with source-grid normalization.
     do is=1,2
        write(path,'(a,i0,a)') 'distribution/test_relativistic.',is,'.array'
        open(newunit=input,file=trim(path),status='old',action='read')
        do i=0,30
           do j=0,60
              read(input,*) p(i,j,1),p(i,j,2),f(i,j)
              f(i,j)=exp(-2*sqrt(1+p(i,j,1)**2+p(i,j,2)**2)-0.15_dp*p(i,j,2)**2)
           enddo
        enddo
        close(input)
        step1=p(1,0,1)-p(0,0,1); step2=p(0,1,2)-p(0,0,2)
        norm=sum(p(:,:,1)*f)*2*acos(-1._dp)*step1*step2
        open(newunit=out,file=trim(path),status='replace')
        do i=0,30
           do j=0,60
              write(out,'(3es26.17e3)') p(i,j,:),f(i,j)/norm
           enddo
        enddo
        close(out)
     enddo
  endif
end program
