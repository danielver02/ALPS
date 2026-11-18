! Derive valid AC_method=1 input from the existing fixture; tests select method 2
! internally only after the normal reader has completed. No input guard bypass.
program make_rel_susceptibility_case
  implicit none
  character(1024) :: repo,arg,line,path
  character(16) :: kind
  integer :: ng,nu,nr,nz,unit,out,ios,is,i,j,species
  double precision :: va,mass(2),rmax,umax,b,r,u,g,norm,dr,du
  double precision, allocatable :: f(:,:),pperp(:,:),ppar(:,:)
  call get_command_argument(1,repo); call get_command_argument(2,kind)
  call get_command_argument(3,arg); read(arg,*) ng
  call get_command_argument(4,arg); read(arg,*) nu
  nr=30; nz=60; rmax=5.d0; umax=5.d0; va=1.d0; mass=1.d0; b=.15d0
  if(kind=='iso') b=0.d0
  if(kind=='coarse') then
     nr=8; nz=16; rmax=1.2d0; umax=1.d0; va=.6d0; mass=[1.2d0,2.d0]
  endif
  open(newunit=unit,file=trim(repo)//'/tests/test_relativistic.in',status='old')
  open(newunit=out,file='case.in',status='new')
  species=0
  do
     read(unit,'(a)',iostat=ios) line
     if(ios<0) exit
     if(ios/=0) error stop 'reading susceptibility fixture'
     line=adjustl(line)
     if(index(line,'&spec_1')==1) species=1
     if(index(line,'&spec_2')==1) species=2
     if(index(line,'n_scan=')==1) line='n_scan=0'
     if(index(line,'profile_startup=')==1) line='profile_startup=F'
     if(index(line,'writeOut=')==1) line='writeOut=F'
     if(index(line,'fit_check=')==1) line='fit_check=F'
     if(index(line,'rel_cache=')==1) line="rel_cache='off'"
     if(index(line,'arrayName=')==1) line="arrayName='susceptibility'"
     if(index(line,'ngamma=')==1) write(line,'(a,i0)') 'ngamma=',ng
     if(index(line,'npparbar=')==1) write(line,'(a,i0)') 'npparbar=',nu
     if(index(line,'nperp=')==1) write(line,'(a,i0)') 'nperp=',nr
     if(index(line,'npar=')==1) write(line,'(a,i0)') 'npar=',nz
     if(index(line,'vA=')==1) write(line,'(a,es24.16)') 'vA=',va
     if(index(line,'mm=')==1) write(line,'(a,es24.16)') 'mm=',mass(species)
     if(index(line,'Bessel_zero')==1) line='Bessel_zero=1.d-12'
     if(index(line,'n_resonance_interval=')==1) line='n_resonance_interval=100'
     if(index(line,'fit_type_in=')==1.and.b/=0.d0) line='fit_type_in=5'
     if(index(line,'fit_2=')==1) write(line,'(a,es24.16)') 'fit_2=',b
     if(index(line,'fit_3=')==1) line='fit_3=0.d0'
     write(out,'(a)') trim(line)
  enddo
  close(unit); close(out)
  allocate(f(0:nr,0:nz),pperp(0:nr,0:nz),ppar(0:nr,0:nz))
  do is=1,2
     do i=0,nr
        do j=0,nz
           r=rmax*i/nr; u=-umax+2.d0*umax*j/nz; g=sqrt(1.d0+r*r+u*u)
           pperp(i,j)=r*mass(is)/va; ppar(i,j)=u*mass(is)/va
           f(i,j)=exp(-2.d0*g-b*u*u)
        enddo
     enddo
     dr=pperp(1,0)-pperp(0,0); du=ppar(0,1)-ppar(0,0)
     norm=2.d0*acos(-1.d0)*dr*du*sum(pperp*f)
     write(path,'(a,i0,a)') 'distribution/susceptibility.',is,'.array'
     open(newunit=out,file=trim(path),status='new')
     do i=0,nr
        do j=0,nz
           write(out,'(3es26.17e3)') pperp(i,j),ppar(i,j),f(i,j)/norm
        enddo
     enddo
     close(out)
  enddo
end program
