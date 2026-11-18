! Geometry oracle fixtures only; no transformation or altered support policy.
program test_rel_source_support
  use alps_var
  use alps_fns_rel, only: report_rel_source_support
  implicit none
  integer :: i,j,test
  double precision :: rmin,rmax,umin,umax,gmax,g,u
  proc0=.true.;ngamma=50;npparbar=200
  allocate(pp(1,0:1,0:1,2),ms(1),ACmethod(1),gamma_rel(1,0:50,0:200),pparbar_rel(1,0:50,0:200))
  ACmethod=2;ms=2.d0;vA=.5d0
  do test=1,2
     rmin=0.d0;rmax=5.d0;umin=-5.d0;umax=5.d0
     if(test==2) then
        rmin=.4d0;rmax=1.2d0;umin=-.7d0;umax=1.3d0
     endif
     do i=0,1
        do j=0,1
           pp(1,i,j,1)=(rmin+(rmax-rmin)*i)*ms(1)/vA
           pp(1,i,j,2)=(umin+(umax-umin)*j)*ms(1)/vA
        enddo
     enddo
     gmax=sqrt(1.d0+rmax*rmax+max(umin*umin,umax*umax))
     do i=0,50
        g=1.d0+(gmax-1.d0)*i/50.d0
        do j=0,200
           u=umin+(umax-umin)*j/200.d0
           gamma_rel(1,i,j)=g;pparbar_rel(1,i,j)=u
        enddo
     enddo
     write(*,'(a,i0)') 'SUPPORT_FIXTURE=',test
     call report_rel_source_support(1,1,.false.)
  enddo
  ! Four points on source boundaries, perturbed less than the stated tolerance.
  deallocate(gamma_rel,pparbar_rel)
  allocate(gamma_rel(1,0:0,0:3),pparbar_rel(1,0:0,0:3))
  do j=0,3
     u=umin
     if(j>=2)u=umax
     g=rmin
     if(mod(j,2)==1)g=rmax
     pparbar_rel(1,0,j)=u+2.d0*epsilon(1.d0)
     gamma_rel(1,0,j)=sqrt(1.d0+g*g+u*u)
  enddo
  write(*,'(a)') 'SUPPORT_FIXTURE=3'
  call report_rel_source_support(1,1,.true.)
end program
