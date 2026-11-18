program test_coordinate_contract
  use relativistic_references
  use ieee_arithmetic, only : ieee_is_finite
  implicit none
  real(dp), parameter :: temperatures(4)=[0.01_dp,0.1_dp,1._dp,10._dp]
  real(dp) :: g,h,lambda,c,beta,theta,a,pp,pz,g2,u2,pp1,pp2,pz1,pz2,jac,expected
  real(dp) :: ecenter(2),eforward(2),eu(2),scale,value_error,stencil_error,coord_error
  complex(dp) :: u,f,fg,fu,fp,fm,dummy1,dummy2,expect_f,exponent,deriv
  integer :: model,it,iu,level,out,ibeta
  value_error=0; stencil_error=0; coord_error=0
  open(newunit=out,file='coordinate_contract.dat',status='new')
  write(out,'(a)') '# model theta beta Re(u) Im(u) centered_error forward_error u_scaled_error; h halves'
  do model=1,3
     do ibeta=1,merge(3,1,model==3)
     do it=1,4
        theta=temperatures(it)
        beta=0.1_dp+0.2_dp*(ibeta-1); c=1.7_dp; g=1.2_dp
        lambda=1/theta
        if(model==2) lambda=2._dp
        if(model==3) lambda=1/(theta*sqrt(1-beta*beta))
        do iu=1,2
           u=cmplx(0.2_dp,0._dp,dp)
           if(iu==2) u=cmplx(0.3_dp,-0.07_dp,dp)
           call evaluate(model,g,u,theta,beta,c,f,fg,fu)
           ! Independent exponential evaluation from real amplitude and phase.
           select case(model)
           case(1)
              exponent=cmplx(-g/theta,0._dp,dp)
           case(2)
              exponent=-2*g-0.15_dp*u*u
           case(3)
              exponent=-(g-beta*u)/(theta*sqrt(1-beta*beta))
           end select
           expect_f=c*exp(real(exponent,dp))*cmplx(cos(aimag(exponent)),sin(aimag(exponent)),dp)
           value_error=max(value_error,abs(f-expect_f)/abs(expect_f))
           call require(abs(f-expect_f)<=3.e-13_dp*abs(expect_f),'analytic complex value')
           do level=1,2
              h=0.1_dp/lambda/2**(level-1)
              call evaluate(model,g+h,u,theta,beta,c,fp,dummy1,dummy2)
              call evaluate(model,g-h,u,theta,beta,c,fm,dummy1,dummy2)
              deriv=(fp-fm)/(2*h)
              ecenter(level)=abs(deriv-fg)/abs(fg)
              expected=sinh(lambda*h)/(lambda*h)
              stencil_error=max(stencil_error,abs(deriv/fg-expected))
              call require(abs(deriv/fg-expected)<2.e-12_dp,'centered gamma stencil at fixed u')
              deriv=(fp-f)/h
              eforward(level)=abs(deriv-fg)/abs(fg)
              expected=(1-exp(-lambda*h))/(lambda*h)
              call require(abs(deriv/fg-expected)<3.e-12_dp,'first-row forward gamma stencil')
              h=0.01_dp/2**(level-1)
              call evaluate(model,g,u+h,theta,beta,c,fp,dummy1,dummy2)
              call evaluate(model,g,u-h,theta,beta,c,fm,dummy1,dummy2)
              eu(level)=abs((fp-fm)/(2*h)-fu)/max(abs(fu),abs(f))
              write(out,'(i2,8es24.15e3)') model,theta,beta,real(u,dp),aimag(u),ecenter(level), &
                   eforward(level),eu(level),h
           enddo
           call require(ecenter(1)/ecenter(2)>3.99_dp.and.ecenter(1)/ecenter(2)<4.01_dp, &
                'second-order centered gamma convergence')
           call require(eforward(1)/eforward(2)>1.96_dp.and.eforward(1)/eforward(2)<2.01_dp, &
                'first-order forward gamma convergence')
           if(model/=1) call require(eu(2)<eu(1)/3.8_dp,'parallel derivative convergence')
           if(model==1) call require(eu(1)==0._dp,'isotropic parallel derivative is zero')
        enddo
     enddo
     enddo
  enddo
  ! Fixed u can leave the neighbor's real cone; the analytic function remains defined.
  g=1.05_dp; h=0.01_dp; u=cmplx(0.3_dp,-0.02_dp,dp)
  call require(real(u,dp)**2<g*g-1.and.real(u,dp)**2>(g-h)**2-1, 'neighbor-interval fixture')
  call anisotropic(g-h,u,2._dp,0.15_dp,1._dp,f,fg,fu)
  call require(ieee_is_finite(real(f,dp)).and.ieee_is_finite(aimag(f)), 'neighbor analytic continuation')
  ! Along a pole u(g)=s*g+t, the total derivative differs from the required partial derivative.
  g=1.2_dp; h=1.e-5_dp; u=cmplx(0.2_dp,-0.05_dp,dp)
  call anisotropic(g,u,2._dp,0.15_dp,1._dp,f,fg,fu)
  call anisotropic(g+h,u+0.4_dp*h,2._dp,0.15_dp,1._dp,fp,dummy1,dummy2)
  call anisotropic(g-h,u-0.4_dp*h,2._dp,0.15_dp,1._dp,fm,dummy1,dummy2)
  call require(abs((fp-fm)/(2*h)-(fg+0.4_dp*fu))<1.e-9_dp*abs(f), 'total derivative identity')
  call require(abs(0.4_dp*fu)>0.01_dp*abs(f), 'partial versus pole-path derivative distinguished')
  ! Non-unit mass and vA catch normalization mistakes hidden by a_s=1 examples.
  a=0.3_dp
  call to_relativistic(0.7_dp,-0.4_dp,0.6_dp,2._dp,g,u2)
  call from_relativistic(g,u2,0.6_dp,2._dp,pp,pz)
  coord_error=max(abs(pp-0.7_dp),abs(pz+0.4_dp))
  call require(coord_error<2.e-14_dp,'momentum round trip')
  g=1.4_dp; u2=0.2_dp; h=1.e-5_dp
  call from_relativistic(g,u2,0.6_dp,2._dp,pp,pz)
  call from_relativistic(g+h,u2,0.6_dp,2._dp,pp1,pz1)
  call from_relativistic(g-h,u2,0.6_dp,2._dp,pp2,pz2)
  jac=pp*(pp1-pp2)/(2*h)/a
  expected=g/a**3
  call require(abs(jac/expected-1)<2.e-10_dp,'cylindrical coordinate Jacobian')
  write(out,'(a,es24.15e3)') '# maximum analytic value relative error ',value_error
  write(out,'(a,es24.15e3)') '# maximum centered-stencil factor error ',stencil_error
  write(out,'(a,es24.15e3)') '# coordinate round-trip absolute error ',coord_error
  write(out,'(a,es24.15e3)') '# finite-difference Jacobian relative error ',abs(jac/expected-1)
  close(out)
  print '(a)', 'PASS: analytic references, coordinate Jacobian, fixed-u derivative contract'
contains
  subroutine evaluate(model,g,u,theta,beta,c,f,fg,fu)
    integer, intent(in) :: model
    real(dp), intent(in) :: g,theta,beta,c
    complex(dp), intent(in) :: u
    complex(dp), intent(out) :: f,fg,fu
    select case(model)
    case(1)
       call juttner(g,u,theta,c,f,fg,fu)
    case(2)
       call anisotropic(g,u,2._dp,0.15_dp,c,f,fg,fu)
    case(3)
       call drifting(g,u,theta,beta,c,f,fg,fu)
    end select
  end subroutine
  subroutine require(ok,message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if(.not.ok) then
       print '(a)', 'FAIL: '//message
       error stop 1
    endif
  end subroutine
end program
