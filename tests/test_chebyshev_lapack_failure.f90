! Test-only replacement of DGELSD verifies errors without production test hooks.
module lapack_fault
  implicit none
  integer :: mode=0
end module

program test_chebyshev_lapack_failure
  use alps_analyt, only : cheb_model,cheb_fit,CHEB_LAPACK_FAILURE,CHEB_RANK_FAILURE,CHEB_NONFINITE_COEFFICIENTS, &
       build_rel_cheb_rows,get_rel_cheb_row,release_rel_cheb_rows,rel_cheb_row_info,REL_CHEB_FIT_FAILURE
  use lapack_fault
  implicit none
  type(cheb_model) :: model
  type(rel_cheb_row_info) :: info
  double precision :: gamma(1,0:0,0:1),u(1,0:0,0:1),f(1,0:0,0:1)
  integer :: status,expected,passed=0
  gamma=1.5d0; u(1,0,:)=[-.2d0,.2d0]; f(1,0,:)=[1.d0,2.d0]
  do mode=1,4
     call cheb_fit([0.d0,1.d0],[1.d0,2.d0],0,model,status)
     select case(mode)
     case(1)
        if(status/=CHEB_LAPACK_FAILURE.or.model%lapack_info/=-12) error stop 'query error propagation'
     case(2)
        if(status/=CHEB_LAPACK_FAILURE.or.model%lapack_info/=1) error stop 'SVD error propagation'
     case(3)
        if(status/=CHEB_RANK_FAILURE) error stop 'degree-zero rank failure'
     case(4)
        if(status/=CHEB_NONFINITE_COEFFICIENTS) error stop 'nonfinite LAPACK result'
     end select
     if(allocated(model%coefficients).or.allocated(model%derivative_coefficients)) error stop 'failed model retained storage'
     if(model%fit_status/=status) error stop 'failure not recorded in model'
     passed=passed+1
     expected=status
     call build_rel_cheb_rows(gamma,u,f,[3],[0],status)
     if(status/=expected) error stop 'row builder did not propagate solver failure'
     call get_rel_cheb_row(1,0,model,info,status)
     if(status/=expected.or.info%state/=REL_CHEB_FIT_FAILURE) error stop 'row numerical failure classification'
     if(allocated(model%coefficients).or.allocated(model%derivative_coefficients)) error stop 'failed row retained storage'
     if(info%species/=3.or.info%gamma_index/=0) error stop 'failed row lost diagnostic identity'
     passed=passed+1
  enddo
  call release_rel_cheb_rows
  print '(a,i0,a)', 'Chebyshev LAPACK fault checks: passed=',passed,', failed=0'
end program

subroutine dgelsd(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,iwork,info)
  use lapack_fault
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, intent(in) :: m,n,nrhs,lda,ldb,lwork
  double precision, intent(inout) :: a(lda,*),b(ldb,*)
  double precision, intent(out) :: s(*),work(*)
  double precision, intent(in) :: rcond
  integer, intent(out) :: rank,iwork(*),info
  info=0; rank=n; work(1)=1.d0; iwork(1)=1
  if(lwork==-1) then
     if(mode==1) info=-12
     return
  endif
  if(mode==2) info=1
  if(mode==3) rank=0
  if(mode==4) b(1,1)=ieee_value(0.d0,ieee_quiet_nan)
  s(1:n)=1.d0
end subroutine
