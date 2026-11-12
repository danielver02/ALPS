! Exercise current construction/cache loading, then the explicit root-side row
! builder. No production guard, continuation, or integrator is bypassed/changed.
program test_rel_chebyshev_transformed
  use alps_var, only : f0,pp,ms,vA,pi,nperp,npar,ngamma,npparbar,nspec,nspec_rel, &
       gamma_rel,pparbar_rel,f0_rel,df0_rel,proc0,writeOut,profile_startup,arrayName, &
       rel_cache,rel_cache_dir,poly_order
  use alps_analyt
  use alps_fns_rel, only : derivative_f0_rel
  use alps_rel_cache, only : rel_cache_identity,cache_identity,cache_load
  use relativistic_references, only : anisotropic
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  integer, parameter :: dp=kind(1.d0)
  type(rel_cache_identity) :: identities(2)
  type(cheb_model), allocatable :: fresh_models(:,:)
  type(cheb_model) :: model
  type(rel_cheb_row_info) :: info
  real(dp), allocatable :: saved_g(:,:,:),saved_u(:,:,:),saved_f(:,:,:),saved_df(:,:,:,:)
  real(dp), allocatable :: shapes(:,:),rows(:)
  real(dp) :: a,b,g,u,r,dg,du,norm,x,value,derivative,abs_err,rel_err,log_res,du_err,interp_err,peak
  real(dp) :: vs,ds,scale,complex_err
  complex(dp) :: z,want,fg,fu,got,dgot
  integer :: s,i,j,status,passed=0,failed=0,min_degree,max_degree,retained,lower_zeros,upper_zeros
  integer :: reduced,sample_reduced,rank_reduced
  logical :: identified,hit
  character(256) :: reason

  ! Current production transformation is tested at non-unit mass and vA.
  ! Profiling disabled means these serial calls require no MPI initialization.
  proc0=.true.; writeOut=.false.; profile_startup=.false.
  nspec=3; nspec_rel=2; nperp=8; npar=16; ngamma=20; npparbar=40
  vA=.6d0; pi=4.d0*atan(1.d0); arrayName='row-layer-fixture'
  allocate(ms(nspec),poly_order(nspec),f0(nspec,0:nperp,0:npar),pp(nspec,0:nperp,0:npar,2))
  ms=[1.2d0,9.d0,2.d0]; poly_order=[0,7,12]
  f0=0.d0; pp=0.d0
  allocate(gamma_rel(2,0:ngamma,0:npparbar),pparbar_rel(2,0:ngamma,0:npparbar), &
       f0_rel(2,0:ngamma,0:npparbar),df0_rel(2,0:ngamma,0:npparbar,2))
  allocate(fresh_models(2,0:ngamma),shapes(0:ngamma,0:npparbar),rows(0:ngamma))
  do s=1,3,2
     a=vA/ms(s); b=0.d0
     if(s==3) b=.15d0
     do i=0,nperp
        do j=0,npar
           r=1.2d0*i/nperp; u=-1.d0+2.d0*j/npar; g=sqrt(1.d0+r*r+u*u)
           pp(s,i,j,1)=r/a; pp(s,i,j,2)=u/a
           f0(s,i,j)=exp(-2.d0*g-b*u*u)
        enddo
     enddo
  enddo
  ! No cache implementation changes: the existing APIs use only this run's
  ! scratch cache directory, never a source-tree or pre-existing cache entry.
  rel_cache='refresh'; rel_cache_dir='cache/'
  do s=1,2
     i=2*s-1
     call cache_identity(pp(i,:,:,:),f0(i,:,:),ms(i),vA,pi,ngamma,npparbar,identities(s),identified,reason)
     if(.not.identified) then
        write(*,'(a)') 'Cache identity unavailable: '//trim(reason)
        error stop 'fixture requires existing Linux/OpenBLAS cache support'
     endif
     call derivative_f0_rel(i,s)
  enddo
  saved_g=gamma_rel; saved_u=pparbar_rel; saved_f=f0_rel; saved_df=df0_rel
  call build_rel_cheb_rows(gamma_rel,pparbar_rel,f0_rel,[1,3],poly_order([1,3]),status)
  call check(status==CHEB_SUCCESS,'completed transformed table row construction')
  call check(all(saved_g==gamma_rel).and.all(saved_u==pparbar_rel).and.all(saved_f==f0_rel).and. &
       all(saved_df==df0_rel),'builder preserves transformed coordinates/values/derivatives')
  do s=1,2
     b=0.d0
     if(s==2) b=.15d0
     dg=gamma_rel(s,1,0)-gamma_rel(s,0,0); du=pparbar_rel(s,0,1)-pparbar_rel(s,0,0)
     shapes=exp(-2.d0*gamma_rel(s,:,:)-b*pparbar_rel(s,:,:)**2)
     ! Independent adjacent-panel normalization of the analytic reference,
     ! using the current target-grid geometry and full Gamma endpoint convention.
     rows=0.d0
     do i=0,ngamma
        do j=0,npparbar-1
           if(f0_rel(s,i,j)==-1.d0.or.f0_rel(s,i,j+1)==-1.d0) cycle
           rows(i)=rows(i)+.5d0*(shapes(i,j)+shapes(i,j+1))*du*gamma_rel(s,i,j)
        enddo
     enddo
     norm=sum(.5d0*(rows(0:ngamma-1)+rows(1:ngamma)))*dg*2.d0*pi*(ms(2*s-1)/vA)**3
     peak=maxval(f0_rel(s,:,:))
     abs_err=0.d0; rel_err=0.d0; log_res=0.d0; du_err=0.d0; interp_err=0.d0; complex_err=0.d0
     min_degree=huge(0); max_degree=0; retained=0; lower_zeros=0; upper_zeros=0
     reduced=0; sample_reduced=0; rank_reduced=0
     do i=0,ngamma
        call get_rel_cheb_row(s,i,model,info,status)
        call check(status==CHEB_SUCCESS.or.status==CHEB_REDUCED,'actual transformed row available')
        fresh_models(s,i)=model
        min_degree=min(min_degree,model%effective_degree); max_degree=max(max_degree,model%effective_degree)
        retained=retained+model%retained_sample_count
        lower_zeros=lower_zeros+model%lower_tail_zeros; upper_zeros=upper_zeros+model%upper_tail_zeros
        if(status==CHEB_REDUCED) reduced=reduced+1
        if(model%sample_limited) sample_reduced=sample_reduced+1
        if(model%rank_reduced) rank_reduced=rank_reduced+1
        log_res=max(log_res,model%residual_norm)
        do j=0,npparbar
           if(f0_rel(s,i,j)==-1.d0) cycle
           call cheb_evaluate(model,pparbar_rel(s,i,j),value,status,derivative)
           call check(status==CHEB_SUCCESS,'actual retained-node evaluation status')
           scale=max(f0_rel(s,i,j),1.d-12*peak)
           abs_err=max(abs_err,abs(value-f0_rel(s,i,j)))
           rel_err=max(rel_err,abs(value-f0_rel(s,i,j))/scale)
           x=shapes(i,j)/norm
           interp_err=max(interp_err,abs(f0_rel(s,i,j)-x)/max(x,1.d-12*peak))
           ! Isotropic derivative is exactly zero. The anisotropic singleton
           ! is an intentionally underdetermined constant and is not an oracle.
           if(s==1.or.model%effective_degree>=2) then
              du_err=max(du_err,abs(derivative+2.d0*b*pparbar_rel(s,i,j)*x)/peak)
           endif
        enddo
        if(s==1.or.(i>1.and.model%effective_degree>=2)) then
           z=cmplx(.15d0,-.025d0,dp)
           call cheb_evaluate(model,z,got,status,dgot)
           call check(status==CHEB_SUCCESS,'actual model complex evaluation status')
           call anisotropic(gamma_rel(s,i,0),z,2.d0,b,1.d0/norm,want,fg,fu)
           complex_err=max(complex_err,abs(got-want)/max(abs(want),1.d-12*peak))
        endif
     enddo
     ! These bounds assess fitting separately from interpolation. Analytic-row
     ! tests impose much tighter errors; this coarse spline fixture permits 1%
     ! row approximation and 3% peak-scaled derivative error for anisotropy.
     if(s==1) then
        call check(rel_err<1.d-10.and.du_err==0.d0.and.complex_err<1.d-10,'affine-log transformed isotropic accuracy')
     else
        call check(rel_err<1.d-2.and.du_err<3.d-2.and.complex_err<3.d-2,'coarse transformed anisotropic error budget')
     endif
     write(*,'(a,i0,a,6es13.4)') 'TRANSFORMED species=',2*s-1, &
          ' abs,rel,log_L2,du_scaled,interp_rel,complex_rel=',abs_err,rel_err,log_res,du_err,interp_err,complex_err
     write(*,'(a,2i4,a,i0,a,2i3,a,3i3)') 'TRANSFORMED degree range=',min_degree,max_degree, &
          ' retained_total=',retained,' tail_zeros=',lower_zeros,upper_zeros, &
          ' reduced/sample/rank=',reduced,sample_reduced,rank_reduced
     ! Exact row-level diagnostics, including rank and counts, for cone/interior/top rows.
     do i=0,ngamma,ngamma/2
        call get_rel_cheb_row(s,i,model,info,status)
        write(*,'(a,3i4,a,4i5,a,2l2,a,es13.4)') 'TRANSFORMED row species/packed/gamma=', &
             info%species,s,i,' degree/rank/unmasked/retained=',model%effective_degree,model%rank, &
             info%geometric_samples,model%retained_sample_count,' sample/rank reduction=', &
             model%sample_limited,model%rank_reduced,' log_L2=',model%residual_norm
     enddo
  enddo

  ! Exercise the actual cache-hit initialization path, then rebuild from its
  ! identical completed payload. This is a new VDF initialization, not a refit
  ! during evaluation, and does not add a polynomial cache.
  rel_cache='auto'
  gamma_rel=-777.d0; pparbar_rel=-777.d0; f0_rel=-777.d0; df0_rel=-777.d0
  do s=1,2
     call derivative_f0_rel(2*s-1,s)
     call cache_load(rel_cache_dir,identities(s),gamma_rel(s,:,:),pparbar_rel(s,:,:), &
          f0_rel(s,:,:),df0_rel(s,:,:,:),x,hit,reason)
     call check(hit,'actual existing relativistic cache hit')
  enddo
  call check(all(saved_g==gamma_rel).and.all(saved_u==pparbar_rel).and.all(saved_f==f0_rel).and. &
       all(saved_df==df0_rel),'cache-restored tables exactly equal constructed tables')
  call build_rel_cheb_rows(gamma_rel,pparbar_rel,f0_rel,[1,3],poly_order([1,3]),status)
  call check(status==CHEB_SUCCESS,'cache-restored row construction')
  do s=1,2
     do i=0,ngamma
        call get_rel_cheb_row(s,i,model,info,status)
        call check(model%effective_degree==fresh_models(s,i)%effective_degree.and. &
             model%rank==fresh_models(s,i)%rank.and.model%fit_status==fresh_models(s,i)%fit_status, &
             'cache/restored model rank and degree identical')
        call check(all(model%coefficients==fresh_models(s,i)%coefficients).and. &
             all(model%derivative_coefficients==fresh_models(s,i)%derivative_coefficients), &
             'cache/restored coefficients and derivatives byte-equivalent values')
        call check(model%center==fresh_models(s,i)%center.and.model%scale==fresh_models(s,i)%scale.and. &
             model%residual_norm==fresh_models(s,i)%residual_norm,'cache/restored mapping and residual identical')
     enddo
  enddo
  call release_rel_cheb_rows
  write(*,'(a,i0,a,i0)') 'Relativistic transformed rows: passed=',passed,', failed=',failed
  if(failed/=0) error stop 1
contains
  subroutine check(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if(condition) then
       passed=passed+1
    else
       failed=failed+1; write(*,'(a)') 'FAIL: '//label
    endif
  end subroutine
end program test_rel_chebyshev_transformed
