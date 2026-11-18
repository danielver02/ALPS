! Test-only observations; never linked into production ALPS.
module rel_cheb_continuation_diagnostics
  use alps_var, only : nspec,ngamma,npparbar,gamma_rel,pparbar_rel,f0_rel,kpar,vA
  use alps_analyt
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  integer, parameter :: dp=kind(1.d0),qp=selected_real_kind(30,300),max_queries=20000
  public :: diagnostic_open,diagnostic_capture,diagnostic_query,diagnostic_stencil,diagnostic_close
  type(cheb_model), allocatable :: table_models(:,:),shape_models(:,:)
  real(dp) :: amplitude(2),quadratic,nan
  integer :: nodes,coeffs,queries,terms,partials,stencils,metadata,seen_count=0,current_degree=-1
  integer :: seen_species(max_queries),seen_rows(max_queries),seen_point(max_queries)
  complex(dp) :: seen_u(max_queries)
contains
  subroutine diagnostic_open(c,b)
    real(dp), intent(in) :: c(2),b
    amplitude=c; quadratic=b; nan=ieee_value(0.d0,ieee_quiet_nan)
    allocate(table_models(nspec,0:ngamma),shape_models(nspec,0:ngamma))
    open(newunit=nodes,file='diagnostic_nodes.csv',status='new')
    open(newunit=coeffs,file='diagnostic_coefficients.csv',status='new')
    open(newunit=queries,file='diagnostic_queries.csv',status='new')
    open(newunit=terms,file='diagnostic_terms.csv',status='new')
    open(newunit=partials,file='diagnostic_truncation.csv',status='new')
    open(newunit=stencils,file='diagnostic_stencils.csv',status='new')
    open(newunit=metadata,file='diagnostic_rows.csv',status='new')
    write(nodes,'(a)') 'requested,effective,species,row,node,u,table,analytic,fit,fit_status,fit_log,analytic_log,table_minus_analytic_log,fit_minus_table_log'
    write(coeffs,'(a)') 'requested,effective,species,row,n,coefficient,abs_coefficient,shape_fit_coefficient,exact_shape_coefficient,table_minus_exact,shape_fit_minus_exact'
    write(queries,'(a)') 'requested,effective,species,row,point,query,u_re,u_im,x_re,x_im,abs_x,real_projection_inside,p_re,p_im,p_status,p_u_re,p_u_im,value_re,value_im,value_status,derivative_re,derivative_im,with_derivative_status,derivative_status,analytic_re,analytic_im,relative_error,analytic_du_re,analytic_du_im,log_magnitude,phase,log_derivative_magnitude,sum_abs_terms,quad_p_re,quad_p_im,clenshaw_quad_abs,shape_fit_relative_error,shape_fit_status'
    write(terms,'(a)') 'requested,effective,species,row,point,query,u_re,u_im,n,coefficient,T_re,T_im,abs_T,term_re,term_im,abs_term,shape_fit_abs_term,exact_shape_abs_term'
    write(partials,'(a)') 'source_requested,source_effective,species,row,point,query,u_re,u_im,partial_degree,p_re,p_im,p_status,value_re,value_im,value_status,derivative_re,derivative_im,with_derivative_status,relative_error,quad_p_re,quad_p_im,clenshaw_quad_abs'
    write(stencils,'(a)') 'requested,species,parent_row,point,harmonic,u_re,u_im,status_gp,status_gm,status_up,status_um,fg_re,fg_im,fu_re,fu_im,numerator_re,numerator_im,analytic_fd_fg_re,analytic_fd_fg_im,analytic_fd_fu_re,analytic_fd_fu_im,analytic_fd_numerator_re,analytic_fd_numerator_im,fg_relative_error,fu_relative_error'
    write(metadata,'(a)') 'requested,effective,species,row,gamma,samples,rank,sample_limited,rank_reduced,center,scale,left,right,svd_tolerance,sigma_max,sigma_min,condition,log_residual,shape_fit_log_residual'
  end subroutine

  function shape(is,g,u) result(f)
    integer, intent(in) :: is,g
    complex(dp), intent(in) :: u
    complex(dp) :: f
    f=amplitude(is)*exp(-2.d0*gamma_rel(is,g,0)-quadratic*u*u)
  end function

  function exact_coeff(model,is,g,n) result(c)
    type(cheb_model), intent(in) :: model
    integer, intent(in) :: is,g,n
    real(dp) :: c
    c=0.d0
    select case(n)
    case(0); c=log10(amplitude(is))-(2.d0*gamma_rel(is,g,0)+quadratic*(model%center**2+.5d0*model%scale**2))/log(10.d0)
    case(1); c=-2.d0*quadratic*model%center*model%scale/log(10.d0)
    case(2); c=-.5d0*quadratic*model%scale**2/log(10.d0)
    end select
  end function

  subroutine diagnostic_capture(is,g,model)
    integer, intent(in) :: is,g
    type(cheb_model), intent(in) :: model
    real(dp), allocatable :: u(:),f(:)
    real(dp) :: value,polynomial,want,c
    integer :: count,j,n,code
    if(model%requested_degree/=current_degree) then
       current_degree=model%requested_degree; seen_count=0
    endif
    table_models(is,g)=model ! Diagnostic copy made during setup only.
    if(g/=2.and.g/=50) return
    count=count_valid(is,g); allocate(u(count),f(count)); count=0
    do j=0,npparbar
       if(f0_rel(is,g,j)==-1.d0) cycle
       count=count+1; u(count)=pparbar_rel(is,g,j)
       f(count)=real(shape(is,g,cmplx(u(count),0.d0,dp)),dp)
    enddo
    call cheb_fit(u,f,model%requested_degree,shape_models(is,g),code,representation=CHEB_LOG10)
    if(code/=CHEB_SUCCESS.and.code/=CHEB_REDUCED) error stop 'analytic-shape control fit unavailable'
    if(shape_models(is,g)%effective_degree/=model%effective_degree) error stop 'control/table ranks differ'
    write(metadata,'(*(g0,:,","))') model%requested_degree,model%effective_degree,is,g,gamma_rel(is,g,0), &
         model%sample_count,model%rank,model%sample_limited,model%rank_reduced,model%center,model%scale, &
         model%retained_min,model%retained_max,model%svd_tolerance,model%singular_value_max,model%singular_value_min, &
         model%singular_value_max/model%singular_value_min,model%residual_norm,shape_models(is,g)%residual_norm
    do n=0,model%effective_degree
       c=exact_coeff(model,is,g,n)
       write(coeffs,'(*(g0,:,","))') model%requested_degree,model%effective_degree,is,g,n,model%coefficients(n), &
            abs(model%coefficients(n)),shape_models(is,g)%coefficients(n),c,model%coefficients(n)-c, &
            shape_models(is,g)%coefficients(n)-c
    enddo
    do j=0,npparbar
       if(f0_rel(is,g,j)==-1.d0) cycle
       call cheb_evaluate(model,pparbar_rel(is,g,j),value,code)
       call cheb_evaluate_polynomial(model,pparbar_rel(is,g,j),polynomial,n)
       want=real(shape(is,g,cmplx(pparbar_rel(is,g,j),0.d0,dp)),dp)
       write(nodes,'(*(g0,:,","))') model%requested_degree,model%effective_degree,is,g,j,pparbar_rel(is,g,j), &
            f0_rel(is,g,j),want,value,code,polynomial,log10(want),log10(f0_rel(is,g,j))-log10(want), &
            polynomial-log10(f0_rel(is,g,j))
    enddo
  end subroutine

  integer function count_valid(is,g)
    integer, intent(in) :: is,g
    count_valid=count(f0_rel(is,g,:)/=-1.d0)
  end function

  subroutine diagnostic_query(point,is,g,u,adapter_value,adapter_status)
    integer, intent(in) :: point,is,g,adapter_status
    complex(dp), intent(in) :: u,adapter_value
    complex(dp) :: x,p,pu,f,fd,d,ref,refdu,shape_f,quad_p
    complex(qp) :: xq,t0,t1,t2,sumq
    complex(dp) :: tn(0:40)
    real(dp) :: relative,shape_relative,logd,sum_abs
    integer :: i,n,ps,vs,ds,overall,shape_status,dummy
    logical :: inside
    type(cheb_model) :: m
    ! Deduplicate identical actual queries within each independently fitted degree.
    do i=1,seen_count
       if(seen_species(i)==is.and.seen_rows(i)==g.and.seen_point(i)==point) then
          if(seen_u(i)==u) return
       endif
    enddo
    seen_count=seen_count+1
    if(seen_count>max_queries) error stop 'diagnostic query capacity exceeded'
    seen_species(seen_count)=is; seen_rows(seen_count)=g; seen_point(seen_count)=point; seen_u(seen_count)=u
    ! ASSOCIATE avoids a new coefficient copy at each diagnostic query.
    associate(model=>table_models(is,g))
    n=model%effective_degree; x=(u-model%center)/model%scale
    call cheb_evaluate_polynomial(model,u,p,ps,pu)
    call cheb_evaluate(model,u,f,vs)
    if(vs/=adapter_status) error stop 'adapter/kernel value status mismatch'
    if(vs==CHEB_SUCCESS.or.vs==CHEB_LOG_UNDERFLOW) then
       if(f/=adapter_value) error stop 'adapter/kernel value mismatch'
    endif
    call cheb_evaluate(model,u,fd,overall,d,derivative_status=ds)
    if(ps/=CHEB_SUCCESS) error stop 'raw polynomial unexpectedly unavailable'
    ref=shape(is,g,u); refdu=-2.d0*quadratic*u*ref
    relative=nan
    if(vs==CHEB_SUCCESS.or.vs==CHEB_LOG_UNDERFLOW) relative=abs(f-ref)/abs(ref)
    inside=real(u,dp)>=model%retained_min.and.real(u,dp)<=model%retained_max
    ! Independent REAL128 forward series diagnoses evaluation arithmetic only.
    xq=cmplx(real(x,dp),aimag(x),qp); t0=1._qp; t1=xq; sumq=real(model%coefficients(0),qp)
    tn(0)=1.d0
    if(n>0) tn(1)=x
    do i=1,n
       if(i==1) then
          sumq=sumq+real(model%coefficients(i),qp)*t1
       else
          t2=2._qp*xq*t1-t0
          sumq=sumq+real(model%coefficients(i),qp)*t2
          t0=t1; t1=t2
          tn(i)=2.d0*x*tn(i-1)-tn(i-2)
       endif
    enddo
    quad_p=cmplx(real(sumq,qp),aimag(sumq),dp)
    sum_abs=sum(abs(model%coefficients*tn(0:n)))
    if(abs(p-quad_p)>1.d-11*max(1.d0,sum_abs)) error stop 'Clenshaw versus REAL128 series mismatch'
    logd=nan
    if(abs(pu)>0.d0) logd=real(p,dp)*log(10.d0)+log(log(10.d0))+log(abs(pu))
    shape_relative=nan; shape_status=CHEB_NOT_REQUESTED
    if(g==2.or.g==50) then
       call cheb_evaluate(shape_models(is,g),u,shape_f,shape_status)
       if(shape_status==CHEB_SUCCESS) shape_relative=abs(shape_f-ref)/abs(ref)
    endif
    write(queries,'(*(g0,:,","))') model%requested_degree,n,is,g,point,seen_count,real(u,dp),aimag(u), &
         real(x,dp),aimag(x),abs(x),inside,real(p,dp),aimag(p),ps,real(pu,dp),aimag(pu), &
         real(f,dp),aimag(f),vs,real(d,dp),aimag(d),overall,ds,real(ref,dp),aimag(ref),relative, &
         real(refdu,dp),aimag(refdu),real(p,dp)*log(10.d0),aimag(p)*log(10.d0),logd, &
         sum_abs,real(quad_p,dp),aimag(quad_p),abs(p-quad_p),shape_relative,shape_status
    if(g/=2.and.g/=50) return
    do i=0,n
       write(terms,'(*(g0,:,","))') model%requested_degree,n,is,g,point,seen_count,real(u,dp),aimag(u), &
            i,model%coefficients(i),real(tn(i),dp),aimag(tn(i)),abs(tn(i)), &
            real(model%coefficients(i)*tn(i),dp),aimag(model%coefficients(i)*tn(i)), &
            abs(model%coefficients(i)*tn(i)),abs(shape_models(is,g)%coefficients(i)*tn(i)), &
            abs(exact_coeff(model,is,g,i)*tn(i))
    enddo
    if(model%requested_degree/=40) return
    ! Partial sums of ONE accepted high-degree fit, never a lower-degree refit.
    m=model
    do i=0,n
       m%coefficients=model%coefficients
       if(i<n) m%coefficients(i+1:n)=0.d0
       ! Direct differentiation identity: d_j=2*sum_{k=j+1,j+3,...} k*c_k;
       ! d_0 has half this weight. Used ONLY to evaluate diagnostic partial sums.
       ! Keep source allocation/rank metadata; zeroed tails specify the partial series.
       call partial_derivative(m)
       call cheb_evaluate_polynomial(m,u,p,ps)
       call cheb_evaluate(m,u,f,vs)
       call cheb_evaluate(m,u,fd,overall,d,derivative_status=ds)
       relative=nan
       if(vs==CHEB_SUCCESS.or.vs==CHEB_LOG_UNDERFLOW) relative=abs(f-ref)/abs(ref)
       xq=cmplx(real(x,dp),aimag(x),qp); t0=1._qp; t1=xq; sumq=real(m%coefficients(0),qp)
       do dummy=1,i
          if(dummy>1) then
             t2=2._qp*xq*t1-t0; t0=t1; t1=t2
          endif
          sumq=sumq+real(m%coefficients(dummy),qp)*t1
       enddo
       quad_p=cmplx(real(sumq,qp),aimag(sumq),dp)
       write(partials,'(*(g0,:,","))') model%requested_degree,n,is,g,point,seen_count,real(u,dp),aimag(u), &
            i,real(p,dp),aimag(p),ps,real(f,dp),aimag(f),vs,real(d,dp),aimag(d),overall,relative, &
            real(quad_p,dp),aimag(quad_p),abs(p-quad_p)
    enddo
    end associate
  end subroutine

  subroutine partial_derivative(model)
    type(cheb_model), intent(inout) :: model
    integer :: j,k,n
    n=model%effective_degree; model%derivative_coefficients=0.d0
    do j=0,n-1
       do k=j+1,n,2
          model%derivative_coefficients(j)=model%derivative_coefficients(j)+2.d0*k*model%coefficients(k)
       enddo
    enddo
    model%derivative_coefficients(0)=.5d0*model%derivative_coefficients(0)
  end subroutine

  subroutine diagnostic_stencil(point,is,g,harmonic,omega,u,du)
    integer, intent(in) :: point,is,g,harmonic
    complex(dp), intent(in) :: omega,u
    real(dp), intent(in) :: du
    complex(dp) :: gp,gm,up,um,fg,fu,numerator,agp,agm,aup,aum,afg,afu,anum
    real(dp) :: dg,er_g,er_u
    integer :: sgp,sgm,sup,sum_code,lo
    dg=gamma_rel(is,2,2)-gamma_rel(is,1,2); lo=g-1
    if(g==1) lo=g
    call cheb_evaluate(table_models(is,g+1),u,gp,sgp)
    call cheb_evaluate(table_models(is,lo),u,gm,sgm)
    call cheb_evaluate(table_models(is,g),u+du,up,sup)
    call cheb_evaluate(table_models(is,g),u-du,um,sum_code)
    agp=shape(is,g+1,u); agm=shape(is,lo,u); aup=shape(is,g,u+du); aum=shape(is,g,u-du)
    afg=(agp-agm)/((g+1-lo)*dg); afu=(aup-aum)/(2.d0*du)
    anum=omega*afg+(kpar/vA)*afu
    fg=cmplx(nan,nan,dp); fu=fg; numerator=fg; er_g=nan; er_u=nan
    if(sgp==CHEB_SUCCESS.and.sgm==CHEB_SUCCESS) then
       fg=(gp-gm)/((g+1-lo)*dg)
       if(abs(afg)>0.d0) er_g=abs(fg-afg)/abs(afg)
    endif
    if(sup==CHEB_SUCCESS.and.sum_code==CHEB_SUCCESS) then
       fu=(up-um)/(2.d0*du)
       if(abs(afu)>0.d0) er_u=abs(fu-afu)/abs(afu)
    endif
    if(all([sgp,sgm,sup,sum_code]==CHEB_SUCCESS)) numerator=omega*fg+(kpar/vA)*fu
    write(stencils,'(*(g0,:,","))') current_degree,is,g,point,harmonic,real(u,dp),aimag(u),sgp,sgm,sup,sum_code, &
         real(fg,dp),aimag(fg),real(fu,dp),aimag(fu),real(numerator,dp),aimag(numerator), &
         real(afg,dp),aimag(afg),real(afu,dp),aimag(afu),real(anum,dp),aimag(anum),er_g,er_u
  end subroutine

  subroutine diagnostic_close
    close(nodes); close(coeffs); close(queries); close(terms); close(partials); close(stencils); close(metadata)
    deallocate(table_models,shape_models)
    print '(a)', 'PASS: coefficient/query/truncation observations; public Clenshaw agrees with independent REAL128 series'
  end subroutine
end module
