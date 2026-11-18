! Task 3F.3 test-only state and output. Never linked into normal production builds.
module spani_diagnostics
 use alps_analyt
 use legacy_glls_reference
 implicit none
 logical::capture=.false.
 complex(kind(1.d0))::previous_e(3)=(0.d0,0.d0)
 integer::iterations=0,query_unit=0,point=0,backend=0
 integer::seen_count(0:200)=0
 complex(kind(1.d0))::seen(0:200,2000)
contains
 subroutine prepare_models(choice)
 use alps_var
 integer,intent(in)::choice
 integer::code,r,j,file,cf,k,diag,strict
 integer,allocatable::meta(:,:)
 double precision,allocatable::mapping(:,:),coefs(:,:),x(:)
 type(cheb_model)::model,fit
 type(nr_cheb_row_info)::info
 backend=choice
 if(choice==3)then
 call initialize_nr_cheb_svd(NR_FIT_NATIVE,code)
 else
 call initialize_historical_nr(code)
 endif
 if(code/=0)error stop 'SPAN-I initialization failed'
 if(choice==2)then
 call nr_cheb_export(meta,mapping,coefs,code)
 allocate(x(0:npar));x=[(-1.d0+2.d0*j/npar,j=0,npar)]
 open(newunit=diag,file='index_svd.csv',status='new')
 write(diag,'(a)')'row,status,rank,degree,residual,rcond,singular_min,singular_max'
 do r=0,nperp
 call cheb_fit(x,log10(f0(1,r,:)),poly_order(1),fit,code,representation=CHEB_DIRECT)
 write(diag,'(*(g0,:,","))')r,code,fit%rank,fit%effective_degree,fit%residual_norm,fit%svd_tolerance,fit%singular_value_min,fit%singular_value_max
 if(code/=0.and.code/=1)error stop 'index SVD failure'
 k=1+r*nspec
 meta(:,k)=[NR_FIT_LEGACY_COMPAT,merge(NR_ROW_FITTED,NR_ROW_REDUCED,code==0),code,fit%requested_degree,fit%effective_degree, &
 CHEB_LOG10,fit%rank,fit%initial_rank,npar+1,npar+1,npar+1,0,0,merge(1,0,fit%sample_limited),merge(1,0,fit%rank_reduced),0,fit%fit_attempts]
 coefs(:,k)=0;coefs(1:size(fit%coefficients),k)=fit%coefficients
 mapping(3:6,k)=[pp(1,r,0,2),pp(1,r,npar,2),pp(1,r,0,2),pp(1,r,npar,2)]
 enddo
 close(diag)
 call nr_cheb_import(nspec,nperp,meta,mapping,coefs,code)
 if(code/=0)error stop 'index control import failure'
 endif
 end subroutine
 subroutine export_fits
 use alps_var
 integer::r,j,file,cf,nf,code
 double precision::v,e,peak,rel,lr,wp,wu,w,mom(4),ref(4),maxabs,base,vs
 complex(kind(1.d0))::z,f,d
 type(cheb_model)::model
 type(nr_cheb_row_info)::info
 peak=maxval(f0(1,:,:));mom=0;ref=0
 open(newunit=file,file='fit_rows.csv',status='new')
 open(newunit=cf,file='coefficients.csv',status='new')
 open(newunit=nf,file='nodes.csv',status='new')
 write(file,'(a)')'row,state,status,requested,effective,rank,original,retained,lower_zeros,upper_zeros,center,scale,fit_residual,diag_status,max_abs,relative_L1,log_RMS,suppressed_nodes'
 write(cf,'(a)')'row,order,coefficient'
 write(nf,'(a)')'row,node,pperp,u,original,reconstructed,status'
 do r=0,nperp
 call get_nr_cheb_row(1,r,model,info,code)
 maxabs=0;rel=0;lr=0;base=0;vs=0
 wp=.5d0*(pp(1,min(r+1,nperp),0,1)-pp(1,max(0,r-1),0,1))
 do j=0,npar
 z=cmplx(pp(1,r,j,2),0.d0,kind(1.d0))
 f=eval_fit(1,r,z,status=code)
 v=real(f);e=abs(v-f0(1,r,j));maxabs=max(maxabs,e);rel=rel+e;base=base+f0(1,r,j)
 if(v>0)then
 lr=lr+(log10(v)-log10(f0(1,r,j)))**2
 else
 vs=vs+1
 endif
 wu=.5d0*(pp(1,r,min(npar,j+1),2)-pp(1,r,max(0,j-1),2))
 w=wp*wu*2*pi*pp(1,r,j,1)
 ref=ref+w*f0(1,r,j)*[1.d0,real(z),real(z)**2,pp(1,r,j,1)**2]
 mom=mom+w*v*[1.d0,real(z),real(z)**2,pp(1,r,j,1)**2]
 write(nf,'(*(g0,:,","))')r,j,pp(1,r,j,1),real(z),f0(1,r,j),v,code
 enddo
 write(file,'(*(g0,:,","))')r,info%state,model%fit_status,model%requested_degree,model%effective_degree,model%rank, &
 model%original_sample_count,model%retained_sample_count,model%lower_tail_zeros,model%upper_tail_zeros,model%center,model%scale, &
 model%residual_norm,model%residual_norm_status,maxabs,rel/base,sqrt(lr/max(1.d0,npar+1-vs)),int(vs)
 do j=0,model%effective_degree
 write(cf,'(*(g0,:,","))')r,j,model%coefficients(j)
 enddo
 enddo
 close(file);close(cf);close(nf)
 open(newunit=file,file='moments.csv',status='new')
 write(file,'(a)')'kind,density,first_parallel,second_parallel,second_perpendicular'
 write(file,'(*(g0,:,","))')'original',ref
 write(file,'(*(g0,:,","))')'reconstructed',mom
 close(file)
 end subroutine
 subroutine record_result(omega,det,label)
 use alps_var
 use alps_fns,only:calc_eigen
 use mpi
 character(*),intent(in)::label
 complex(kind(1.d0)),intent(in)::omega,det
 integer::file,s,i,j
 complex(kind(1.d0))::electric(3),magnetic(3),vm(3,nspec),dn(nspec)
 double precision::ps(nspec),pss(4,nspec),energy,scaled,alignment
 if(.not.proc0)return
 open(newunit=file,file='tensors.csv',position='append',status='unknown')
 do s=1,nspec
 do i=1,3
 do j=i,3
 write(file,'(*(g0,:,","))')label,point,kpar,real(omega),aimag(omega),s,i,j,real(chi0(s,i,j)),aimag(chi0(s,i,j))
 enddo
 enddo
 enddo
 do i=1,3
 do j=1,3
 write(file,'(*(g0,:,","))')label,point,kpar,real(omega),aimag(omega),0,i,j,real(wave(i,j)),aimag(wave(i,j))
 enddo
 enddo
 close(file)
 call calc_eigen(omega,electric,magnetic,vm,dn,ps,pss,energy,.false.,.false.)
 scaled=abs(det)/max(maxval(abs(wave))**3,tiny(1.d0))
 open(newunit=file,file='results.csv',position='append',status='unknown')
 write(file,'(*(g0,:,","))')label,point,kpar,real(omega),aimag(omega),iterations,real(det),aimag(det),scaled, &
 real(electric),aimag(electric),real(magnetic),aimag(magnetic)
 close(file)
 if(label=='root')previous_e=electric
 if(label=='scan')then
 alignment=abs(sum(conjg(previous_e)*electric))/sqrt(sum(abs(previous_e)**2)*sum(abs(electric)**2))
 if(scaled>1.d-7.or.alignment<.95d0.or.iterations>=49)then
 write(*,*)'SPAN_BRANCH_STOP',point,scaled,alignment,iterations
 call mpi_abort(MPI_COMM_WORLD,97,ierror)
 endif
 previous_e=electric
 endif
 end subroutine
end module
subroutine span_query(s,r,u,wants_derivative,code)
 use spani_diagnostics,only:capture,query_unit,seen,seen_count
 use alps_var,only:iproc
 implicit none
 integer,intent(in)::s,r,code
 logical,intent(in)::wants_derivative
 double complex,intent(in)::u
 integer::j
 character(64)::name
 if(.not.capture)return
 if(r<0.or.r>200)error stop 'SPAN-I instrumentation row limit'
 if(code==0)then
 do j=1,seen_count(r)
 if(seen(r,j)==u)return
 enddo
 endif
 if(query_unit==0)then
 write(name,'(a,i0,a)')'queries_rank',iproc,'.csv'
 open(newunit=query_unit,file=trim(name),status='new')
 write(query_unit,'(a)')'species,row,u_re,u_im,derivative_requested,status'
 endif
 write(query_unit,'(*(g0,:,","))')s,r,real(u),aimag(u),wants_derivative,code
 flush(query_unit)
 if(seen_count(r)<2000)then
 seen_count(r)=seen_count(r)+1;seen(r,seen_count(r))=u
 else
 error stop 'SPAN-I query dedup capacity exceeded'
 endif
end subroutine
subroutine span_step(j)
 use spani_diagnostics,only:iterations
 integer,intent(in)::j
 iterations=j
end subroutine
subroutine span_scan_result(omega,det)
 use spani_diagnostics,only:record_result,point
 double complex,intent(in)::omega,det
 point=point+1
 call record_result(omega,det,'scan')
end subroutine
