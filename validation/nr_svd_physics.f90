! Test-only hooks used in a scratch copy of the established NR baseline driver.
subroutine nr_svd_prepare
 use alps_analyt
 use alps_var,only:proc0,nspec,nperp
 use legacy_glls_reference
 implicit none
 character(32)::backend
 integer::code,s,r,file,cfile,j
 type(cheb_model)::model
 type(nr_cheb_row_info)::info
 call get_environment_variable('ALPS_NR_SVD_BACKEND',backend)
 if(trim(backend)=='svd_compat')then
 call initialize_nr_cheb_svd(NR_FIT_LEGACY_COMPAT,code)
 elseif(trim(backend)=='svd_native')then
 call initialize_nr_cheb_svd(NR_FIT_NATIVE,code)
 else
 return
 endif
 if(code/=CHEB_SUCCESS)error stop 'NR SVD construction'
 open(newunit=file,file='svd_rows.csv',status='new')
 open(newunit=cfile,file='svd_coefficients.csv',status='new')
 write(file,'(a)')'species,row,policy,state,status,requested,effective,original,retained,rank,lower_zeros,upper_zeros,center,scale,original_min,original_max,retained_min,retained_max,residual,relative_residual'
 write(cfile,'(a)')'species,row,order,coefficient'
 do s=1,nspec
 do r=0,nperp
 call get_nr_cheb_row(s,r,model,info,code)
 write(file,'(*(g0,:,","))')s,r,info%policy,info%state,code,model%requested_degree,model%effective_degree,model%original_sample_count, &
 model%retained_sample_count,model%rank,model%lower_tail_zeros,model%upper_tail_zeros,model%center,model%scale, &
 info%original_min,info%original_max,model%retained_min,model%retained_max,model%residual_norm,model%relative_residual
 if(allocated(model%coefficients))then
 do j=0,model%effective_degree
 write(cfile,'(*(g0,:,","))')s,r,j,model%coefficients(j)
 enddo
 endif
 enddo
 enddo
 close(file);close(cfile)
end subroutine
subroutine nr_svd_activate
 use legacy_glls_reference,only:initialize_historical_nr,select_nr_test_evaluation
 use alps_analyt
 use alps_com,only:pass_nr_cheb_models
 implicit none
 character(32)::backend
 integer::code
 call get_environment_variable('ALPS_NR_SVD_BACKEND',backend)
 select case(trim(backend))
 case('shared')
 call initialize_historical_nr(code)
 if(code/=CHEB_SUCCESS)error stop 'NR shared construction'
 call select_nr_test_evaluation(.true.)
 case('svd_compat','svd_native')
 call pass_nr_cheb_models(code)
 if(code/=CHEB_SUCCESS)error stop 'NR SVD transfer'
 call select_nr_test_evaluation(.true.)
 end select
end subroutine
subroutine nr_svd_query(row,z,value,kind)
 use alps_var,only:iproc,sproc
 implicit none
 integer,intent(in)::row,kind
 double complex,intent(in)::z,value
 integer,save::file=0
 character(64)::name
 if(file==0)then
 write(name,'(a,i0,a)')'resonant_queries_rank',iproc,'.csv'
 open(newunit=file,file=trim(name),status='new')
 write(file,'(a)')'species,row,kind,u_re,u_im,value_re,value_im'
 endif
 write(file,'(*(g0,:,","))')sproc,row,kind,real(z),aimag(z),real(value),aimag(value)
end subroutine
