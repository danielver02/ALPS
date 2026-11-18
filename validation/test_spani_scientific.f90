program test_spani
 use mpi
 use alps_var
 use alps_io,only:init_param,read_f0,alps_error_init
 use alps_com,only:pass_instructions,pass_distribution,pass_nr_cheb_models,nr_cheb_distribution_info,complete_nr_cheb_initialization
 use alps_analyt
 use alps_fns,only:derivative_f0,determine_nmax,split_processes,determine_bessel_array,disp,span_serial_disp,rtsec,om_scan
 use spani_diagnostics
 use, intrinsic::iso_fortran_env,only:int64
 use legacy_glls_reference
 implicit none
 character(32)::arg,phase
 integer::choice,degree,selector,iw,code,unit,s,r,nb,nf,ns0,rows,evals,setup,tx,nc,oldb,oldf,oldt,oldc,flag,ios,pscode,mscode,vscode,dscode
 double precision::probe,ur,ui,xn
 double complex::omega,det,root,value,df,fd,plus,minus,poly,pd
 integer,allocatable::check_meta(:,:)
 double precision,allocatable::check_maps(:,:),check_coefs(:,:)
 integer(int64),allocatable::bits(:),bmin(:),bmax(:)
 integer::at,nbits
 type(cheb_model)::model
 type(nr_cheb_row_info)::info
 call mpi_init(ierror);call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror);call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)
 proc0=iproc==0;sproc=0
 call get_command_argument(2,arg);read(arg,*)choice
 call get_command_argument(3,arg);read(arg,*)degree
 call get_command_argument(4,arg);read(arg,*)selector
 call get_command_argument(5,phase)
 call get_command_argument(6,arg);read(arg,*)iw
 call get_command_argument(7,arg);read(arg,*)probe
 if(proc0)then
 call alps_error_init;call init_param
 endif
 call pass_instructions
 if(nproc>1.and.nproc<nspec+1)then
 if(proc0)print *,'UNSUPPORTED: mixed-species production assembly needs at least nspec+1 ranks'
 call mpi_abort(MPI_COMM_WORLD,90,ierror)
 endif
 poly_order(1)=degree;poly_derivative_type(1)=selector;kpar=probe
 if(allocated(historical_coefficients))deallocate(historical_coefficients)
 if(choice<3)then
 allocate(historical_coefficients(nspec,0:nperp,0:degree));historical_coefficients=0
 endif
 if(proc0)then
 allocate(f0(nspec,0:nperp,0:npar));call read_f0;call derivative_f0
 if(choice==4)then
 call initialize_nr_cheb_production(code)
 if(code/=0)error stop 'production initialization'
 call determine_param_fit
 else
 if(choice<2)then
 call legacy_determine_param_fit
 write(*,*)'SPAN_GLLS_FIT_ONCE',iproc,nperp+1
 endif
 call prepare_models(choice)
 endif
 if(choice/=4)call select_nr_test_evaluation(choice/=0)
 if(phase=='fit')call export_fits
 if(phase=='fit'.and.choice==3)then
 call initialize_nr_cheb_svd(NR_FIT_LEGACY_COMPAT,code)
 open(newunit=unit,file='strict_compatibility.csv',status='new')
 do r=0,nperp
 call get_nr_cheb_row(1,r,model,info,code);write(unit,'(*(g0,:,","))')r,code
 enddo
 close(unit)
 call prepare_models(choice)
 endif
 deallocate(f0)
 endif
 history_transfer_enabled=choice<3
 call pass_distribution
 if(choice==4)then
 call complete_nr_cheb_initialization(code)
 if(code/=0)error stop 'production transfer'
 elseif(choice==1.or.choice==0)then
 if(.not.proc0)call initialize_historical_nr(code)
 else
 call pass_nr_cheb_models(code)
 endif
 if(choice/=4)call select_nr_test_evaluation(choice/=0)
 call nr_cheb_export(check_meta,check_maps,check_coefs,code)
 nbits=size(check_meta)+size(check_maps)+size(check_coefs)
 allocate(bits(nbits),bmin(nbits),bmax(nbits));at=size(check_meta)
 bits(:at)=int(reshape(check_meta,[at]),int64)
 bits(at+1:at+size(check_maps))=transfer(check_maps,[0_int64],size(check_maps));at=at+size(check_maps)
 bits(at+1:)=transfer(check_coefs,[0_int64],size(check_coefs))
 call mpi_allreduce(bits,bmin,nbits,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,ierror)
 call mpi_allreduce(bits,bmax,nbits,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,ierror)
 if(any(bmin/=bmax))error stop 'SPAN-I model transfer mismatch'
 if(proc0)print *,'SPAN_MODELS_BIT_IDENTICAL',nbits
 deallocate(bits,bmin,bmax,check_meta,check_maps,check_coefs)
 if(phase=='replay')then
 open(newunit=unit,file='query_union.csv',status='old')
 open(74,file='replay.csv',status='new')
 write(74,'(a)')'species,row,u_re,u_im,x_re,x_im,outside_real,status,value_status,derivative_status,plus_status,minus_status,value_re,value_im,derivative_re,derivative_im,FD_re,FD_im,rawP_re,rawP_im'
 do
 read(unit,*,iostat=ios)s,r,ur,ui
 if(ios/=0)exit
 omega=cmplx(ur,ui,kind(1.d0));call get_nr_cheb_row(s,r,model,info,code)
 value=eval_fit(s,r,omega,dfpar=df,status=code,value_status=vscode,derivative_status=dscode)
 xn=abs(pp(s,r,2,2)-pp(s,r,1,2))
 plus=eval_fit(s,r,omega+xn,status=pscode);minus=eval_fit(s,r,omega-xn,status=mscode);fd=(plus-minus)/(2*xn)
 call cheb_evaluate_polynomial(model,omega,poly,ios,pd)
 write(74,'(*(g0,:,","))')s,r,ur,ui,(ur-model%center)/model%scale,ui/model%scale, &
 abs(ur-model%center)>model%scale,code,vscode,dscode,pscode,mscode,real(value),aimag(value),real(df),aimag(df),real(fd),aimag(fd),real(poly),aimag(poly)
 enddo
 close(unit);close(74)
 endif
 if(phase=='fit'.or.phase=='replay')then
 call mpi_finalize(ierror);stop
 endif
 call determine_nmax
 if(nproc==1)then
 sproc=1;nlim=[0,nmax(1)];call determine_bessel_array;sproc=0
 else
 call split_processes
 if(sproc/=0)call determine_bessel_array
 endif
 capture=.true.
 omega=cmplx(-.067722d0,-1.2364d-7,kind(1.d0))
 if(iw==4)omega=wroots(1) ! Generic existing-fixture seed for migration controls.
 if(iw==2)omega=cmplx(real(omega),-1.d-4,kind(1.d0))
 if(iw==3)omega=cmplx(real(omega),-1.d-3,kind(1.d0))
 call nr_cheb_collection_info(ns0,rows,oldb,oldf,evals,setup);call nr_cheb_distribution_info(oldt,oldc)
 det=frequency(omega);call record_result(omega,det,'preflight')
 if(phase=='root'.or.phase=='scan')then
 iterations=0;flag=0;root=rtsec(frequency,omega,1,flag)
 det=frequency(root);call record_result(root,det,'root')
 if(phase=='scan')then
 if(nproc==1)error stop 'scan uses production MPI assembly'
 wroots(1)=root;nroots=1
 if(allocated(scan))deallocate(scan)
 allocate(scan(1));scan(1)%type_s=4;scan(1)%n_out=6;scan(1)%n_res=1
 scan(1)%log_scan=.true.;scan(1)%diff=(log10(.09d0)-log10(.075d0))/6
 scan(1)%heat_s=.false.;scan(1)%eigen_s=.false.
 call om_scan(1)
 endif
 endif
 call nr_cheb_collection_info(ns0,rows,nb,nf,evals,setup);call nr_cheb_distribution_info(tx,nc)
 if(nb/=oldb.or.nf/=oldf.or.tx/=oldt.or.nc/=oldc)error stop 'initialization repeated during evaluation'
 write(*,*)'COUNTERS',iproc,nb,nf,tx,nc,evals
 call mpi_finalize(ierror)
contains
 function frequency(w)result(d)
 double complex,intent(in)::w
 double complex::d
 if(nproc==1)then
 d=span_serial_disp(w)
 else
 d=disp(w)
 endif
 end function
end program
