! Read-only diagnostics linked only into a copied production main.
subroutine r1_dump
 use alps_var
 use alps_analyt,only:get_rel_cheb_row,cheb_model,rel_cheb_row_info,cheb_evaluate
 use, intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 type(cheb_model)::model
 type(rel_cheb_row_info)::info
 integer::s,packed,g,j,file,rows,bounds,code
 double precision::radius,lo,hi,umin,umax,perp0,perp1,up,rr
 double complex::fitted,dfitted
 integer::evalcode
 logical::inside
 open(newunit=file,file='r1_nodes.csv',status='replace')
 open(newunit=rows,file='r1_rows.csv',status='replace')
 open(newunit=bounds,file='r1_bounds.csv',status='replace')
 write(file,'(a)')'species,packed,g,j,gamma,u,r,pperp,ppar,inside_source,table_value,fit_status,fit_value,fit_derivative'
 write(rows,'(a)')'species,packed,g,state,status,requested,effective,rank,representation,geometric,positive,center,scale'
 write(bounds,'(a)')'species,packed,mass,va,pperp_min,pperp_max,ppar_min,ppar_max,gamma_min,gamma_max,ngamma,npparbar'
 packed=0
 do s=1,nspec
 if(.not.relativistic(s))cycle
 packed=packed+1
 if(useBM(s))cycle
 perp0=minval(pp(s,:,:,1));perp1=maxval(pp(s,:,:,1));umin=minval(pp(s,:,:,2));umax=maxval(pp(s,:,:,2))
 write(bounds,'(*(g0,:,","))')s,packed,ms(s),vA,perp0,perp1,umin,umax,minval(gamma_rel(packed,:,:)),maxval(gamma_rel(packed,:,:)),ngamma,npparbar
 do g=0,ngamma
 call get_rel_cheb_row(packed,g,model,info,code)
 write(rows,'(*(g0,:,","))')s,packed,g,info%state,code,model%requested_degree,model%effective_degree,model%rank,model%representation,info%geometric_samples,info%positive_samples,model%center,model%scale
 do j=0,npparbar
 if(f0_rel(packed,g,j)==-1.d0)cycle
 up=pparbar_rel(packed,g,j);radius=sqrt(max(0.d0,gamma_rel(packed,g,j)**2-1.d0-up*up))
 rr=radius*ms(s)/vA;lo=up*ms(s)/vA
 inside=rr>=perp0-1.d-12.and.rr<=perp1+1.d-12.and.lo>=umin-1.d-12.and.lo<=umax+1.d-12
 fitted=cmplx(ieee_value(0.d0,ieee_quiet_nan),0.d0,kind(1.d0));dfitted=fitted;evalcode=code
 if(info%state==1.or.info%state==2)call cheb_evaluate(model,cmplx(up,0.d0,kind(1.d0)),fitted,evalcode,dfitted)
 write(file,'(*(g0,:,","))')s,packed,g,j,gamma_rel(packed,g,j),up,radius,rr,lo,inside,f0_rel(packed,g,j),evalcode,real(fitted),real(dfitted)
 enddo
 enddo
 enddo
 close(file);close(rows);close(bounds)
end subroutine
subroutine r1_stats(phase)
 use alps_var,only:iproc
 use alps_analyt,only:rel_cheb_collection_info
 use alps_com,only:rel_cheb_distribution_info
 use iso_fortran_env,only:int64
 implicit none
 character(*),intent(in)::phase
 integer::ns,ng,fits,builds,tx,collect
 integer(int64)::bytes,evals,wire
 call rel_cheb_collection_info(ns,ng,fits,bytes,builds,evals)
 call rel_cheb_distribution_info(tx,collect,wire)
 write(*,'(3a,8i8)')'R1_STATS ',phase,' ',iproc,ns,ng,builds,fits,tx,collect,evals
 if(phase=='initial'.and.ns>0)call r2_export_initial
end subroutine

subroutine r2_tensor(om,det)
 use alps_var,only:chi0,wave,nspec
 implicit none
 double complex,intent(in)::om,det
 integer::s,i,j,unit
 logical,save::first=.true.
 open(newunit=unit,file='r2_tensor.csv',position='append')
 if(first)write(unit,'(a)')'omega_re,omega_im,species,i,j,re,im'
 first=.false.
 do s=1,nspec+2
 do i=1,3
 do j=1,3
 if(s<=nspec)then
 write(unit,'(*(g0,:,","))')real(om),aimag(om),s,i,j,real(chi0(s,i,j)),aimag(chi0(s,i,j))
 elseif(s==nspec+1)then
 write(unit,'(*(g0,:,","))')real(om),aimag(om),0,i,j,real(sum(chi0(:,i,j))),aimag(sum(chi0(:,i,j)))
 else
 write(unit,'(*(g0,:,","))')real(om),aimag(om),-1,i,j,real(wave(i,j)),aimag(wave(i,j))
 endif
 enddo
 enddo
 enddo
 write(unit,'(*(g0,:,","))')real(om),aimag(om),-2,0,0,real(det),aimag(det)
 close(unit)
end subroutine

subroutine r2_root_record(x,f,iteration)
 use alps_var,only:wave,proc0
 implicit none
 double complex,intent(in)::x,f
 integer,intent(in)::iteration
 double complex::cross(3),best(3)
 double precision::norm,bestnorm,ez,ey,helicity,eta
 integer::i,j,unit
 if(.not.proc0)return
 bestnorm=0.d0;best=(0.d0,0.d0)
 do i=1,2
 do j=i+1,3
 cross=[wave(i,2)*wave(j,3)-wave(i,3)*wave(j,2), &
 wave(i,3)*wave(j,1)-wave(i,1)*wave(j,3),wave(i,1)*wave(j,2)-wave(i,2)*wave(j,1)]
 norm=sum(abs(cross)**2)
 if(norm>bestnorm)then
 best=cross;bestnorm=norm
 endif
 enddo
 enddo
 ez=0.d0;ey=0.d0;helicity=0.d0
 if(bestnorm>0.d0)then
 ez=abs(best(3))**2/bestnorm;ey=abs(best(2))**2/bestnorm
 endif
 norm=abs(best(1))**2+abs(best(2))**2
 if(norm>0.d0)helicity=2.d0*aimag(conjg(best(1))*best(2))/norm
 eta=abs(f)/max(maxval(abs(wave))**3,tiny(1.d0))
 open(newunit=unit,file='r2_root.csv',status='replace')
 write(unit,'(a)')'omega_re,omega_im,iterations,det_re,det_im,scaled_det,ez_fraction,ey_fraction,helicity'
 write(unit,'(*(g0,:,","))')real(x),aimag(x),iteration,real(f),aimag(f),eta,ez,ey,helicity
 close(unit)
end subroutine

subroutine r2_preflight
 use alps_fns,only:disp
 implicit none
 integer::unit,ios
 double precision::re,im
 double complex::det
 logical::exists
 inquire(file='preflight.txt',exist=exists)
 if(.not.exists)return
 open(newunit=unit,file='preflight.txt',status='old')
 read(unit,*,iostat=ios)re,im
 close(unit)
 if(ios/=0)error stop 'malformed production comparison frequency'
 det=disp(cmplx(re,im,kind(1.d0)))
end subroutine

subroutine r2_export_initial
 use alps_var,only:iproc
 use alps_analyt,only:rel_cheb_export,CHEB_SUCCESS
 implicit none
 integer,allocatable::metadata(:,:)
 double precision,allocatable::mapping(:,:),coefficients(:,:)
 integer::code,unit
 character(80)::name
 call rel_cheb_export(metadata,mapping,coefficients,code)
 if(code/=CHEB_SUCCESS)error stop 'diagnostic initial-model export failed'
 write(name,'(a,i0,a)')'models_rank',iproc,'.bin'
 open(newunit=unit,file=trim(name),status='replace',access='stream',form='unformatted')
 write(unit)metadata,mapping,coefficients
 close(unit)
end subroutine
