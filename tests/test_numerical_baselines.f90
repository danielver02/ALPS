! Calls production modules without changing the production main or algorithms.
program test_numerical_baselines
  use mpi
  use alps_var
  use alps_io, only : init_param,read_f0,alps_error_init
  use alps_com, only : pass_instructions,pass_distribution
  use alps_analyt, only : determine_param_fit,eval_fit
  use alps_fns, only : derivative_f0,determine_nmax,split_processes,determine_bessel_array,disp,rtsec,determine_resonances
  use relativistic_references, only : dp,anisotropic
  use ieee_arithmetic, only : ieee_is_finite
  implicit none
  integer :: out,iw,is,i,j,k,row,iflag,root_unit
  real(dp) :: dg,du,x,normalization,fit_error,interp_error,deriv_error,fref,peak,b,weight,inside_error
  complex(dp) :: omega,d,z,f,fg,fu,fminus,fplus,df,root,expected
  character(32) :: mode
  character(128) :: key
  call mpi_init(ierror)
  call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierror)
  call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)
  proc0=iproc==0
  sproc=0
  call get_command_argument(2,mode)
  if(nproc<=2.or.mod(nproc,2)/=0) call mpi_abort(MPI_COMM_WORLD,99,ierror)
  if(proc0) then
     call alps_error_init
     call init_param
  endif
  call pass_instructions
  if(proc0) then
     allocate(f0(nspec,0:nperp,0:npar)); f0=0
     call read_f0
     call derivative_f0
     call determine_param_fit
     open(newunit=out,file='observations.dat',status='new')
     do is=1,nspec
        call record_real('source_density',is,0,density_int(is))
        if(relativistic(is)) then
           b=0
           if(trim(mode)=='rel5') b=0.15_dp
           dg=gamma_rel(is,2,2)-gamma_rel(is,1,2)
           du=pparbar_rel(is,2,2)-pparbar_rel(is,2,1)
           weight=2*pi*dg*du*(ms(is)/vA)**3
           normalization=0
           do i=0,ngamma
              do j=0,npparbar
                 if(f0_rel(is,i,j)<0) cycle
                 normalization=normalization+gamma_rel(is,i,j)* &
                      exp(-2*gamma_rel(is,i,j)-b*pparbar_rel(is,i,j)**2)*weight
              enddo
           enddo
           fit_error=0; interp_error=0; deriv_error=0; inside_error=0; x=0
           peak=maxval(f0_rel(is,:,:))
           do i=0,ngamma
              do j=0,npparbar
                 if(f0_rel(is,i,j)<0) cycle
                 z=cmplx(pparbar_rel(is,i,j),0._dp,dp)
                 f=eval_fit(is,i,z)
                 fref=exp(-2*gamma_rel(is,i,j)-b*real(z,dp)**2)/normalization
                 interp_error=max(interp_error,abs(f0_rel(is,i,j)-fref)/max(fref,1.e-12_dp*peak))
                 if((ms(is)/vA)**2*(gamma_rel(is,i,j)**2-1-real(z,dp)**2)<=maxval(pp(is,:,:,1))**2) &
                      inside_error=max(inside_error,abs(f0_rel(is,i,j)-fref)/max(fref,1.e-12_dp*peak))
                 fit_error=max(fit_error,abs(f-f0_rel(is,i,j))/max(f0_rel(is,i,j),1.e-12_dp*peak))
                 x=x+gamma_rel(is,i,j)*f0_rel(is,i,j)*weight
                 if(i<1.or.i>ngamma-1) cycle
                 ! The outer parallel edges are intentionally zero placeholders in production.
                 if(j==0.or.j==npparbar) cycle
                 if(f0_rel(is,i-1,j)<=0.or.f0_rel(is,i+1,j)<=0) cycle
                 expected=-2*fref*sinh(2*dg)/(2*dg)
                 deriv_error=max(deriv_error,abs(df0_rel(is,i,j,1)-expected)/(2*peak))
              enddo
           enddo
           if(abs(x-1)>2.e-12_dp) error stop 'completed relativistic normalization'
           call record_real('completed_density',is,0,x)
           call record_real('interpolation_max_relative',is,0,interp_error)
           call record_real('interpolation_in_domain_relative',is,0,inside_error)
           call record_real('fit_max_relative',is,0,fit_error)
           call record_real('grid_dgamma_scaled_error',is,0,deriv_error)
           if(trim(mode)=='rel4') then
              if(interp_error>1.e-10_dp.or.fit_error>1.e-10_dp.or.deriv_error>1.e-10_dp) &
                   error stop 'isotropic affine-log interpolation/fit/stencil contract'
           endif
           row=max(2,ngamma/8)
           z=cmplx(0.2_dp,-0.03_dp,dp)
           f=eval_fit(is,row,z)
           fg=(eval_fit(is,row+1,z)-eval_fit(is,row-1,z))/(2*dg)
           fu=(eval_fit(is,row,z+du)-eval_fit(is,row,z-du))/(2*du)
           call record('fit_complex',is,row,f)
           call record('fit_dgamma_fixed_u',is,row,fg)
           call record('fit_du_fd',is,row,fu)
           call anisotropic(gamma_rel(is,row,1),z,2._dp,b,1/normalization,f,fg,fu)
           call record('analytic_complex',is,row,f)
           call record('analytic_dgamma',is,row,fg)
           call record('analytic_du',is,row,fu)
        else
           fit_error=0
           peak=maxval(f0(is,:,:))
           do i=0,nperp
              do j=0,npar
                 z=cmplx(pp(is,i,j,2),0._dp,dp)
                 f=eval_fit(is,i,z)
                 fit_error=max(fit_error,abs(f-f0(is,i,j))/max(f0(is,i,j),peak*1.e-12_dp))
              enddo
           enddo
           call record_real('fit_max_relative',is,0,fit_error)
           row=nperp/10
           du=pp(is,row,2,2)-pp(is,row,1,2)
           z=cmplx(0.2_dp,-0.03_dp,dp)*sqrt(ms(is))
           f=eval_fit(is,row,z,dfpar=df)
           fplus=eval_fit(is,row,z+du); fminus=eval_fit(is,row,z-du)
           call record('fit_complex',is,row,f)
           call record('fit_du_analytic',is,row,df)
           call record('fit_du_fd',is,row,(fplus-fminus)/(2*du))
           expected=exp(-pp(is,row,0,1)**2/ms(is)-z*z/ms(is))/(pi*ms(is))**1.5_dp
           call record('analytic_complex',is,row,expected)
           call record('analytic_du',is,row,-2*z/ms(is)*expected)
           if(abs(f-expected)>2.e-10_dp*abs(expected)) error stop 'Maxwellian polynomial reconstruction'
           if(abs(df+2*z/ms(is)*expected)>2.e-9_dp*max(abs(expected)/sqrt(ms(is)),abs(df))) &
                error stop 'Maxwellian analytic polynomial derivative'
        endif
     enddo
     if(any(relativistic)) call sparse_rows
     deallocate(f0)
  endif
  call pass_distribution
  call determine_nmax
  call split_processes
  if(sproc/=0) call determine_bessel_array
  do iw=1,3
     if(any(relativistic)) then
        omega=cmplx(0.062713_dp,-1.e-4_dp,dp)
        if(iw==2) omega=conjg(omega)
        if(iw==3) omega=cmplx(1._dp,-1.e-4_dp,dp)
     else
        omega=cmplx(0.03_dp,-1.e-5_dp,dp)
        if(iw==2) omega=conjg(omega)
        if(iw==3) omega=cmplx(0.045_dp,-0.019_dp,dp)
     endif
     d=disp(omega)
     if(proc0) then
        call record('omega',0,iw,omega)
        call record('determinant',0,iw,d,1.e-24_dp,1.e-10_dp)
        do is=1,nspec
           do i=1,3
              do j=i,3
                 write(key,'(a,i0,i0)') 'chi',i,j
                 call record(trim(key),is,iw,chi0(is,i,j))
              enddo
           enddo
        enddo
     endif
  enddo
  ! Secondary baseline only; convergence is checked independently of the solver's exit message.
  root=wroots(1); iflag=0
  root=rtsec(disp,root,1,iflag)
  d=disp(root)
  if(proc0) then
     x=abs(d)/max(maxval(abs(wave))**3,tiny(1._dp))
     open(newunit=root_unit,file='root_status.dat',status='new')
     write(root_unit,'(a,i0,a,es24.15e3)') 'iterations=',iflag,' scaled_residual=',x
     write(root_unit,'(a,2es26.17e3)') 'root=',root
     if(iflag>=numiter-1.or.x>1.e-7_dp) then
        write(root_unit,'(a)') 'UNRESOLVED: root is not accepted as a numerical reference'
     else
        write(root_unit,'(a)') 'ALGEBRAICALLY_CONVERGED: branch identity and grid convergence require separate qualification'
     endif
     close(root_unit)
     close(out); close(unit_error)
     print '(a)', 'PASS: finite numerical observations captured from freshly built production modules'
  endif
  call mpi_finalize(ierror)
contains
  subroutine record(label,species,index,value,abs_tol,rel_tol)
    character(*), intent(in) :: label
    integer, intent(in) :: species,index
    complex(dp), intent(in) :: value
    real(dp), intent(in), optional :: abs_tol,rel_tol
    real(dp) :: at,rt
    character(128) :: name
    at=1.e-12_dp; rt=1.e-10_dp
    if(present(abs_tol)) at=abs_tol
    if(present(rel_tol)) rt=rel_tol
    if(.not.ieee_is_finite(real(value,dp)).or..not.ieee_is_finite(aimag(value))) then
       print *, 'nonfinite observation: ',label,species,index
       call mpi_abort(MPI_COMM_WORLD,98,ierror)
    endif
    write(name,'(a,a,i0,a,i0)') label,'_s',species,'_i',index
    write(out,'(a,4es26.17e3)') trim(name)//' ',real(value,dp),aimag(value),at,rt
  end subroutine
  subroutine record_real(label,species,index,value)
    character(*), intent(in) :: label
    integer, intent(in) :: species,index
    real(dp), intent(in) :: value
    call record(label,species,index,cmplx(value,0._dp,dp))
  end subroutine
  subroutine sparse_rows
    integer :: file,summary,gidx,pidx,n,lo,hi,count,zero,one,few,used,cross,neighbor
    real(dp) :: ur,lower,upper
    logical :: active,needs_outside,gate(-1:1),unused_gate
    sproc=1
    call determine_resonances(cmplx(0.062713_dp,-0.0001_dp,dp),0,gate(0),unused_gate)
    call determine_resonances(cmplx(0.062713_dp,-0.0001_dp,dp),1,gate(1),gate(-1))
    sproc=0
    open(newunit=file,file='sparse_rows.dat',status='new')
    write(file,'(a)') '# row Gamma count lower_u upper_u active_residue neighbor_extrapolation; species 1'
    zero=0; one=0; few=0; used=0; cross=0
    do gidx=0,ngamma
       lo=-1; hi=-1; count=0
       do pidx=0,npparbar
          if(f0_rel(1,gidx,pidx)<0) cycle
          count=count+1
          if(lo==-1) lo=pidx
          hi=pidx
       enddo
       lower=0; upper=0
       if(count>0) then
          lower=pparbar_rel(1,gidx,lo); upper=pparbar_rel(1,gidx,hi)
       endif
       if(count==0) zero=zero+1
       if(count==1) one=one+1
       if(count<31) few=few+1
       active=.false.; needs_outside=.false.
       if(gidx>0.and.gidx<ngamma) then
       do n=-1,1
          if(.not.gate(n)) cycle
          ur=(gamma_rel(1,gidx,0)*0.062713_dp-real(n,dp)*qs(1)/ms(1))*vA/kpar
          if(ur**2>gamma_rel(1,gidx,0)**2-1) cycle
          active=.true.
          do neighbor=max(1,gidx-1),gidx+1
             if(neighbor==gidx) cycle
             lo=-1; hi=-1
             do pidx=0,npparbar
                if(f0_rel(1,neighbor,pidx)<0) cycle
                if(lo==-1) lo=pidx
                hi=pidx
             enddo
             if(lo==-1) then
                needs_outside=.true.
             else if(ur<pparbar_rel(1,neighbor,lo).or.ur>pparbar_rel(1,neighbor,hi)) then
                needs_outside=.true.
             endif
          enddo
       enddo
       endif
       write(file,'(i6,es25.16e3,i8,2es25.16e3,2l3)') gidx,gamma_rel(1,gidx,0),count,lower,upper,active,needs_outside
       if(active) used=used+1
       if(needs_outside) cross=cross+1
    enddo
    close(file)
    open(newunit=summary,file='sparse_summary.dat',status='new')
    write(summary,'(a,2i8)') 'ngamma,npparbar ',ngamma,npparbar
    write(summary,'(a,3i8)') 'zero,singleton,fewer_than_31 ',zero,one,few
    write(summary,'(a,2i8)') 'active_residue_rows,neighbor_interval_extrapolation_rows ',used,cross
    write(summary,'(a,3l3)') 'production found_res for n=-1,0,1 ',gate
    write(summary,'(a)') 'Rows 1:ngamma-1 are visited by landau_integrate_rel; row 0 is not used by its gamma stencil.'
    write(summary,'(a)') 'Active counts use omega=0.062713-0.0001i, n=-1,0,1 and the production real-projection cone gate.'
    write(summary,'(a)') 'Sample-interval endpoints differ from exact cone endpoints; complex u is never real grid data.'
    close(summary)
  end subroutine
end program
