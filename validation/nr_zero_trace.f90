module nr_zero_diagnostics
 implicit none
 integer::case_id=-1,selector=-1,file=0
 logical::b5_recording=.false.,b5_inject=.false.
 integer::b5_row_current=-1,b5_degree=0,b5_selector=0,b5_counts(3,2,2)=0,b5_file=0
 character(16)::b5_variant=""
contains
 subroutine b5_begin(degree,selected,variant,inject)
 use alps_var,only:iproc
 integer,intent(in)::degree,selected
 character(*),intent(in)::variant
 logical,intent(in)::inject
 character(80)::name
 b5_recording=.true.;b5_inject=inject;b5_degree=degree;b5_selector=selected;b5_variant=variant;b5_counts=0
 if(b5_file==0)then
 write(name,'(a,i0,a)')'query_trace_rank',iproc,'.csv'
 open(newunit=b5_file,file=trim(name),status='new')
 write(b5_file,'(a)')'degree,selector,variant,contour_row,query_row,kind,endpoint,derivative,u_re,u_im,injected_status'
 endif
 end subroutine
 subroutine context(c,s)
 use alps_var,only:iproc
 integer,intent(in)::c,s
 character(64)::name
 case_id=c;selector=s
 if(file==0)then
 write(name,'(a,i0,a)')'integrands_rank',iproc,'.csv'
 open(newunit=file,file=trim(name),status='new')
 write(file,'(a)')'case,selector,variant,row,u_re,u_im,dfpar_re,dfpar_im,dfperp_re,dfperp_im,accumulator_re,accumulator_im'
 endif
 end subroutine
end module
subroutine zero_trace(variant,row,u,dfpar,dfperp,accumulator)
 use nr_zero_diagnostics
 implicit none
 character(*),intent(in)::variant
 integer,intent(in)::row
 double complex,intent(in)::u,dfpar,dfperp,accumulator
 write(file,'(*(g0,:,","))')case_id,selector,variant,row,real(u),aimag(u),real(dfpar),aimag(dfpar), &
 real(dfperp),aimag(dfperp),real(accumulator),aimag(accumulator)
end subroutine
subroutine zero_gate(variant,row,a,b,c,d)
 use nr_zero_diagnostics
 implicit none
 character(*),intent(in)::variant
 integer,intent(in)::row
 double complex,intent(in)::a,b,c,d
 write(*,*)'ZERO_GATE',case_id,selector,variant,row,a,b,c,d
end subroutine

subroutine b5_row(row)
 use nr_zero_diagnostics,only:b5_row_current
 implicit none
 integer,intent(in)::row
 b5_row_current=row
end subroutine
subroutine b5_query(species,row,u,derivative,inject)
 use nr_zero_diagnostics
 use alps_var,only:nperp
 implicit none
 integer,intent(in)::species,row
 double complex,intent(in)::u
 logical,intent(in)::derivative
 logical,intent(out)::inject
 integer::category,endpoint,df,code
 inject=.false.
 if(.not.b5_recording)return
 category=1
 if(row/=b5_row_current)category=2
 if(real(u)/=0.d0)category=3
 endpoint=1
 if(b5_row_current==0.or.b5_row_current==nperp)endpoint=2
 df=1
 if(derivative)df=2
 b5_counts(category,endpoint,df)=b5_counts(category,endpoint,df)+1
 ! Independently designated query: first interior row, positive parallel shift.
 inject=b5_inject.and.species==1.and.row==1.and.u==cmplx(.5d0,-.125d0,kind(1.d0))
 code=0
 if(inject)code=18
 write(b5_file,'(*(g0,:,","))')b5_degree,b5_selector,trim(b5_variant),b5_row_current,row,category,endpoint-1,derivative,real(u),aimag(u),code
 if(inject)flush(b5_file)
end subroutine
