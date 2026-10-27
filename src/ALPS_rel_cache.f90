! Native relativistic cache support. No interpolation, MPI, or fitting dependencies.
module alps_cache_sha256
  !! Incremental SHA-256, using unsigned 32-bit values held in int64 containers.
  !! All additions fit int64; masking implements modulo 2**32 without signed overflow.
  use iso_fortran_env, only : int64
  implicit none
  private
  integer(int64), parameter :: mask32=int(z'FFFFFFFF',int64)
  integer(int64), parameter :: constants(64)=[ &
       int(z'428A2F98',int64),int(z'71374491',int64),int(z'B5C0FBCF',int64),int(z'E9B5DBA5',int64), &
       int(z'3956C25B',int64),int(z'59F111F1',int64),int(z'923F82A4',int64),int(z'AB1C5ED5',int64), &
       int(z'D807AA98',int64),int(z'12835B01',int64),int(z'243185BE',int64),int(z'550C7DC3',int64), &
       int(z'72BE5D74',int64),int(z'80DEB1FE',int64),int(z'9BDC06A7',int64),int(z'C19BF174',int64), &
       int(z'E49B69C1',int64),int(z'EFBE4786',int64),int(z'0FC19DC6',int64),int(z'240CA1CC',int64), &
       int(z'2DE92C6F',int64),int(z'4A7484AA',int64),int(z'5CB0A9DC',int64),int(z'76F988DA',int64), &
       int(z'983E5152',int64),int(z'A831C66D',int64),int(z'B00327C8',int64),int(z'BF597FC7',int64), &
       int(z'C6E00BF3',int64),int(z'D5A79147',int64),int(z'06CA6351',int64),int(z'14292967',int64), &
       int(z'27B70A85',int64),int(z'2E1B2138',int64),int(z'4D2C6DFC',int64),int(z'53380D13',int64), &
       int(z'650A7354',int64),int(z'766A0ABB',int64),int(z'81C2C92E',int64),int(z'92722C85',int64), &
       int(z'A2BFE8A1',int64),int(z'A81A664B',int64),int(z'C24B8B70',int64),int(z'C76C51A3',int64), &
       int(z'D192E819',int64),int(z'D6990624',int64),int(z'F40E3585',int64),int(z'106AA070',int64), &
       int(z'19A4C116',int64),int(z'1E376C08',int64),int(z'2748774C',int64),int(z'34B0BCB5',int64), &
       int(z'391C0CB3',int64),int(z'4ED8AA4A',int64),int(z'5B9CCA4F',int64),int(z'682E6FF3',int64), &
       int(z'748F82EE',int64),int(z'78A5636F',int64),int(z'84C87814',int64),int(z'8CC70208',int64), &
       int(z'90BEFFFA',int64),int(z'A4506CEB',int64),int(z'BEF9A3F7',int64),int(z'C67178F2',int64)]
  type, public :: sha256_state
     private
     integer(int64) :: h(8),bytes=0
     integer :: used=0
     character(64) :: buffer
   contains
     procedure :: init=>sha_init
     procedure :: update=>sha_update
     procedure :: finish=>sha_finish
  end type
  public :: sha256,sha256_file
contains
  subroutine sha_init(self)
    class(sha256_state), intent(out) :: self
    self%h=[int(z'6A09E667',int64),int(z'BB67AE85',int64),int(z'3C6EF372',int64),int(z'A54FF53A',int64), &
         int(z'510E527F',int64),int(z'9B05688C',int64),int(z'1F83D9AB',int64),int(z'5BE0CD19',int64)]
    self%bytes=0; self%used=0; self%buffer=repeat(achar(0),64)
  end subroutine
  pure integer(int64) function rotate(x,n) result(r)
    integer(int64), intent(in) :: x
    integer, intent(in) :: n
    r=iand(ior(shiftr(x,n),shiftl(x,32-n)),mask32)
  end function
  subroutine compress(self,data)
    class(sha256_state), intent(inout) :: self
    character(64), intent(in) :: data
    integer(int64) :: w(64),v(8),a,b,t1,t2
    integer :: i,j
    do i=1,16
       w(i)=0
       do j=1,4
          w(i)=ior(shiftl(w(i),8),int(iachar(data(4*(i-1)+j:4*(i-1)+j)),int64))
       enddo
    enddo
    do i=17,64
       a=ieor(ieor(rotate(w(i-15),7),rotate(w(i-15),18)),shiftr(w(i-15),3))
       b=ieor(ieor(rotate(w(i-2),17),rotate(w(i-2),19)),shiftr(w(i-2),10))
       w(i)=iand(w(i-16)+a+w(i-7)+b,mask32)
    enddo
    v=self%h
    do i=1,64
       a=ieor(ieor(rotate(v(5),6),rotate(v(5),11)),rotate(v(5),25))
       b=ieor(iand(v(5),v(6)),iand(ieor(v(5),mask32),v(7)))
       t1=iand(v(8)+a+b+constants(i)+w(i),mask32)
       a=ieor(ieor(rotate(v(1),2),rotate(v(1),13)),rotate(v(1),22))
       b=ieor(ieor(iand(v(1),v(2)),iand(v(1),v(3))),iand(v(2),v(3)))
       t2=iand(a+b,mask32)
       v=[iand(t1+t2,mask32),v(1:3),iand(v(4)+t1,mask32),v(5:7)]
    enddo
    self%h=iand(self%h+v,mask32)
  end subroutine
  subroutine sha_update(self,data)
    class(sha256_state), intent(inout) :: self
    character(*), intent(in) :: data
    integer :: p,n,take
    n=len(data); p=1; self%bytes=self%bytes+int(n,int64)
    if(self%used>0) then
       take=min(64-self%used,n)
       self%buffer(self%used+1:self%used+take)=data(1:take)
       self%used=self%used+take; p=p+take
       if(self%used==64) then
          call compress(self,self%buffer); self%used=0
       endif
    endif
    do while(p+63<=n)
       call compress(self,data(p:p+63)); p=p+64
    enddo
    if(p<=n) then
       self%buffer(1:n-p+1)=data(p:n); self%used=n-p+1
    endif
  end subroutine
  function sha_finish(self) result(hex)
    class(sha256_state), intent(inout) :: self
    character(64) :: hex
    character(8) :: length_bytes
    integer(int64) :: bits
    integer :: i,padding
    bits=self%bytes*8_int64
    padding=modulo(55-self%used,64)
    call self%update(achar(128)//repeat(achar(0),padding))
    do i=1,8
       length_bytes(i:i)=achar(int(iand(shiftr(bits,8*(8-i)),255_int64)))
    enddo
    call self%update(length_bytes)
    write(hex,'(8(z8.8))') self%h
    do i=1,64
       if(hex(i:i)>='A'.and.hex(i:i)<='F') hex(i:i)=achar(iachar(hex(i:i))+32)
    enddo
  end function
  function sha256(data) result(hex)
    character(*), intent(in) :: data
    character(64) :: hex
    type(sha256_state) :: state
    call state%init(); call state%update(data); hex=state%finish()
  end function
  subroutine sha256_file(path,hex,ok)
    character(*), intent(in) :: path
    character(64), intent(out) :: hex
    logical, intent(out) :: ok
    type(sha256_state) :: state
    character(65536) :: buffer
    integer :: unit,ios,cs,n
    integer(int64) :: length,p
    ok=.false.; hex=''
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=ios)
    if(ios/=0) return
    inquire(unit=unit,size=length,iostat=ios)
    call state%init(); p=0
    do while(ios==0.and.p<length)
       n=int(min(65536_int64,length-p))
       read(unit,iostat=ios) buffer(:n)
       if(ios/=0) exit
       call state%update(buffer(:n)); p=p+n
    enddo
    close(unit,iostat=cs)
    if(ios/=0.or.cs/=0) return
    hex=state%finish(); ok=.true.
  end subroutine
end module

module alps_rel_cache
  !! Rank-0 caller only. Cache errors are recoverable; live arrays change only on a valid hit.
  !! Native IEEE binary64 stream format; Linux build discovery and POSIX atomic publication.
  use iso_fortran_env, only : int32,int64,real64,compiler_version,compiler_options
  use iso_c_binding, only : c_int,c_char,c_null_char
  use, intrinsic :: ieee_arithmetic
  use alps_cache_sha256
  implicit none
  private
  integer(int64), parameter :: format_version=1,algorithm_revision=1
  integer, parameter :: header_size=312,trailer_size=80
  character(16), parameter :: magic='ALPS_REL_CACHE',complete='ALPS_REL_DONE'
  logical, save :: build_attempted=.false.,build_ok=.false.
  character(64), save :: build_hash=''
  character(256), save :: build_reason=''
  type, public :: rel_cache_identity
     !! Validated construction identity produced by [[cache_identity]], never read from a filename.
     character(64) :: key='',build=''
     !! Lowercase SHA-256 input key and strict executable/library/environment compatibility digest.
     integer :: dims(4)=0
     !! Source nperp/npar and target ngamma/npparbar interval counts, excluding species index.
     real(real64) :: settings(4)=0 ! mass, vA, pi, smoothing
     !! Exact construction parameters; smoothing is currently zero.
  end type
  public :: cache_identity,cache_load,cache_save
  interface
     function posix_mkstemp(template) bind(C,name='mkstemp') result(fd)
       import c_int,c_char
       character(c_char), intent(inout) :: template(*)
       integer(c_int) :: fd
     end function
     function posix_close(fd) bind(C,name='close') result(status)
       import c_int
       integer(c_int), value :: fd
       integer(c_int) :: status
     end function
     function posix_rename(old,new) bind(C,name='rename') result(status)
       import c_int,c_char
       character(c_char), intent(in) :: old(*),new(*)
       integer(c_int) :: status
     end function
     function posix_unlink(path) bind(C,name='unlink') result(status)
       import c_int,c_char
       character(c_char), intent(in) :: path(*)
       integer(c_int) :: status
     end function
  end interface
contains
  subroutine discover_build
    !! Fail closed if executable/library/CPU identity cannot be established.
    type(sha256_state) :: state
    type(ieee_round_type) :: rounding
    character(4096) :: line,path,paths(512),env
    character(64) :: hashes(513),swap
    character(32), parameter :: names(*)=[character(32) :: 'OPENBLAS_NUM_THREADS','OMP_NUM_THREADS', &
         'MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_DYNAMIC','MKL_DYNAMIC', &
         'OPENBLAS_CORETYPE','GOTO_NUM_THREADS','OMP_PROC_BIND','OMP_PLACES','MKL_CBWR','OMP_SCHEDULE']
    character(4) :: elf
    integer :: unit,ios,i,j,p,n,k,fu,cs,fields,length
    logical :: ok,known_backend,flags_found
    if(build_attempted) return
    build_attempted=.true.; build_reason='cannot establish Linux executable/library identity'
    if(storage_size(0.d0)/=64.or.storage_size(0_real64)/=64.or.radix(0.d0)/=2.or.digits(0.d0)/=53) return
    if(storage_size(0_int64)/=64.or.storage_size(0_int32)/=32.or..not.ieee_support_datatype(0.d0)) return
    call sha256_file('/proc/self/exe',hashes(1),ok)
    if(.not.ok) return
    n=0; known_backend=.false.
    open(newunit=unit,file='/proc/self/maps',status='old',action='read',iostat=ios)
    if(ios/=0) return
    do
       read(unit,'(a)',iostat=ios) line
       if(ios/=0) exit
       if(len_trim(line)==len(line)) exit
       p=1
       do fields=1,5
          do while(p<=len_trim(line))
             if(line(p:p)/=' ') exit
             p=p+1
          enddo
          do while(p<=len_trim(line))
             if(line(p:p)==' ') exit
             p=p+1
          enddo
       enddo
       path=adjustl(line(p:))
       if(path(1:1)/='/') cycle
       ! Executable mappings identify the executable and loaded ELF objects, not MPI shared-memory files.
       if(index(line(1:min(p-1,len_trim(line))),'x')==0) cycle
       if(index(path,' (deleted)')/=0) then
          ios=1; exit
       endif
       if(any(paths(1:n)==path)) cycle
       open(newunit=fu,file=trim(path),access='stream',form='unformatted',status='old',action='read',iostat=ios)
       if(ios/=0) exit
       read(fu,iostat=ios) elf
       close(fu,iostat=cs)
       if(ios/=0.or.cs/=0) exit
       if(elf/=achar(127)//'ELF') cycle
       if(index(path,'flexiblas')/=0.or.index(path,'mkl_rt')/=0) then
          ios=1; exit
       endif
       if(n==512) then
          ios=1; exit
       endif
       n=n+1; paths(n)=path
       if(index(path,'libopenblas')/=0) known_backend=.true.
    enddo
    close(unit,iostat=cs)
    if(ios>=0.or.cs/=0.or..not.known_backend) return
    do i=1,n
       call sha256_file(trim(paths(i)),hashes(i+1),ok)
       if(.not.ok) return
    enddo
    ! Sort content digests, excluding relocation addresses and file names.
    do i=2,n+1
       swap=hashes(i); j=i-1
       do while(j>=1)
          if(hashes(j)<=swap) exit
          hashes(j+1)=hashes(j); j=j-1
       enddo
       hashes(j+1)=swap
    enddo
    call state%init(); call state%update('ALPS Linux native cache build v1')
    do i=1,n+1
       if(i>1) then
          if(hashes(i)==hashes(i-1)) cycle
       endif
       call state%update(hashes(i))
    enddo
    call framed(state,compiler_version()); call framed(state,compiler_options())
    build_reason='cannot establish stable Linux CPU identity'
    open(newunit=unit,file='/proc/cpuinfo',status='old',action='read',iostat=ios)
    if(ios/=0) return
    flags_found=.false.
    do
       read(unit,'(a)',iostat=ios) line
       if(ios/=0.or.len_trim(line)==0) exit
       do j=1,len_trim(line)
          if(iachar(line(j:j))==9) line(j:j)=' '
       enddo
       k=index(line,':'); if(k==0) cycle
       path=trim(adjustl(line(:k-1)))
       select case(trim(path))
       case('vendor_id','cpu family','model','model name','stepping','microcode','flags')
          call framed(state,trim(path)//':'//trim(adjustl(line(k+1:))))
          if(trim(path)=='flags') flags_found=.true.
       end select
    enddo
    close(unit,iostat=cs)
    if(.not.flags_found.or.cs/=0) return
    call ieee_get_rounding_mode(rounding)
    build_reason='cannot identify floating-point rounding mode'
    if(rounding==ieee_nearest) then
       call state%update('nearest')
    elseif(rounding==ieee_up) then
       call state%update('upward')
    elseif(rounding==ieee_down) then
       call state%update('downward')
    elseif(rounding==ieee_to_zero) then
       call state%update('tozero')
    else
       return
    endif
    do i=1,size(names)
       build_reason='cannot read numerical environment: '//trim(names(i))
       call get_environment_variable(trim(names(i)),env,length=length,status=ios)
       if(ios==-1.or.length>len(env)) return
       if(ios==1) env='<unset>'
       if(ios/=0.and.ios/=1) return
       call framed(state,trim(names(i))//'='//trim(env))
    enddo
    build_hash=state%finish(); build_ok=.true.; build_reason=''
  end subroutine
  subroutine framed(state,text)
    type(sha256_state), intent(inout) :: state
    character(*), intent(in) :: text
    call state%update(transfer(int(len(text),int64),repeat(' ',8)))
    call state%update(text)
  end subroutine
  subroutine cache_identity(coordinates,source,mass,va,pi,ngamma,npparbar,id,ok,reason)
    !! Hash all consumed source data and construction settings, never species index or filename.
    real(real64), intent(in) :: coordinates(0:,0:,:),source(0:,0:),mass,va,pi
    integer, intent(in) :: ngamma,npparbar
    type(rel_cache_identity), intent(out) :: id
    logical, intent(out) :: ok
    character(*), intent(out) :: reason
    type(sha256_state) :: state
    integer :: j,c
    ok=.false.; reason='unsupported cache dimensions'
    if(ngamma<2.or.npparbar<4.or.ngamma>huge(0)/8-1) return
    if(int(ngamma,int64)+1>(huge(0_int64)-392)/40/(int(npparbar,int64)+1)) return
    call discover_build
    ok=build_ok; reason=build_reason
    if(.not.ok) return
    id%build=build_hash
    id%dims=[size(source,1)-1,size(source,2)-1,ngamma,npparbar]
    id%settings=[mass,va,pi,0.d0]
    if(any(shape(coordinates)/=[size(source,1),size(source,2),2])) then
       ok=.false.; reason='source shape mismatch'; return
    endif
    if(.not.all(ieee_is_finite(source)).or..not.all(ieee_is_finite(coordinates)).or. &
         .not.all(ieee_is_finite(id%settings))) then
       ok=.false.; reason='nonfinite construction input'; return
    endif
    call state%init(); call state%update('ALPS relativistic identity v1')
    call state%update(transfer([format_version,algorithm_revision],repeat(' ',16)))
    call state%update(id%build)
    call state%update(transfer(int(id%dims,int64),repeat(' ',32)))
    call state%update(transfer(id%settings,repeat(' ',32)))
    do c=1,2
       do j=0,id%dims(2)
          call state%update(transfer(coordinates(:,j,c),repeat(' ',8*size(source,1))))
       enddo
    enddo
    do j=0,id%dims(2)
       call state%update(transfer(source(:,j),repeat(' ',8*size(source,1))))
    enddo
    id%key=state%finish()
  end subroutine
  function cache_path(directory,id) result(path)
    character(*), intent(in) :: directory
    type(rel_cache_identity), intent(in) :: id
    character(:), allocatable :: path
    path=trim(directory)//'/alps_rel_'//id%key//'.cache'
  end function
  function metadata(id) result(items)
    type(rel_cache_identity), intent(in) :: id
    integer(int64) :: items(16)
    items=[format_version,algorithm_revision,int(id%dims,int64),0_int64,0_int64,2_int64, &
         64_int64,2_int64,53_int64,32_int64,int(z'01020304',int64), &
         40_int64*(int(id%dims(3),int64)+1)*(int(id%dims(4),int64)+1),1_int64]
  end function
  function header(id,normalization) result(bytes)
    type(rel_cache_identity), intent(in) :: id
    real(real64), intent(in) :: normalization
    character(header_size) :: bytes
    bytes=magic//transfer(metadata(id),repeat(' ',128))// &
         transfer([id%settings,normalization],repeat(' ',40))//id%key//id%build
  end function
  logical function valid_values(values,normalization) result(ok)
    real(real64), intent(in) :: values(0:,0:,:),normalization
    integer :: g,p,ng,np
    real(real64) :: limit
    ok=.false.; ng=ubound(values,1); np=ubound(values,2)
    if(.not.all(ieee_is_finite(values)).or..not.ieee_is_finite(normalization)) return
    if(normalization<=0.d0.or.ng<2.or.np<4) return
    limit=sqrt(huge(1.d0))*0.25d0
    if(any(values(:,:,1)<1.d0).or.any(abs(values(:,:,1:2))>limit)) return
    if(any(values(1:,0,1)<=values(:ng-1,0,1)).or.any(values(0,1:,2)<=values(0,:np-1,2))) return
    do p=0,np
       if(any(values(:,p,1)/=values(:,0,1))) return
       if(any(values(:,p,2)/=values(0,p,2))) return
       do g=0,ng
          if(values(g,p,1)**2-1.d0<values(g,p,2)**2) then
             if(values(g,p,3)/=-1.d0) return
          else
             if(values(g,p,3)<0.d0) return
          endif
       enddo
    enddo
    if(any(values(0,:,4:5)/=0.d0).or.any(values(ng,:,4:5)/=0.d0).or. &
         any(values(:,0,4:5)/=0.d0).or.any(values(:,np,4:5)/=0.d0)) return
    ok=.true.
  end function
  subroutine read_entry(path,id,values,normalization,ok,reason)
    character(*), intent(in) :: path
    type(rel_cache_identity), intent(in) :: id
    real(real64), allocatable, intent(out) :: values(:,:,:)
    real(real64), intent(out) :: normalization
    logical, intent(out) :: ok
    character(*), intent(out) :: reason
    character(header_size) :: bytes,expected
    character(64) :: checksum
    character(16) :: ending
    character(:), allocatable :: row
    type(sha256_state) :: state
    integer(int64) :: length,items(16),wanted(16)
    integer :: unit,ios,cs,p,c,ng,np
    real(real64) :: settings(5)
    ok=.false.; normalization=0; reason='missing or unreadable entry'
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=ios)
    if(ios/=0) return
    wanted=metadata(id); ng=id%dims(3); np=id%dims(4)
    inquire(unit=unit,size=length,iostat=ios)
    if(ios==0) then
       if(length/=int(header_size+trailer_size,int64)+wanted(15)) ios=1
    endif
    if(ios==0) read(unit,iostat=ios) bytes
    reason='truncated entry or incompatible header'
    if(ios==0) then
       items=transfer(bytes(17:144),items)
       if(any(items/=wanted)) ios=1
       if(bytes(1:16)/=magic) ios=1
    endif
    if(ios==0) then
       settings=transfer(bytes(145:184),settings)
       normalization=settings(5); expected=header(id,normalization)
       if(bytes/=expected) ios=1
    endif
    if(ios==0) then
       allocate(values(0:ng,0:np,5),stat=ios)
       if(ios==0) then
          allocate(character(8*(ng+1)) :: row,stat=ios)
       endif
       call state%init(); call state%update(bytes)
       do c=1,5
          do p=0,np
             if(ios/=0) exit
             read(unit,iostat=ios) row
             if(ios/=0) exit
             call state%update(row)
             values(:,p,c)=transfer(row,values(:,p,c))
          enddo
          if(ios/=0) exit
       enddo
       if(ios==0) read(unit,iostat=ios) checksum,ending
       reason='truncated payload or failed integrity check'
       if(ios==0) then
          if(checksum/=state%finish().or.ending/=complete) ios=1
       endif
       if(ios==0) then
          reason='invalid numerical representation'
          if(.not.valid_values(values,normalization)) ios=1
       endif
    endif
    close(unit,iostat=cs)
    ok=ios==0.and.cs==0
    if(ok) then
       reason='hit'
    else
       if(allocated(values)) deallocate(values)
    endif
  end subroutine
  subroutine cache_load(directory,id,gamma,parallel,distribution,derivatives,normalization,hit,reason)
    !! Transactional load: rejected entries never modify caller arrays.
    !! Arrays retain zero-based target-grid bounds; derivative component bounds are 1:2.
    !! On hit, normalization is the saved pre-division integral, not a new integration.
    !! The caller must use an identity from [[cache_identity]] and matching destination shapes.
    character(*), intent(in) :: directory
    type(rel_cache_identity), intent(in) :: id
    real(real64), intent(inout) :: gamma(0:,0:),parallel(0:,0:),distribution(0:,0:),derivatives(0:,0:,:)
    real(real64), intent(out) :: normalization
    logical, intent(out) :: hit
    character(*), intent(out) :: reason
    real(real64), allocatable :: values(:,:,:)
    call read_entry(cache_path(directory,id),id,values,normalization,hit,reason)
    if(.not.hit) return
    gamma=values(:,:,1); parallel=values(:,:,2); distribution=values(:,:,3)
    derivatives(:,:,1)=values(:,:,4); derivatives(:,:,2)=values(:,:,5)
  end subroutine
  subroutine cache_save(directory,id,gamma,parallel,distribution,derivatives,normalization,ok,reason)
    !! Write/close/verify a unique same-directory temporary entry before atomic POSIX rename.
    !! Inputs are read-only. Failure returns ok=false and leaves an existing entry untouched.
    !! Uses libc via ISO_C_BINDING; no C source or additional library is required on Linux.
    character(*), intent(in) :: directory
    type(rel_cache_identity), intent(in) :: id
    real(real64), intent(in) :: gamma(0:,0:),parallel(0:,0:),distribution(0:,0:),derivatives(0:,0:,:)
    real(real64), intent(in) :: normalization
    logical, intent(out) :: ok
    character(*), intent(out) :: reason
    character(:), allocatable :: template,temporary,destination,row
    character(header_size) :: bytes
    character(64) :: checksum
    type(sha256_state) :: state
    real(real64), allocatable :: verify_values(:,:,:)
    real(real64) :: verify_normalization
    integer(c_int) :: fd,status
    integer :: unit,ios,cs,p,c,ng,np
    logical :: verified
    ok=.false.; reason='cannot create cache temporary file'
    template=trim(directory)//'/.alps-rel-XXXXXX'//c_null_char
    fd=posix_mkstemp(template)
    if(fd<0) return
    temporary=template(:len(template)-1); status=posix_close(fd)
    ios=1
    if(status==0) then
       open(newunit=unit,file=temporary,access='stream',form='unformatted',status='old',action='write',iostat=ios)
    endif
    if(ios==0) then
       ng=id%dims(3); np=id%dims(4)
       bytes=header(id,normalization)
       call state%init(); call state%update(bytes)
       write(unit,iostat=ios) bytes
       do c=1,5
          do p=0,np
             if(ios/=0) exit
             select case(c)
             case(1)
                row=transfer(gamma(:,p),repeat(' ',8*(ng+1)))
             case(2)
                row=transfer(parallel(:,p),repeat(' ',8*(ng+1)))
             case(3)
                row=transfer(distribution(:,p),repeat(' ',8*(ng+1)))
             case(4:5)
                row=transfer(derivatives(:,p,c-3),repeat(' ',8*(ng+1)))
             end select
             call state%update(row)
             write(unit,iostat=ios) row
          enddo
          if(ios/=0) exit
       enddo
       checksum=state%finish()
       if(ios==0) write(unit,iostat=ios) checksum,complete
       if(ios==0) flush(unit,iostat=ios)
       close(unit,iostat=cs)
       if(cs/=0) ios=cs
    endif
    reason='cache write or close failed'
    if(ios==0) then
       call read_entry(temporary,id,verify_values,verify_normalization,verified,reason)
       if(verified) then
          destination=cache_path(directory,id)//c_null_char
          status=posix_rename(temporary//c_null_char,destination)
          ok=status==0
          if(ok) then
             reason='saved'
          else
             reason='atomic cache publication failed'
          endif
       endif
    endif
    if(.not.ok) status=posix_unlink(temporary//c_null_char)
  end subroutine
end module
