"""Compare the actual Fortran routine with its former full-scan search.

Run directly with Python, or with pytest. Requires gfortran; builds only in a
temporary directory. Tensor/Bessel functions are taken unchanged from source.
"""

from pathlib import Path
import shutil
import subprocess
import tempfile


LEGACY_SEARCH = """
  ipparbar_close=-2
  do ipparbar=0,npparbar-1
     if (pparbar_rel(sproc_rel,igamma,ipparbar+1)>pparbar.and. &
         pparbar_rel(sproc_rel,igamma,ipparbar)<=pparbar) ipparbar_close=ipparbar
  enddo
  bracket=ipparbar_close
  if (f0_rel(sproc_rel,igamma,ipparbar_close+1)<=-1.d0) ipparbar_close=ipparbar_close-1
  if (f0_rel(sproc_rel,igamma,ipparbar_close-1)<=-1.d0) ipparbar_close=ipparbar_close+1
  if (pparbar==pparbar_rel(sproc_rel,igamma,npparbar)) ipparbar_close=npparbar-2
  if (ipparbar_close>=npparbar-1) ipparbar_close=npparbar-2
  if (ipparbar_close<=1) ipparbar_close=2
"""

VARIABLES = """
module alps_var
  implicit none
  integer :: sproc=1, npparbar
  real(8) :: ms(1)=1.d0, qs(1)=1.d0, kpar=0.3d0, kperp=0.2d0, vA=0.5d0
  real(8) :: pi=4.d0*atan(1.d0)
  logical :: kperp_norm=.true.
  real(8), allocatable :: pparbar_rel(:,:,:), gamma_rel(:,:,:), f0_rel(:,:,:), df0_rel(:,:,:,:)
end module
"""

DRIVER = """
program compare_search
  use alps_var
  use search_test
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer :: sizes(3)=[4,200,500], sz,g,j,r,m,nn,sgn,b,c,old_b,old_c
  integer :: indices=0, values=0, unsafe=0, masked=0, seed_size
  integer, allocatable :: seed(:)
  real(8) :: p,h,u,limit,imaginary
  complex(8) :: actual,expected,omega
  character(16) :: argument
  call random_seed(size=seed_size)
  allocate(seed(seed_size)); seed=7351; call random_seed(put=seed)
  call get_command_argument(1,argument)
  do sz=1,size(sizes)
     npparbar=sizes(sz)
     allocate(pparbar_rel(1,0:4,0:npparbar),gamma_rel(1,0:4,0:npparbar), &
              f0_rel(1,0:4,0:npparbar),df0_rel(1,0:4,0:npparbar,2))
     do g=0,4
        ! Narrow support, both mask tails, and a fully valid slice.
        limit=0.05d0+0.7d0*g
        do j=0,npparbar
           p=-2.1d0+4.2d0*j/npparbar
           pparbar_rel(1,g,j)=p
           gamma_rel(1,g,j)=sqrt(1.d0+limit**2)
           f0_rel(1,g,j)=exp(-p*p)
           if (gamma_rel(1,g,j)**2-1.d0<p*p) f0_rel(1,g,j)=-1.d0
           df0_rel(1,g,j,1)=cos(p)+0.2d0*g
           df0_rel(1,g,j,2)=sin(p)-0.1d0*g
        enddo
     enddo
     h=pparbar_rel(1,2,2)-pparbar_rel(1,2,1)
     if (len_trim(argument)>0) then
        p=pparbar_rel(1,4,0)-h
        if (argument=='above') p=pparbar_rel(1,4,npparbar)+h
        if (argument=='nan') p=ieee_value(0.d0,ieee_quiet_nan)
        actual=funct_g_rel(1,p,4,cmplx(1.d0,0.1d0,8),1,1)
        error stop 'expected out-of-domain diagnostic'
     endif
     do g=0,4
        do j=0,npparbar
           p=pparbar_rel(1,g,j); call compare(p)
           if (j>0) call compare(nearest(p,-1.d0))
           if (j<npparbar) call compare(nearest(p,1.d0))
           if (j<npparbar) call compare(0.5d0*(p+pparbar_rel(1,g,j+1)))
        enddo
        do r=1,1000
           call random_number(u)
           call compare(pparbar_rel(1,g,0)+u*(pparbar_rel(1,g,npparbar)-pparbar_rel(1,g,0)))
        enddo
     enddo
     deallocate(pparbar_rel,gamma_rel,f0_rel,df0_rel)
  enddo
  if (indices==0.or.values==0.or.unsafe==0.or.masked==0) error stop 'missing coverage'
  print *, 'PASS: exact bracket/center comparisons:',indices,'; exact complex comparisons:',values
  print *, 'Legacy-unsafe inputs checked safely:',unsafe,'; masked stencils diagnosed:',masked
contains
  subroutine compare(x)
    real(8), intent(in) :: x
    logical :: defined
    ! Determine whether the old sequential neighbor accesses are in bounds.
    old_b=-2
    do b=0,npparbar-1
       if (pparbar_rel(1,g,b)<=x.and.x<pparbar_rel(1,g,b+1)) old_b=b
    enddo
    old_c=old_b; defined=.false.
    if (old_c>=0) then
       if (f0_rel(1,g,old_c+1)<=-1.d0) old_c=old_c-1
       if (old_c>=1) then
          if (f0_rel(1,g,old_c-1)<=-1.d0) old_c=old_c+1
          defined=.true.
       endif
    endif
    actual=funct_g_rel(1,x,g,cmplx(1.d0,0.1d0,8),1,1)
    c=center
    if (x==pparbar_rel(1,g,0).and.c/=2) error stop 'lower endpoint center'
    if (x==pparbar_rel(1,g,npparbar).and.c/=npparbar-2) error stop 'upper endpoint center'
    if (c<2.or.c>npparbar-2) error stop 'invalid stencil bounds'
    if (.not.defined) then
       unsafe=unsafe+1
       return
    endif
    b=bracket
    expected=legacy_g(1,x,g,cmplx(1.d0,0.1d0,8),1,1)
    if (b/=bracket.or.b/=old_b.or.c/=center) error stop 'index mismatch'
    indices=indices+1
    if (any(f0_rel(1,g,c-1:c+1)<=-1.d0)) then
       masked=masked+1
       return
    endif
    do m=1,6
       do nn=-1,1
          do sgn=-1,1
             imaginary=sgn*0.07d0
             omega=cmplx(0.9d0,imaginary,8)
             actual=funct_g_rel(1,x,g,omega,nn,m)
             expected=legacy_g(1,x,g,omega,nn,m)
             if (.not.ieee_is_finite(real(actual)).or..not.ieee_is_finite(aimag(actual))) &
                error stop 'nonfinite complex result'
             if (actual/=expected) error stop 'complex mismatch'
             values=values+1
          enddo
       enddo
    enddo
  end subroutine
end program
"""


def test_funct_g_rel_search():
    compiler = shutil.which("gfortran")
    if compiler is None:
        raise RuntimeError("gfortran is required for the Fortran search regression")
    source = (Path(__file__).resolve().parents[1] / "src/ALPS_fns_rel.f90").read_text()
    start = source.index("double complex function funct_g_rel(")
    end = source.index("end function funct_g_rel", start) + len("end function funct_g_rel")
    candidate = source[start:end]
    search_start = candidate.index("  if (.not.ieee_is_finite(pparbar)")
    search_end = candidate.index("\t! calculate the function on the grid")
    legacy = candidate[:search_start] + LEGACY_SEARCH + candidate[search_end:]
    legacy = legacy.replace("funct_g_rel", "legacy_g")
    candidate = candidate.replace(
        "     if (ipparbar_close.EQ.-2) then",
        "     bracket=ipparbar_close\n     if (ipparbar_close.EQ.-2) then",
    )
    # Instrument copies only: expose selected indices and avoid known invalid
    # square roots for masked stencils. All unmasked calls use the real tensor.
    def instrument(routine, name):
        marker = "\t! calculate the function on the grid"
        return routine.replace(marker, f"""
  center=ipparbar_close
  if (any(f0_rel(sproc_rel,igamma,center-1:center+1)<=-1.d0)) then
     {name}=cmplx(0.d0,0.d0,8)
     return
  endif
{marker}""")

    tensor_start = source.index("double complex function int_T_rel(")
    tensor_end = source.index("end function int_T_rel", tensor_start) + len("end function int_T_rel")
    bessel_start = source.index("DOUBLE PRECISION FUNCTION BESSJ (")
    helpers = source[tensor_start:tensor_end] + "\n" + source[bessel_start:source.rindex("end module")]
    module = (VARIABLES + "\nmodule search_test\ninteger :: bracket,center\ncontains\n"
              + instrument(candidate, "funct_g_rel") + "\n"
              + instrument(legacy, "legacy_g") + "\n" + helpers + "\nend module\n" + DRIVER)
    with tempfile.TemporaryDirectory(prefix="alps-search-test-") as tmp:
        tmp = Path(tmp)
        (tmp / "test.f90").write_text(module)
        subprocess.run([compiler, "-O2", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
                        "-ffree-line-length-none", "test.f90", "-o", "test"], cwd=tmp, check=True)
        subprocess.run([str(tmp / "test")], cwd=tmp, check=True)
        for argument in ("below", "above", "nan"):
            result = subprocess.run([str(tmp / "test"), argument], cwd=tmp, capture_output=True, text=True)
            assert result.returncode != 0
            diagnostic = "nonfinite parallel momentum" if argument == "nan" else "momentum outside parallel grid"
            assert diagnostic in result.stderr


if __name__ == "__main__":
    test_funct_g_rel_search()
