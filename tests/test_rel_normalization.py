"""Exercise the actual normalization/division block with independent panel sums.

Run with python3 or pytest. Requires NumPy and gfortran; builds in /tmp only.
"""

from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np


def test_rel_normalization():
    source = (Path(__file__).resolve().parents[1] / "src/ALPS_fns_rel.f90").read_text()
    routine = source.split("subroutine construct_f0_rel(is,is_rel,integrate)", 1)[1]
    start = routine.index("      integrate = 0.d0")
    end = routine.index("   if(writeOut) write", start)
    block = routine[start:end]
    driver = """
program normalization_test
  implicit none
  integer :: ngamma,npparbar,is=1,is_rel=1,igamma,ipparbar,ipparbar_lower,ipparbar_upper
  real(8) :: integrate,dgamma_rel,dpparbar,slice_integral,gamma_weight
  real(8) :: pi=4.d0*atan(1.d0),ms(1)=2.d0,vA=0.7d0
  real(8), allocatable :: gamma_rel(:,:,:),pparbar_rel(:,:,:),f0_rel(:,:,:)
  read(*,*) ngamma,npparbar
  allocate(gamma_rel(1,0:ngamma,0:npparbar),pparbar_rel(1,0:ngamma,0:npparbar), &
           f0_rel(1,0:ngamma,0:npparbar))
  do igamma=0,ngamma
     do ipparbar=0,npparbar
        read(*,*) gamma_rel(1,igamma,ipparbar),pparbar_rel(1,igamma,ipparbar),f0_rel(1,igamma,ipparbar)
     enddo
  enddo
""" + block + """
  write(*,'(es26.17)') integrate
  do igamma=0,ngamma
     do ipparbar=0,npparbar
        write(*,'(es26.17)') f0_rel(1,igamma,ipparbar)
     enddo
  enddo
end program
"""
    compiler = shutil.which("gfortran")
    if compiler is None:
        raise RuntimeError("gfortran is required for the normalization regression")
    with tempfile.TemporaryDirectory(prefix="alps-normalization-test-") as tmp:
        tmp = Path(tmp)
        (tmp / "test.f90").write_text(driver)
        subprocess.run([compiler, "-O2", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
                        "-ffree-line-length-none", "test.f90", "-o", "test"], cwd=tmp, check=True)
        for counts in ([5, 5, 5], [0, 1, 2, 5, 3], [1, 0, 4, 1, 0], [2, 2, 2]):
            gamma = np.linspace(1.25, 3.25, len(counts))
            momentum = np.linspace(-2, 2, 5)
            fields = np.full((len(counts), 5), -1.0)
            for i, count in enumerate(counts):
                lower = (5 - count) // 2
                fields[i, lower:lower + count] = 1 + 0.2 * gamma[i] + 0.1 * momentum[lower:lower + count]
            # Independent interval panels: neither isolated points nor masked
            # neighbors create width, and no extrapolation to a cone is used.
            def quadrature(values):
                f = gamma[:, None] * values
                panels = 0.5 * (f[:, 1:] + f[:, :-1]) * np.diff(momentum)
                panels[(values[:, 1:] <= -1) | (values[:, :-1] <= -1)] = 0
                rows = panels.sum(axis=1)
                return np.sum(0.5 * (rows[1:] + rows[:-1]) * np.diff(gamma)) * 2 * np.pi * (2 / 0.7)**3

            records = [f"{len(gamma)-1} 4"]
            records += [f"{g:.17g} {p:.17g} {fields[i,j]:.17g}"
                        for i, g in enumerate(gamma) for j, p in enumerate(momentum)]
            result = subprocess.run([str(tmp / "test")], input="\n".join(records),
                                    text=True, capture_output=True, check=True)
            values = np.fromstring(result.stdout, sep=" ")
            normalized = values[1:].reshape(fields.shape)
            np.testing.assert_allclose(values[0], quadrature(fields), rtol=2e-15)
            np.testing.assert_allclose(quadrature(normalized), 1, rtol=2e-15)
            np.testing.assert_array_equal(normalized[fields == -1], -1)
            np.testing.assert_allclose(normalized[fields > -1], fields[fields > -1] / values[0], rtol=2e-15)
    print("PASS: full, masked, empty, singleton, two-point slices; endpoint weights, prefactor, unit integral, masks")


if __name__ == "__main__":
    test_rel_normalization()
