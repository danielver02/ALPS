! Adversarial files for the numerical regression reader, generated only in scratch.
program test_comparison_fixtures
  implicit none
  integer :: file
  open(newunit=file,file='reader_reference.dat',status='new')
  write(file,'(a)') 'known_value 1.0 0.0 1e-12 1e-10'
  close(file)
  open(newunit=file,file='reader_changed.dat',status='new')
  write(file,'(a)') 'known_value 1.01 0.0 1e-12 1e-10'
  close(file)
  open(newunit=file,file='reader_nonfinite.dat',status='new')
  write(file,'(a)') 'known_value NaN 0.0 1e-12 1e-10'
  close(file)
  open(newunit=file,file='reader_empty.dat',status='new')
  close(file)
end program
