!===============================================================================
!I                                                                             I
!I                              A  L  P  S                                     I
!I                     Arbitrary Linear Plasma Solver                          I
!I                                                                             I
!I                              Version 1.0                                    I
!I                                                                             I
!I  Kristopher Klein   (kgklein@arizona.edu)                                   I
!I  Daniel Verscharen  (d.verscharen@ucl.ac.uk)                                I
!I                                                                             I
!===============================================================================

module alps_check
	implicit none

	public  :: check_parameters
	!private ::

contains


subroutine check_parameters
  implicit none

  write (*,'(a)') "Running parameter checks..."


  ! THIS IS THE POINT WHERE WE WILL ADD THE PARAMETER CHECKS

  write (*,'(a)') " Parameter checks completed."
  write (*,'(a)') "-=-=-=-=-=-=-=-=-"

end subroutine


end module
