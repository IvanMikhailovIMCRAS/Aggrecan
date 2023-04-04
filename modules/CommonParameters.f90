!##############################################################################
module CommonParameters
!##############################################################################
	implicit none
	real(8), public, parameter :: tolerance = 1e-7 ! calculation accuracy
	integer(4), public, parameter :: max_iter = 90000000 ! limit of iterations
!!!!!!!!! descriptors for input files !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer(4), public, parameter :: n_input = 11
	character(5), public, parameter :: name_input = 'INPUT'
!!!!!!!!! descriptors fot output files !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer(4), public, parameter :: n_error = 21
	character(5), public, parameter :: name_error = 'ERROR'
	integer(4), public, parameter :: n_energy = 22
	character(10), public, parameter :: name_energy = 'energy.out'
	integer(4), public, parameter :: n_pro = 23
	integer(4), public, parameter :: n_bp = 24
	integer(4), public, parameter :: n_pX = 25
	character(6), public, parameter :: name_pX = 'pX.out'
	integer(4), public, parameter :: n_FX = 26
	character(10), public, parameter :: name_FX = 'FXdata.out'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
Contains

!##############################################################################
end module CommonParameters
!##############################################################################
