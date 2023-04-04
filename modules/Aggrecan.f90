!##############################################################################
program TWOCOMB
!##############################################################################
	use CommonParameters
	use ErrorList
	use InputData
	use OutputData
	use OnePoint
	use Walking

	implicit none

	integer(4) Nb, n1, m1, P1, n2, m2, P2, D, n_free
	integer(4) prn_inf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8) eta, sigma	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8), parameter :: initial_eta = 0.01 ! the step size of convergence in first time
		
	!!!! reading of input parameters
	call ReadInputData(Nb, n1, m1, P1, n2, m2, P2, D, sigma, n_free)
	!!!! Initial parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Nb = m1*P1 + m2*P2 
	!!!!	
	prn_inf = 1

	eta = initial_eta
				
	call CalcOnePoint(Nb, n1, m1, P1, n2, m2, P2, D, sigma, n_free, prn_inf, eta)	
	
	stop
!##############################################################################
end program TWOCOMB
!##############################################################################



