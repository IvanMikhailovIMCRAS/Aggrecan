!##############################################################################
module InputData
!##############################################################################
	use CommonParameters
	use ErrorList
	implicit none

Contains

!******************************************************************************
subroutine ReadInputData(Nb, n1, m1, P1, n2, m2, P2, D, sigma, n_free)
!******************************************************************************
	integer(4),intent(out) :: Nb, n1, m1, P1, n2, m2, P2, D, n_free
	real(8), intent(out) :: sigma
	integer(4) ioer
	
	open(n_input,file=name_input,iostat=ioer,status='old')
	if (ioer.ne.0) call ERROR(1) 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	read(n_input,*,iostat=ioer) P1  ! number of side chains
	if (ioer.ne.0) call ERROR(2)  
	read(n_input,*,iostat=ioer) n1  ! length of arm
	if (ioer.ne.0) call ERROR(2)  
	read(n_input,*,iostat=ioer) m1  ! spacer length of backbone
	if (ioer.ne.0) call ERROR(2)  
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	read(n_input,*,iostat=ioer) P2  ! number of side chains
	if (ioer.ne.0) call ERROR(2) 
	read(n_input,*,iostat=ioer) n2  ! length of arm
	if (ioer.ne.0) call ERROR(2)  
	read(n_input,*,iostat=ioer) m2  ! spacer length of backbone
	if (ioer.ne.0) call ERROR(2)  
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	read(n_input,*,iostat=ioer) D  ! distance between walls
	if (ioer.ne.0) call ERROR(2)
	read(n_input,*,iostat=ioer) sigma  ! grafting density
	if (ioer.ne.0) call ERROR(2)  
	read(n_input,*,iostat=ioer) n_free  ! число холостых шагов сходимости
	if (ioer.ne.0) call ERROR(2) 
	close(n_input)
	
		
	if (n1.lt.0.or.n2.lt.0)				call ERROR(5)
	if (m1.le.0.or.m2.lt.0)				call ERROR(6)
	
	Nb = m1*P1 + m2*P2
		
	if ((Nb + n1*P1 + n2*P2)*sigma .gt. Nb)	then
		write(*,*) 'Error! M*sigma > backbone length! You should cange parameters ;)'
		call ERROR(10)
	
	endif
		
	n1 = n1 + 1   !!!!!! plus one segment of backbone
	n2 = n2 + 1   !!!!!! plus one segment of backbone
	
	
	return
!******************************************************************************
end subroutine ReadInputData
!******************************************************************************

!##############################################################################
end module InputData
!##############################################################################
