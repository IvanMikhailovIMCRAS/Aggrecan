!##############################################################################
module ErrorList
!##############################################################################
	use CommonParameters
	implicit none
	
Contains

!******************************************************************************
subroutine ERROR(n)
!******************************************************************************
	implicit none
	integer(4),intent(in) :: n

	open(n_error,file=name_error)

	select case (n)
		case(1)
			write(n_error,'(a,I0,a)') 'ERROR ',n,'. Input file is not found!'
		case(2)
			write(n_error,'(a,I0,a)') 'ERROR ',n,'. Input file is not correct!'
		case(3)
			write(n_error,'(a,I0,a)') 'ERROR ',n,'. Input: Xstart > Xfinish !'
		case(4)
			write(n_error,'(a,I0,a)') 'ERROR ',n,'. Input: dX <= 0.0 !'
		case(5)
			write(n_error,'(a,I0,a)') 'ERROR ',n,'. Input: n < 0 !'
		case(6)
			write(n_error,'(a,I0,a)') 'ERROR ',n,'. Input: m <= 0 !'
		case(7)
			write(n_error,'(a,I0,a)') 'ERROR ',n,'. Input: mod(Np,m) is not 1!'
		case(8)
			write(n_error,'(a,I0,a)') 'ERROR ',n,'. Input: Xstart < 0.0 !'
		case(9)
			write(n_error,'(a,I0,a)') 'ERROR ',n,'. Input: Xfinish > 1.0 !'
		case(10)
			write(n_error,'(a,I0,a)') 'ERROR ',n,'. M*sigma > backbone length!'
		case(11)
			write(n_error,'(a,I0,a)') 'ERROR ',n,'.	New eta < tolerance !'
		case default
			write(n_error,'(a,I0,a)') 'ERROR ',n,'. Unknown error!'
	end select
	
	close(n_error)
	stop
!******************************************************************************
end subroutine ERROR
!******************************************************************************

!##############################################################################
end module ErrorList
!##############################################################################
