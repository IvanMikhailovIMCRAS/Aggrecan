!##############################################################################
module OutputData
!##############################################################################
	use CommonParameters
	implicit none

Contains

!******************************************************************************
subroutine print_free_energy_title()
!******************************************************************************
	implicit none
	open(n_energy,file=name_energy)
	write(n_energy,'(12x,a,19x,a,18x,a,14x,a,14x,a,16x,a)') &
	&              'X','F','Fmix','Fbridge','Floop','Fpol'
	close(n_energy)
	return
!******************************************************************************
end subroutine print_free_energy_title
!******************************************************************************

!******************************************************************************
subroutine print_free_energy(X,F,Fbridge,Floop)
!******************************************************************************
	implicit none
	real(8),intent(in) :: X, F, Fbridge, Floop
	real(8) Fmix
	if (X.eq.0) then
		Fmix = 0.0
	else if (X.eq.1.0) then
		Fmix = 0.0
	else
		Fmix = (1.0 - X)*log(1.0 - X) + X*log(X)
	endif
	open(n_energy,file=name_energy,access='APPEND')
	write(n_energy,'(E20.10e3,E20.10e3,E20.10e3,E20.10e3,E20.10e3,E20.10e3)') & 
			&		X, F, Fmix, Fbridge, Floop, Fbridge + Floop
	close(n_energy)
	return
!******************************************************************************
end subroutine print_free_energy
!******************************************************************************

!******************************************************************************
subroutine print_profile(D,phi0_br,phi0_ll,phi0_rl, &
	&					 phi_br,phi_ll,phi_rl,alpha,phi_m,X)
!******************************************************************************
	implicit none
	integer(4),intent(in) :: D
	real(8),intent(inout) :: X
	real(8),intent(in) :: phi_br(0:D+1), phi_ll(0:D+1), phi_rl(0:D+1)
	real(8),intent(in) :: phi0_br(0:D+1), phi0_ll(0:D+1), phi0_rl(0:D+1)
	real(8),intent(in) :: alpha(0:D+1),phi_m(0:D+1) 
	character(128) name_pro
	integer(4) z
	!!!!
	
	call create_name('profile',X,'pro',name_pro)
	open(n_pro,file=name_pro) 
	write(n_pro,'(9(1x,a,1x))') &
&	'z', 'phi0_BR', 'phi0_RL', 'phi0_LL', 'phi_BR', 'phi_RL', 'phi_LL', 'alpha', 'end_m'
	do z = 1, D
		write(n_pro,'(I9,1x,8(E20.10e3))') &
		z, phi0_br(z), phi0_rl(z), phi0_ll(z), &
		phi_br(z), phi_rl(z), phi_ll(z), alpha(z)-alpha((D+1)/2), phi_m(z)
	enddo
		
	close(n_pro)
		
	return
!******************************************************************************
end subroutine print_profile
!******************************************************************************

!******************************************************************************
subroutine print_branching_points(D,Nb,m,n_br,n_ll,n_rl,X)
!******************************************************************************
	implicit none
	integer(4),intent(in) :: D, Nb, m
	real(8),intent(inout) :: X
	real(8),intent(in) :: n_br(0:D+1,1:Nb/m)
	real(8),intent(in) :: n_rl(0:D+1,1:Nb/m)
	real(8),intent(in) :: n_ll(0:D+1,1:Nb/m)
	character(128) name_bp
	integer(4) i, k
	integer(4) z
	
	do i = 1, 3
		select case (i)
			case(1)
				call create_name('bpBridge',X,'bp',name_bp)
				open(n_bp,file=name_bp)
			case(2)
				call create_name('bpLloop',X,'bp',name_bp)
				open(n_bp,file=name_bp)
			case(3)
				call create_name('bpRloop',X,'bp',name_bp)
				open(n_bp,file=name_bp)
		end select
		
		write(n_bp,'(1x,a,1x)',ADVANCE='no') 'z'
		
		do k = 1, Nb/m
			write(n_bp,'(1x,a,I0,1x)',ADVANCE='no') 'n', k
		enddo
		
		write(n_bp,'(a)',ADVANCE='yes') ' '
		
		do z = 1, D
			write(n_bp,'(I9,1x)',ADVANCE='no') z
			do k = 1, Nb/m
				if (i.eq.1) write(n_bp,'(E20.10e3)',ADVANCE='no') n_br(z,k)
				if (i.eq.2) write(n_bp,'(E20.10e3)',ADVANCE='no') n_ll(z,k)
				if (i.eq.3) write(n_bp,'(E20.10e3)',ADVANCE='no') n_rl(z,k)
			enddo
			write(n_bp,'(a)',ADVANCE='yes') ' '
		enddo		
		close(n_bp)
	enddo
	return
!******************************************************************************
end subroutine print_branching_points
!******************************************************************************

!******************************************************************************
subroutine print_end_points(D,Nb,m,n_br,n_ll,n_rl,X)
!******************************************************************************
	implicit none
	integer(4),intent(in) :: D, Nb, m
	real(8),intent(inout) :: X
	real(8),intent(in) :: n_br(0:D+1,1:Nb/m)
	real(8),intent(in) :: n_rl(0:D+1,1:Nb/m)
	real(8),intent(in) :: n_ll(0:D+1,1:Nb/m)
	character(128) name_bp
	integer(4) i, k
	integer(4) z
	
	do i = 1, 3
		select case (i)
			case(1)
				call create_name('endBridge',X,'bp',name_bp)
				open(n_bp,file=name_bp)
			case(2)
				call create_name('endLloop',X,'bp',name_bp)
				open(n_bp,file=name_bp)
			case(3)
				call create_name('endRloop',X,'bp',name_bp)
				open(n_bp,file=name_bp)
		end select
		
		write(n_bp,'(1x,a,1x)',ADVANCE='no') 'z'
		
		do k = 1, Nb/m-1
			write(n_bp,'(1x,a,I0,1x)',ADVANCE='no') 'n', k
		enddo
		
		write(n_bp,'(a)',ADVANCE='yes') ' '
		
		do z = 1, D
			write(n_bp,'(I9,1x)',ADVANCE='no') z
			do k = 1, Nb/m-1
				if (i.eq.1) write(n_bp,'(E20.10e3)',ADVANCE='no') n_br(z,k)
				if (i.eq.2) write(n_bp,'(E20.10e3)',ADVANCE='no') n_ll(z,k)
				if (i.eq.3) write(n_bp,'(E20.10e3)',ADVANCE='no') n_rl(z,k)
			enddo
			write(n_bp,'(a)',ADVANCE='yes') ' '
		enddo		
		close(n_bp)
	enddo
	return
!******************************************************************************
end subroutine print_end_points
!******************************************************************************

!******************************************************************************
subroutine create_name(prefix,X,ext,name_file)
!******************************************************************************
	implicit none
	character(*),intent(in) :: prefix, ext
	real(8),intent(in) :: X
	character(128),intent(out) :: name_file
	character(128) format_name
	integer(4) i
	
	write(format_name,'(a,I0,a,I0,a)') '(a,F', 8, '.', 6, ')'
	write(name_file,format_name) trim(prefix), X
	do i = len(trim(name_file)), 1, -1
		if (name_file(i:i).eq.'.') then
			name_file(i:i) = '_'
			exit
		endif
	enddo
	write(name_file,'(a,a,a)') trim(name_file),'.',trim(ext)
	
	return
!******************************************************************************
end subroutine create_name
!******************************************************************************


!##############################################################################
end module OutputData
!##############################################################################
