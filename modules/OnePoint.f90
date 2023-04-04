!##############################################################################
module OnePoint
!##############################################################################
	use CommonParameters
	use Walking
	use ErrorList
	use OutputData
	implicit none

Contains

!******************************************************************************
subroutine CalcOnePoint(Nb, n1, m1, P1, n2, m2, P2, D, n_free, prn_inf, eta)
!******************************************************************************
	implicit none
	integer(4),intent(in) :: Nb, n1, m1, P1, n2, m2, P2, D, n_free
	integer(4),intent(in) :: prn_inf
	real(8),intent(inout) :: eta ! the step size of convergence
	real(8) alpha(0:D+1) !Lagrange field (potential)	 
	real(8) F ! free energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer(4) s, z 	! s - number of segment, z - coordinate
	integer(4) i, j, k 	! iterators
	integer(4) iter 	! main iterator
	!!!! Green's functions (back (_b) and forward (_f) propagators)
	real(8) Grl_b(0:D+1,1:Nb), Grl_f(0:D+1,1:Nb)	    ! for (r)ight loop backbone
	real(8) Gll_b(0:D+1,1:Nb), Gll_f(0:D+1,1:Nb)	    ! for (l)eft loop backbone
	real(8) G1arm_b(0:D+1,1:n1), G2arm_b(0:D+1,1:n2)	! back propagator for arms
	!!!! forward propagators for arms (:,:,P)
	real(8) G1Rarm_f(0:D+1,1:n1,1:P1), G2Rarm_f(0:D+1,1:n2,1:P2)     ! (R)ight loop 
	real(8) G1Larm_f(0:D+1,1:n1,1:P1), G2Larm_f(0:D+1,1:n2,1:P2)     ! (L)eft loop
	
	
	
	real(8) phi_rl(0:D+1), phi_ll(0:D+1), WB(0:D+1)
	!!!! volume fraction profile of backbones
	real(8) phi0_rl(0:D+1), phi0_ll(0:D+1), phi_m_rl(0:D+1), phi_m_ll(0:D+1)
	!!!! distributions of branching points
	real(8) sigma
	!!!! partition functions
	real(8) Qll, Qrl
	!!!!
	real(8) Fmix, AL, AR
	!!!! 
	real(8) deviation, deviation_old  
	!!!!
	integer(4) counter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!! grafting densities:
	sigma = 0.5*dble(D)/dble(Nb+P1*(n1-1)+P2*(n2-1))   ! grafting density 
	!!!!
	iter = 0
	deviation = dble(D)*1000000.0
	counter = 0
	
	!!!!!!!!!!!!!!!!!!!!! zero-initialization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	Grl_b(:,:) = 0.0; Grl_f(:,:) = 0.0
	Gll_b(:,:) = 0.0; Gll_f(:,:) = 0.0	    
	G1arm_b(:,:) = 0.0; G2arm_b(:,:) = 0.0	
	G1Rarm_f(:,:,:) = 0.0; G1Rarm_f(:,:,:) = 0.0   
	G2Larm_f(:,:,:) = 0.0; G2Larm_f(:,:,:) = 0.0  
	phi_rl(:) = 0.0; phi_ll(:) = 0.0; WB(:) = 0.0
	phi0_rl(:) = 0.0; phi0_ll(:) = 0.0; phi_m_rl(:) = 0.0; phi_m_ll(:) = 0.0
	
	alpha(:) = 0.0
	
!!!!!!!!!!!!!!!!!!!! ENGINE OF PROGRAM  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do while (deviation.gt.tolerance)
		iter = iter + 1
		if (iter.ge.max_iter) then
			write(n_error,'(a,I0)') 'Iteration steps > ', max_iter
			exit
		endif
		
		!alpha(1:D) = alpha(1:D) - alpha((D+1)/2) 
				
		!!!! express Boltzmann weight through the Lagrange field
		WB(1:D) = exp(-alpha(1:D))   
		!!!! initial conditions
       	! Left loop backbone :
		Gll_f(1,1) = WB(1)	! first segment in first layer
		Gll_b(1:D,Nb) = WB(1:D) ! and last segment  in first layer
		! Right loop backbone :
		Grl_f(D,1) = WB(D)	! first segment in last layer
		Grl_b(1:D,Nb) = WB(1:D) ! and last segment  in last layer 
        ! Arms
		G1arm_b(1:D,n1) = WB(1:D) ! end of arms is free   
		G2arm_b(1:D,n2) = WB(1:D) ! end of arms is free   
				
		!!!!!!!!!! MIRROR BOUNDARY CONDITIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		WB(0) = WB(1); WB(D+1) = WB(D)
		
		Gll_f(0,:) = Gll_f(1,:); Gll_f(D+1,:) = Gll_f(D,:)
		Gll_b(0,:) = Gll_b(1,:); Gll_b(D+1,:) = Gll_b(D,:)
		
		Grl_f(0,:) = Grl_f(1,:); Grl_f(D+1,:) = Grl_f(D,:)
		Grl_b(0,:) = Grl_b(1,:); Grl_b(D+1,:) = Grl_b(D,:)	
		
		G1arm_b(0,:) = G1arm_b(1,:); G1arm_b(D+1,:) = G1arm_b(D,:)
		G2arm_b(0,:) = G2arm_b(1,:); G2arm_b(D+1,:) = G2arm_b(D,:)
		
									
		!!!! Now we iterate over all Green's functions (propagators) 
		call back_prop(G1arm_b,D,n1,WB) ! back walking of arms
		call back_prop(G2arm_b,D,n2,WB) ! back walking of arms
		
					
		! walking of left loop
		call Propagator(Nb,n1,m1,P1,n2,m2,P2,D,G1arm_b,G2arm_b,WB,Gll_b,Gll_f,G1Larm_f,G2Larm_f)
		! walking of right loop
		call Propagator(Nb,n1,m1,P1,n2,m2,P2,D,G1arm_b,G2arm_b,WB,Grl_b,Grl_f,G1Rarm_f,G2Rarm_f)
		
			
		
		G1Rarm_f(0,:,:) = G1Rarm_f(1,:,:); G1Rarm_f(D+1,:,:) = G1Rarm_f(D,:,:) 
		G1Larm_f(0,:,:) = G1Larm_f(1,:,:); G1Larm_f(D+1,:,:) = G1Larm_f(D,:,:) 
		
		G2Rarm_f(0,:,:) = G2Rarm_f(1,:,:); G2Rarm_f(D+1,:,:) = G2Rarm_f(D,:,:) 
		G2Larm_f(0,:,:) = G2Larm_f(1,:,:); G2Larm_f(D+1,:,:) = G2Larm_f(D,:,:) 
		
		
		!!!!! Calculation phi(z)  according to the Compositions law 
		phi_rl(:) = 0.0; phi_ll(:) = 0.0
		!!!! segment-by-segment summation (within a backbone)
				
		do s = 1, P1*m1
			!!!! branching points not accepted
			if (modulo(s+m1/2+1,m1).ne.1.or.s.eq.1.or.s.eq.P1*m1) then
				phi_rl(:) = phi_rl(:) + Grl_f(:,s)*Grl_b(:,s)/WB(:)
				phi_ll(:) = phi_ll(:) + Gll_f(:,s)*Gll_b(:,s)/WB(:)
			endif
		enddo
		
		do s = P1*m1+1, Nb
			!!!! branching points not accepted
			if (modulo(s+m2/2+1-P1*m1,m2).ne.1.or.s.eq.1.or.s.eq.Nb) then
				phi_rl(:) = phi_rl(:) + Grl_f(:,s)*Grl_b(:,s)/WB(:)
				phi_ll(:) = phi_ll(:) + Gll_f(:,s)*Gll_b(:,s)/WB(:)
			endif
		enddo
		
		!!!! segment-by-segment summation (within each k-th arm)
    	do k = 1, P1
    		do s = 1, n1
				phi_ll(:) = phi_ll(:) + G1arm_b(:,s)*G1Larm_f(:,s,k)/WB(:)
				phi_rl(:) = phi_rl(:) + G1arm_b(:,s)*G1Rarm_f(:,s,k)/WB(:)
			enddo  
		enddo
		
		do k = 1, P2
    		do s = 1, n2
				phi_ll(:) = phi_ll(:) + G2arm_b(:,s)*G2Larm_f(:,s,k)/WB(:)
				phi_rl(:) = phi_rl(:) + G2arm_b(:,s)*G2Rarm_f(:,s,k)/WB(:)
			enddo  
		enddo
		!!!! partition functions
		Qll = sum(phi_ll(1:D))
		Qrl = sum(phi_rl(1:D))
		
		!!!! Normalization of volume fractions (and freedom conditions)
		phi_ll(:) = phi_ll(:) * 0.5 * dble(D) / Qll
		phi_rl(:) = phi_rl(:) * 0.5 * dble(D) / Qrl
		
		!!!! Gradient descent
		do z = 1, D 
			alpha(z) = alpha(z) + eta*(phi_ll(z) + phi_rl(z) - 1.0) 
		enddo
		
				
		deviation_old = deviation
		deviation = 0.0
       	do z = 1, D
       		deviation = deviation + (phi_ll(z) + phi_rl(z) - 1.0)**2 
        enddo
        
        deviation = sqrt(deviation)
        
       	!!!! automatic selection of step size $eta$ 
		if ((deviation_old-deviation).lt.tolerance**2 &
		&				.or.deviation.ne.deviation) then
			if (eta.lt.tolerance) then
				call ERROR(11)
			endif
			if (counter.lt.n_free.and.deviation.eq.deviation) then
				counter = counter + 1
			else
				counter = 0
				alpha(:) = 0.0
				
				Grl_b(:,:) = 0.0; Grl_f(:,:) = 0.0
				Gll_b(:,:) = 0.0; Gll_f(:,:) = 0.0	    
				G1arm_b(:,:) = 0.0; G2arm_b(:,:) = 0.0	
				G1Rarm_f(:,:,:) = 0.0; G2Rarm_f(:,:,:) = 0.0   
				G1Larm_f(:,:,:) = 0.0; G2Larm_f(:,:,:) = 0.0  
				WB(:) = 0.0
				
			    deviation = dble(D)*1000000.0
			    eta = 2.0*eta/(sqrt(5.0)+1.0) ! golden ratio
			    write(*,'(a,I0,5x,a,E16.10,5x,a,E16.10)') &
	&		 'iter = ', iter, 'dev = ', deviation, 'eta = ', eta	
			endif
		
		endif
		
		!!!!
		!if (prn_inf.ne.0) write(*,'(a,I0,5x,a,E16.10,5x,a,E16.10)') &
		!&		 'iter = ', iter, 'dev = ', deviation, 'eta = ', eta	
		
		
						
!!!!!!!!!!!!!!!!!!!! END OF ENGINE PROGRAM CYCLE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	enddo
	
	
	if (prn_inf.ne.0) then
	phi0_ll(:) = 0.0
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do s = 1, P1*m1
		!!!! branching points not accepted
		if (modulo(s+m1/2+1,m1).ne.1.or.s.eq.1.or.s.eq.P1*m1) then
			phi0_ll(:) = phi0_ll(:) + Gll_f(:,s)*Gll_b(:,s)/WB(:)
		endif
	enddo
		
	do s = P1*m1+1, Nb
		!!!! branching points not accepted
		if (modulo(s+m2/2+1-P1*m1,m2).ne.1.or.s.eq.1.or.s.eq.Nb) then
			phi0_ll(:) = phi0_ll(:) + Gll_f(:,s)*Gll_b(:,s)/WB(:)
		endif
	enddo
		
	!!!! segment-by-segment summation (within each k-th arm)
    do k = 1, P1
    	phi0_ll(:) = phi0_ll(:) + G1arm_b(:,1)*G1Larm_f(:,1,k)/WB(:)
	enddo
		
	do k = 1, P2
    	phi0_ll(:) = phi0_ll(:) + G2arm_b(:,1)*G2Larm_f(:,1,k)/WB(:)
	enddo
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	phi0_ll(:) = phi0_ll(:) * 0.5 * dble(D) / Qll ! * sigma / Gll_b(1,1)
	
	
	phi_m_ll(:) = Gll_f(:,m1*P1)*Gll_b(:,m1*P1)/WB(:)
	phi_m_rl(:) = Gll_f(:,Nb)*Gll_b(:,Nb)/WB(:)
	
	phi_m_ll(:) = phi_m_ll(:) / Gll_b(1,1) !/ Qll * dble(Nb+n1*P1+n2*P2)
	phi_m_rl(:) = phi_m_rl(:) / Gll_b(1,1) !/ Qll * dble(Nb+n1*P1+n2*P2)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	open(112,file='profile.pro')
	do z = 1, D
		write(112,*) z, phi0_ll(z), phi_ll(z), alpha(z)-alpha(D/2+1), phi_m_ll(z), phi_m_rl(z)
	enddo
	close(112)
	
	endif	
	!if (prn_inf.eq.0) write(*,'(a,I0,5x,a,E16.10,5x,a,E16.10)') &
	!&		 'iter = ', iter, 'dev = ', deviation, 'eta = ', eta	
	
	write(*,'(a,I0,5x,a,E16.10,5x,a,E16.10)') &
	&		 'iter = ', iter, 'dev = ', deviation, 'eta = ', eta	
	
!!!!!!!!!!!!!!!!! calculation of sum(U*phi) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	AL = 0.0
	AR = 0.0
	do z = 1, D
		AL = AL + alpha(z)*phi_ll(z)
		AR = AR + alpha(z)*phi_rl(z)
	enddo
!!!!!!!!!!!!!!!!!! caclulation of free energy and energy parts !!!!!!!!!!!!!!!!
	
		
	F = log(dble(D)) - log(Gll_b(1,1))  &
	& - (AL + AR)*dble(Nb+P1*(n1-1)+P2*(n2-1))/dble(D) &
	& - log(dble(Nb+P1*(n1-1)+P2*(n2-1))) 
	
	
	open(78,file='F.txt',access='append')
		write(78,'(I6,F16.8)') D, F
	close(78)
	!!!!!!*********!!!!!!!!!!*********!!!!!!!!!*********!!!!!!!!!********!!!***
	

	
			
	return
!******************************************************************************
end subroutine CalcOnePoint
!******************************************************************************

!##############################################################################
end module OnePoint
!##############################################################################
