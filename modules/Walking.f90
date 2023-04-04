!##############################################################################
module Walking
!##############################################################################
	implicit none

Contains

!******************************************************************************
subroutine Propagator(Nb,n1,m1,P1,n2,m2,P2,D,G1arm_b,G2arm_b,WB,GB_b,GB_f,G1arm_f,G2arm_f)
!******************************************************************************
	implicit none
	integer(4),intent(in) :: Nb,n1,m1,P1,n2,m2,P2,D
	real(8),intent(inout) :: WB(0:D+1)
	real(8),intent(inout) :: G1arm_b(0:D+1,1:n1), G1arm_f(0:D+1,1:n1,P1)
	real(8),intent(inout) :: G2arm_b(0:D+1,1:n2), G2arm_f(0:D+1,1:n2,P2)
	real(8),intent(inout) :: GB_b(0:D+1,1:Nb),GB_f(0:D+1,1:Nb)
	integer(4) k, z ! iterator
	
	!!!! forward walking for backbone
	call forward_prop(GB_f(0:D+1,1:m1/2+1),D,m1/2+1,WB)
	GB_f(1:D,m1/2+1) = GB_f(1:D,m1/2+1)*G1arm_b(1:D,1)/WB(1:D)
		
	do k = 2, P1
		call forward_prop(GB_f(0:D+1,(k-1)*m1-m1/2:k*m1-m1/2),D,m1+1,WB)
		GB_f(1:D,k*m1-m1/2) = GB_f(1:D,k*m1-m1/2)*G1arm_b(1:D,1)/WB(1:D)
	enddo
	
	call forward_prop(GB_f(0:D+1,P1*m1-m1/2:P1*m1+1),D,m1/2+2,WB)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	call forward_prop(GB_f(0:D+1,P1*m1+1:P1*m1+m2/2+1),D,m2/2+1,WB)
	GB_f(1:D,P1*m1+m2/2+1) = GB_f(1:D,P1*m1+m2/2+1)*G2arm_b(1:D,1)/WB(1:D)
	
	do k = 2, P2
		call forward_prop(GB_f(0:D+1,P1*m1+(k-1)*m2-m2/2:P1*m1+k*m2-m2/2),D,m2+1,WB)
		GB_f(1:D,P1*m1+k*m2-m2/2) = GB_f(1:D,P1*m1+k*m2-m2/2)*G2arm_b(1:D,1)/WB(1:D)
	enddo
	
	call forward_prop(GB_f(0:D+1,Nb-m2/2:Nb),D,m2/2+1,WB)
	
		
	!!!! back walking for backbone 
	
	call back_prop(GB_b(0:D+1,Nb-m2/2:Nb),D,m2/2+1,WB)
	GB_b(1:D,Nb-m2/2) = GB_b(1:D,Nb-m2/2)*G2arm_b(1:D,1)/WB(1:D)

	do k = P2, 2, -1
		
		do z = 1, D
				G2arm_f(z,1,k) = GB_f(z,P1*m1+k*m2-m2/2)*GB_b(z,P1*m1+k*m2-m2/2) &
			&					/G2arm_b(z,1)/G2arm_b(z,1)*WB(z)
			!write(*,*) GB_f(z,P1*m1+k*m2-m2/2), GB_b(z,P1*m1+k*m2-m2/2)
		enddo
		
		!stop
		
		call back_prop(GB_b(0:D+1,P1*m1+(k-1)*m2-m2/2:P1*m1+k*m2-m2/2),D,m2+1,WB)
		GB_b(1:D,P1*m1+(k-1)*m2-m2/2) = GB_b(1:D,P1*m1+(k-1)*m2-m2/2)*G2arm_b(1:D,1)/WB(1:D)	
	enddo
	
		
	do z = 1, D
		G2arm_f(z,1,1) = GB_f(z,P1*m1+m2-m2/2)*GB_b(z,P1*m1+m2-m2/2) &
			&					/G2arm_b(z,1)/G2arm_b(z,1)*WB(z)
	enddo
	
	call back_prop(GB_b(0:D+1,P1*m1:P1*m1+m2/2+1),D,m2/2+2,WB)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call back_prop(GB_b(0:D+1,Nb-m1/2-P2*m2:Nb-P2*m2),D,m1/2+1,WB)
	GB_b(1:D,Nb-m1/2-P2*m2) = GB_b(1:D,Nb-m1/2-P2*m2)*G1arm_b(1:D,1)/WB(1:D)
				
	do k = P1, 2, -1
		
		do z = 1, D
				G1arm_f(z,1,k) = GB_f(z,k*m1-m1/2)*GB_b(z,k*m1-m1/2) &
			&					/G1arm_b(z,1)/G1arm_b(z,1)*WB(z)
		enddo
		
		call back_prop(GB_b(0:D+1,(k-1)*m1-m1/2:k*m1-m1/2),D,m1+1,WB)
		GB_b(1:D,(k-1)*m1-m1/2) = GB_b(1:D,(k-1)*m1-m1/2)*G1arm_b(1:D,1)/WB(1:D)
	enddo
	
		
	do z = 1, D
		G1arm_f(z,1,1) = GB_f(z,m1-m1/2)*GB_b(z,m1-m1/2) &
			&					/G1arm_b(z,1)/G1arm_b(z,1)*WB(z)  
			
	enddo
	
	call back_prop(GB_b(0:D+1,1:m1/2+1),D,m1/2+1,WB)
	
	
	
	
	!!!! forward walking for arms	
	do k = 1, P1
		call forward_prop(G1arm_f(0:D+1,1:n1,k),D,n1,WB)
	enddo 
	
	do k = 1, P2
		call forward_prop(G2arm_f(0:D+1,1:n2,k),D,n2,WB)
	enddo  

	return
!******************************************************************************
end subroutine Propagator
!******************************************************************************


!******************************************************************************
subroutine forward_prop(G,D,N,W)
!******************************************************************************
implicit none
integer(4), intent(in) :: D, N
real(8), intent(in) :: W(0:D+1)
real(8), intent(inout) :: G(0:D+1,1:N)
integer(4) z, s

do s = 2, N
	do z = 1, D
		G(z,s) = W(z)/6.0*(G(z-1,s-1) + 4.0*G(z,s-1) + G(z+1,s-1))
	enddo
enddo

return
end subroutine forward_prop
!******************************************************************************

!******************************************************************************
subroutine back_prop(G,D,N,W)
!******************************************************************************
implicit none
integer(4), intent(in) :: D, N
real(8), intent(in) :: W(0:D+1)
real(8), intent(inout) :: G(0:D+1,1:N)
integer(4) z, s


do s = 2, N
   do z = 1, D
		G(z,N-s+1) = W(z)/6.0*(G(z-1,N-s+2) + 4.0*G(z,N-s+2) + G(z+1,N-s+2))
   enddo
enddo


return
end subroutine back_prop
!****************************************************************************** 

!##############################################################################
end module Walking
!##############################################################################
