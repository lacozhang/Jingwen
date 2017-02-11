program eignvaluetest
	
	implicit none 
	real(kind=8),dimension(3,3):: a
	real(kind=8),dimension(18):: work
	real(kind=8),dimension(3)::w
	integer::info,lwork,LDA
	LDA = 3
	lwork=18


	a(1,:)=0
	a(1,2)=1
	a(2,:)=0
	a(2,3)=sqrt(0.5)
	a(3,:)=0

	print*, 'test start'
	call dsyev('V','U',3, a,LDA,w,work,lwork,info)

end program