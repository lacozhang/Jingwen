program test
	implicit none
	real(kind=4):: A(3,3), b(3),pivot(3)
	integer :: ok, i
	A(1,1)=3.1
	A(1,2)=1.3
	A(1,3)=-5.7
	A(2,1)=1.0
	A(2,2)=-6.9
	A(2,3)=5.8
	A(3,1)=3.4
	A(3,2)=7.2
	A(3,3)=-8.8
	b(1)=-1.3
	b(2)=-0.1
	b(3)=1.8
	call sgesv(3, 1, A, 3, pivot, b, 3, ok)
	print*,b
end program