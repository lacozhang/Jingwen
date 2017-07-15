program random
	use delay
	implicit none
	real(kind = 8) :: t
	integer :: n
	real(kind = 8), dimension(:),allocatable :: w
	t =1
	n = 500
	allocate(w(0:n))
	call brownian (t,n,w)
	print*,"w=",w

end program random