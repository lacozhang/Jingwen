program testgh
	use gausshermite
	implicit none
	integer:: n
	real(kind=8), dimension(:), allocatable:: w
	real(kind=8),dimension(:), allocatable:: x
	n = 5
	allocate(w(n))
	allocate(x(n))
	call gauss_hermite(n, w, x)
	
	print*, "w", w 
	print*, "x", x
	deallocate(w)
	deallocate(x)
end program