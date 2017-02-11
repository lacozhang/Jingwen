program test
	use diag
	implicit none
	real(kind=8), dimension(3):: a
	real(kind=8), dimension(:,:), allocatable:: b
	integer:: c
	a = (/1, 2, 3/)
	c = -1
	allocate(b(size(a) + abs(c),size(a)+abs(c)))
	call diagin(a, b, c)
	print*, b 
	deallocate(b)
end program