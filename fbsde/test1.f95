program test1
	use main
	
	implicit none
	real(kind=8)::stepsize
	integer:: xint
	stepsize=1.0/8.0
	call thelengthofx(10, 8, stepsize, xint)
	print*, xint
end program test1
