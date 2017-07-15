program test
	use example
	implicit none
	real(kind=8), dimension(5)::error
	integer, dimension(5)::timestep
	real(kind=8)::order
	timestep = (/4,8,16,32,64/)
	error = (/1.517e-02,7.418e-03,3.660e-03,1.817e-03,9.055e-04/)
	call leastsquare(timestep, error, order)
	print*,"order",order
end program 

