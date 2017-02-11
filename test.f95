program test
	use leastsq
	implicit none
	real(kind=8), dimension(5)::timesteplength, error
	real(kind=8)::order
	timesteplength = (/0.25, 1.0/8, 1.0/16, 1.0/32, (1.0/64)/)
	print*,"timesteplength",timesteplength
	error = (/1.218e-2, 5.703e-3, 2.772e-3, 1.335e-3, 6.642e-4/)
	print*,"print",error
	call leastsquare(timesteplength,error,order)
	print*,"order",order
end program