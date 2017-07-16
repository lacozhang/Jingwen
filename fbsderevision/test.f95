program test
	use main 
	implicit none
	real(kind=8):: x0
	integer :: k, timestep
	real(kind=8):: xrange
	x0 = 1.0d0
	timestep = 4
	xrange = 10
	k = 8
	call thelengthofx(x0,k,timestep,xrange)
	print*, "over"
	print*,"xrange=",xrange 
end program
