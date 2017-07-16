program example34
	use main
	implicit none
	integer,dimension(4)::timestep
	real(kind=8),dimension(4)::yerror
	real(kind=8)::ordery, theta1
	integer::i, xrange
	print*,"theta1 = "
	read*,theta1
	print*,"xrange = "
	read*,xrange
	
	timestep = (/32,64,126,252/)
	do i = 1, 4
		call li34y(timestep(i),xrange,yerror(i), theta1)
	end do 
	print*,"yerror=",yerror
	
	call leastsquare(timestep, yerror, ordery)
	print*,"ordery=",ordery
	
end program