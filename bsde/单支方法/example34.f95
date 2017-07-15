program example34
	use main
	implicit none
	integer,dimension(5)::timestep
	real(kind=8),dimension(5)::yerror
	real(kind=8)::ordery, theta
	integer::i,xrange
	print*,"xrange="
	read*, xrange
	print*,"theta="
	read*, theta
	timestep = (/4, 8, 16, 32, 64/)
	do i = 1, 5
		call li34y(timestep(i),xrange,yerror(i), theta)
	end do
	print*,"yerror=",yerror
	
	call leastsquare(timestep, yerror, ordery)
	print*,"ordery=",ordery
	
end program