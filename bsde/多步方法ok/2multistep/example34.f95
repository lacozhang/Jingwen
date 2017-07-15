program example34
	use main
	implicit none
	integer,dimension(7)::timestep
	real(kind=8),dimension(7)::yerror
	real(kind=8)::ordery
	integer::i,xrange
	print*,"xrange="
	read*, xrange
	timestep = (/4,8,16,32,64,128,256/)
	do i = 1, 7
		call li34y(timestep(i),xrange,yerror(i))
	end do 
	print*,"yerror=",yerror
	
	call leastsquare(timestep, yerror, ordery)
	print*,"ordery=",ordery
	
end program