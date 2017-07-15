program example34
	use main
	implicit none
	integer,dimension(2)::timestep
	real(kind=8),dimension(2)::yerror
	real(kind=8)::ordery
	integer::i
	timestep = (/ 8, 16/)
	do i = 1, 2
		call li34y(timestep(i),yerror(i))
	end do 
	print*,"yerror=",yerror
	
	call leastsquare(timestep, yerror, ordery)
	print*,"ordery=",ordery
	
end program