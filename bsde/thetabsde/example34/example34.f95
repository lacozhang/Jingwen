program example34
	use main

	implicit none
	integer,dimension(5)::timestep
	real(kind=8),dimension(5)::yerror, zerror
	real(kind=8)::ordery, orderz, theta
	integer::i
	print*,"theta = "
	read*,theta 
	timestep = (/4, 8, 16, 32, 64/)
	do i = 1, 5
		call li34(timestep(i),yerror(i),zerror(i), theta)
	end do 
	print*,"yerror=",yerror
	print*,"zerror=",zerror
	call leastsquare(timestep, yerror, ordery)
	print*,"ordery=",ordery
	call leastsquare(timestep, zerror, orderz)
	print*,"orderz=",orderz

end program