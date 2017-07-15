program example51
	use main

	implicit none
	integer,dimension(5)::timestep
	real(kind=8),dimension(5)::yerror, zerror
	real(kind=8)::ordery, orderz, theta
	integer::i

!	read*,theta 
	timestep = [32, 64, 128, 256, 512]
	do i = 1, 5
		call li51(timestep(i), yerror(i), zerror(i))
		print*, "y error", yerror(i)
		print*, "z error", zerror(i)
	end do 
	print*,"yerror=",yerror
	print*,"zerror=",zerror
	call leastsquare(timestep, yerror, ordery)
	print*,"ordery=",ordery
	call leastsquare(timestep, zerror, orderz)
	print*,"orderz=",orderz

end program