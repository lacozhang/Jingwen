program example31
	use main

	implicit none
	integer,dimension(4)::timestep
	real(kind=8),dimension(4)::yerror, zerror
	real(kind=8)::ordery, orderz, theta
	integer::i ,xrange
	print*,"theta = "
	read*,theta
	!theta = 1.0d0/3.0d0
	print*,"xrange = "
	read*,xrange
	timestep = [ 8, 16, 32, 64]
	do i = 1, 4
		print*, "steps ", timestep(i)
		call li31(timestep(i),xrange,yerror(i),zerror(i), theta)
	end do 
	print*,"yerror=",yerror
	print*,"zerror=",zerror
	call leastsquare(timestep, yerror, ordery)
	print*,"ordery=",ordery
	call leastsquare(timestep, zerror, orderz)
	print*,"orderz=",orderz

end program
