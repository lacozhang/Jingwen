program example34
	use main
	implicit none
	integer,dimension(4)::timestep
	real(kind=8),dimension(4)::zerror
	real(kind=8)::ordery, orderz, theta2
	integer::i, xrange
	print*,"xrange = "
	read*,xrange
	print*,"theta2 = "
	read*,theta2
	
	timestep = (/32,64,128, 256/)
	do i = 1, 4
		call li34z(timestep(i),xrange,zerror(i), theta2)
	end do 

	print*,"zerror=",zerror
	call leastsquare(timestep, zerror, orderz)
	print*,"orderz=",orderz

end program