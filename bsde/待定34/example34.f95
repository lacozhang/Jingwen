program example34
	use main
	implicit none
	integer,dimension(4)::timestep
	real(kind=8),dimension(4)::yerror, zerror
	real(kind=8)::ordery, orderz, theta1,theta2
	integer::i,xrange
	print*,"xrange = "
	read*,xrange
	print*,"theta1 = "
	read*,theta1
	print*,"theta2 = "
	read*,theta2
	timestep = (/32,64,128, 256/)
	do i = 1, 4
		call li34(timestep(i),xrange,yerror(i),zerror(i), theta1, theta2)
	end do 
	print*,"yerror=",yerror
	print*,"zerror=",zerror
	call leastsquare(timestep, yerror, ordery)
	print*,"ordery=",ordery
	call leastsquare(timestep, zerror, orderz)
	print*,"orderz=",orderz

end program