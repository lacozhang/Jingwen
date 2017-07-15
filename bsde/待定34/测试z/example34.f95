program example34
	use main
	implicit none
	integer,dimension(5)::timestep
	real(kind=8),dimension(5)::zerror
	real(kind=8)::ordery, orderz, theta2
	integer::i
	print*,"theta2 = "
	read*,theta2
	
	timestep = (/4,8,16, 32, 64/)
	do i = 1, 5
		call li34z(timestep(i),zerror(i), theta2)
	end do 

	print*,"zerror=",zerror
	call leastsquare(timestep, zerror, orderz)
	print*,"orderz=",orderz

end program