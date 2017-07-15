program test
	use main
	implicit none
	
	integer,dimension(3)::timestep
	real(kind=8),dimension(3)::yerror, zerror
	real(kind=8)::ordery, orderz, theta
	integer::i

!	read*,theta 
	timestep = [4,8,16]
	do i = 1,3
		call li52(timestep(i), yerror(i), zerror(i))
		print*, "y error", yerror(i)
		print*, "z error", zerror(i)
	end do 
	print*,"yerror=",yerror
	print*,"zerror=",zerror
	call leastsquare(timestep, yerror, ordery)
	print*,"ordery=",ordery
	call leastsquare(timestep, zerror, orderz)
	print*,"orderz=",orderz
end program test