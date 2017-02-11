program example31
	use example
	use leastsq
	implicit none
	integer,dimension(5)::timestep
	real(kind=8),dimension(5)::yerror, zerror
	real(kind=8)::ordery, orderz
	integer::i
	timestep = (/4, 8, 16, 32, 64/)
	do i = 1, 5
		call li31(timestep(i),yerror(i),zerror(i))
	end do 
	print*,"yerror=",yerror
	print*,"zerror=",zerror
	call leastsquare(timestep, yerror, ordery)
	print*,"ordery=",ordery
	call leastsquare(timestep, zerror, orderz)
	print*,"orderz=",orderz

end program