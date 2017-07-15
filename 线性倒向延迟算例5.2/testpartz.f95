program test
	use main
	implicit none
	
	integer,dimension(5)::timestep
	real(kind=8),dimension(5)::yerror, zerror
	real(kind=8)::ordery, orderz,temp
	integer::i
	temp = 0
!	read*,theta 
	timestep = [ 8,16,32,64,128]
	do i = 1,5
		call li52(timestep(i), yerror(i), zerror(i))
		print*, "y error", yerror(i)
		print*, "z error", zerror(i)
	end do 
	print*,"yerror=",yerror
	print*,"zerror=",zerror
	print*,"0.0= ",temp
	call leastsquare(timestep, zerror, orderz)
	print*,"orderz=",orderz
end program test