program test
	use thetaderivate
	implicit none
	real(kind=8)::theta, order
	integer,dimension(2)::timestep
	real(kind=8),dimension(:),allocatable::y_aim, error
	integer::i
	allocate(y_aim(2))
	allocate(error(2))
	timestep = (/10,20/)
	print*,"theta="
	read*,theta
	do i = 1, 2
		print*,"theta=",theta 
		call thetaderi(timestep(i), theta,  y_aim(i),error(i))
		print*,"theta=",theta 
		print*,"the timestep",timestep(i)
		print*,"y_aim=",y_aim(i)
		print*,"error",error(i)
			
	end do 
	print*,"error",error
	call leastsquare(timestep, error, order)
	print*,"order",order
	deallocate(y_aim)
	deallocate(error)

end program