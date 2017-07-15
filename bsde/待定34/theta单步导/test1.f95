program test1
	use thetaderivate
	implicit none
	real(kind=8)::theta,  error
	integer::timestep
	real(kind=8)::y_aim

	timestep = 10
	print*,"theta="
	read*,theta
	
	call thetaderi(timestep, theta, y_aim, error)
	print*,"the timestep",timestep
	print*,"y_aim=",y_aim
	print*,"error",error
			

end program