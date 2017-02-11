program testpicard

	use nolinear_equation
	implicit none
	real(kind=8),dimension(2)::x,y
	real::t1,t2	
	interface
		pure subroutine function_value(x,y)
			real(kind=8),dimension(:),intent(in)::x
	   		real(kind=8),dimension(:),intent(out)::y
	   	end subroutine
	end interface
	call cpu_time(t1)
	x=(/0,0/)
	call picard_iteration(function_value,x)
	call cpu_time(t2)
	print*,"cpu time",t2-t1
end program

pure subroutine function_value(x,y)
	implicit none
	real(kind=8),dimension(:),intent(in)::x
	real(kind=8),dimension(:),intent(out)::y
	y(1)= 0.5*(cos(x(1))-sin(x(2)))
	y(2)= 0.5*(sin(x(1))+cos(x(2)))
end subroutine