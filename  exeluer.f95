program exeluer
	implicit none
	interface
	pure subroutine function_value(t,y,f)
		real(kind=8),intent(inout)::t
		real(kind=8),intent(inout)::y
		real(kind=8),intent(out)::f
	end subroutine
	end interface
	real(kind=8),dimension(11)::y,t
	real(kind=8)::h,f
	integer ::i
	h=0.1
	y(1)=0
	t(1)=0
	do i=2,11,h
		t(i)=(i-1)*h
		call function_value(t(i-1),y(i-1),f)
		y(i)=y(i-1)+h*f
		print*,"t",t(i),"y",y(i)
	end do
end program
pure subroutine function_value(t,y,f)
	real(kind=8),intent(inout)::t
	real(kind=8),intent(inout)::y
	real(kind=8),intent(out)::f
	f=(2/(t+1))*y+((t+1)**2)*exp(t)
end subroutine
