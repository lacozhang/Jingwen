program exeluer
	implicit none
	interface
		pure subroutine function_value(t,y,f)
			real(kind=8),intent(inout)::t
			real(kind=8),intent(inout)::y
			real(kind=8),intent(out)::f
		end subroutine
		pure subroutine refunction_value(t,y)
			real(kind=8),intent(in)::t
			real(kind=8),intent(out)::y
		end subroutine
	end interface
	integer, parameter:: num_steps=11
	real(kind=8),dimension(:),allocatable::y,t,diff
	real(kind=8)::h,f,ref
	integer ::i
	allocate(y(num_steps))
	allocate(t(num_steps))
	allocate(diff(num_steps))
	h=0.1
	y(1)=1
	t(1)=0
	diff(1)=0
	do i=2,num_steps
		t(i)=(i-1)*h
		call function_value(t(i-1),y(i-1),f)
		y(i)=y(i-1)+h*f
		call refunction_value(t(i),ref)
		diff(i)=abs(y(i)-ref)
		print"(a,f3.1)","t=",t(i)
		print"(a, f6.4)","y=",y(i)
		print*,"error=",diff(i)
	end do
	deallocate(y)
	deallocate(t)
	deallocate(diff)
end program
pure subroutine function_value(t,y,f)
	real(kind=8),intent(inout)::t
	real(kind=8),intent(inout)::y
	real(kind=8),intent(out)::f
	f=(2/(t+1))*y+((t+1)**2)*exp(t)
end subroutine
pure subroutine refunction_value(t,y)
	real(kind=8),intent(in)::t
	real(kind=8),intent(out)::y
	y=((t+1)**2)*exp(t)
end subroutine