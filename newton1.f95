program newton1
	use newton
	implicit none
	real(kind=8),dimension(2):: x
	interface
		pure subroutine f(x,fy)
			real(kind=8),dimension(:),intent(in):: x
			real(kind=8),dimension(:),intent(out):: fy
		end subroutine f
		pure subroutine f1(x,f1y)
			real(kind=8),dimension(:,:),intent(in)::x
			real(kind=8),dimension(:), intent(out)::f1y
		end subroutine f1
	end interface
	x= (/0.0,0.0/)
	call newton_iteration(f, f1, x)
	
end program
pure function f(x)	
	implicit none
	real(kind=8),dimension(:), allocatable::f
	real(kind=8),dimension(:),intent(in) ::x
	f(1)= -0.5*(2*x(1)-cos(x(1))+sin(x(2)))
	f(2)= -0.5*(2*x(2)-sin(x(1))-cos(x(2)))
end function f
pure function f1(x)
	implicit none
	real(kind=8),dimension(:,:), allocatable::f1
	real(kind=8),dimension(:), intent(in)::x
	f1(1,1)= 0.5*(2+sin(x(1)))
	f1(1,2)= 0.5*(cos(x(2)))
	f1(2,1)= 0.5*(-cos(x(1)))
	f1(2,2)= 0.5*(2+sin(x(2)))
end function f1