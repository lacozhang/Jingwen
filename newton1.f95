program newton1
	use newton
	implicit none
	real(kind=8),dimension(3):: x

	interface
		pure subroutine f(x,fy)
			real(kind=8),dimension(:),intent(in):: x
			real(kind=8),dimension(:),intent(out):: fy
		end subroutine f
		pure subroutine f1(x,f1y)
			real(kind=8),dimension(:),intent(in)::x
			real(kind=8),dimension(:,:), intent(out)::f1y
		end subroutine f1
	end interface
	x= (/0.0,0.0,0.0/)
	call newton_iteration(f, f1, x)
	
end program
pure subroutine f(x, fy)	
	implicit none
	real(kind=8),dimension(:),intent(in)::x
	real(kind=8),dimension(:),intent(out)::fy
	fy(1)= -1*(5*x(1)+2*x(2)+x(3)+12)
	fy(2)= -1*(-x(1)+4*x(2)+2*x(3)-20)
	fy(3)= -1*(2*x(1)-3*x(2)+10*x(3)-3)
end subroutine f
pure subroutine f1(x, f1y)
	implicit none
	real(kind=8),dimension(:), intent(in)::x
	real(kind=8),dimension(:,:), intent(out)::f1y
	f1y(1,1:3)= (/5,2,1/)
	f1y(2,1:3)= (/-1,4,2/)
	f1y(3,1:3)= (/2,-3,10/)
end subroutine f1