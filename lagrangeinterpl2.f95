module interpolation
	implicit none
contains
	subroutine lagrangeinterpl2(x,y,x_aim)
		real(kind=8), dimension(:), intent(in):: x
		real(kind=8), dimension(:), intent(in):: y
		real(kind=8), intent(inout):: x_aim
		x_aim = y(1)*(x_aim-x(2))/(x(1)-x(2))+y(2)*(x_aim-x(1))/(x(2)-x(1))
	end subroutine
	subroutine newtoninterl(x, y, x_aim)
		real(kind=8), dimension(:), intent(in):: x
		real(kind=8), dimension(:), intent(in):: y
		real(kind=8), intent(inout):: x_aim
		real(kind=8), dimension(:,:):: z
		
