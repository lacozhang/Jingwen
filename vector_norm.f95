module norm
	implicit none
	private
	public::vector_norm
contains
function vector_norm(x)
	real(kind=8), dimension(:):: x
	integer:: sizex 
	real(kind=8):: vector_norm 
	real(kind=8):: vector_sum
	integer:: i
	sizex = size(x)
	vector_sum = 0
	do i = 1, sizex
		vector_sum = x(i)**2+vector_sum
	end do
	vector_norm = sqrt(vector_sum)
end function
end module norm