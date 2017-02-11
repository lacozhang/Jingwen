module kronecker
	implicit none
contains
	subroutine kron(a,b,c)
		real(kind = 8), dimension(:,:), intent(inout):: a, b
		real(kind = 8), dimension(:,:), intent(out):: c
		integer, parameter:: sizea, sizeb 
		integer:: i, j, start_point, end_point
		sizea = size(a,1)
		sizeb = size(b,1)
		i=0
		j=0
		do 
			


	end subroutine
end module