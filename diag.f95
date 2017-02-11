module diag
	implicit none
	private
	public :: diagin
contains
	subroutine diagin(a,b,c)
		real(kind=8), dimension(:), intent(in):: a
		real(kind=8), dimension(:,:), intent(inout):: b
		integer, intent(in):: c 
		integer:: i
		integer:: sizea,sizeanew
		sizea = size(a)
		select case (c)
		case (0)
			do i = 1,sizea
				b(i, 1:sizea) = 0
				b(i, i) = a(i)
			end do
		case (1:)
			sizeanew = sizea + c
			do i = 1,sizea
				b(i, 1:sizeanew) = 0
				b(i, i+c) = a(i)
			end do
			b(sizea:sizeanew,:) = 0 
		case (:-1)
			sizeanew = size(a) - c
			b(1:-c,1:sizea) = 0
			do i = -c, sizea
			    b(i, 1:sizeanew) = 0
			    b(i, i+c) = 0 
			end do 
		end select
	end subroutine
end module