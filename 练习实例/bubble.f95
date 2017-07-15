module bubblesort
	implicit none
	private
	public::bubble_sort
contains
	recursive subroutine bubble_sort(array, low, high)
	!low : high means the range of array needed resort
	!array returns the ascending array
		real(kind=8), dimension(:), intent(inout):: array
		integer , intent(in):: low, high
		real(kind=8) :: temp
		integer :: j, i
		if (low < high) then
			do i= low, high-1
				if (array(i) > array(i+1)) then
					temp= array(i)
					array(i)= array(i+1)
					array(i+1)= temp
				end if
			end do
			call bubble_sort(array, low, high-1)
		end if 
	
	end subroutine bubble_sort 

end module bubblesort