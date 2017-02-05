program bubble
	implicit none
	real, dimension(10):: array
	integer :: i, n
	interface 
		subroutine bubble_sort(array, low, high)
			real, dimension(10), intent(inout):: array
			integer , intent(in):: low, high
		end subroutine bubble_sort
	end interface
	do i= 1, 10
		read*, n
		array(i)=n
	end do 
	print*,array
	call bubble_sort(array, 1, 10)
	print*, array
end program

recursive subroutine bubble_sort(array, low, high)
	implicit none
	real, dimension(10), intent(inout):: array
	integer , intent(in):: low, high
	real :: temp
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