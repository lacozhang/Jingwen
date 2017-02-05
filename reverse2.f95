program reverse2
	implicit none
	real, dimension(10)::array
	integer:: low, high
	array= (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 0 /)
	print*,"the old array is:",array
	print*,"please input 2 numbers:"
	print*, "the lower:"
	read*, low
	print*, "the higher:"
	read*, high
	array(low:high)=array(high:low:-1)
	print*, "the new array is:",array
end program