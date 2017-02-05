program re
	implicit none
	integer, dimension(10):: array
	integer:: i, low, high, differece, temp, h, o
	do i=1,10
		array(i)=i**2
	end do
	print*,array
	read*, low, high
	differece = low - high
	if (differece == 0) then
		print*,"no changed"
	else if (differece > 0) then
		temp = low
		low = high
		high = temp
	end if
	print*,"low:", low, "high:", high
	do
		if (low >= high) then 
			exit
		end if 
		temp= array(low)
		array(low)= array(high)
		array(high)= temp
		low = low + 1
		high = high - 1
	end do 
	print*,array
end program