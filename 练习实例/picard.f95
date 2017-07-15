program picard
	implicit none
	real, dimension(2):: x, xnext
	real:: distance, times
	interface 
		function vector_norm(x)
		real, dimension(:):: x
		real:: vector_norm 
		end function
	end interface
	xnext=(/0, 0 /)
	times=0
	do
		times=1+times
		x=xnext
		xnext(1)=0.5*(cos(x(1))-sin(x(2)))
		xnext(2)=0.5*(sin(x(1))+cos(x(2)))
		distance= vector_norm(x-xnext)
		print*, distance
		if(distance < 1e-8) then
			exit
		end if
	end do	
	print*, x, times
end program
function vector_norm(x)
	implicit none
	real, dimension(:):: x
	integer:: sizex 
	real:: vector_norm 
	real:: vector_sum
	integer:: i
	sizex = size(x)
	vector_sum = 0
	do i = 1, sizex
		vector_sum = x(i)**2+vector_sum
	end do
	vector_norm = sqrt(vector_sum)
end function