program ff
	implicit none
	integer:: limit, sum1, sum2,sum3, i
	read*, limit
	sum1=1
	sum2=1
	
	select case (limit)
	case (:0) 
		print*, "sorry!byebye!*_*"
	case (1) 
		print*,1
	case (2) 
		print*,1
	case default
		do i=1,limit-2
			sum3=sum1+sum2
			sum1=sum2
			sum2=sum3			
		end do
		print*, sum3
	end select
end program