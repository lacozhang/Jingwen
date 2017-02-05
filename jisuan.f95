program jisuan
	implicit none
	integer:: m, u
	interface
		recursive function number1(n) result(need)
			integer:: n, need
		end function number1
	end interface
	read*,m
	u= number1(m)
	print*, u
end program

recursive function number1(n) result(need)
	implicit none
	integer:: n, need
	select case (n)
		case (1)
			need= 1
		case (2)
			need= 1
		case default
			need= number1(n-2)+number1(n-1)
	end select
			
end function number1
