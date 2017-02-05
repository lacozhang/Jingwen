program mj
	implicit none
	real :: mianji,b,zhijing 
	real, parameter::pi=3.14159268
	read*, zhijing
	b=zhijing/2
	print*, b, pi
	mianji=pi*(b**2)
	print*, mianji
end program