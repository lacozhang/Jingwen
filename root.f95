program root
	implicit none
	real:: a, b, c, d
	complex:: x1, x2
	read*, a, b, c
	d= b**2-4*a*c
	if (d > 0) then
		x1=(-b-sqrt(d))/(2*a)
		x2=(-b+sqrt(d))/(2*a)
		print*,"x1:", x1, "x2:", x2
	else if (d == 0) then
		x1= -b/(2*a)
		print*,"x1:", x1, "x2:", x1
	else
		x1= cmplx(-b/(2*a), sqrt(-d)/(2*a))
		x2= cmplx(-b/(2*a), -sqrt(-d)/(2*a))
		print*,"x1:", x1, "x2:", x2
	end if
	
	
end program