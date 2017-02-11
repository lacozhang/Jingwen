subroutine thelengthofx(k, step, timesteplength,x)

	!x returns result
	!need bubble.f95 gausshermite.f95
	use gausshermite
	use bubblesort
	implicit none
	integer, intent(in):: k, step
	real(kind=8), intent(out)::x
	real(kind=8), dimension(:),allocatable::w, x0
	real(kind=8), dimension(:),allocatable::a
	real(kind=8), dimension(:),allocatable::xnext
	real(kind=8)::x_max, x_min,timesteplength
	integer::i
	allocate(w(k))
	allocate(a(k))
	allocate(xnext(k))
	allocate(x0(k))
	x0(:) = 1
	x_max = 0
	call gauss_hermite(k, w, a)
	print*,"a",a
	do i = 1, step
		xnext = x_max * x0 + (sqrt(2.0 * timesteplength) * a)
		print*,"xnext=",xnext
		call bubble_sort(xnext,1,k)
		x_min = abs(xnext(1))
		x_max = abs(xnext(k))
		print*,"x_min=", x_min, "x_max", x_max
		if (x_max < x_min) then
			x_max = x_min
		end if
		x = x_max
	end do 
		x = x + 5
		print*,"x", x
	deallocate(w)
	deallocate(a)
	deallocate(xnext)
	deallocate(x0)
end subroutine