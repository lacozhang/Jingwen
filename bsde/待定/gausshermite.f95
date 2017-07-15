module gausshermite
	!need -llapack -lblas	
	
	implicit none
	private
	public::gauss_hermite, thelengthofx, bubble_sort
contains
	subroutine gauss_hermite(n,w,x)
	!w is weight 
		integer,intent(in):: n
		real(kind=8),dimension(:),intent(out):: w
		real(kind=8),dimension(:),intent(out):: x
		real(kind=8),dimension(:),allocatable:: a 
		real(kind=8),dimension(:,:),allocatable:: cm
		real(kind=8),dimension(10*n)::work
		real(kind=8),dimension(n)::l
		integer:: temp
		integer:: lwork, info
		allocate(a(n-1))
		allocate(cm(n,n))

		do temp = 1,n-1
			a(temp) = sqrt(temp*0.5)

		end do 
		
		do temp = 1,n-1
			cm(temp,:) = 0
			cm(temp,temp+1) = a(temp)
			
		end do
		lwork = 10*n
		call dsyev('V','U', n, cm, n, l, work, lwork, info)
		x = l

		w = sqrt(4*atan(1.0))*cm(1,:)**2
		deallocate(a)
		deallocate(cm)	
	end subroutine
	
	
	
	
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
	
	
	
	
	subroutine thelengthofx(k, step, timesteplength,xint)
	!xint returns result which is the length of x in positive
		integer, intent(in):: k, step
		integer, intent(out):: xint
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
		!print*,"a",a
		do i = 1, step
			xnext = x_max * x0 + (sqrt(2.0 * timesteplength) * a)
			!print*,"xnext=",xnext
			call bubble_sort(xnext,1,k)
			x_min = abs(xnext(1))
			x_max = abs(xnext(k))
			!print*,"x_min=", x_min, "x_max", x_max
			if (x_max < x_min) then
				x_max = x_min
			end if
		end do 
		xint = floor(x_max) + 3
		deallocate(w)
		deallocate(a)
		deallocate(xnext)
		deallocate(x0)
	end subroutine
end module