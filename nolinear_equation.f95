module nolinear_equation
	implicit none
	private
	public::newton_iteration, picard_iteration
contains
	subroutine newton_iteration(function_value, function_jacobi, x)
		use norm
		interface
			pure subroutine  function_value(x,fv)
				real(kind=8), dimension(:), intent(in):: x
				real(kind=8), dimension(:), intent(out):: fv
		    end subroutine
			pure  subroutine function_jacobi(x,fj)
				real(kind=8), dimension(:), intent(in):: x
				real(kind=8), dimension(:,:), intent(out):: fj 
			end subroutine
		end interface
		real(kind=8), dimension(:):: x
		real(kind=8), dimension(:), allocatable:: f, pivot, deltax, xnext, fv
		real(kind=8), dimension(:,:), allocatable::f1, fj
		integer :: sizex
		real(kind=8) :: diff, times, ok
		sizex=size(x)
		allocate( f(sizex) )
		allocate(pivot(sizex))
		allocate(deltax(sizex))
		allocate(f1(sizex,sizex))
		allocate(xnext(sizex))
		allocate(fv(sizex))
		allocate(fj(sizex,sizex))
		times = 0
		xnext = x		
		do 
			times = 1 + times
			x = xnext
			call function_value(x, fv)			
			f = fv		
			call function_jacobi(x, fj)			
			f1= fj			
			call DGESV(sizex, 1, f1, sizex, pivot, f, sizex, ok)
			deltax= f			
			diff= vector_norm(deltax)			
			xnext=deltax+x
			print*, "times", times,"x",xnext
			if (diff <= 1e-8) then
				exit
			end if
		end do
		x= xnext
		print*,"x",x
		deallocate(f)
		deallocate(pivot)
		deallocate(deltax)
		deallocate(f1)
		deallocate(xnext)
		deallocate(fv)
		deallocate(fj)
 	end subroutine 
	
	subroutine picard_iteration(function_value,x)
		use norm
		real(kind=8),dimension(:), intent(inout)::x
		interface 
			pure subroutine function_value(x,y)
				real(kind=8), dimension(:),intent(in)::x
				real(kind=8), dimension(:),intent(out)::y
			end subroutine
		end interface
		real(kind=8),dimension(:),allocatable::xnext,ynext
		integer::times,sizex
		real(kind=8),dimension(:),allocatable::diff
		real(kind=8)::deltax
		sizex=size(x)
		allocate(xnext(sizex))
		allocate(ynext(sizex))
		allocate(diff(sizex))
		times=0
		print*,"x0",x
		do
			times=1+times
			call function_value(x,xnext)
			diff= x-xnext
			deltax= vector_norm(diff)
			x=xnext
			print*,"times",times
			print*,"x",x
			if (deltax <= 1e-8)then
				exit
			end if 
		end do
		deallocate(xnext)
		deallocate(ynext)
		deallocate(diff)
	end subroutine
end module