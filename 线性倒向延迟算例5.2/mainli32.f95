module main
!-llapck -lblas
	implicit none
	private
	public::true_zvalue,gauss_hermite, thelengthofx, bubble_sort, newtoninterpl, leastsquare,  li52

contains
	pure subroutine true_zvalue(t,x,z)
		real(kind=8),intent(inout)::x,t
		real(kind=8),intent(out)::z
		z =2.0/(exp(1.0-t)+exp(t-1.0))
	end subroutine
	
	subroutine gauss_hermite(n,w,x)
	!-llapck -lblas
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
		
		cm(:,:) = 0
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
	
	subroutine thelengthofx(k, step, timesteplength, xint)
	!different equation have different thelengthofx!!!
	!k is gaussherimite's points
	!step is the scheme's step
	!xint returns result which is the length of x in positive
		integer, intent(in):: k, step
		integer, intent(out):: xint
		real(kind=8),intent(inout):: timesteplength
		real(kind=8), dimension(:),allocatable::w
		real(kind=8), dimension(:),allocatable::a
		real(kind=8), dimension(:),allocatable::xnext
		real(kind=8)::x_max, x_temp, x_prev
		integer::i, j
		allocate(w(k))
		allocate(a(k))
		x_max = 0 !(x0 = 0)
		call gauss_hermite(k, w, a)
		do i = 1, step
			!print*,"i=",i
			x_prev = x_max
			x_max = 0
			do j=1,k
				x_temp = x_prev + sqrt(2*timesteplength)*a(j)
				if(abs(x_temp) > x_max) then
					x_max = x_temp
				end if
			end do
		end do 
		xint = floor(x_max) + 3
		!print*, xint
		deallocate(w)
		deallocate(a)
	end subroutine
	
	subroutine newtoninterpl(x, y, x_aim, y_aim)
		!finall y_aim which is our need
		real(kind=8), dimension(:), intent(inout):: x
		real(kind=8), dimension(:), intent(in):: y
		real(kind=8), intent(inout):: x_aim
		real(kind=8), intent(out):: y_aim 
		real(kind=8), dimension(:,:), allocatable:: z
		integer:: sizex, temp1, temp2,  w
		real::n
		sizex = size(x)
		allocate(z(sizex,sizex))
		z(:,:) = 0
		n = 0
		z(:,1) = y
		do temp2 = 2,sizex
			do temp1 = temp2, sizex
				z(temp1,temp2) = (z(temp1-1,temp2-1)-z(temp1,temp2-1))/(x(temp1-(temp2-1))-x(temp1))
			end do 
		end do
		x(2:sizex) = x(1:(sizex-1))
		x(1) = x_aim -1
		x = x_aim - x
		do temp1 = 2, sizex
			x(temp1) = x(temp1-1)*x(temp1)
		end do
		do temp1 = 1, sizex
			n = n+z(temp1,temp1)*x(temp1)
		end do 
		y_aim = n 
		deallocate(z)
	end subroutine
	
	subroutine leastsquare(timestep, error, order)
	!sizeerror >= 2
		real(kind=8), dimension(:),intent(inout)::error
		integer, dimension(:),intent(inout)::timestep
		real(kind=8), dimension(:), allocatable::timesteplength
		real(kind=8), intent(out)::order
		real(kind=8), dimension(:,:), allocatable:: a
		real(kind=8), dimension(2)::s
		real(kind=8), dimension(:), allocatable:: work
		integer::sizeerror, info, rank, lwork
		sizeerror = size(error)
		lwork = 6 + max(4, sizeerror)
		allocate(timesteplength(sizeerror))
		allocate(a(sizeerror,2))
		allocate(work(lwork))
		timesteplength = 1.0/(timestep)
		error = log(error)
		timesteplength = log(timesteplength)
		a = 1
		a(:,2) = timesteplength
		call dgels('N', sizeerror, 2, 1, a, sizeerror, error, sizeerror, work, lwork, info)
		order = error(2)
		deallocate(timesteplength)
		deallocate(a)
		deallocate(work)
	end subroutine
	
	
	subroutine li52(timestep, yerror, zerror)
	!cranck nicolson method
	!-llapck -lblas
		
		integer, intent(inout):: timestep
		real(kind=8), intent(out)::yerror, zerror


		real(kind=8), dimension(:,:), allocatable::y, z, prev_y
		real(kind=8), dimension(:),allocatable:: t, x
		real(kind=8):: timesteplength, xsteplength, yb
		real(kind=8):: prev_y1, z_int,z0
! 	
		integer::  i, timeindex, xindex, xstep, xrange, k, x_point1, x_point2, j, cycletimes
		integer:: overflow,times
		real(kind=8):: b, tau, f,  expectyw, expecty, pi, diff,  x_point , yture, zture
		real(kind=8),dimension(:),allocatable:: w, a
		real(kind=8), dimension(2)::yinterpl,  xinterpl 
		logical:: keep, firstcycle

		pi = 4*atan(1.0)

		timesteplength = (1.0)/timestep
		print*, "time step size  : ", timesteplength
		
		
	! GH
		k = 8
		print*, "number of steps : ", timestep
		print*, "gauss hermite   : ", k

		allocate(w(k))
		allocate(a(k))
		call gauss_hermite(k, w, a)
		!print*,"w=:",w
		!print*,"a=:",a
	! end GH
		
		call thelengthofx(k, timestep, timesteplength, xrange)
		xsteplength = timesteplength
		xstep = ceiling((xrange+xsteplength)/xsteplength)
		print*,"xrange:",xrange
		print*,"xstep:",xstep
		allocate(t(0:timestep))
		allocate(x(-xstep:xstep))

	!	t sapce

		do timeindex = 0, timestep
			t(timeindex) = timeindex 
		end do 
		t = t * timesteplength
	! x space 
		x(0) = 0
		do xindex = -xstep,xstep
			x(xindex) = x(0)+xindex*xsteplength
		end do
		
		allocate(y(0:timestep,-xstep:xstep))
		allocate(z(0:timestep,-xstep:xstep))
		allocate(prev_y(0:timestep, -xstep:xstep))

		y(:,:)=0
		z(:,:)=0
		prev_y(:,:)=0
		
		! last line of Y matrix (timestep)
		do xindex = -xstep, xstep
			y(timestep,xindex)= x(xindex)
		end do 
		!print*,"last Y",y(timestep,:)
		
	
		keep = .true.
		firstcycle = .true.
		cycletimes = 0
		do
			!cycletimes = cycletimes + 1
			!print*,"cycletimes:",cycletimes
			do timeindex = timestep-1, 0, -1
				do xindex = -xstep, xstep
					expectyw = 0
					expecty = 0 
					do i = 1, k
					
						x_point= x(xindex)+sqrt(2*timesteplength)*a(i)
						x_point1 = floor(x_point/xsteplength)
						x_point2 = x_point1+1
						overflow = 0

						if (x_point1 < -xstep) then
							x_point1 = -xstep
							x_point2 = -xstep + 1
							overflow = 1
						end if
						
						if (x_point2 > xstep ) then 
							x_point2 = xstep
							x_point1 = xstep - 1
							overflow = 1
						end if

						if ((x(x_point1) > x_point) .and. (overflow == 0)) then
							print*, "compute error"
							print*, "low val", x(x_point1)
							print*, "point1", x(x_point1)
							print*, "val ", x_point
						end if

						if ((x(x_point2) < x_point) .and. (overflow == 0)) then
							print*, "compute error"
							print*, "high val", x(x_point2)
							print*, "point1", x(x_point1)
							print*, "val ", x_point
						end if
						!start to interpolate
						yinterpl = (/y(timeindex+1,x_point1),y(timeindex+1,x_point2)/)
						xinterpl = (/x(x_point1),x(x_point2)/)
						call newtoninterpl(xinterpl,yinterpl, x_point, yb)

						!call true_yvalue(t(timeindex+1), x_point, yb)
						
						expectyw = expectyw + w(i)*yb*sqrt(2*timesteplength)*a(i)
						expecty = expecty + w(i)*yb	
					end do 

					expectyw = expectyw/sqrt(pi)
					expecty = expecty/sqrt(pi)
					z(timeindex,xindex) = (expectyw)/(timesteplength)
					!call true_zvalue(t(timeindex), x(xindex),z(timeindex, xindex))
					y(timeindex,xindex) = expecty
					do j=1, timeindex-1
						y(timeindex, xindex) = y(timeindex, xindex) - 0.5*timesteplength*(prev_y(j, xindex) + prev_y(j+1, xindex))*timesteplength
					end do
					!call true_yvalue(t(timeindex), x(xindex), y(timeindex, xindex))
				end do 
			end do
			!print*,"based on z ", z0
			!print*,"z(0,0)=",z(0,0)
			!print*,"y(0,0)=",y(0,0)
			if (firstcycle) then
				prev_y1 = y(0,0)
				firstcycle = .false.
			else if(abs(prev_y1 - y(0,0)) < 1e-8) then
				print*, "exit criteria met"
				exit
			end if
			prev_y1 = y(0,0)

			prev_y = y
		end do
		yture = 0
		call true_zvalue(t(0),x(0),zture)
		print*,"yture=:",yture,"zture=:",zture
		print*, "y(0,0) =:", y(0,0), "z(0,0) =:", z(0,0)
	
		zerror = abs(zture - z(0,0))
		yerror = abs(yture - y(0,0))
		!print*,"yerror", yerror
		!print*, "zerror",zerror
		deallocate(w)
		deallocate(a)
		deallocate(y)
		deallocate(z)
		deallocate(t)
		deallocate(x)
		deallocate(prev_y)
	end subroutine
end module main