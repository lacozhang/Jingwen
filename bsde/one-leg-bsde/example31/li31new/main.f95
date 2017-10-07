module main
	implicit none
	private
	public:: gauss_hermite, thelengthofx, bubble_sort, function3, newtoninterpl, leastsquare, li31


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
	subroutine thelengthofx(k, step, timesteplength,xint)
	!xint returns result which is the length of x in positive
		integer, intent(in):: k, step
		integer, intent(out):: xint
		real(kind=8),intent(inout):: timesteplength
		real(kind=8), dimension(:),allocatable::w, x0
		real(kind=8), dimension(:),allocatable::a
		real(kind=8), dimension(:),allocatable::xnext
		real(kind=8)::x_max, x_min
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
	pure subroutine function3(y,f)
		real(kind=8),intent(inout)::y
		real(kind=8),intent(out)::f
		f = -y**3+2.5*(y**2)-1.5*y	
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
		integer, dimension(:)::timestep
		real(kind=8), dimension(:), allocatable::timesteplength
		real(kind=8), intent(out)::order
		real(kind=8), dimension(:,:), allocatable:: a
		real(kind=8), dimension(2)::s
		real(kind=8), dimension(:), allocatable:: work
		integer::sizeerror, info, rank, lwork
		sizeerror = size(error)
		lwork = 10*sizeerror
		allocate(timesteplength(sizeerror))
		allocate(a(sizeerror,2))
		allocate(work(lwork))
		timesteplength = 1.0/(timestep)
		error = log(error)
		timesteplength = log(timesteplength)
		a(:,:)= 1
		a(:,2) = timesteplength
		!print*, "size of error", sizeerror
		!print*,"a = ",a
		!print*,"error",error
!		call dgelss(sizeerror, 2, 1, a, sizeerror, error, sizeerror, s, 0.01, rank, work, lwork, info)
		call dgels('N', sizeerror, 2, 1, a, sizeerror, error, sizeerror, work, lwork, info)
		!print*, "error", error
		order = error(2)
		deallocate(timesteplength)
		deallocate(a)
		deallocate(work)
	end subroutine
	subroutine li31(timestep, yerror, zerror, theta2)
	!moudle interpolation gaussherimite function31
	!-llapck -lblas
		
		integer, intent(inout):: timestep
		real(kind=8), intent(out)::yerror, zerror
		real(kind=8), intent(inout)::theta2

		real(kind=8), dimension(:,:), allocatable:: y, z
		real(kind=8),dimension(:),allocatable:: t, x
		real(kind=8):: timesteplength, xsteplength, yb, zb

! 	
		integer::  i, timeindex, xindex, xstep, k, x_pointlow, x_pointhigh
		integer:: overflow
		real(kind=8):: f, expecty, expectf, expectyw,expectfw, expectz, pi, diff, y0, y1, x_point 
		real(kind=8),dimension(:),allocatable:: w, a
		real(kind=8), dimension(2)::yinterpl, xinterpl,zinterpl
		real(kind=8), parameter:: yture=0.5, zture=0.25
	

		pi = 4*atan(1.0)

		timesteplength = (1.0)/timestep
		print*, "time step size  : ", timesteplength
		
		
	! GH
		k = 10
		print*, "number of steps : ", timestep
		print*, "gauss hermite   : ", k

		allocate(w(k))
		allocate(a(k))
		call gauss_hermite(k, w, a)
	! end GH

		call thelengthofx(k, timestep, timesteplength, xstep)
		xsteplength = sqrt(timesteplength**3)
		print*, "x radius        : ", xstep
		xstep = xstep/xsteplength

		allocate(t(0:timestep))
		allocate(x(-xstep:xstep))
 
	!
		t(0) = 0
		do timeindex = 1, timestep
			t(timeindex) = t(timeindex-1) + timesteplength
		end do 

	! x space range
		x(0) = 0
		do xindex = -xstep, xstep
			x(xindex) = xindex
		end do
		x = x * xsteplength

		allocate(y(0:timestep,-xstep:xstep))
		allocate(z(0:timestep,-xstep:xstep))

	! last line of Y matrix
		do xindex = -xstep, xstep
			y(timestep,xindex) = exp(x(xindex)+t(timestep))/(exp(x(xindex)+t(timestep))+1)
		end do 
	! last but one line of Y matrix
		do xindex = -xstep, xstep
			expecty = 0
			expectf = 0
			do i = 1, k 
				yb = exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))&
				/(exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))+1)
				call function3(yb,f)
				expecty = expecty+w(i)*yb 
				expectf = expectf+w(i)*f
			end do 
			expecty = expecty/sqrt(pi)
			expectf = expectf/sqrt(pi)
			y0 = y(timestep,xindex)
			do 
				call function3(y0,f)
				y1 = expecty + 0.5*timesteplength*expectf + 0.5 *timesteplength*f
				diff = abs(y1-y0)
				if (diff <= 1e-10) then 
					exit
				end if
				y0 = y1
			end do 
			y(timestep-1,xindex) = y1
		end do 
		!print*,"y(timestep-1,:)", y(timestep-1,:)

	!last  line of Z matrix
		do xindex = -xstep, xstep
			z(timestep,xindex) = exp(x(xindex)+t(timestep))/((exp(x(xindex)+t(timestep))+1)**2)
		end do 
		
	! last but one line of Z matrix	
		do xindex = -xstep, xstep
			expectyw = 0
			expectfw = 0
			expectz = 0
			do i = 1, k
				yb = exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))/(exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))+1)
				zb = exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))/((exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))+1)**2)
				call function3(yb,f)
				expectyw = expectyw + w(i)*yb*sqrt(2*timesteplength)*a(i)
				expectfw = expectfw + w(i)*f*sqrt(2*timesteplength)*a(i)
				expectz = expectz + w(i)*zb 
			end do 
			expectyw = expectyw/sqrt(pi)
			expectfw = expectfw/sqrt(pi)
			expectz = expectz/sqrt(pi)
			z(timestep-1,xindex) = (expectyw+theta2*timesteplength*expectfw-(1-theta2)*timesteplength*expectz)/(theta2*timesteplength)
		end do 
		!print*,"z(timestep-1,:)=:",z(timestep-1,:)
! when 0<= n <=N-2

		do timeindex = timestep-2, 0, -1
			do xindex = -xstep, xstep
				expecty = 0
				expectf = 0
				expectyw = 0
				expectfw = 0
				expectz = 0				
				!find the interpolation point
				do i = 1, k
					x_point = x(xindex)+sqrt(2*timesteplength)*a(i)

					x_pointlow = floor(x_point/xsteplength)
					x_pointhigh = x_pointlow+1
					overflow = 0

					if (x_pointhigh <= -xstep) then
						x_pointhigh = -xstep + 1
						x_pointlow = -xstep
						overflow = 1
					end if
					if (x_pointlow >= xstep ) then 
						x_pointlow = xstep - 1
						x_pointhigh = xstep 
						overflow = 1
					end if

					if ((x(x_pointlow) > x_point) .and. (overflow == 0)) then
						print*, "compute error"
						print*, "low val", x(x_pointlow)
						print*, "val ", x_point
					end if

					if ((x(x_pointhigh) < x_point) .and. (overflow == 0)) then
						print*, "compute error"
						print*, "high val", x(x_pointhigh)
						print*, "val ", x_point
					end if
					!start to interpolate
					yinterpl = (/y(timeindex+1,x_pointlow),y(timeindex+1,x_pointhigh)/)
					xinterpl = (/x(x_pointlow),x(x_pointhigh)/)
					call newtoninterpl(xinterpl,yinterpl, x_point, yb)	
					x_point = x(xindex)+sqrt(2*timesteplength)*a(i)		
					zinterpl = (/z(timeindex+1,x_pointlow),z(timeindex+1,x_pointhigh)/)
					xinterpl = (/x(x_pointlow),x(x_pointhigh)/)					
					call newtoninterpl(xinterpl, zinterpl, x_point, zb)
				
					call function3(yb,f)
					expecty = expecty+w(i)*yb 
					expectf = expectf+w(i)*f
					expectyw = expectyw+w(i)*yb*sqrt(2*timesteplength)*a(i)
					expectfw = expectfw +w(i)*f*sqrt(2*timesteplength)*a(i)
					expectz = expectz + w(i)*zb 
				end do 
				expecty = expecty/sqrt(pi)
				expectf = expectf/sqrt(pi)
				expectyw = expectyw/sqrt(pi)
				expectfw = expectfw/sqrt(pi)
				expectz = expectz/sqrt(pi)
	
				z(timeindex,xindex) = (expectyw+theta2*timesteplength*expectfw-(1-theta2)*timesteplength*expectz)/(theta2*timesteplength)

				y0 = y(timeindex+1,xindex)
				do 
					call function3(y0,f)
					y1 = expecty + 0.5*timesteplength*expectf + 0.5*timesteplength*f
					diff=abs(y0-y1)
					if (diff <= 1e-10) then 
						exit
					end if
					y0 = y1
				end do
				y(timeindex,xindex) = y1
			end do 
		end do 
		!print*,"z(0,0)=",z(0,0)
		!print*,"y(0,0)=",y(0,0)
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
	end subroutine
end module main
