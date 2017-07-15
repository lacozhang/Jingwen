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
	pure subroutine function3(y,f,f1,f2)
		real(kind=8),intent(inout)::y
		real(kind=8),intent(out)::f,f1,f2
		f = -y**3+2.5*y**2-1.5*y
		f1 = -3*y**2+5*y-1.5
		f2 = -6*y+5
		
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
	subroutine li31(timestep, yerror, zerror, theta1, theta2)
		integer, intent(inout):: timestep
		real(kind=8), intent(out)::yerror, zerror
		real(kind=8), intent(inout)::theta1, theta2

		real(kind=8), dimension(:,:), allocatable:: y, z
		real(kind=8), dimension(:),allocatable:: t, x
		real(kind=8):: timesteplength, xsteplength, yb, zb

! 	
		integer::  i, timeindex, xindex, xrange, xstep, k, x_point1, x_point2,x_point3,x_point4,x_point5 
		integer:: overflow
		real(kind=8):: f, f1, f2
		real(kind=8):: expecty, expectf, expectf1f, expectf2zz, expectyw, expectz, expectfw,expectf1fw,expectf2zzw
		real(kind=8):: pi, diff, y0, y1, x_point 
		real(kind=8), dimension(:),allocatable:: w, a
		real(kind=8), dimension(5)::yinterpl, xinterpl,zinterpl
		real(kind=8), parameter:: yture=0.5, zture=0.25
		integer:: cycleindex
		
	

		pi = 4*atan(1.0)
		!print*,"pi=",pi
		print*, "number of steps : ", timestep
		timesteplength = (1.0)/timestep
		print*, "timesteplength  : ", timesteplength
		
		
	! GH
		k = 2
		
		!print*, "gauss hermite number  : ", k
		

		allocate(w(k))
		allocate(a(k))
		!print*,"w",w,"a",a
		call gauss_hermite(k, w, a)
		!print*,"w=",w,"a=",a
	! end GH
		!print*,"k=", k, "timestep=", timestep, "timesteplength=", timesteplength
		call thelengthofx(k, timestep, timesteplength, xrange)
		
		!print*,"xrange= ", xrange 
		xsteplength = timesteplength**(4.0/5)
		! be detemined by the order of x space and t space
		!print*, "xsteplength   = : ", xsteplength
		xstep = xrange/xsteplength
		!print*,"xstep=",xstep 
		
		
		allocate(t(0:timestep))
		allocate(x(-xstep:xstep))
 
	!
		t(0) = 0
		do timeindex = 1, timestep
			t(timeindex) = t(timeindex-1) + timesteplength
		end do 
		!print*,"time space:",t
		
	! x space range
		x(0) = 0
		do xindex = -xstep, xstep
			x(xindex) = xindex
		end do
		x = x * xsteplength
		!print*,"x space:",x
		
		
		allocate(y(0:timestep,-xstep:xstep))
		allocate(z(0:timestep,-xstep:xstep))
		
	!last  line of Z matrix
		do xindex = -xstep, xstep
			z(timestep,xindex) = exp(x(xindex)+t(timestep))/((exp(x(xindex)+t(timestep))+1)**2)
		end do 
		!print*,"z(timestep,:)=", z(timestep,:)
		
		
	! last but one line of z matrix	
		do xindex = -xstep ,xstep
			expectyw = 0
			expectz = 0
			expectfw = 0
			expectf1fw = 0
			expectf2zzw = 0
		!caculate the expect
			do i = 1,k
				yb = exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))/(exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))+1)
				zb = exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))/((exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))+1)**2)
				!print*,"yb =",yb,"zb =",zb 
				call function3(yb,f,f1,f2)
				expectyw = expectyw + w(i)*yb*sqrt(2*timesteplength)*a(i)
				expectz = expectz + w(i)*zb
				expectfw = expectfw + w(i)*f*sqrt(2*timesteplength)*a(i)
				expectf1fw = expectf1fw + w(i)*f1*f*sqrt(2*timesteplength)*a(i)
				expectf2zzw = expectf2zzw + w(i)*f2*(zb**2)*sqrt(2*timesteplength)*a(i)
			
			end do
			expectyw = expectyw/sqrt(pi)
			expectz = expectz/sqrt(pi)
			expectfw = expectfw/sqrt(pi)
			expectf1fw = expectf1fw/sqrt(pi)
			expectf2zzw = expectf2zzw/sqrt(pi)
			!print*,"expectyw:=",expectyw, "expectz=", expectz,"expectfw=",expectfw,"expectf1fw=",expectf1fw,"expectf2zzw=" ,expectf2zzw
			
		
			z(timestep-1,xindex) = (expectyw - (1-theta2)*timesteplength*expectz + theta2*timesteplength*expectfw & 
								&	-theta2*0.5*(timesteplength**2)*(-expectf1fw &
									+0.5*expectf2zzw))/(theta2*timesteplength)
									
			!print*,"z(timestep-1,xindex)=",z(timestep-1,xindex)
		end do 
	!print*,"z(timestep-1,:)=",z(timestep-1,:)
	
! last line of Y matrix
		do xindex = -xstep, xstep
			y(timestep,xindex) = exp(x(xindex)+t(timestep))/(exp(x(xindex)+t(timestep))+1)
		end do 
		!print*,"y(timestep,:)=",y(timestep,:)
! last but one line of y matrix 
	do xindex = -xstep ,xstep
		expecty = 0
		expectf = 0
		expectf1f = 0
		expectf2zz = 0
		do i = 1,k
			yb = exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))/(exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))+1)
			zb = exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))/((exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))+1)**2)
			call function3(yb,f,f1,f2)
			expecty = expecty +w(i)*yb
			expectf = expectf + w(i)*f
			expectf1f= expectf1f +w(i)*f1*f
			expectf2zz = expectf2zz +w(i)*f2*(zb**2)
		
		end do
		expecty = expecty/sqrt(pi)
		expectf = expectf/sqrt(pi)
		expectf1f = expectf1f/sqrt(pi)
		expectf2zz = expectf2zz/sqrt(pi)
		!print*,"expecty= ",expecty ,"expectf=",expectf,"expectf1f=",expectf1f,"expectf2zz=",expectf2zz
		y0 = y(timestep,xindex)
		!print*,"y0=",y0
		do 
			call function3(y0,f,f1,f2)
			y1 = expecty + theta1*timesteplength*f + (1-theta1)*timesteplength*expectf + & 
				& ((3*theta1-1)/6.0)*(timesteplength**2)*(-f1*f+0.5*f2*(z(timestep-1,xindex)**2))+&
				((3*theta1-2)/6.0)*(timesteplength**2)*(-expectf1f +0.5 * expectf2zz)
			diff = abs(y1-y0)
			if (diff <= 1e-12) then
				exit
			end if
			y0 = y1
		end do 
		y(timestep-1,xindex) = y0 
		!print*,"y(timestep-1,xindex)=",y(timestep-1,xindex)
	end do 
	!print*,"y(timestep-1,:)=",y(timestep-1,:)
	
! when 0<= n <=N-2

		do timeindex = timestep-2, 0, -1
		!time fixed x change ,Z
			do xindex = -xstep, xstep
				expectyw = 0
				expectz = 0	
				expectfw = 0
				expectf1fw = 0
				expectf2zzw = 0				
				!find the interpolation point
				do i = 1, k
					x_point = x(xindex)+sqrt(2*timesteplength)*a(i)
					
					x_point3 = floor(x_point/xsteplength)
					x_point2 = x_point3-1
					x_point1 = x_point3-2
					x_point4 = x_point3+1
					x_point5 = x_point3+2
					overflow = 0

					if (x_point1 < -xstep) then
						x_point1 = -xstep
						x_point2 = -xstep + 1
						x_point3 = -xstep + 2
						x_point4 = -xstep + 3
						x_point5 = -xstep + 4
						overflow = 1
					end if
					if (x_point5 > xstep ) then 
						x_point5 = xstep
						x_point4 = xstep - 1
						x_point3 = xstep - 2
						x_point2 = xstep - 3
						x_point1 = xstep - 4
						overflow = 1
					end if

					if ((x(x_point1) > x_point) .and. (overflow == 0)) then
						print*, "compute error"
						print*, "low val", x(x_point1)
						print*, "point3", x(x_point3)
						print*, "val ", x_point
					end if

					if ((x(x_point5) < x_point) .and. (overflow == 0)) then
						print*, "compute error"
						print*, "high val", x(x_point5)
						print*, "point3", x(x_point3)
						print*, "val ", x_point
					end if
					!start to interpolate
					yinterpl = (/y(timeindex+1,x_point1),y(timeindex+1,x_point2),y(timeindex+1,x_point3), &
									y(timeindex+1,x_point4),y(timeindex+1,x_point5)/)
					xinterpl = (/x(x_point1),x(x_point2),x(x_point3),x(x_point4),x(x_point5)/)
					!print*,"xinterpl=",xinterpl
					!print*,"yinterpl=",yinterpl
					!print*,"x_point=",x_point
					call newtoninterpl(xinterpl,yinterpl, x_point, yb)	
					!print*,"yb",yb
					x_point = x(xindex)+sqrt(2*timesteplength)*a(i)		
					zinterpl = (/z(timeindex+1,x_point1),z(timeindex+1,x_point2),z(timeindex+1,x_point3), &
									z(timeindex+1,x_point4),z(timeindex+1,x_point5)/)
					xinterpl =  (/x(x_point1),x(x_point2),x(x_point3),x(x_point4),x(x_point5)/)
					!print*,"xinterpl=",xinterpl
					!print*,"zinterpl=",zinterpl
					!print*,"x_point=",x_point
					call newtoninterpl(xinterpl, zinterpl, x_point, zb)
					!print*,"zb",zb
					call function3(yb,f,f1,f2)
					expectyw = expectyw +w(i)*yb*sqrt(2*timesteplength)*a(i)
					expectz = expectz + w(i)*zb 
					expectfw = expectfw +w(i)*f*sqrt(2*timesteplength)*a(i)
					expectf1fw = expectf1fw +w(i)*f1*f*sqrt(2*timesteplength)*a(i)
					expectf2zzw = expectf2zzw +w(i)*f2*(zb**2)*sqrt(2*timesteplength)*a(i)
					!print*, 'do k'
				end do
				!print*, 'end k'
				expectyw = expectyw/sqrt(pi)
				expectz =expectz/sqrt(pi)
				expectfw = expectfw/sqrt(pi)
				expectf1fw = expectf1fw/sqrt(pi)
				expectf2zzw = expectf2zzw/sqrt(pi)
				!print*,"expectyw:=",expectyw, "expectz=", expectz,"expectfw=",expectfw,"expectf1fw=",expectf1fw,"expectf2zzw=" ,expectf2zzw
				z(timeindex,xindex) = (expectyw - (1-theta2)*timesteplength*expectz + theta2*timesteplength*expectfw & 
								&	-theta2*0.5*(timesteplength**2)*(-expectf1fw &
									+0.5*expectf2zzw))/(theta2*timesteplength)
				!print*,"z(timeindex,xindex)=",z(timeindex,xindex)
			end do 
				!print*,"z(timeindex,:)=",z(timeindex,:)
	! time fixed x change ,Y 	
		!print*, 'calculate y'
		do xindex = -xstep ,xstep
			expecty = 0
			expectf = 0
			expectf1f = 0
			expectf2zz = 0
			do i = 1,k
				x_point = x(xindex)+sqrt(2*timesteplength)*a(i)
				x_point3 = floor(x_point/xsteplength)
				x_point2 = x_point3-1
				x_point1 = x_point3-2
				x_point4 = x_point3+1
				x_point5 = x_point3+2
				overflow = 0

				if (x_point1 < -xstep) then
					x_point1 = -xstep
					x_point2 = -xstep + 1
					x_point3 = -xstep + 2
					x_point4 = -xstep + 3
					x_point5 = -xstep + 4
					overflow = 1
				end if
				if (x_point5 > xstep ) then 
					x_point5 = xstep
					x_point4 = xstep - 1
					x_point3 = xstep - 2
					x_point2 = xstep - 3
					x_point1 = xstep - 4
					overflow = 1
				end if

				if ((x(x_point1) > x_point) .and. (overflow == 0)) then
					print*, "compute error"
					print*, "low val", x(x_point1)
					print*, "val ", x_point
				end if

				if ((x(x_point5) < x_point) .and. (overflow == 0)) then
					print*, "compute error"
					print*, "high val", x(x_point5)
					print*, "val ", x_point
				end if
				!start to interpolate
				yinterpl = (/y(timeindex+1,x_point1),y(timeindex+1,x_point2),y(timeindex+1,x_point3),&
								y(timeindex+1,x_point4),y(timeindex+1,x_point5)/)
				xinterpl = (/x(x_point1),x(x_point2),x(x_point3),x(x_point4),x(x_point5)/)
				!print*,"xinterpl=",xinterpl
				!print*,"yinterpl=",yinterpl
				!print*,"x_point=",x_point
				call newtoninterpl(xinterpl,yinterpl, x_point, yb)	
				!print*,"yb",yb
				x_point = x(xindex)+sqrt(2*timesteplength)*a(i)		
				zinterpl = (/z(timeindex+1,x_point1),z(timeindex+1,x_point2),z(timeindex+1,x_point3),&
								z(timeindex+1,x_point4),z(timeindex+1,x_point5)/)
				xinterpl =  (/x(x_point1),x(x_point2),x(x_point3),x(x_point4),x(x_point5)/)
				call newtoninterpl(xinterpl, zinterpl, x_point, zb)
				!print*,"zb",zb
				call function3(yb,f,f1,f2)
				expecty = expecty +w(i)*yb
				expectf = expectf +w(i)*f
				expectf1f= expectf1f +w(i)*f1*f
				expectf2zz = expectf2zz +w(i)*f2*(zb**2)
		
			end do

			expecty = expecty/sqrt(pi)
			expectf = expectf/sqrt(pi)
			expectf1f = expectf1f/sqrt(pi)
			expectf2zz = expectf2zz/sqrt(pi)
			!print*,"expecty= ",expecty ,"expectf=",expectf,"expectf1f=",expectf1f,"expectf2zz=",expectf2zz
			
			y0 = y(timeindex+1,xindex)
			!print*,"y0=",y0
			cycleindex = 0
			do 
				cycleindex = cycleindex + 1
				!print*, 'circle ', cycleindex
				call function3(y0,f,f1,f2)
				y1 = expecty + theta1*timesteplength*f + (1-theta1)*timesteplength*expectf + & 
				& ((3*theta1-1)/6.0)*(timesteplength**2)*(-f1*f+0.5*f2*(z(timeindex,xindex)**2))+&
				((3*theta1-2)/6.0)*(timesteplength**2)*(-expectf1f +0.5 * expectf2zz)
				diff = abs(y1-y0)
				if (diff <= 1e-12) then
					EXIT
				end if
				y0 = y1
			end do 
			y(timeindex,xindex) = y0 
			!print*,"y(timeindex,xindex)=",y(timeindex,xindex)		
		end do 
		!print*, 'end compute y'
		!print*,"y(timeindex,:)=",y(timeindex,:)
	end do
		!print*,"z=",z,"y=",y
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
