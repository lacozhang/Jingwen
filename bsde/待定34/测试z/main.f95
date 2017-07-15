module main
	implicit none
	private
	public:: gauss_hermite, thelengthofx, bubble_sort, function34,finalvaluez,valuey, newtoninterpl, leastsquare, li34z


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
	pure subroutine function34(y,f,f1,f2)
		real(kind=8),intent(inout)::y
		real(kind=8),intent(out)::f,f1,f2
		! f = -y*(1-y)*(0.75- y)
		! f1 = -3*(y**2) +3.5*y-0.75
		! f2 = -6*y +3.5
		f = - (y**3) + 2.5* (y ** 2) - 1.5*y
		f1 = - 3 * (y**2) + 5*y - 1.5
		f2 = -6*y+5
	end subroutine
	pure subroutine finalvaluez(t,x,z)
		real(kind=8),intent(inout)::t,x
		real(kind=8),intent(out)::z
		!z = (exp(-x-0.25*t))/((1+exp(-x-0.25*t))**2)
		z = exp(x+t)/((exp(x+t)+1)**2)
	end subroutine
	pure subroutine valuey(t,x,y)
		real(kind=8),intent(inout)::t,x
		real(kind=8),intent(out)::y
		!y = (1.0)/(1+exp(-x-0.25*t))
		y = exp(t+x)/(exp(t+x)+1)
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
	subroutine li34z(timestep ,zerror, theta2)
		integer, intent(inout):: timestep
		real(kind=8), intent(out):: zerror
		real(kind=8), intent(inout):: theta2

		real(kind=8), dimension(:,:), allocatable::  z
		real(kind=8), dimension(:),allocatable:: t, x
		real(kind=8):: timesteplength, xsteplength, yb, zb, z_true_value

! 	
		integer::  i, timeindex, xindex, xrange, xstep, k, x_point1, x_point2,x_point3,x_point4
		integer:: overflow
		real(kind=8):: f, f1, f2
		real(kind=8):: expectyw, expectz, expectfw,expectf1fw,expectf2zzw
		real(kind=8):: pi, diff, y0, y1, x_point 
		real(kind=8), dimension(:),allocatable:: w, a
		real(kind=8), dimension(4):: xinterpl,zinterpl
		real(kind=8), parameter:: zture=0.25
	
		pi = 4*atan(1.0)
		!print*,"pi=",pi
		print*, "number of steps : ", timestep
		timesteplength = (1.0)/timestep
		print*, "timesteplength  : ", timesteplength
		
		
	! GH
		k = 10
		
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
		xsteplength = timesteplength**(0.75)
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
		
		
		allocate(z(0:timestep,-xstep:xstep))
		
	!last  line of Z matrix
		do xindex = -xstep, xstep
			call finalvaluez(t(timestep),x(xindex),z(timestep,xindex))
		end do 
		
		
	
	


! when 0<= n <=N-2

		do timeindex = timestep-1, 0, -1
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
			
					overflow = 0

					if (x_point1 < -xstep) then
						x_point1 = -xstep
						x_point2 = -xstep + 1
						x_point3 = -xstep + 2
						x_point4 = -xstep + 3
					
						overflow = 1
					end if
					if (x_point4 > xstep ) then 
						x_point4 = xstep
						x_point3 = xstep - 1
						x_point2 = xstep - 2
						x_point1 = xstep - 3
	
						overflow = 1
					end if

					if ((x(x_point1) > x_point) .and. (overflow == 0)) then
						print*, "compute error"
						print*, "low val", x(x_point1)
						print*, "point3", x(x_point3)
						print*, "val ", x_point
					end if

					if ((x(x_point4) < x_point) .and. (overflow == 0)) then
						print*, "compute error"
						print*, "high val", x(x_point4)
						print*, "point3", x(x_point3)
						print*, "val ", x_point
					end if
					!start to interpolate	
					zinterpl = (/z(timeindex+1,x_point1),z(timeindex+1,x_point2),z(timeindex+1,x_point3),&
									z(timeindex+1,x_point4)/)
					xinterpl =  (/x(x_point1),x(x_point2),x(x_point3),x(x_point4)/)
					!print*,"xinterpl=",xinterpl
					!print*,"zinterpl=",zinterpl
					!print*,"x_point=",x_point
					call newtoninterpl(xinterpl, zinterpl, x_point, zb)
					!print*,"zb",zb
					call valuey(t(timeindex+1),x_point,yb)
					call function34(yb,f,f1,f2)
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
				! z(timeindex,xindex) = (expectyw - (1-theta2)*timesteplength*expectz+ theta2*timesteplength*expectfw)&
				! /(theta2*timesteplength)					
				!print*,"z(timeindex,xindex)=",z(timeindex,xindex)
				 if (timeindex > 0) then
					 call finalvaluez(t(timeindex), x(xindex), z_true_value)
					 z(timeindex, xindex) = z_true_value
				 end if
			end do 
				!print*,"z(timeindex,:)=",z(timeindex,:)
		
				
		end do 
	

	
		zerror = abs(zture - z(0,0))
	
		deallocate(w)
		deallocate(a)
		deallocate(z)
		deallocate(t)
		deallocate(x)
	end subroutine
end module main
