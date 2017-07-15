module main
!-llapck -lblas
	implicit none
	private
	public::function51sde, function51, turevalue, &
			gauss_hermite, thelengthofx, bubble_sort, newtoninterpl, leastsquare, li51

contains


	pure subroutine function51sde(t,x,b,tau)
		real(kind=8),intent(inout)::t, x
		real(kind=8),intent(out)::b, tau
		b = cos(x)
		tau = (sin(x)+3.0)*cos(x)
	end subroutine function51sde
	
	pure subroutine function51(t,x,y,z,f)
		real(kind=8),intent(inout)::t, x, y, z
		real(kind=8),intent(out)::f
		f = - y - z/(sin(x)+3.0) + 0.5*z*(1+3.0*sin(x))
	end subroutine function51
	
	pure subroutine turevalue(t,x,y,z)
		real(kind=8),intent(inout)::t, x
		real(kind=8),intent(out)::y,z
		y = exp(t)*log(sin(x)+3.0)
		z = exp(t)*((cos(x))**2)
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
	
	
	recursive subroutine getmaxlength(coeffs, k, tidx, step, timesteplength, x_current, x_max)
	
		integer, intent(in)::k, tidx, step
		integer:: i
		real(kind=8), intent(inout)::x_current, x_max
		real(kind=8), dimension(k)::coeffs
		real(kind=8):: x_temp, x_temp_max, t, timesteplength, b, tau
		
		if(tidx >= step) then
			if(abs(x_current) > x_max) then	
				x_max = x_current
			end if
			return
		end if
	
		t = tidx * timesteplength
		call function51sde(t, x_current, b, tau)
	
		do i=1,k
			x_temp = x_current + b*timesteplength + tau*(sqrt(2.0*timesteplength)*coeffs(i))
			x_temp_max = 0.0
			call getmaxlength(coeffs, k, tidx+1, step, timesteplength, x_temp, x_temp_max)
			if (abs(x_temp_max)>x_max) then
				x_max = x_temp_max
			end if
		end do
	
	end subroutine
	
	subroutine thelengthofx(k, step, timesteplength, xint)
	!k is gaussherimite's points
	!step is the scheme's step
	!xint returns result which is the length of x in positive
		integer, intent(in):: k, step
		integer, intent(out):: xint
		real(kind=8),intent(inout):: timesteplength
		real(kind=8), dimension(:),allocatable::w
		real(kind=8), dimension(:),allocatable::a
		real(kind=8)::x_max, x_prev, x_temp
		real(kind=8)::b, tau , t
		integer::i, j
		allocate(w(k))
		allocate(a(k))
		x_max = 1.0
		x_temp = 1.0
		call gauss_hermite(k, w, a)
		
		
		do i=1,step
			x_prev = x_max
			t = i*timesteplength
			call function51sde(t, x_prev, b, tau)
			x_max = 0.0
			do j=1,k
				x_temp = x_prev + timesteplength + 4.0*sqrt(2.0*timesteplength)*a(j)
				if(abs(x_temp) > x_max) then
					x_max = abs(x_temp)
				end if
			end do
		end do
		!call getmaxlength(a, k, 1, step, timesteplength, x_temp, x_max)
		xint = floor(x_max) + 3
		!print*,floor(x_max)
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
	
	subroutine li51(timestep, yerror, zerror)
	!cranck nicolson method
	!-llapck -lblas
		
		integer, intent(inout):: timestep
		real(kind=8), intent(out)::yerror, zerror


		real(kind=8), dimension(:,:), allocatable:: y, z
		real(kind=8), dimension(:),allocatable:: t, x
		real(kind=8):: x0,timesteplength, xsteplength, yb, zb

! 	
		integer::i, timeindex, xindex, xstep, xrange, k, x_point1, x_point2, x_point3, x_point4, max_try
		integer:: overflow
		real(kind=8):: b, tau, f, expecty, expectf, expectyw, expectfw, expectz, pi, diff, y0, y1, x_point , yture, zture
		real(kind=8),dimension(:),allocatable:: w, a
		real(kind=8), dimension(4)::yinterpl,  xinterpl , zinterpl

		print*, "start program"
		pi = 4*atan(1.0)

		timesteplength = (1.0)/timestep
		print*, "time step size  : ", timesteplength
		
		
	! GH
		k = 8
		!print*, "number of steps : ", timestep
		!print*, "gauss hermite   : ", k

		allocate(w(k))
		allocate(a(k))
		call gauss_hermite(k, w, a)
	! end GH
	
		call thelengthofx(k, timestep, timesteplength, xrange)
		xsteplength = timesteplength**(3.0/4.0)
		xstep = ceiling(xrange*1.0/xsteplength)
		print*,"xrange:",xrange
		print*,"xstep:",xstep
		allocate(t(0:timestep))
		allocate(x(-xstep:xstep))

	!	t sapce

		do timeindex = 0, timestep
			t(timeindex) = timeindex 
		end do 
		t = t * timesteplength
!		print*,"t=:",t
	! x space 
		x(0) = 1
		do xindex = -xstep,xstep
			x(xindex) = x(0)+xindex*xsteplength
		end do
	
		!print*,"x=:",x
		
		allocate(y(0:timestep,-xstep:xstep))
		allocate(z(0:timestep,-xstep:xstep))
		y(:,:)=0
		z(:,:)=0

	! last line of Y matrix and Z marix(timestep)
		do xindex = -xstep, xstep
			call turevalue(t(timestep),x(xindex),y(timestep,xindex),z(timestep,xindex))
		end do 
		!print*,"last Y:",y(timestep,:),"last Z:",z(timestep,:)	
		
		
	! last but one line of Z matrix	and Y marix(timestep-1)
		do xindex = -xstep, xstep
			expectyw = 0
			expectfw = 0
			expectz = 0
			expecty = 0
			expectf = 0
			do i = 1, k
				call function51sde(t(timestep-1),x(xindex),b,tau)
				!print*,"t=:",t(timestep-1),"x=:",x(xindex),"b=:",b,"tau=:",tau
				x_point = x(xindex)+b*timesteplength+tau*sqrt(2*timesteplength)*a(i)
				call turevalue(t(timestep),x_point,yb,zb)
				call function51(t(timestep),x_point,yb,zb,f)
				expectyw = expectyw + w(i)*yb*sqrt(2*timesteplength)*a(i)
				expectfw = expectfw + w(i)*f*sqrt(2*timesteplength)*a(i)
				expectz = expectz + w(i)*zb 
				expecty = expecty + w(i)*yb 
				expectf = expectf + w(i)*f
			end do 
			expectyw = expectyw/sqrt(pi)
			expectfw = expectfw/sqrt(pi)
			expectz = expectz/sqrt(pi)
			expecty = expecty/sqrt(pi)
			expectf = expectf/sqrt(pi)
			z(timestep-1,xindex) = (expectyw+0.5*timesteplength*expectfw-0.5*timesteplength*expectz)/(0.5*timesteplength)
			y0 = y(timestep,xindex)
			do 
				call function51(t(timestep-1),x(xindex),y0,z(timestep-1,xindex),f)
				y1 = expecty + 0.5*timesteplength*expectf + 0.5*timesteplength*f
				diff = abs(y1-y0)
				if (diff <= 1e-10) then 
					exit
				end if
				y0 = y1
			end do 
			y(timestep-1,xindex) = y1
		end do 
		!print*,"z(timestep-1,:)=:",z(timestep-1,:)	
		!print*,"y(timestep-1,:)=", y(timestep-1,:)


		
	
! when 0<= n <=N-2

		do timeindex = timestep-1, 0, -1
			!print*,"times:",timeindex
			do xindex = -xstep, xstep
				expectyw = 0
				expectfw = 0
				expectz = 0	
				expecty = 0
				expectf = 0	
				
				call function51sde(t(timeindex),x(xindex),b,tau)
				do i = 1, k
					!print*,"t=:",t(timeindex),"x=:",x(xindex),"b=:",b,"tau=:",tau
					x_point = x(xindex)+b*timesteplength+tau*sqrt(2*timesteplength)*a(i)
					!find the interpolation point
					x_point3 = floor((x_point-1)/xsteplength)
					!print*, "x point3", x_point3
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
					yinterpl = (/y(timeindex+1,x_point1),y(timeindex+1,x_point2),y(timeindex+1,x_point3),y(timeindex+1,x_point4)/)
					xinterpl = (/x(x_point1),x(x_point2),x(x_point3),x(x_point4)/)
					call newtoninterpl(xinterpl,yinterpl, x_point, yb)	
	
					!print*, "yb ", yb
					xinterpl = (/x(x_point1),x(x_point2),x(x_point3),x(x_point4)/)
					zinterpl = (/z(timeindex+1,x_point1),z(timeindex+1,x_point2),z(timeindex+1,x_point3),z(timeindex+1,x_point4)/)
					!print*, "z inter", zinterpl	
					!print*, "x inter", xinterpl
					!print*, "x point", x_point
					call newtoninterpl(xinterpl, zinterpl, x_point, zb)
					!print*, "zb ", zb
					call function51(t(timeindex+1),x_point,yb,zb,f)
				
					expectyw = expectyw + w(i)*yb*sqrt(2*timesteplength)*a(i)
					expectfw = expectfw + w(i)*f*sqrt(2*timesteplength)*a(i)
					expectz = expectz + w(i)*zb 
					expecty = expecty + w(i)*yb 
					expectf = expectf + w(i)*f
				end do 

				expectyw = expectyw/sqrt(pi)
				expectfw = expectfw/sqrt(pi)
				expectz = expectz/sqrt(pi)
				expecty = expecty/sqrt(pi)
				expectf = expectf/sqrt(pi)		
				z(timeindex,xindex) = (expectyw+0.5*timesteplength*expectfw-0.5*timesteplength*expectz)/(0.5*timesteplength)

				y0 = y(timeindex+1,xindex)
				max_try = 0
				do 
					if(max_try > 10000) then
						print*, "quit cycle"
						exit
					end if
					call function51(t(timeindex),x(xindex),y0,z(timeindex,xindex),f)
					y1 = expecty + 0.5*timesteplength*expectf + 0.5*timesteplength*f
					diff=abs(y0-y1)
					if (diff <= 1e-10) then 
						exit
					end if
					y0 = y1
					max_try = max_try + 1
				end do
				y(timeindex,xindex) = y1
				!print*, "y final value", y1
			end do
			!print*, "y value", y(timeindex, :)
			!print*, "z value", z(timeindex, :)
		end do 
		!print*,"z(0,0)=",z(0,0)
		!print*,"y(0,0)=",y(0,0)
		call turevalue(t(0),x(0),yture,zture)
		print*,"yture",yture,"zture",zture
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