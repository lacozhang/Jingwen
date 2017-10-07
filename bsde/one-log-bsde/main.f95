module main
	implicit none
	private
	public:: gauss_hermite,  interpolateIndex, function34,finalvaluey,&
			valuez, thelengthofx,  newtoninterpl, leastsquare, li31,func34derivative


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
			a(temp) = dsqrt(dfloat(temp)*0.5d0)
		end do 
		
		cm(:,:) = 0.0d0
		do temp = 1,n-1
			cm(temp,:) = 0.0d0
			cm(temp,temp+1) = a(temp)
		end do

		lwork = 10*n
		call dsyev('V','U', n, cm, n, l, work, lwork, info)
		x = l
		w = dsqrt(4.0d0*dATAN(1.0d0))*cm(1,:)**2
		deallocate(a)
		deallocate(cm)
	end subroutine		

	pure subroutine function34(y,f)
		real(kind=8),intent(in)::y
		real(kind=8),intent(out)::f
		!f = -y*(1-y)*(0.75- y)
		f = - (y**3) + 2.5d0 * (y ** 2) - 1.5d0 * y
	end subroutine

	subroutine func34derivative(y, grad_y)
		real(kind=8),intent(in)::y
		real(kind=8),intent(out)::grad_y
		grad_y = -3.0d0 * y ** 2 + 5.0d0 * y - 1.5d0 
	end subroutine
	
	pure subroutine finalvaluey(t,x,y)
		real(kind=8),intent(inout)::t,x
		real(kind=8),intent(out)::y
		!y = (1.0)/(1+dexp(-x-0.25*t))		
		y = dexp(t+x)/(dexp(t+x)+1.0d0)
	end subroutine

	pure subroutine valuez(t,x,z)
		real(kind=8),intent(inout)::t,x
		real(kind=8),intent(out)::z
		!z = (dexp(-x-0.25*t))/((1.0d0+dexp(-x-0.25*t))**2)
		z = dexp(x+t)/((dexp(x+t)+1.0d0)**2)
	end subroutine

	subroutine thelengthofx(x0,k,timestep,xrange)
		real(kind=8), intent(in):: x0
		integer, intent(in):: k, timestep
		real(kind=8),intent(inout):: xrange
		real(kind=8),dimension(k):: w, a
		real(kind=8):: timesteplength, xsteplength,  x_max, x_min
		real(kind=8),dimension(:),allocatable:: x, x_need
		real(kind=8),dimension(0:timestep)::t
		integer:: xstep, xindex, timeindex, xneedindex
		allocate(x_need(k))
		w(:)=0.0d0
		a(:)=0.0d0
		call gauss_hermite(k, w, a)
		!print*,"w=",w
		!print*,"a=",a
		timesteplength = 1.0d0/dfloat(timestep)
		!print*,"timesteplength=",timesteplength
		xsteplength = timesteplength**(0.75)
		!print*,"xsteplength=",xsteplength
		!xsteplength depended on the order of x and the order of t
		xstep = floor(xrange/xsteplength)
		allocate(x(-xstep:xstep)) 
		do timeindex = 0, timestep
			t(timeindex) = dfloat(timeindex) 
		end do 
		t = t * timesteplength
		x(0) = x0
		do xindex = -xstep,xstep
			x(xindex) = x(0)+dfloat(xindex)*xsteplength
		end do
		xneedindex = 0
		do timeindex = 0,timestep
			x_need = x(xneedindex)+dsqrt(2.0d0*timesteplength)*a
			x_max = maxval(dabs(x_need))
			!print*, "x_max= ", x_max
			xneedindex = floor((x_max-x0)/xsteplength)+1
			
			if (xneedindex > xstep)then
				print*,"xneedindex= ", xneedindex
				print*,"xstep=", xstep
				print*, "xmax=", x_max
				print*, "x0= ", x0
				print*, "xmax - x0 = ", x_max-x0
				print*, "xsteplength= ", xsteplength
				print*,"please,extend the xrange"
				exit
			end if 
			
		end do 
		xrange = x0 + xneedindex*xsteplength
		!print*,"xrange= ", xrange
		deallocate(x)
		deallocate(x_need)
	end subroutine




subroutine newtoninterpl(x, y, x_aim, y_aim)
		!finall y_aim which is our need
		real(kind=8), dimension(:), intent(in):: x
		real(kind=8), dimension(:), intent(in):: y
		real(kind=8), intent(inout):: x_aim
		real(kind=8), intent(out):: y_aim 
		real(kind=8), dimension(:,:), allocatable:: z
		real(kind=8), dimension(:),allocatable::x_
		integer:: sizex, temp1, temp2,  w
		real(kind=8)::n
		sizex = size(x)
		allocate(z(sizex,sizex))
		allocate(x_(sizex))
		x_ = x 
		z(:,:) = 0.0d0
		n = 0.0d0
		z(:,1) = y
		do temp2 = 2,sizex
			do temp1 = temp2, sizex
				z(temp1,temp2) = (z(temp1-1,temp2-1)-z(temp1,temp2-1))/(x_(temp1-(temp2-1))-x_(temp1))
			end do
		end do
		x_(2:sizex) = x_(1:(sizex-1))
		x_(1) = x_aim - 1.0d0
		x_ = x_aim - x_
		do temp1 = 2, sizex
			x_(temp1) = x_(temp1-1)*x_(temp1)
		end do
		do temp1 = 1, sizex
			n = n+z(temp1,temp1)*x_(temp1)
		end do
		y_aim = n
		deallocate(z)
		deallocate(x_)
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
		timesteplength = 1.0d0 / dfloat(timestep)
		error = dlog(error)
		timesteplength = dlog(timesteplength)
		a(:,:) = 1.0d0
		a(:,2) = timesteplength
		call dgels('N', sizeerror, 2, 1, a, sizeerror, error, sizeerror, work, lwork, info)
		order = error(2)
		deallocate(timesteplength)
		deallocate(a)
		deallocate(work)
	end subroutine

	!using 4 points
	subroutine interpolateIndex(xvalues, yvalues, xpoint, xsteplength, timesteps, xsteps, timeindex, xintvalues, yintvalues)
	!xvalues:is the all values of x   xintvalues: we need (xpoint in them) 
		real(kind=8), dimension(-xsteps:xsteps), intent(in):: xvalues
		real(kind=8), dimension(0:timesteps,-xsteps:xsteps), intent(in):: yvalues
		real(kind=8), intent(in):: xpoint
		real(kind=8), intent(in):: xsteplength
		integer, intent(in):: xsteps, timeindex, timesteps
		real(kind=8), dimension(4), intent(inout)::xintvalues, yintvalues
		integer, dimension(4):: xindex
		integer:: overflow

		overflow = 0
		xindex(3) = floor(xpoint/xsteplength)
		xindex(4) = xindex(3)+1
		xindex(2) = xindex(3)-1
		xindex(1) = xindex(3)-2
		if(xindex(1) < -xsteps) then
			xindex(1) = -xsteps
			xindex(2) = xindex(1)+1
			xindex(3) = xindex(2)+1
			xindex(4) = xindex(3)+1
			overflow = 1
		end if

		if(xindex(4) > xsteps) then
			xindex(4) = xsteps
			xindex(3) = xindex(4)-1
			xindex(2) = xindex(3)-1
			xindex(1) = xindex(2)-1
			overflow = 1
		end if

		xintvalues(1) = xvalues(xindex(1))
		xintvalues(2) = xvalues(xindex(2))
		xintvalues(3) = xvalues(xindex(3))
		xintvalues(4) = xvalues(xindex(4))

		yintvalues(1) = yvalues(timeindex, xindex(1))
		yintvalues(2) = yvalues(timeindex, xindex(2))
		yintvalues(3) = yvalues(timeindex, xindex(3))
		yintvalues(4) = yvalues(timeindex, xindex(4))


		if ((xvalues(xindex(1)) > xpoint) .and. (overflow == 0)) then
			print*, "compute error"
			print*, "low val", xvalues(xindex(1))
			print*, "val ", xpoint
		end if

		if ((xvalues(xindex(4)) < xpoint) .and. (overflow == 0)) then
			print*, "compute error"
			print*, "high val", xvalues(xindex(4))
			print*, "val ", xpoint
		end if
	end subroutine

	subroutine li31(timestep,xrange, yerror, zerror, theta)
	!moudle interpolation gaussherimite function31
	!-llapck -lblas

		integer, intent(inout):: timestep
		real(kind=8), intent(out)::yerror, zerror
		real(kind=8), intent(inout)::theta
		integer, intent(in)::xrange 
		real(kind=8), dimension(:,:), allocatable:: y, z
		real(kind=8),dimension(:),allocatable:: t, x
		real(kind=8):: timesteplength, xsteplength, yb, zb, actual_range
		integer:: cycleindex
		integer:: i, timeindex, xindex, xstep, k,tindex
		integer:: overflow,signal 
		real(kind=8):: f, fy, expecty, expectf, expectfyz, expectz, pi, y0, y1, y_i, z_i, x_point 
		real(kind=8):: diff_y, diff_z
		real(kind=8),dimension(:),allocatable:: w, a
		real(kind=8), dimension(4)::yinterpl, xinterpl,zinterpl
		real(kind=8), parameter:: yture=0.5d0, zture=0.25d0, x0= 0.0d0

		pi = 4.0d0*dATAN(1.0d0)
		!print*,"pi=",pi
		print*, "number of steps : ", timestep
		timesteplength = (1.0d0)/dfloat(timestep)
		print*, "timesteplength  : ", timesteplength
		
		
	! GH
		k = 10
		print*, "number of steps : ", timestep
		print*, "gauss hermite   : ", k

		allocate(w(k))
		allocate(a(k))
		call gauss_hermite(k, w, a)
	! end GH
     
		if(dabs(theta -0.5d0) < 1d-15) then
			xsteplength = (timesteplength) ** (0.75)
			print*, "using 0.75"
		else
			xsteplength = (timesteplength) ** (0.5)
			print*, "using 0.5"
		end if
		actual_range = xrange
		call thelengthofx(x0,k,timestep,actual_range)
		xstep = floor(actual_range/xsteplength)
		print*,"actual_range:",actual_range
		

		allocate(t(0:timestep))
		allocate(x(-xstep:xstep))

	!
		t(0) = 0.0d0
		do timeindex = 1, timestep
			t(timeindex) = t(timeindex-1) + timesteplength
		end do 

	! x space range
		x(0) = 0.0d0
		do xindex = -xstep, xstep
			x(xindex) = dfloat(xindex)
		end do
		x = x * xsteplength

		allocate(y(0:timestep,-xstep:xstep))
		allocate(z(0:timestep,-xstep:xstep))
		
		y(:,:)=0.0d0
		z(:,:)=0.0d0

	! last line of Y matrix
		do xindex = -xstep, xstep
			call finalvaluey(t(timestep),x(xindex),y(timestep,xindex))
		end do 

	!last  line of Z matrix
		do xindex = -xstep, xstep
			call valuez(t(timestep),x(xindex),z(timestep,xindex))
		end do 
		
! when 0<= n <=N-1

		do timeindex = timestep-1, 0, -1
			do xindex = -xstep, xstep

				y_i = y(timeindex+1, xindex)
				z_i = z(timeindex+1, xindex)

				do
					expecty = 0.0d0
					expectf = 0.0d0
					expectfyz = 0.0d0
					expectz = 0.0d0			
					!find the interpolation point
					do i = 1, k
						x_point = x(xindex)+dsqrt(2.0d0*timesteplength)*a(i)
						tindex = timeindex+1
						call interpolateIndex(x, y, x_point, xsteplength, timestep, xstep, tindex,xinterpl, yinterpl)
						call interpolateIndex(x, z, x_point, xsteplength, timestep, xstep, tindex,xinterpl, zinterpl)
						call newtoninterpl(xinterpl, yinterpl, x_point, yb)
						call newtoninterpl(xinterpl, zinterpl, x_point, zb)
				
						expecty = expecty+w(i)*yb
						expectz = expectz + w(i)*zb

						call function34(yb*(1.0d0-theta)+theta*y_i,f)
						expectf = expectf+w(i)*f
						call func34derivative(yb*(1.0d0-theta)+theta*y_i, fy)
						expectfyz = expectfyz + w(i)*fy*((1.0d0-theta)*zb+z_i)
					end do
!					print*, expectz 
					expecty = expecty/dsqrt(pi)
					expectf = expectf/dsqrt(pi)
					expectfyz = exptctfyz/dsqrt(pi)
					expectz = expectz/dsqrt(pi)

					y(timeindex, xindex) = expecty + timesteplength * expectf
					z(timeindex, xindex) = expectz + timesteplength * expectfyz
					
					diff_y = dabs(y(timeindex, xindex) - y_i)
					diff_z = dabs(z(timeindex, xindex) - z_i)
					
					if (diff_y <= 10-12) .and. (diff_z <= 10-12) then
						exit
					end if

				end do
			end do 
			
		end do 
		!print*,"z(0,0)=",z(0,0)
		!print*,"y(0,0)=",y(0,0)
		zerror = dabs(zture - z(0,0))
		yerror = dabs(yture - y(0,0))
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
