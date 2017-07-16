module main
	implicit none
	private
	public:: gauss_hermite, thelengthofx,interpolateIndex,&
		function34,finalvaluez,valuey, newtoninterpl, leastsquare, li34z


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
	pure subroutine function34(y,f,f1,f2)
		real(kind=8),intent(in)::y
		real(kind=8),intent(out)::f,f1,f2
		!f = -y*(1-y)*(0.75- y)
		!f1 = -3*(y**2) +3.5*y-0.75
		!f2 = -6*y +3.5
		f = - (y**3) + 2.5d0* (y ** 2) - 1.5d0*y
		f1 = - 3.0d0 * (y**2) + 5.0d0*y - 1.5d0
		f2 = -6.0d0*y+5.0d0
		
	end subroutine
	pure subroutine finalvaluez(t,x,z)
		real(kind=8),intent(inout)::t,x
		real(kind=8),intent(out)::z
		!z = (exp(-x-0.25*t))/((1+exp(-x-0.25*t))**2)
		z = dexp(x+t)/((dexp(x+t)+1.0d0)**2)
	end subroutine
	pure subroutine valuey(t,x,y)
		real(kind=8),intent(inout)::t,x
		real(kind=8),intent(out)::y
		!y = (1.0)/(1+exp(-x-0.25*t))		
		y = dexp(t+x)/(dexp(t+x)+1.0d0)
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
	subroutine li34z(timestep ,xrange, zerror, theta2)
		integer, intent(inout):: timestep
		real(kind=8), intent(out):: zerror
		real(kind=8), intent(inout):: theta2
		integer, intent(in)::xrange 
		real(kind=8), dimension(:,:), allocatable::  z
		real(kind=8), dimension(:),allocatable:: t, x
		real(kind=8):: timesteplength, xsteplength, yb, zb, z_true_value,actual_range

! 	
		integer::  i, timeindex,tindex, xindex, xstep, k
		integer:: overflow
		real(kind=8):: f, f1, f2
		real(kind=8):: expectyw, expectz, expectfw,expectf1fw,expectf2zzw
		real(kind=8):: pi, diff, y0, y1, x_point 
		real(kind=8), dimension(:),allocatable:: w, a
		real(kind=8), dimension(4):: xinterpl,zinterpl
		real(kind=8), parameter:: zture=0.25d0, x0=0.0d0
	
		pi = 4.0d0*datan(1.0d0)
		!print*,"pi=",pi
		print*, "number of steps : ", timestep
		timesteplength = (1.0d0)/dfloat(timestep)
		print*, "timesteplength  : ", timesteplength
		
		
	! GH
		k = 10
		

		allocate(w(k))
		allocate(a(k))
		!print*,"w",w,"a",a
		call gauss_hermite(k, w, a)
		!print*,"w=",w,"a=",a
	! end GH
		xsteplength = timesteplength
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
		!print*,"time space:",t
		
	! x space range
		x(0) = 0.0d0
		do xindex = -xstep, xstep
			x(xindex) = dfloat(xindex)
		end do
		x = x * xsteplength
		!print*,"x space:",x
		
		allocate(z(0:timestep,-xstep:xstep))
		
	!last  line of Z matrix
		do xindex = -xstep, xstep
			call finalvaluez(t(timestep),x(xindex),z(timestep,xindex))
		end do 
		
! when 0<= n <=N-1

		do timeindex = timestep-1, 0, -1
		!time fixed x change ,Z
			do xindex = -xstep, xstep
				expectyw = 0.0d0
				expectz = 0.0d0	
				expectfw = 0.0d0
				expectf1fw = 0.0d0
				expectf2zzw = 0.0d0				
				!find the interpolation point
				do i = 1, k
					x_point = x(xindex)+dsqrt(2.0d0*timesteplength)*a(i)
					tindex = timeindex+1
					call interpolateIndex(x, z, x_point, xsteplength, timestep, xstep, tindex, xinterpl, zinterpl)
					call newtoninterpl(xinterpl, zinterpl, x_point, zb)
					!print*,"zb",zb
					call valuey(t(tindex),x_point,yb)
					call function34(yb,f,f1,f2)
					expectyw = expectyw +w(i)*yb*dsqrt(2.0d0*timesteplength)*a(i)
					expectz = expectz + w(i)*zb 
					expectfw = expectfw +w(i)*f*dsqrt(2.0d0*timesteplength)*a(i)
					expectf1fw = expectf1fw +w(i)*f1*f*dsqrt(2.0d0*timesteplength)*a(i)
					expectf2zzw = expectf2zzw +w(i)*f2*(zb**2)*dsqrt(2.0d0*timesteplength)*a(i)
					!print*, 'do k'
				end do
				!print*, 'end k'
				expectyw = expectyw/dsqrt(pi)
				expectz =expectz/dsqrt(pi)
				expectfw = expectfw/dsqrt(pi)
				expectf1fw = expectf1fw/dsqrt(pi)
				expectf2zzw = expectf2zzw/dsqrt(pi)
				!print*,"expectyw:=",expectyw, "expectz=", expectz,"expectfw=",expectfw,"expectf1fw=",expectf1fw,"expectf2zzw=" ,expectf2zzw
				z(timeindex,xindex) = (expectyw - (1.0d0-theta2)*timesteplength*expectz + theta2*timesteplength*expectfw & 
				&-theta2*0.5d0*(timesteplength**2)*(-expectf1fw &
				+0.5d0*expectf2zzw))/(theta2*timesteplength)
				! z(timeindex,xindex) = (expectyw - (1-theta2)*timesteplength*expectz+ theta2*timesteplength*expectfw)&
				! /(theta2*timesteplength)					
				!print*,"z(timeindex,xindex)=",z(timeindex,xindex)
				 !if (timeindex > 0) then
					 !call finalvaluez(t(timeindex), x(xindex), z_true_value)
					 !z(timeindex, xindex) = z_true_value
				 !end if
			end do 
				!print*,"z(timeindex,:)=",z(timeindex,:)
		
				
		end do 
	

	
		zerror = dabs(zture - z(0,0))
	
		deallocate(w)
		deallocate(a)
		deallocate(z)
		deallocate(t)
		deallocate(x)
	end subroutine
end module main
