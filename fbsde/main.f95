module main
!-llapck -lblas
	implicit none
	private
	public::function51sde, function51, turevalue, &
			thelengthofx, gauss_hermite,  newtoninterpl, leastsquare, li51

contains


	pure subroutine function51sde(t,x,b,tau)
		real(kind=8),intent(inout)::t, x
		real(kind=8),intent(out)::b, tau
		b = dcos(x)
		tau = (dsin(x)+3.0d0)*dcos(x)
	end subroutine function51sde
	
	pure subroutine function51(t,x,y,z,f)
		real(kind=8),intent(inout)::t, x, y, z
		real(kind=8),intent(out)::f
		f = - y - z/(dsin(x)+3.0d0) + 0.5d0*z*(1.0d0+3.0d0*dsin(x))
	end subroutine function51
	
	pure subroutine turevalue(t,x,y,z)
		real(kind=8),intent(inout)::t, x
		real(kind=8),intent(out)::y,z
		y = dexp(t)*dlog(dsin(x)+3.0d0)
		z = dexp(t)*((dcos(x))**2)
	end subroutine
		
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

	subroutine thelengthofx(x0,k,timestep,xrange)
		real(kind=8), intent(in):: x0
		integer, intent(in):: k, timestep
		real(kind=8),intent(inout):: xrange
		real(kind=8),dimension(k):: w, a
		real(kind=8):: timesteplength, xsteplength, b, tau, x_max, x_min
		real(kind=8),dimension(:),allocatable:: x, x_need
		real(kind=8),dimension(0:timestep)::t
		integer:: xstep, xindex, timeindex, xneedindex
		allocate(x_need(k))
		w(:)=0.0d0
		a(:)=0.0d0
		call gauss_hermite(k, w, a)
		print*,"w=",w
		print*,"a=",a
		timesteplength = 1.0d0/dfloat(timestep)
		print*,"timesteplength=",timesteplength
		xsteplength = timesteplength**(0.75)
		print*,"xsteplength=",xsteplength
		!xsteplength depended on the order of x and the order of t
		xstep = floor(xrange/xsteplength)
		print*,"xstep=",xstep
		allocate(x(-xstep:xstep)) 
		do timeindex = 0, timestep
			t(timeindex) = dfloat(timeindex) 
		end do 
		t = t * timesteplength
		print*,"t=",t
		x(0) = x0
		do xindex = -xstep,xstep
			x(xindex) = x(0)+dfloat(xindex)*xsteplength
		end do
		print*,"x=",x
		xneedindex = 0
		do timeindex = 0,timestep
			print*, "timeindex=", timeindex
			call function51sde(t(timeindex),x(xneedindex),b,tau)
			print*, "b=", b, "tau = ", tau 
			x_need = x(xneedindex)+b*timesteplength+tau*dsqrt(2.0d0*timesteplength)*a
			print*, "x_need=", x_need
			x_max = dabs(maxval(x_need))
			print*, "x_max=", x_max
			x_min = dabs(minval(x_need))
			print*, "x_min=", x_min
			if (x_min > x_max) then 
				x_max = x_min
			end if 
			print*, "x_max= ", x_max
			xneedindex = floor((x_max-x0)/xsteplength)+1
			print*,"xneedindex= ", xneedindex 
		end do 
		xrange = x0 + xneedindex*xsteplength
		print*,"xrange= ", xrange
		deallocate(x)
		deallocate(x_need)
	end subroutine
subroutine newtoninterpl(x, y, x_aim, y_aim)
		!finall y_aim which is our need
		real(kind=8), dimension(:), intent(inout):: x
		real(kind=8), dimension(:), intent(in):: y
		real(kind=8), intent(inout):: x_aim
		real(kind=8), intent(out):: y_aim 
		real(kind=8), dimension(:,:), allocatable:: z
		integer:: sizex, temp1, temp2,  w
		real(kind=8)::n
		sizex = size(x)
		allocate(z(sizex,sizex))
		z(:,:) = 0.0d0
		n = 0.0d0
		z(:,1) = y
		do temp2 = 2,sizex
			do temp1 = temp2, sizex
				z(temp1,temp2) = (z(temp1-1,temp2-1)-z(temp1,temp2-1))/(x(temp1-(temp2-1))-x(temp1))
			end do
		end do
		x(2:sizex) = x(1:(sizex-1))
		x(1) = x_aim - 1.0d0
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
	subroutine interpolateIndex(xvalues, yvalues, xpoint,x0, xsteplength, timesteps, xsteps, timeindex, xintvalues, yintvalues)
	!xvalues:is the all values of x   xintvalues: we need (xpoint in them) 
		real(kind=8), dimension(-xsteps:xsteps), intent(in):: xvalues
		real(kind=8), dimension(0:timesteps,-xsteps:xsteps), intent(in):: yvalues
		real(kind=8), intent(in):: xpoint,x0
		real(kind=8), intent(in):: xsteplength
		integer, intent(in):: xsteps, timeindex, timesteps
		real(kind=8), dimension(4), intent(inout)::xintvalues, yintvalues
		integer, dimension(4):: xindex
		integer:: overflow

		overflow = 0
		xindex(3) = floor((xpoint-x0)/xsteplength)
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
	subroutine li51(timestep, xrange,yerror, zerror)
	!cranck nicolson method
	!-llapck -lblas
		
		integer, intent(inout):: timestep,xrange
		real(kind=8), intent(out)::yerror, zerror


		real(kind=8), dimension(:,:), allocatable:: y, z
		real(kind=8), dimension(:),allocatable:: t, x
		real(kind=8):: x0,timesteplength, xsteplength, yb, zb

! 	
		integer::i, timeindex, xindex, xstep, k,  max_try
		integer:: overflow,tindex
		real(kind=8):: b, tau, f, expecty, expectf, expectyw, expectfw, expectz, pi, diff, y0, y1, x_point , yture, zture
		real(kind=8),dimension(:),allocatable:: w, a
		real(kind=8), dimension(4)::yinterpl,  xinterpl , zinterpl

		!print*, "start program"
		pi = 4.0d0*dATAN(1.0d0)
		x0 = 1.0d0
		timesteplength = (1.0d0)/dfloat(timestep)
		print*, "time step size  : ", timesteplength
		
		
	! GH
		k = 10
		print*, "number of steps : ", timestep
		print*, "gauss hermite   : ", k

		allocate(w(k))
		allocate(a(k))
		call gauss_hermite(k, w, a)
	! end GH
	
		xsteplength = timesteplength**(0.75)
		xstep = floor(dfloat(xrange)/xsteplength)
		!print*,"xrange:",xrange
		!print*,"xstep:",xstep
		allocate(t(0:timestep))
		allocate(x(-xstep:xstep))

	!	t sapce

		do timeindex = 0, timestep
			t(timeindex) = dfloat(timeindex) 
		end do 
		t = t * timesteplength
!		print*,"t=:",t
	! x space 
		x(0) = x0
		do xindex = -xstep,xstep
			x(xindex) = x(0)+dfloat(xindex)*xsteplength
		end do
	
		!print*,"x=:",x
		
		allocate(y(0:timestep,-xstep:xstep))
		allocate(z(0:timestep,-xstep:xstep))
		y(:,:)=0.0d0
		z(:,:)=0.0d0

	! last line of Y matrix and Z marix(timestep)
		do xindex = -xstep, xstep
			call turevalue(t(timestep),x(xindex),y(timestep,xindex),z(timestep,xindex))
		end do 
		!print*,"last Y:",y(timestep,:),"last Z:",z(timestep,:)	
		
		
! when 0<= n <=N-1

		do timeindex = timestep-1, 0, -1
			!print*,"times:",timeindex
			do xindex = -xstep, xstep
				expectyw = 0.0d0
				expectfw = 0.0d0
				expectz = 0.0d0
				expecty = 0.0d0
				expectf = 0.0d0
				
				call function51sde(t(timeindex),x(xindex),b,tau)
				do i = 1, k
					!print*,"t=:",t(timeindex),"x=:",x(xindex),"b=:",b,"tau=:",tau
					x_point = x(xindex)+b*timesteplength+tau*dsqrt(2.0d0*timesteplength)*a(i)
					!find the interpolation point
					tindex = timeindex + 1
					call interpolateIndex(x, y, x_point,x0, xsteplength, timestep, xstep, tindex,xinterpl, yinterpl)
					call newtoninterpl(xinterpl,yinterpl, x_point, yb)
					call interpolateIndex(x, z, x_point,x0, xsteplength, timestep, xstep, tindex,xinterpl, zinterpl)
					call newtoninterpl(xinterpl, zinterpl, x_point, zb)
					!print*, "zb ", zb
					call function51(t(timeindex+1),x_point,yb,zb,f)
				
					expectyw = expectyw + w(i)*yb*dsqrt(2.0d0*timesteplength)*a(i)
					expectfw = expectfw + w(i)*f*dsqrt(2.0d0*timesteplength)*a(i)
					expectz = expectz + w(i)*zb 
					expecty = expecty + w(i)*yb 
					expectf = expectf + w(i)*f
				end do 

				expectyw = expectyw/dsqrt(pi)
				expectfw = expectfw/dsqrt(pi)
				expectz = expectz/dsqrt(pi)
				expecty = expecty/dsqrt(pi)
				expectf = expectf/dsqrt(pi)		
				z(timeindex,xindex) = (expectyw+0.5d0*timesteplength*expectfw-0.5d0*timesteplength*expectz)/(0.5d0*timesteplength)

				y0 = y(timeindex+1,xindex)
				max_try = 0
				do 
					if(max_try > 10000) then
						print*, "quit cycle"
						exit
					end if
					call function51(t(timeindex),x(xindex),y0,z(timeindex,xindex),f)
					y1 = expecty + 0.5d0*timesteplength*expectf + 0.5d0*timesteplength*f
					diff=dabs(y0-y1)
					if (diff <= 1d-10) then 
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
		!print*,"yture",yture,"zture",zture
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