module main
	implicit none
	private
	public:: gauss_hermite, function34,valuey,valuez, newtoninterpl, leastsquare, li34y

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
		f = - (y**3) + 2.5* (y ** 2) - 1.5*y

		
	end subroutine
	pure subroutine valuey(t,x,y)
		real(kind=8),intent(inout)::t,x
		real(kind=8),intent(out)::y
		!y = (1.0)/(1+dexp(-x-0.25*t))		
		y = dexp(t+x)/(dexp(t+x)+1.0d0)
	end subroutine
	pure subroutine valuez(t,x,z)
		real(kind=8),intent(inout)::t,x
		real(kind=8),intent(out)::z
		!z = (dexp(-x-0.25*t))/((1.0d0+dexp(-x-0.25*t))**2)
		z = dexp(x+t)/((dexp(x+t)+1.0d0)**2.0d0)
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
		x(1) = x_aim -1.0d0
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
		timesteplength = 1.0d0/(timestep)
		error = dlog(error)
		timesteplength = dlog(timesteplength)
		a = 1.0d0
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

	subroutine li34y(timestep,xrange, yerror)
		integer, intent(inout):: timestep
		real(kind=8), intent(out)::yerror

		real(kind=8), dimension(:,:), allocatable:: y
		real(kind=8), dimension(:),allocatable:: t, x
		real(kind=8):: timesteplength, xsteplength, yb
		integer, intent(in)::xrange 
! 	
		integer:: i, timeindex, xindex, xstep, k, x_point1, x_point2,x_point3,x_point4
	
		real(kind=8):: f
		real(kind=8):: expecty3, expecty1, expectfy1, expecty2, expectfy2, expectfy3
		real(kind=8):: pi, diff, y0, y1, x_point, y_true_value 
		real(kind=8), dimension(1):: temp_point, temp_y
		real(kind=8), dimension(:),allocatable:: w, a
		real(kind=8), dimension(4)::yinterpl, xinterpl
		real(kind=8):: yture
		integer:: cycleindex
		
		pi = 4.0d0*datan(1.0d0)
		!print*,"pi=",pi
		print*, "number of steps : ", timestep
		timesteplength = (1.0d0)/timestep
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
		
		xsteplength = (timesteplength)**(1.25)
		!!!!!!!!be detemined by the order of x space and t space
		!print*, "xsteplength   = : ", xsteplength
		xstep = (xrange+ xsteplength) /xsteplength
		!print*,"xstep=",xstep 
		
		
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
			x(xindex) = xindex
		end do
		x = x * xsteplength
		!print*,"x space:",x
		
		
		allocate(y(0:timestep,-xstep:xstep))
		
		
	!the starting  line of  Y matrix 
		do timeindex = timestep, timestep-2, -1
			do xindex = -xstep, xstep
				call valuey(t(timeindex),x(xindex),y(timeindex,xindex))
			end do 
		end do 
! when 0<= n <=N-3

		do timeindex = timestep-3, 0, -1
		! time fixed x change ,Y 	
			do xindex = -xstep ,xstep
				expecty3 = 0.0d0
				expectfy1 = 0.0d0
				expectfy2 = 0.0d0
				expectfy3 = 0.0d0
				
				
				
				do i = 1,k
					!expectfy1
!					print*, "timeindex ", timeindex, " xindex ", xindex
					x_point = x(xindex)+dsqrt(2.0d0*timesteplength)*a(i)
					call interpolateIndex(x, y, x_point, xsteplength, timestep, xstep, timeindex+1, xinterpl, yinterpl)
					call newtoninterpl(xinterpl,yinterpl, x_point, yb)				
					call function34(yb,f)				
					expectfy1 = expectfy1 +w(i)*f
					
					!expectfy2
					x_point = x(xindex)+dsqrt(4.0d0*timesteplength)*a(i)
					call interpolateIndex(x, y, x_point, xsteplength, timestep, xstep, timeindex+2, xinterpl, yinterpl)
					call newtoninterpl(xinterpl,yinterpl, x_point, yb)	
					call function34(yb,f)				
					expectfy2 = expectfy2 +w(i)*f
					
					
!
					!expecty3.expectfy3
					x_point = x(xindex)+dsqrt(6.0d0*timesteplength)*a(i)
					call interpolateIndex(x, y, x_point, xsteplength, timestep,  xstep, timeindex+3, xinterpl, yinterpl)
					call newtoninterpl(xinterpl,yinterpl, x_point, yb)	
					call function34(yb,f)
					expecty3 = expecty3 +w(i)*yb
					expectfy3 = expectfy3 +w(i)*f					
				end do

				expecty3 = expecty3/dsqrt(pi)
				expectfy1 = expectfy1/dsqrt(pi)
				expectfy2 = expectfy2/dsqrt(pi)
				expectfy3 = expectfy3/sqrt(pi)
			    !print*,"expecty= ",expecty ,"expectf=",expectf,"expectf1f=",expectf1f,"expectf2zz=",expectf2zz
				
				y0 = y(timeindex+1,xindex)
				!print*,"y0=",y0
				! cycleindex = 0
				do 
					! cycleindex = cycleindex + 1
					!print*, 'circle ', cycleindex
					call function34(y0,f)
					
					y1 = expecty3 + 3.0d0*timesteplength*((1.0d0/8.0d0)*f+(3.0d0/8.0d0)*expectfy1 &
					+(3.0d0/8.0d0)*expectfy2+(1.0d0/8.0d0)*expectfy3)										
					
					diff = dabs(y1-y0)

					if (diff <= 1d-15) then
						EXIT
					end if					
					y0 = y1

					! if (cycleindex > 10000) then
						! EXIT
					! end if
				end do 
				!print*,"cycleindex:",cycleindex
				y(timeindex,xindex) = y1 
				call valuey(t(timeindex), x(xindex), y_true_value)
				! if (timeindex > 0) then
					! y(timeindex, xindex) = y_true_value
				! end if
				!print*,"y(timeindex,xindex)=",y(timeindex,xindex)		
			end do 
		!print*, 'end compute y'
		!print*,"y(timeindex,:)=",y(timeindex,:)
		end do

		yture = 0.5
		yerror = dabs(yture - y(0,0))	
		
		deallocate(w)
		deallocate(a)
		deallocate(y)
		deallocate(t)
		deallocate(x)
	end subroutine
end module main
