module thetaderivate
	implicit none
	private
	public::functionf, functionfture, leastsquare, thetaderi
	
	
contains
pure subroutine functionf(f,ft,fy,t,y)
	real(kind=8),intent(in)::t,y
	real(kind=8),intent(out)::f,ft,fy
	f = (2.0/(t+1))*y+(t+1)**2*exp(t)
	ft = (-2.0)/((t+1)**2)*y+2*(t+1)*exp(t)+exp(t)*(t+1)**2
	fy = (2.0)/(t+1)
end subroutine
pure subroutine functionfture(t,y)
	real(kind=8),intent(inout)::t
	real(kind=8),intent(out)::y
	y = (t+1)**2*exp(t)
end subroutine
subroutine leastsquare(timestep, error, order)
	!sizeerror >= 2
		real(kind=8), dimension(:),intent(inout)::error
		integer, dimension(:),intent(inout)::timestep
		real(kind=8), intent(out)::order
		real(kind=8), dimension(:), allocatable::timesteplength
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
subroutine thetaderi( timestep, theta,  y_aim, error)
	integer,intent(inout):: timestep
	real(kind=8), intent(inout):: theta
	real(kind=8), intent(out)::y_aim, error
	real(kind=8), dimension(:), allocatable::t, y, yture
	
	integer::i
	real(kind=8)::f1,f1t,f1y,f2,f2t,f2y,ystart,y1,h, diff
	h = 1.0/timestep
	print*,"h=",h
	allocate(y(0:timestep))
	allocate(t(0:timestep))
	allocate(yture(0:timestep))
	
	
	
	do i=0,timestep
		t(i) = i	
	end do 
	t = t*h
	print*,"t=",t
	
	do i=0,timestep
		call functionfture(t(i),yture(i))
	end do 
	print*,"yture=",yture
	
	call functionfture(t(0),y(0)) 
	
	print*,"y0",y(0)
	print*,"before cycle the theta is:",theta 
	do i=1,timestep
		call functionf(f1,f1t,f1y,t(i-1),y(i-1))		
		ystart=y(i-1)+h*f1
		
		do
			print*,"the beginning value:",ystart
			call functionf(f2,f2t,f2y,t(i),ystart)
			y1 = y(i-1)+h*theta*f1+h*(1-theta)*f2+((3*theta-1)/6.0) &
					*(h**2)*(f1t+f1*f1y)+((3*theta-2)/6.0)*(h**2)*(f2t+f2*f2y)
			print*,"y1",y1
			diff = abs(y1-ystart)
			print*,"diff",diff
			if(diff <= 1e-10) then
				exit
			end if
			ystart = y1
			
		end do 
		y(i) = ystart
	end do 
	print*,"y=",y
	
	y_aim = y(timestep)
	error =abs(y_aim-yture(timestep))
	deallocate(y)
	deallocate(t)
	deallocate(yture)

end subroutine

end module