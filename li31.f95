module example
	implicit none
	private
	public::li31


contains

subroutine li31(timestep, yerror, zerror)
	!moudle interpolation gaussherimite
	!bublle
	!-llapck -lblas
	use interpolation
	use gausshermite



	integer, intent(inout):: timestep
	real(kind=8), intent(out)::yerror, zerror
	
	real(kind=8), dimension(:,:), allocatable:: y, z
	real(kind=8),dimension(:),allocatable:: t, x
	real(kind=8):: timesteplength, xsteplength, yb

! 	
	integer::  i, timeindex, xindex, xstep, k, x_pointlow, x_pointhigh
	integer:: overflow
	real(kind=8):: expecty, expectyw, pi, diff, y0, y1, x_point
	real(kind=8),dimension(:),allocatable:: w, a
	real(kind=8), dimension(2)::yinterpl, xinterpl
	real(kind=8), parameter:: yture=0.5, zture=0.25
	

	pi = 4*atan(1.0)

	timestep = 8
	timesteplength = (1.0)/timestep
	print*, "time step size  : ", timesteplength
! GH
	k = 10
	print*, "number of steps : ", timestep
	print*, "gauss hermite   : ", k

	allocate(w(k))
	allocate(a(k))
	call gauss_hermite(k, w, a)
	print*, "w=", w,"a=", a
! end GH

	call thelengthofx(k, timestep, timesteplength,xstep)
	print*, "x radius        : ", xsteplength
	xsteplength = timesteplength
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
	allocate(z(0:timestep-1,-xstep:xstep))

! last line of Y matrix
	do xindex = -xstep, xstep
		y(timestep,xindex) = exp(x(xindex)+t(timestep))/(exp(x(xindex)+t(timestep))+1)
	end do 
! last but one line of Y matrix
	do xindex = -xstep, xstep
		expecty=0
		do i = 1, k 
			expecty = expecty+w(i)*(exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))&
			/(exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))+1))
		end do 
		expecty = expecty/sqrt(pi)
		y0 = y(timestep,xindex)
		do 
			y1 = expecty + timesteplength*(-y0**3+2.5*y0**2-1.5*y0)
			diff = abs(y1-y0)
			if (diff <= 1e-8) then 
				exit
			end if
			y0 = y1
		end do 
		y(timestep-1,xindex) = y1
	end do 

!last  line of Z matrix
	do xindex = -xstep, xstep
		expectyw = 0
		do i = 1, k
			expectyw = expectyw+w(i)*(exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))&
			/(exp(x(xindex)+sqrt(2*timesteplength)*a(i)+t(timestep))+1))*sqrt(2*timesteplength)*a(i)
		end do 
		expectyw = expectyw/sqrt(pi)
		z(timestep-1,xindex) = expectyw/timesteplength
	end do 


	do timeindex = timestep-2, 0, -1
		do xindex = -xstep, xstep
			expectyw = 0
			expecty = 0
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

				yinterpl = (/y(timeindex+1,x_pointlow),y(timeindex+1,x_pointhigh)/)
				xinterpl = (/x(x_pointlow),x(x_pointhigh)/)
				call newtoninterpl(xinterpl,yinterpl, x_point)
				yb=x_point

				expecty = w(i)*yb+expecty
				expectyw = w(i)*yb*sqrt(2*timesteplength)*a(i)+expectyw
			end do 
			expecty = expecty/sqrt(pi)
			expectyw = expectyw/sqrt(pi)

			z(timeindex,xindex) = expectyw/timesteplength

			y0 = y(timeindex+1,xindex)
			do 
				y1 = expecty + timesteplength*(-y0**3+2.5*(y0**2)-1.5*y0)
				diff = abs(y1-y0)
				if (diff <= 1e-8) then 
					exit
				end if
				y0 = y1
			end do
			y(timeindex,xindex) = y1
		end do 
	end do 
	print*,"z(0,0)=",z(0,0)
	print*,"y(0,0)=",y(0,0)
	zerror = abs(zture - z(0,0))
	yerror = abs(yture - y(0,0))
	print*,"yerror", yerror
	print*, "zerror",zerror
	deallocate(w)
	deallocate(a)
	deallocate(y)
	deallocate(z)
	deallocate(t)
	deallocate(x)
end module example