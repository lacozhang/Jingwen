module delay
	implicit none
	private
	public::randmnuniform, randn, brownian, numericalx
	
contains
	subroutine randmnuniform(x)
		real(kind = 8), dimension(:),intent(inout) :: x
        call random_seed () ! 系统根据日期和时间随机地提供种子
		call random_number (x) ! 每次的随机数就都不一样
	end subroutine randmnuniform

	subroutine randn(x1,x2)
		real(kind = 8),intent(out) :: x1, x2 
		real(kind = 8), dimension(2) :: x
		real(kind = 8) :: pi, u1, u2
		pi = 4*atan(1.0)
		x = (/0, 1/)
		!print*,"early x =", x
		call  randmnuniform(x)
		!print*,"after x = ", x
		u1 = x(1)
		u2 = x(2)
		x1 = sqrt(-2*log(u1)) * cos(2*pi*u2);
		x2 = sqrt(-2*log(u1)) * sin(2*pi*u2);
	end subroutine randn
	
	subroutine brownian (t,n,w)
		real(kind = 8), intent(inout):: t
		integer, intent(inout)::  n
		real(kind = 8), dimension(0:n), intent(inout):: w
		real(kind = 8), dimension(1:n):: dw
		real(kind = 8):: dt, x1, x2 
		integer::i
		dt = t/n
		print*, "dt =",dt
		w(0) = 0
		call randn(x1,x2)
		print*, "x1 =", x1
		do i = 1,n
			dw(i) = sqrt(dt)*x1
			w(i) = w(i-1) + dw(i)
		end do 
		
		
	end subroutine brownian
	
	subroutine numericalx(x0,n,x)
		real(kind = 8), intent(in)::x0
		integer, intent(in)::n
		real(kind = 8), dimension(0:n), intent(inout):: x
		x(0) = x0
	end subroutine  numericalx
	
end module