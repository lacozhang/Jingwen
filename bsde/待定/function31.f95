module function31
	implicit none
	private
	public::function3
contains
	
	pure subroutine function3(y,f,f1,f2)
		real(kind=8),intent(inout)::y
		real(kind=8),intent(out)::f,f1,f2
		f = -y**3+2.5*y**2-1.5*y
		f1 = -3*y**2+5*y-1.5
		f2 = -6*y+5
		
	end subroutine
end module function31