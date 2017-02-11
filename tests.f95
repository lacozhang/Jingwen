program tests
	implicit none
	interface
		subroutine the_length_ofx(k, step, timesteplength,x)
			integer, intent(in):: k, step,timesteplength
			real(kind=8), intent(out)::x
		end subroutine
	end interface
	integer::k, step, timesteplength
	real(kind=8)::x
	k = 2
	step = 3
	timesteplength = 0.5
	call the_length_ofx(k, step, timesteplength, x)
	print*,"x=",x
end program