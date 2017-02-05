program pool
	implicit none
	interface 
		recursive subroutine remove(n, start, dest, help)
			integer, intent(in):: n, start, dest, help
		end subroutine
	end interface
	call remove(3, 2, 3, 1) 
end program
recursive subroutine remove(n, start, dest, help)
	implicit none
	integer, intent(in):: n, start, dest, help
	select case (n)
		case (1)
		print*, "move ", n, " from ", start, " to ", dest
		case default
		call remove(n-1, start, help, dest)
		print*, "move ", n," from ", start, " to ", dest
		call remove(n-1, help, dest, start)
	end select
end subroutine