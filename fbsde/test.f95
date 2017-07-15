program test
	use main
	
	implicit none
	integer:: timestep
	real(kind=8)::yerror, zerror
	timestep = 4
	call li51(timestep, yerror, zerror)
	print*,"yerror", yerror
	print*, "zerror",zerror 
end program test