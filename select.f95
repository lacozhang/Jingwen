program selector
	implicit none
	integer :: readnumber, oddnumber
	
	interface
		function test(input)
			integer:: test
			integer:: input
		end function test
	end interface
	
	oddnumber = 0
	print*,"please,input integer number"
	do
		read*, readnumber
		if (readnumber < 0) then
			exit
		end if
		if (mod(readnumber, 2) == 1) then
			oddnumber=oddnumber + 1
		end if
		print*, oddnumber
	end do
	readnumber = test(20)
	print*, readnumber
end program

function test(input)
	implicit none
	integer:: test
	integer:: input
	test = 2 * input
end function test