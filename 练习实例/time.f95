program timew
	implicit none
	integer:: time, hour, minute, second
	read*, time
	hour = time/3600
	time = mod(time, 3600)
	print*, time 
	minute = time/60
	second = mod(time, 60)
	print*, "hour:",hour , "minute:",minute, "second:",second
end program	
	
	