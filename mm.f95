program mm
	implicit none
	real,dimension(5,5):: m1, m2,m3
	integer:: i, j, l, he,k
	do i=1,5
		do j=1,5
			m1(i,j)=i+j
			m2(i,j)=i*j
		end do 
	end do
	
	do i=1,5
		do j=1,5
			he=0
			do l=1,5
		     	he=he+m1(i,l)*m2(l,j)
			end do 
		m3(i, j)=he
		end do 
	end do 
	print*, m1, m2, m3
	
	iloop: do i = 1, 3
		print*, "i: ", i
		jloop: do j = 1, 3
			print*, "j: ", j
			kloop: do k = 1, 3
				print*, "k: ", k
				if (k==2) then
				exit jloop
				end if
			end do kloop
		end do jloop
	end do iloop
	
end program
		
	
