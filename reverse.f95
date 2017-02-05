program reverse
	implicit none
	integer, dimension(10):: array
	integer:: i, j, k, l, d, h, o
	do i=1,10
		array(i)=i**2
	end do
	print*,array
	read*, j, k
	d=j-k
	print*,d
	if (d==0) then
		print*,"wuxiao"
	else if  (d<0) then
		d=j
		j=k
		k=d
	end if
	d= j-k
	print*,"genggaihoudeshunxv:", k, j
	l=mod(d,2)
	print*,"mod",l
	if (l==0)then
		d=d/2-1
		do i=0,d
		o=array(k+i)
		array(k+i)=array(j-i)	
		array(j-i)=o
		end do
		print*,array(k:j)
	else
	    d=(d-1)/2-1
		do i=0,d
		o=array(k+i)
		array(k+i)=array(j-i)	
		array(j-i)=o
		end do
		print*,array(k:j)
	end if
	print*,array
end program 