program student_list
	use student_class
	implicit none
	real :: mark
	call create_student (1, "Peter Peterson", 8.5)
	call create_student (2, "John Johnson", 6.3)
	call get_mark ("Peter Peterson", mark)
	call delete_student ("John Johnson")
	print*, "Peter Peterson:", mark
end program student_list