module student_class
	implicit none
	private
	public :: create_student, get_mark, delete_student
	type student_data
		character (len=50) :: name
		real :: mark
	end type student_data
	type (student_data), dimension (100) :: student
contains
	subroutine create_student (student_n, name, mark)
	! here some code to set the name and mark of a student
		integer , intent(in):: student_n
		character(len= *), intent(in):: name
		real, intent(in):: mark
		student(student_n)% name = name
		student(student_n)% mark = mark
	end subroutine create_student
	subroutine get_mark (name, mark)
	! dummy arguments
		character (len=*), intent (in) :: name
		real, intent (out) :: mark
	! local variables
		integer :: i
		do i = 1,100
			if (student(i)%name == name) then
				mark = student(i)%mark
			end if
		end do
	end subroutine get_mark
	subroutine delete_student (student_name)
		character (len= *), intent (in) :: student_name
		integer :: i
		do i= 1, 100
		    if (student(i)%name == student_name) then
				student(i)%name = ""
				student(i)%mark = 0
				print*, "gaoding!"
			end if
		end do
	end subroutine delete_student
		
end module student_class