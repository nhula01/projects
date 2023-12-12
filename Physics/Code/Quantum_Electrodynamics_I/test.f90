program StringConcatenation
    implicit none
    integer :: i
    character(len=5) :: str 
    character(len=3) :: str2 = "2_0"

    ! Concatenate the strings
    character(len=8) :: result

    ! Print the result
    do i=300,1000,100
        write(str, '(I5.5)') i ! convert to string
        result = str // str2
        write(*,*) result
    enddo 

end program StringConcatenation
