implicit none

integer, parameter :: row = 3, column = 3
integer :: A(row, column), B(row, column), m, i, j
integer :: C_1(row,column), C_2(row,column), result
A = reshape([1,2,3,&
             3,2,1,&
             1,1,1],[row, column])
B = reshape([1,1,1,&
             1,2,3,&
             3,2,1],[row, column])
! using builtin functions

C_1 = matmul(A,B) ! matrix multiplication
C_2 = matmul(A,B)-matmul(B,A) ! commutator
write(*,*) C_1(1,1), C_1(1,2), C_1(1,3)
write(*,*) C_1(2,1), C_1(2,2), C_1(2,3)
write(*,*) C_1(3,1), C_1(3,2), C_1(3,3)

!using logical

do i=1, row
    do j=1, column
        do m = 1, 3
            result = result + A(j,m) * B(m,i)
        enddo
        C_1(i,j) = result
        result = 0
    enddo
enddo

write(*,*) C_1(1,1), C_1(1,2), C_1(1,3)
write(*,*) C_1(2,1), C_1(2,2), C_1(2,3)
write(*,*) C_1(3,1), C_1(3,2), C_1(3,3)

end