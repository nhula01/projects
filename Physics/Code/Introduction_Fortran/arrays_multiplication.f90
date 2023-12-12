implicit none
integer, parameter :: array1_size1=3,array1_size2=3 !dimension for A
double precision A(array1_size1,array1_size2)
integer, parameter :: array2_size1=3,array2_size2=3 !dimension for B
double precision :: B(array2_size1,array2_size2), T(array1_size1,array2_size2), D(array1_size1,array2_size2)
double precision :: C(array1_size1,array2_size2) !if we calculate A1xA2 => dimensions of C are [array1_size1 x array2_size2]
                                               !note that array1_size2=array2_size1 [otherwise you cannot calculate A1xA2]
double precision tmp
integer n,n1,n2

!quick check to ensure everything is correct

if(array1_size2.ne.array2_size1)then
  write(*,*) 'check A and B dimensions!'
  stop !incorrect dimensions, program is terminated
endif

!define the matrix A and B
A = reshape([1,2,3,&
             3,2,1,&
             1,1,1],[array1_size1, array1_size2])
B = reshape([1,1,1,&
             1,2,3,&
             3,2,1],[array2_size1, array2_size2])

do n1=1,array1_size1
  do n2=1,array2_size2
    tmp=0.0
    do n=1,array1_size2 !or array2_size1 'cause they must be the same
      tmp=tmp+A(n1,n)*B(n,n2) !
    enddo
    C(n1,n2)=tmp
  enddo
enddo

!try the builtin function
T=matmul(A,B)
D=matmul(B,A)
D=T-D !the commutator 

open(file='/Users/phihung/Documents/PHY/first/homework/hw_01/matrix_multiplication.dat',unit=10)
write(10,*) "Matrix multiplication of A and B gives"
write(10,*) "using logic and loops"
write(10,*) C(3,3), C(2,3), C(1,3)
write(10,*) C(3,2), C(2,2), C(1,2)
write(10,*) C(3,1), C(2,1), C(1,1)
write(10,*) "using built-in functions"
write(10,*) T(3,3), T(2,3), T(1,3)
write(10,*) T(3,2), T(2,2), T(1,2)
write(10,*) T(3,1), T(2,1), T(1,1)
!do n1=1,array1_size1
 ! write(10,*) (C(n1,n2),n2=1,array2_size2)
!enddo
write(10,*) "The commutator of A and B gives"
write(10,*) D(3,3), D(2,3), D(1,3)
write(10,*) D(3,2), D(2,2), D(1,2)
write(10,*) D(3,1), D(2,1), D(1,1)
close(unit=10)


end
