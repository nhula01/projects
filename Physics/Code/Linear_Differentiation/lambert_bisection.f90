implicit none
double precision, parameter :: tol=1.0d-6 ! maximum relative error
double precision, parameter :: x_left=-10.0,x_right=-1.0 ! initial interval to check for intersection
! function must have different signs on the left and right
double precision :: f,x_left_n,x_right_n,x_middle_n,z
double precision :: x_iteration_old,x_iteration,relative_error
double precision :: tmp,tmp1,tmp2
integer, parameter :: Nz=1000 ! creating a thoundsand points
double precision, parameter :: z1=-0.3678794,z2=-.01 ! going from -e^(-1) to -.01
double precision :: W1(Nz) ! create the W0 array of 1000 points
integer n,k,iterations(Nz) ! See each iteration

do k=1,Nz ! going through a thoundsand values
  z=z1+(z2-z1)*float(k-1)/float(Nz-1) ! create each value of z

  ! find the value of f=W0*exp(W0)-z where W0=x_left/x_right
  call lambert(z,x_left,tmp1)
  call lambert(z,x_right,tmp2)

  ! check if interval is correct
  tmp=tmp1*tmp2 
  if(tmp>0.0)then
    write(*,*) 'function must have different signs in x_left and x_right'
    write(*,*) 'exit with error'
    stop
  endif

  n=0 !zeroth iteration
  x_left_n=x_left     !initial point on the left
  x_right_n=x_right   !initial point on the right
  x_middle_n=(x_left_n+x_right_n)/2.0 !finding the middle point
  relative_error=1.0 !set inital error to any value higher than tol
  x_iteration_old=x_middle_n !zeroth iteration

  do while (relative_error>tol) !iterations untill relative error is less than tol
    n=n+1 ! update the iteration values

    !calculate the value f(w) again on both ends
    call lambert(z,x_left_n,tmp1)
    call lambert(z,x_middle_n,tmp2)

    !Update the boundary
    if(tmp1*tmp2<=0.0)then 
        x_right_n=x_middle_n !<--- the root is in left subinterval
      else
        x_left_n=x_middle_n  !<--- the root is in right subinterval
    endif

    !update data
    x_middle_n=(x_left_n+x_right_n)/2.0 !finding the middle point
    x_iteration=x_middle_n 
    relative_error=abs((x_iteration-x_iteration_old)/x_iteration) !find the differences btw two values
    x_iteration_old=x_iteration 
  enddo

  !update the value of each z for W1 and the iteration array
  W1(k)=x_iteration
  iterations(k)=n
enddo

open(file='/Users/phihung/Documents/PHY/first/homework/hw_02/lambert_function_W1.dat',unit=32) !file to record Lambert function W0
do k=1,Nz
  z=z1+(z2-z1)*float(k-1)/float(Nz-1)
  write(32,*) z,W1(k),iterations(k)
enddo
close(unit=32)

end

subroutine lambert(z,W1,f)
  implicit none
  double precision z,W1,f
  f=W1*exp(W1)-z
end subroutine lambert
