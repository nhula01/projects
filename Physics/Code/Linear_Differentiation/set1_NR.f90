implicit none
double precision, parameter :: tol=1.0d-6 !maximum relative error
double precision, parameter :: x_initial=-2.51 !initial interval !-3->-2
!function must have different signs on the left and right
double precision f,x_left_n,x_right_n,x_middle_n,z,dfdx1
double precision x_iteration_old,x_iteration,relative_error
double precision tmp,tmp1,tmp2
integer, parameter :: Nz=10000
double precision, parameter :: z1=-0.3678794,z2=-.02
double precision W1(Nz)
integer n,k,iterations(Nz)

do k=1,Nz
  z=z1+(z2-z1)*float(k-1)/float(Nz-1) ! find z value 

  n=0 !zeroth iteration
  call lambert(z,x_initial,tmp1,dfdx1) ! fidn the value of function and the derivative
  relative_error=1.0 !set inital error to any value higher than tol
  x_iteration_old=x_initial-tmp1/dfdx1 !zeroth iteration
  do while (relative_error>tol) !iterations untill relative error is less than tol
    n=n+1
    call lambert(z,x_iteration_old,tmp1,dfdx1)
    x_iteration=x_iteration_old-tmp1/dfdx1
    relative_error=abs((x_iteration-x_iteration_old)/x_iteration)
    x_iteration_old=x_iteration
  enddo
  W1(k)=x_iteration
  iterations(k)=n
enddo

open(file='/Users/phihung/Documents/PHY/first/homework/hw_02/lambert_function_W1_NR.dat',unit=32) !file to record Lambert function W0
do k=1,Nz
  z=z1+(z2-z1)*float(k-1)/float(Nz-1)
  write(32,*) z,W1(k),iterations(k)
enddo
close(unit=32)

end

subroutine lambert(z,W1,f,df_dx)
  implicit none
  double precision z,W1,f, df_dx
  f=W1*exp(W1)-z
  df_dx=W1*exp(W1)+exp(W1)
end subroutine lambert