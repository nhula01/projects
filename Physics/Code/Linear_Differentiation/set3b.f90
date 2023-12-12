implicit none
double precision, parameter :: tol=1.0d-6 !maximum relative error
double precision, parameter :: x_initial=-10. !initial guess
double precision f,df_dx
double precision x_iteration_old,x_iteration,relative_error
double precision tmp
integer n

n=0 !zeroth iteration
call f_and_deriv(x_initial,f,df_dx)
x_iteration_old=x_initial-f/df_dx
relative_error=1.0 !set inital error to any value higher than tol

!newton method
open(file='/Users/phihung/Documents/PHY/first/homework/hw_02/Newton_Raphson.dat',unit=11) !file to record convergence
  do while (relative_error>tol) !iterations untill relative error is less than tol
    n=n+1
    call f_and_deriv(x_iteration_old,f,df_dx)
    x_iteration=x_iteration_old-f/df_dx
    relative_error=abs((x_iteration-x_iteration_old)/x_iteration)
    x_iteration_old=x_iteration
    write(11,*) n,x_iteration,f
  enddo
close(unit=11)   
end

!define function
subroutine f_and_deriv(x,f,df_dx)
  implicit none
  double precision x,f,df_dx, pi
  pi = 4.0*atan(1.0) !define pi
  f=(sin(x**2-pi/2))**2-x !function
  df_dx=4*x*sin(x**2-pi/2)*cos(x**2-pi/2)-1 !the derivative
end subroutine f_and_deriv
