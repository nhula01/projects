implicit none
double precision, parameter :: tol=1.0d-6 !maximum relative error
double precision, parameter :: x_left=-0.0,x_right=1000 !initial interval
                                  !function must have different signs on the left and right
double precision f,x_left_n,x_right_n,x_middle_n
double precision x_iteration_old,x_iteration,relative_error
double precision tmp
integer n

tmp=f(x_left)*f(x_right) !check if interval is correct
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

open(file='/Users/phihung/Documents/PHY/first/homework/hw_02/bisectional_-0_1000.dat',unit=11) !file to record convergence

  do while (relative_error>tol) !iterations untill relative error is less than tol
    n=n+1
    if(f(x_left_n)*f(x_middle_n)<=0.0)then
        x_right_n=x_middle_n !<--- the root is in left subinterval
      else
        x_left_n=x_middle_n  !<--- the root is in right subinterval
    endif
    x_middle_n=(x_left_n+x_right_n)/2.0 !finding the middle point
    x_iteration=x_middle_n
    relative_error=abs((x_iteration-x_iteration_old)/x_iteration)
    x_iteration_old=x_iteration
    write(11,*) n,x_iteration,f(x_iteration)
  enddo
close(unit=11)

end

double precision function f(x)
  implicit none
  double precision pi
  double precision x
  pi = 4.0*atan(1.0)
  f=(sin(x**2-pi/2))**2-x
end function f