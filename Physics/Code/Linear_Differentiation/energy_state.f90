implicit none
double precision, parameter :: tol=1.0d-6 !maximum relative error
double precision x_left_n,x_right_n,x_middle_n,f1, f2
double precision x_iteration_old,x_iteration,relative_error
double precision :: tmp_odd,tmp_even, double_check
integer n
double precision :: r=20.0, dx=0.1, bound_1=0.00, bound_2=0.00, even_value=0.0, odd_value=0.0

open(file='/Users/phihung/Documents/PHY/first/homework/hw_02/energies.dat', unit=12) !file to record convergence
    !find energy of the tangent expression
    do while(bound_2<r) !when the bound given is not searched
        bound_2 = bound_2 + dx !update second bound
        tmp_even = f1(bound_1,r)*f1(bound_2,r) ! check if the two bounds have the same sign
        if(tmp_even>0.0)then !if yes
            continue  ! continue to search for a good range 
        else !initiate the zeroth iteration
            n=0 !zeroth iteration
            x_left_n=bound_1     !initial point on the left
            x_right_n=bound_2   !initial point on the right
            x_middle_n=(x_left_n+x_right_n)/2.0 !finding the middle point
            relative_error=1.0 !set inital error to any value higher than tol
            x_iteration_old=x_middle_n !zeroth iteration
        endif

        do while (relative_error>tol) !iterations untill relative error is less than tol
            n=n+1
            if(f1(x_left_n,r)*f1(x_middle_n,r)<=0.0)then
                x_right_n=x_middle_n !<--- the root is in left subinterval
                else
                x_left_n=x_middle_n  !<--- the root is in right subinterval
            endif
            x_middle_n=(x_left_n+x_right_n)/2.0 !finding the middle point
            x_iteration=x_middle_n
            relative_error=abs((x_iteration-x_iteration_old)/x_iteration)
            x_iteration_old=x_iteration
        enddo

        double_check = x_iteration*tan(x_iteration) !check if the value is positive
        if(even_value<x_iteration .AND. double_check>0.0)then !if satisfies the condition
            write(12,*) x_iteration  !write
            even_value = x_iteration !update
        endif
        bound_1 = bound_2 ! skip that range to initiate left bound as the old right bound
    enddo
!reset bounds
bound_1=0.00
bound_2=0.00
    !find energy of the cotangent expression
    do while(bound_2<r)
        bound_2 = bound_2 + dx !update second bound
        !tmp=f1(x_left)*f1(x_right) !check if interval is correct
        tmp_odd = f2(bound_1,r)*f2(bound_2,r) ! check if the two bounds have the same sign
        if(tmp_odd>0.0)then !if yes
            continue  ! continue until the range is good  
        else !initiate the zeroth iteration
            n=0 !zeroth iteration
            x_left_n=bound_1     !initial point on the left
            x_right_n=bound_2   !initial point on the right
            x_middle_n=(x_left_n+x_right_n)/2.0 !finding the middle point
            relative_error=1.0 !set inital error to any value higher than tol
            x_iteration_old=x_middle_n !zeroth iteration
        endif
        do while (relative_error>tol) !iterations untill relative error is less than tol
            n=n+1
            if(f2(x_left_n,r)*f2(x_middle_n,r)<=0.0)then
                x_right_n=x_middle_n !<--- the root is in left subinterval
                else
                x_left_n=x_middle_n  !<--- the root is in right subinterval
            endif
            x_middle_n=(x_left_n+x_right_n)/2.0 !finding the middle point
            x_iteration=x_middle_n
            relative_error=abs((x_iteration-x_iteration_old)/x_iteration)
            x_iteration_old=x_iteration
        enddo
        double_check = -x_iteration/tan(x_iteration)
        if(odd_value<x_iteration .AND. (-x_iteration/tan(x_iteration))>0.0)then
            write(12,*) x_iteration 
            odd_value = x_iteration
        endif

        bound_1 = bound_2 ! skip that range to initiate left bound as the old right bound
    enddo
close(unit=12)
end

double precision function f1(x,z_0)
  implicit none
  double precision x
  double precision z_0
  f1=(x*tan(x))**2+x**2-z_0**2
end function f1

double precision function f2(x,z_0)
  implicit none
  double precision x
  double precision z_0
  f2=(-x/tan(x))**2+x**2-z_0**2
end function f2
