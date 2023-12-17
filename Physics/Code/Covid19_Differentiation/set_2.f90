subroutine exact_derivative(x,df_dx) !a
    implicit none 
    double precision, INTENT(IN) :: x
    double precision, INTENT(OUT) :: df_dx
    df_dx = sinh(x)/(cosh(x))
end subroutine exact_derivative

subroutine generate_exact(h, N_max) !a
    implicit none
    integer :: n 
    integer, INTENT(IN) :: N_max
    double precision :: x_min=-3, x_max=3 !define the range of x
    double precision :: x_new, h, df_dx
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_exact_problem2.dat',unit=11) !file to record convergence
        do n=1, N_max
            x_new=x_min+(x_max-x_min)*float(n-1)/float(N_max-1) ! generate data
            call exact_derivative(x_new, df_dx)
            write(11,*) x_new, df_dx
        enddo
    close(unit=11)
end subroutine generate_exact

subroutine function(x,f)
    implicit none 
    double precision, INTENT(IN) :: x
    double precision, INTENT(OUT) :: f
    f = LOG(cosh(x))
end subroutine function

subroutine fourthorder(h,N_max) ! calculate central !a
    implicit none
    double precision :: x_min=-3, x_max=3 !define the range of x
    double precision :: pi, tmp1, tmp2, df_dx, x_new, f, h
    integer :: N_max,n 
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_fourthorder_problem2.dat',unit=11) !file to record convergence
        do n=2,N_max !start at two because we subtract -2h -> out of range
             ! find new value of x
            x_new=x_min+(x_max-x_min)*float(n-1)/float(N_max-1)
            ! fourth order derivative
            call function(x_new-2*h, tmp1)     !function at x_point+dx
            df_dx=tmp1
            call function(x_new-h, tmp2)        !function at x_new
            df_dx=df_dx-8*tmp2              !forward difference
            call function(x_new+h, tmp2)
            df_dx=df_dx+8*tmp2 
            call function(x_new+2*h, tmp2)
            df_dx=(df_dx-tmp2)/(12*h)
            !record the value x and df_dx
            write(11,*) x_new, df_dx 
        enddo
    close(unit=11)
end subroutine fourthorder

subroutine df_dx_and_d2f_dx2(x,df_dx,d2f_dx2) !b
    implicit none 
    double precision, INTENT(IN) :: x
    double precision, INTENT(OUT) :: df_dx, d2f_dx2
    !ln(coshx)
    df_dx = sinh(x)/(cosh(x))
    d2f_dx2 = 1- sinh(x)/(cosh(x)**2)
end subroutine df_dx_and_d2f_dx2

subroutine root_nr !b
    implicit none
    double precision, parameter :: tol=1.0d-6 !maximum relative error
    double precision, parameter :: x_initial=10
    double precision :: x_iteration,relative_error, x_iteration_old
    double precision :: df_dx, d2f_dx2, pi
    double precision :: x1=0,x2
    integer :: n=0

    call df_dx_and_d2f_dx2(x_initial,df_dx,d2f_dx2)
    relative_error=1.0 !set inital error to any value higher than tol
    x_iteration_old=x_initial-df_dx/d2f_dx2 !zeroth iteration
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_iteration_2.dat',unit=11)
        do while (relative_error>tol) !iterations untill relative error is less than tol
            n=n+1
            call df_dx_and_d2f_dx2(x_iteration_old,df_dx,d2f_dx2)
            x_iteration=x_iteration_old-df_dx/d2f_dx2 
            relative_error=abs((x_iteration-x_iteration_old)/x_iteration)
            x_iteration_old=x_iteration
            write(11,*) n, x_iteration
        enddo
    close(unit=11)
end subroutine root_nr

program hw3_2
    implicit none
    double precision :: h=3
    integer :: N_max = 5000
    call fourthorder(h, N_max)
    call generate_exact(h, N_max)
    call root_nr
end program hw3_2