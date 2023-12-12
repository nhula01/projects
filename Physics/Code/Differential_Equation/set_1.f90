double precision function f(y,x)
    implicit none
    double precision :: y,x
    f = exp(-y)
end function f

subroutine euler(file_name, dx) ! done
    implicit none
    double precision, parameter :: x_min=0.0d0,x_max=1.0d0 !interval of integration of ODE
    double precision, parameter :: y0=0 !initial condition
    double precision :: x,y,f
    character(len=100), intent(in) :: file_name
    double precision, intent(in) :: dx
    open(file=file_name, unit=10)
        ! initialize the x-value and initial condistion
        x=x_min
        y=y0
        ! while not over 1 yet 
        do while(x<=x_max)
            y=y+dx*f(y,x) !forward Euler method
            x=x+dx !update the x value
            write(10,*) x,y ! recorrd the value to the file
        enddo
    close(unit=10)
end subroutine euler

subroutine runge_kutta_2nd(file_name, dx) ! done
    implicit none
    double precision, parameter :: x_min=0.0d0,x_max=1.0d0 !interval of integration of ODE
    double precision, parameter :: y0=0 !initial condition
    double precision :: x,y,f,k1,k2,tmp1,tmp2
    character(len=100), intent(in) :: file_name
    double precision, intent(in) :: dx
    open(file=file_name, unit=10)
        ! the initial step and initial condition
        x=x_min
        y=y0
        ! while less than 1
        ! operation: write x,y on a file
        do while(x<=x_max)
            k1=f(y,x) !calculate k1
            !update y and x for k2
            k2=f(y+.5*k1*dx,x+.5*dx) !calculate k2
            ! update values
            y=y+k2*dx 
            x=x+dx
            write(10,*) x,y
        enddo
    close(unit=10)
end subroutine runge_kutta_2nd

subroutine runge_kutta_4th(file_name, dx)
    implicit none
    double precision, parameter :: x_min=0.0d0,x_max=1.0d0 !interval of integration of ODE    double precision, parameter :: y0=0 !initial condition
    double precision :: x,y,f,k1,k2,k3,k4,tmp1,tmp2, y0=0
    character(len=100), intent(in) :: file_name
    double precision, intent(in) :: dx
    open(file=file_name, unit=10)
        ! initial value
        x=x_min
        y=y0
        ! while less than 1
        ! operation : write x,y on a file
        do while(x<=x_max)
            ! find k1
            k1=f(y,x) 
            ! find k2
            k2=f(y+.5*k1*dx,x+.5*dx)
            ! find k3
            k3=f(y+.5*k2*dx,x+.5*dx)
            ! find k4
            k4=f(y+k3*dx,x+dx)
            ! update the values
            y=y+dx*(k1+2*k2+2*k3+k4)/6
            x=x+dx
            write(10,*) x,y
        enddo
    close(unit=10)
end subroutine runge_kutta_4th


program hw5_1
    implicit none 
    character(len=100) :: file_name = '/Users/phihung/Documents/PHY/first/homework/hw_05/' 
    character(len=100) :: name_a_1, name_b_1, name_c_1 
    character(len=100) :: name_a_2, name_b_2, name_c_2 
    character(len=100) :: name_a_3, name_b_3, name_c_3 
    double precision :: dx1 = 1d-2, dx2=5d-2, dx3=.1
    name_a_1 = "/Users/phihung/NumMethod/first/homework/hw_05/euler_1.dat" 
    name_b_1 = "/Users/phihung/NumMethod/first/homework/hw_05/runge_kutta_2nd_1.dat" 
    name_c_1 = "/Users/phihung/NumMethod/first/homework/hw_05/runge_kutta_4th_1.dat"
    name_a_2 = "/Users/phihung/NumMethod/first/homework/hw_05/euler_2.dat"
    name_b_2 = "/Users/phihung/NumMethod/first/homework/hw_05/runge_kutta_2nd_2.dat"
    name_c_2 = "/Users/phihung/NumMethod/first/homework/hw_05/runge_kutta_4th_2.dat"
    name_a_3 = "/Users/phihung/NumMethod/first/homework/hw_05/euler_3.dat"
    name_b_3 = "/Users/phihung/NumMethod/first/homework/hw_05/runge_kutta_2nd_3.dat"
    name_c_3 = "/Users/phihung/NumMethod/first/homework/hw_05/runge_kutta_4th_3.dat"
    
    ! generate first set
    call euler(name_a_1, dx1)
    call euler(name_a_2, dx2)
    call euler(name_a_3, dx3)

    ! generate second set
    call runge_kutta_2nd(name_b_1, dx1)
    call runge_kutta_2nd(name_b_2, dx2)
    call runge_kutta_2nd(name_b_3, dx3)

    !third_set
    call runge_kutta_4th(name_c_1, dx1)
    call runge_kutta_4th(name_c_2, dx2)
    call runge_kutta_4th(name_c_3, dx3)

end program hw5_1