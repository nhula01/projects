double precision function f(y,x)
    implicit none
    double precision :: y,x
    f = 1+2*x*y
end function

subroutine euler(file_name, dx)
    implicit none
    double precision, parameter :: x_min=0.0d0,x_max=1.0d0 !interval of integration of ODE
    double precision, parameter :: y0=0,i=0 !initial condition
    double precision :: x,y,f
    character(len=100), intent(in) :: file_name
    double precision, intent(in) :: dx
    character(len=100) :: str,file_name_2
    integer::j=0, k 
    j = j + 1     ! record multiple files
    k = 1 ! flag for finding the value at .45
    write(str, '(I5.5)') j ! convert to string
    file_name_2 = "/Users/phihung/NumMethod/first/homework/hw_05/"//str !files recording values at .45 (0001-4)
    open(file=file_name, unit=10)
    open(file=file_name_2, unit=11)
        ! initial value
        x=x_min
        y=y0
        ! while less than 1
        ! operations: record data (x,y) for a dx in one file
        ! record (dx,y(.45)) in one file
        do while(x<=x_max)
            !forward Euler method
            y=y+dx*f(y,x)
            x=x+dx
            write(10,*) x,y
            !record value at .45 for that in (0001-4)
            if(x>.45 .and. k .eq. 1)then
                write(11,*) x,y
                k=0 !set flag to false
            endif
        enddo
    close(unit=11)
    close(unit=10)
end subroutine euler

subroutine runge_kutta_4th(file_name, dx)
    implicit none
    double precision, parameter :: x_min=0.0d0,x_max=1.0d0 !interval of integration of ODE
    double precision, parameter :: y0=0 !initial condition
    double precision :: x,y,f,k1,k2,k3,k4,tmp1,tmp2
    character(len=100), intent(in) :: file_name
    double precision, intent(in) :: dx
    integer :: j=4, k 
    character(len=100) :: file_name_2, str
    j = j + 1
    k = 1
    write(str, '(I5.5)') j ! convert integer to string
    file_name_2 = "/Users/phihung/NumMethod/first/homework/hw_05/"//str ! files to record values at .45 (0005-8)
    open(file=file_name, unit=10)
    open(file=file_name_2, unit=11)
        x=x_min
        y=y0
        ! operations:
        ! record (x,y) for different delta steps
        ! record (dx,y(.45)) 
        do while(x<=x_max)
            k1=f(y,x) 
            k2=f(y+.5*k1*dx,x+.5*dx)
            k3=f(y+.5*k2*dx,x+.5*dx)
            k4=f(y+k3*dx,x+dx)
            y=y+dx*(k1+2*k2+2*k3+k4)/6
            x=x+dx
            write(10,*) x,y
            ! record value at .45
            if(x>.45 .and. k .eq. 1)then
                write(11,*)x,y
                k=0
            endif
        enddo
    close(unit=11)
    close(unit=10)
end subroutine runge_kutta_4th

program hw5_2
    implicit none
    character(len=100) :: file_name
    double precision :: dx1=1d-2,dx2=1d-1,dx3=.15, dx4=.3 !dx1=1d-2,dx2=1d-3,dx3=.15, dx4=1d-4
    file_name = "/Users/phihung/NumMethod/first/homework/hw_05/p2_dx1.dat"
    call runge_kutta_4th(file_name,dx1)
    file_name = "/Users/phihung/NumMethod/first/homework/hw_05/p2_euler_dx1.dat"
    call euler(file_name, dx1)
    file_name = "/Users/phihung/NumMethod/first/homework/hw_05/p2_dx2.dat"
    call runge_kutta_4th(file_name,dx2)
    file_name = "/Users/phihung/NumMethod/first/homework/hw_05/p2_euler_dx2.dat"
    call euler(file_name, dx2)
    file_name = "/Users/phihung/NumMethod/first/homework/hw_05/p2_dx3.dat"
    call runge_kutta_4th(file_name,dx3)
    file_name = "/Users/phihung/NumMethod/first/homework/hw_05/p2_euler_dx3.dat"
    call euler(file_name, dx3)
    file_name = "/Users/phihung/NumMethod/first/homework/hw_05/p2_dx4.dat"
    call runge_kutta_4th(file_name,dx4)
    file_name = "/Users/phihung/NumMethod/first/homework/hw_05/p2_euler_dx4.dat"
    call euler(file_name, dx4)
    ! and 2b
end program hw5_2