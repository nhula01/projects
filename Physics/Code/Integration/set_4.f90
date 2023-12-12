double precision function f(x,y)
    implicit none
    double precision, intent(in) :: x, y
    f=(sin((x**3)*(y**2)))**2
end function f

subroutine montecarlo_through_y(x_large, y_large,tmp)
    implicit none
    !character(len=100), intent(in) :: name
    character(len=100) :: file_name
    integer :: y_large, x_large !max number of random numbers
    double precision, parameter :: a_y=-.5,b_y=1.0 !limits of integration in y
    double precision :: tmp1, tmp2,f_av, y1, y2,dy, tmp
    integer n,m
    !file_name = '/Users/phihung/Documents/PHY/first/homework/hw_04/' // name
    !open(file=file_name,unit=32)
        dy = (b_y-a_y)/y_large
        y1=a_y
        tmp=0.0
        do n=1,y_large !loop through all the points
            y2=y1+dy
            call monetecarlo_at_one_y(x_large, y1,tmp1)
            call monetecarlo_at_one_y(x_large, y2,tmp2)
            f_av=(tmp2+tmp1)/2.0
            tmp=tmp+f_av*dy
            !write(32,*) y1, tmp
            y1=y1+dy !<--- next step
        enddo
    !close(unit=32)
end subroutine montecarlo_through_y

subroutine monetecarlo_at_one_y(x_large, y, tmp) !find the integration at a constant y ! treat this as f(y) 
    implicit none
    integer :: x_large
    double precision, allocatable :: numbers_allocate(:)
    double precision, parameter :: a_x=-1.0,b_x=1.0 !limits of integration in x
    double precision :: x,u,tmp,f, y
    integer :: n
    allocate(numbers_allocate(x_large)) ! create a dict to store all random x values at certain y
    numbers_allocate=0.0
    call random_number(numbers_allocate) !generates random numbers betwee 0 and 1
    !sum in x direction
    tmp=0.0
    do n=1,x_large !find the value of f(y)
        u=numbers_allocate(n)
        x=u*b_x+a_x
        tmp=tmp+f(x,y)
    enddo
    tmp=tmp*(b_x-a_x)/float(x_large) !return the value f(y)
    deallocate(numbers_allocate)
end subroutine monetecarlo_at_one_y 


program hw4_4
    implicit none 
    double precision :: tmp
    integer :: n, x_large, y_large
    character(len=100) :: name = "double_integration_3d.dat", file_name
    file_name = '/Users/phihung/Documents/PHY/first/homework/hw_04/' // name
    open(file=file_name, unit=10) 
        do n=1, 10000
            call montecarlo_through_y(n, n, tmp)
            write(10,*) n, n, tmp
        enddo
    close(unit=10)
end program hw4_4