double precision function f1(x)
    implicit none
    double precision, intent(in) :: x
    f1=2*exp(-4*(x**2))/sqrt(3-2*x) !??? 
end function f1

double precision function f2(x)
    implicit none
    double precision, intent(in) :: x
    f2=1/(1+x**2)
end function f2

subroutine montecarlo(name, func)
    implicit none
    integer, intent(in) :: func ! for question a
    character(len=100), intent(in) :: name
    character(len=100) :: file_name
    integer, parameter :: n_large=10000  !max number of random numbers
    double precision, allocatable :: numbers_allocate(:)
    double precision, parameter :: a=0.0,b=1.0 !limits of integration
    double precision x,u,tmp, f1, f2
    integer n,m
    file_name = '/Users/phihung/Documents/PHY/first/homework/hw_04/' // name
    open(file=file_name,unit=32)
        do n=1,n_large
            allocate(numbers_allocate(n)) ! allocate random points
            numbers_allocate=0.0
            call random_number(numbers_allocate) !generates random numbers betwee 0 and 1
            tmp=0.0
            do m=1,n
                u=numbers_allocate(m)
                x=u*b+a
                if(func==0)then 
                    tmp=tmp+f1(x)
                else
                    tmp=tmp+f2(x)
                endif
            enddo
            tmp=tmp*b/float(n)
            write(32,*) m,tmp
            deallocate(numbers_allocate)
        enddo
    close(unit=32)
end subroutine montecarlo 
program hw4_3
    implicit none 
    character(len=100) :: name_a = 'exp.dat', name_b = '1_x.dat'
    double precision :: f1, f2
    call montecarlo(name_a, 0)
    call montecarlo(name_b, 1)
end program hw4_3