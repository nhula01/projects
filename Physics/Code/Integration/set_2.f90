double precision function f(x)
!definee function x
  implicit none
  double precision, intent(in) :: x
  f=(x**2)*exp(-x)
end function f

subroutine simpson(b, N_max, tmp)
! take in the number of iterations
! rreturn the value of the integration
    implicit none
    double precision :: a=0, dx, b , x=0, f
    integer :: m 
    integer, intent(in) :: N_max ! number of iterations
    double precision, intent(out) :: tmp
    tmp = f(a) + f(b)
    dx = (b-a)/N_max
    !N_max = int((b-a)*dx)
    !write(*,*) tmp
    do m=2, N_max-2
        if(mod(m,2)==0)then
            tmp = tmp + 2*f(a+m*dx)
        else
            tmp = tmp + 4*f(a+m*dx)
        endif
    enddo
    tmp = tmp*dx/3.0
end subroutine simpson

program hw4_2
    implicit none 
    double precision :: a, b, tmp, k=1000000
    integer :: N_max, n
    N_max = 5000000 ! number of time steps
    !4900002
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_04/N_points_vs_value.dat', unit=11)
        do n=2, N_max, 10000
            call simpson(k,n,tmp)
            write(11,*) n, tmp
        enddo 
    close(unit=11)

    open(file='/Users/phihung/Documents/PHY/first/homework/hw_04/upper_limit.dat', unit=11)
        do n=10, 1500000, 10000
            k = float(n)
            call simpson(k, N_max, tmp)
            write(11,*) n, tmp
        enddo 
    close(unit=11)

end program hw4_2