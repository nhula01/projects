subroutine trapezoidal(dx, z, name, K_value) 
! takes in the steps/value of z/nameoffile/ return K value
! the subroutine also records a file
    implicit none
    double precision, intent(in) :: dx, z !step !z
    double precision, parameter :: pi=3.141592653589793
    double precision :: x1 = 0, x2, tmp, f_av, f
    integer ::  n_points, n
    character(len=100), intent(in) :: name 
    double precision, intent(out) :: K_value !area
    character(len=100) :: dir_new
    n_points = int((pi/2)/dx) ! the number of points to take the integration
    dir_new =  '/Users/phihung/Documents/PHY/first/homework/hw_04/' // name 
    open(file=dir_new, unit=10)
      x1=0
      do n=1,n_points !loop through all the points
          x2=x1+dx
          f_av=(f(x1,z)+f(x2,z))/2.0
          tmp=tmp+f_av*dx
          write(10,*) x1,tmp ! record the convergence to the value
          x1=x1+dx !<--- next step
      enddo
      K_value = tmp 
    close(unit=10)
end subroutine

subroutine trapezoidal_auto(dx, z, K_value) 
! take in the steps/ value of z/ return the K_value
! automate K_value
    implicit none
    double precision, intent(in) :: dx, z !step !z
    double precision, parameter :: pi=3.141592653589793
    double precision :: x1 = 0, x2, tmp, f_av, f
    integer ::  n_points, n
    double precision, intent(out) :: K_value !area
    character(len=100) :: dir_new
    n_points = int((pi/2)/dx) ! the number of points to take the integration
    do n=1,n_points !loop through all the points
        x2=x1+dx
        f_av=(f(x1,z)+f(x2,z))/2.0
        tmp=tmp+f_av*dx
        x1=x1+dx !<--- next step
    enddo
    K_value = tmp 
end subroutine

double precision function f(x, z)
!define the function
  implicit none
  double precision, intent(in) :: x, z
  f=1/sqrt(1-(z**2)*(sin(x)**2))
end function f

program hw3_1
! record three files of convergence for a,b, c 
  implicit none
  integer :: m, n=500
  double precision :: K_value, dx=.01, T, dn
  double precision, parameter :: pi=3.141592653589793
  double precision :: za, zb, zc, z, theta, d_theta
  character(len=100) :: name_a = 'integral_record_set1a.dat', name_b = 'integral_record_set1b.dat'
  character(len=100) :: name_c = 'integral_record_set1c.dat' 
  !fundamental constants
  double precision, parameter :: g = 9.81, L = .25  
  ! write files with four data values for a,b,c and the one being compared with
  open(file='/Users/phihung/Documents/PHY/first/homework/hw_04/specific_comparison_1.dat', unit=11)
    za = sin(pi/6) !a !at steps dx
    call trapezoidal(dx,za ,name_a, K_value) !in dx,za and name !out files and K value
    T = 4 * (sqrt(L/g)) * K_value
    write(11,*) za, T
    ! find out the convergence of the file
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_04/a_convergence.dat', unit=10)
      do m=1, n
        dn=3*float(m)/float(n)
        call trapezoidal_auto(dn, za, K_value)
        write(10,*) dn, K_value
      enddo
    close(unit=10)

    zb = sin(pi/8) !b
    call trapezoidal(dx,zb ,name_b, K_value) !in dx, zb and name !out files and K value
    T = 4 * (sqrt(L/g)) * K_value
    write(11,*) zb, T
    ! find out the convergence of the file
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_04/b_convergence.dat', unit=10)
      do m=1, n
        dn=3*float(m)/float(n)
        call trapezoidal_auto(dn, zb, K_value)
        write(10,*) dn, K_value
      enddo
    close(unit=10)

    zc = sin(pi/80) !c
    call trapezoidal(dx,zc ,name_c, K_value) !in dx, zc and name !out files and K value
    T = 4 * (sqrt(L/g)) * K_value
    write(11,*) zc, T
    !find out the convergence of the file
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_04/c_convergence.dat', unit=10)
      do m=1, n
        dn=3*float(m)/float(n)
        call trapezoidal_auto(dn, zc, K_value)
        write(10,*) dn, K_value
      enddo
    close(unit=10)
    
   ! T = 4 * (sqrt(L/g)) 
    !write(11,*) z, T
  close(unit=11)

  !investigate how the periods differ through the whole 0->pi
  d_theta = pi/10000
  theta = 0
  open(file='/Users/phihung/Documents/PHY/first/homework/hw_04/theta_vs_T.dat', unit=10)
    do m=1, 10000
      z = sin(theta)
      call trapezoidal_auto(dx,z, K_value) !in dx, zc and name !out files and K value
      T = 4 * (sqrt(L/g)) * K_value
      write(10,*) theta, T
      theta = theta + d_theta
    enddo
  close(unit=10)

end program hw3_1