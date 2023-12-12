implicit none
  integer, parameter ::  n_points= 1000
  real, parameter :: a = 2.0
  real:: t, x, y
  integer:: i
  real :: pi
  double precision :: bound1, bound2, step, current_position

  pi = 3.14159265
  open(unit = 1, file='/Users/phihung/Documents/PHY/first/homework/hw_01/foliumdatx.dat')
  open(unit=2, file='/Users/phihung/Documents/PHY/first/homework/hw_01/foliumdaty.dat')
  ! boundary
  bound1 = -2*pi/9
  bound2= 2*pi/.5
  ! range and step
  step = (bound2-bound1)/n_points
  do i = 1, n_points
     x = (3.0*a*current_position)/(1+current_position**3.0)
     y = (3.0*a*current_position**2.0)/(1+current_position**3.0)
     current_position = current_position + step ! update position
     write(1, *)x
     write(2,*)y
  enddo
  close(unit=1)
end program