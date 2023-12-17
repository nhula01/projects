implicit none
double precision :: t_final   !total time of propagation
double precision, parameter :: dt=.05           !time step in seconds
integer :: Nt!=int(t_final/dt)            !total number of time steps
double precision, parameter :: g=9.8                !free-fall accelleration in m/s**2
double precision, parameter :: x0=0.0,y0=0.0        !initial position
double precision :: v0=35.0,angle=68.0,theta   !initial speed [in m/s] and angle [in degrees]
double precision t,x,y,v0x,v0y,vy,pi
double precision x_ex,y_ex,vy_ex,t_ex
integer n

pi=4.0*atan(1.0)  !calculate pi
theta = angle * pi / 180.0 !convert theta to radian

v0x=v0*cos(theta) !initial velocity along x --- doesn't change since we neglect air resistance
v0y=v0*sin(theta) !initial velocity along y

! find the time of flight
! reformulate will give us t = 2v_0sin(theta) /g analytically
t_final = (2 * v0y) / g

Nt=int(t_final/dt)  !total steps

!open files
open(file='/Users/phihung/Documents/PHY/first/homework/hw_01/x_t.dat',unit=10)
open(file='/Users/phihung/Documents/PHY/first/homework/hw_01/y_t.dat',unit=11)
open(file='/Users/phihung/Documents/PHY/first/homework/hw_01/vy_t.dat',unit=12)

n=0.0
! set initial conditions for exact calculation as well
x=x0
y=x0
vy=v0y

do while (x>=0.0 .and. y>=0.0) !propagate until the projectile hits the ground

  !numerical
  x=x+v0x*dt !update the x position
  vy=vy-g*dt !update velocity in y dir
  y=y+vy*dt !update the y position

  !exact
  n=n+1 !update the iteration
  t=dt*float(n) !find time passed
  x_ex=x0+v0x*t !find the exact x location
  vy_ex=v0y-g*t !update the exact y speed
  y_ex=y0+v0y*t-0.5*g*(t**2) !update the exact y location

  !write the files
  write(10,*) t,x,x_ex
  write(11,*) t,y,y_ex
  write(12,*) t,vy,vy_ex

enddo

close(unit=10)
close(unit=11)
close(unit=12)

write(*,*) 'exact time of flight', t_final
write(*,*) 'calculated time of flight',t
end
