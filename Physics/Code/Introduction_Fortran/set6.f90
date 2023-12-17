implicit none

double precision, parameter :: dt = .05 !timestep
double precision :: x = 0, y = 0 !initialposition
double precision :: time_flight, g = 9.8, angle = 68
double precision :: v_x, v_y, v_o = 35 !speed
double precision :: t, theta
! find the time of flight
! reformulate will give us t = 2v_0sin(theta) /g analytically
! convert theta to radian
theta = angle * 3.14 / 180

v_x = v_o * cos(theta)
v_y = v_o * sin(theta)

time_flight = (2 * v_y) / g

write(*,*) time_flight

open(file='/Users/phihung/Documents/PHY/first/homework/hw_01/projectile.dat', unit=10)
t = 0 
    do while (y >= 0)
        
        write(10,*) t, x, y !record data
        t = t + dt !update time
        v_y = v_y - g*t !update y speed?
        !x speed stays the same
        dx = x + v_x * dt ! change in x
        dy = y + v_y * dt - g * t ! change in y
        write(*,*) y
    enddo

close(unit=10)

end