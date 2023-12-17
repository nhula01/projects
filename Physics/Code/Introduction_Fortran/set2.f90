! generate two columns for folium descartes
implicit none !so it doesnt assign random variables to be integers
! function to be generated is x^3 + y^3 -3axy = 0
! converting function to polar coordinate 
! r = 3acos(theta)sin(theta)/[cos^3(theta)+sin^3(theta)] 
double precision :: step, start 
integer :: n !for looping
integer, parameter :: a1= 2, a2=3, a3=4 !a-values
real(kind=8) :: pi 
integer, parameter :: n_max = 1000 ! number of theta points to be generated
double precision :: theta(n_max), r(n_max) ! array for storing data points

! set default values (theta, r)
theta = 0.0 !angle 
r = 0.0 !radius

! find n_max points from -pi/6 to pi/1.5
pi=4.D0*DATAN(1.D0) !calculate pi values
start = -pi/6 !start from -pi/6
step = ((pi/1.5)+(pi/6)) / n_max !the value of one step

! assign the values (theta, r)
do n=1,n_max ! loop through all coordinates
    theta(n) = start ! take theta value
    start = start + step ! update theta value
    r(n) = (3*cos(theta(n))*sin(theta(n))) / ((cos(theta(n)))**3 + (sin(theta(n)))**3)  ! take radius value
enddo

!record the data for polar values (theta, r*2,r*3,r*4) 
open(file='/Users/phihung/Documents/PHY/first/homework/hw_01/polar_coordinates.dat', unit=10)
do n=1, n_max
    write(10,*) theta(n), r(n)*2, r(n)*3, r(n)*4
enddo
close(unit=10)

!record the data for cartesian values (x1,y1)
open(file='/Users/phihung/Documents/PHY/first/homework/hw_01/cartesian_coordinates_1.dat', unit=11)
do n=1, n_max
    write(11,*) r(n)*2*cos(theta(n)), r(n)*2*sin(theta(n))! x = rcos(theta), y=rsin(theta) a=2
enddo
close(unit=11)
!record the data for cartesian values (x2,y2)
open(file='/Users/phihung/Documents/PHY/first/homework/hw_01/cartesian_coordinates_2.dat', unit=11)
do n=1, n_max
    write(11,*) r(n)*3*cos(theta(n)), r(n)*3*sin(theta(n))! x = rcos(theta), y=rsin(theta) a=3
enddo
close(unit=11)
!record the data for cartesian values (x3,y3)
open(file='/Users/phihung/Documents/PHY/first/homework/hw_01/cartesian_coordinates_3.dat', unit=11)
do n=1, n_max
    write(11,*) r(n)*4*cos(theta(n)), r(n)*4*sin(theta(n))! x = rcos(theta), y=rsin(theta) a=4
enddo
close(unit=11)
end