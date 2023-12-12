implicit  none
double precision :: pi
integer :: Nt=1000,N=1000 ! Nt:number of x values, N: iterations for sigma
double precision :: dx ! x_step
integer :: n_1, n_2
double precision :: x, y ! functions
double precision :: function, x_current


pi=4.D0*DATAN(1.D0)! create pi

dx = 20*pi / float(Nt) !find dx between 0 to 20*pi

!default value
function = 0.0
x_current = 0.0

open(file='/Users/phihung/Documents/PHY/first/homework/hw_01/convergence_1000.dat',unit=10)
do n_1 = 1, Nt ! loop through to calculate all points
    x = x_current ! assign the x value
    function = 0.0 ! reset value for y

    do n_2 = 1, N !use looping thr sigma to find a single point f(x)
        function = function + sin(x/n_2)/(n_2**2) !find the f(x) value
    end do 

    x_current = x_current + dx ! update the x-value 
    y = function ! paste to y
    write(10,*) x , y ! record data
end do
close(unit=10)
end
