implicit none

integer, parameter :: n_max = 500  ! N values
integer :: n_2, n_1
double precision :: value ! the pi values recorded

! open the file to record
open(file='/Users/phihung/Documents/PHY/first/homework/hw_01/gregory_convergence.dat', unit=11)
do n_1=1, n_max !calculate convergence from N = 1 to 500
    value = 0.0 ! reset value each loop
    do n_2=1, n_1 !calculate the series for each N value
        value = value + (((-1)**(float(n_2)+1)) / (2*float(n_2) - 1)) !gregory series
    end do 
    value = value * 4.0 ! times 4 from outside
    write(11,*) n_1, value ! write the N and the pi value recorded
end do
close(unit=11)

end