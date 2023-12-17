implicit none
integer, parameter :: n_max_terms=5000 !maximum number of terms in the series
integer n_terms
double precision x_degrees,x_radians,pi
double precision func
integer n1,n2

!calculate pi
pi=4.0*atan(1.0)

open(file='/Users/phihung/Documents/PHY/first/homework/hw_01/series_convergence.dat',unit=10)

  x_degrees=10.0 !at which x we are investigating convergence
  x_radians=x_degrees*pi/180.0 !convert to radians

do n1=1,n_max_terms
  func=0.0
  do n2=1,n1
    func=func+cos(x_radians/float(n2)**2)/float(n2)
  enddo
  write(10,*) n1,func
enddo

close(unit=10)

end
