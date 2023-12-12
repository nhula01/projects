implicit none
double precision, parameter :: pi=3.141592653589793
double complex, parameter :: Im=(0.0d0,1.0d0) !imaginary unit

integer, parameter :: NN=201 !total number of points
integer, parameter :: a=(NN-1)/2
double precision, parameter :: x0=-5.0 !lower limit, assuming upper limit is -x0 (symmetric interval)
double precision, parameter :: dx=-2.0*x0/(NN-1) !step dx
double precision psi,x,k,dpsi,psi_array(NN)
double complex tmp
integer n,m,l

do n=0,NN-1
  x=x0+dx*n
  psi_array(n)=psi(x)
enddo

open(file='derivative_psi_201points.dat',unit=11)
dpsi=0.0
do l=0,NN-1
  tmp=(0.0,0.0)
  do n=0,NN-1
    do m=0,NN-1
      tmp=tmp+psi_array(n+1)*(m-a)*exp(2.0*pi*Im*(m-a)*(l-n)/NN)
    enddo
  enddo
  dpsi=tmp*2.0*pi*Im/(dx*NN**2)
  write(11,*) x0+dx*(l+1),dpsi
enddo

end

double precision function psi(x)
  double precision x
  psi=exp(-(x**2)/2.0)*sin(4.0*x)
end function psi
