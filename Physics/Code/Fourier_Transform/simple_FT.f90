implicit none
double precision, parameter :: pi=3.141592653589793
double complex, parameter :: Im=(0.0d0,1.0d0) !imaginary unit

integer, parameter :: NN=201 !total number of points
integer, parameter :: a=(NN-1)/2
double precision, parameter :: x0=-5.0 !lower limit, assuming upper limit is -x0 (symmetric interval)
double precision, parameter :: dx=-2.0*x0/(NN-1) !step dx
double precision, parameter :: k0=-pi/dx !lower limit for k0 due to Nyquist condition
double precision, parameter :: dk=2.0*pi/(NN*dx)
double precision psi,x,k
double complex phi !Fourier transform is complex
double complex tmp
integer n,m

open(file='ReIm_phi_201points.dat',unit=11)
open(file='abs_phi_201points.dat',unit=12)

do m=0,NN-1 !k
  k=k0+dk*m
  tmp=(0.0,0.0)
  do n=0,NN-1
    x=x0+dx*n
    tmp=tmp+psi(x)*exp(-2.0*pi*Im*(m-a)*(n-a)/NN)
  enddo

  phi=tmp*dx/sqrt(2.0*pi)

  write(11,*) k,dreal(phi),aimag(phi) !real and imaginary parts
  write(12,*) k,abs(phi)**2       !|phi|^2
enddo

close(unit=11)
close(unit=12)
end

double precision function psi(x)
  double precision x
  psi=exp(-(x**2)/2.0)*sin(4.0*x)
end function psi
