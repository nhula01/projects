implicit none
integer, parameter :: n_size=3
double precision current(n_size),A(n_size,n_size),B(n_size),tmp
integer indx(n_size)

!define A and B
A(1,1)= 1.0 !first row
A(1,2)=-1.0
A(1,3)= 1.0

A(2,1)=-3.0 ! second row
A(2,2)=-6.0
A(2,3)= 0.0

A(3,1)= 0.0 ! third row
A(3,2)=-6.0
A(3,3)=-4.0
! Ax = B 
B(1)= 0.0 
B(2)= -5.0
B(3)= -7.0

call ludcmp(A,n_size,n_size,indx,tmp) !this call performs LU decomposition of A
                                      !note that A is replaced by LU decomposition matrix
call lubksb(A,n_size,n_size,indx,B)   !B is replaced by the solution of A*X=B

write(*,*) 'current I1=',B(1)
write(*,*) 'current I2=',B(2)
write(*,*) 'current I3=',B(3)

end


subroutine lubksb(a,n,np,indx,b)
  implicit none
  integer n,np,indx(n)
  double precision a(np,np),b(n)
  integer i,ii,j,ll
  double precision sum_tmp
  ii=0
  do 12 i=1,n
    ll=indx(i)
    sum_tmp=b(ll)
    b(ll)=b(i)
    if (ii.ne.0)then
      do 11 j=ii,i-1
        sum_tmp=sum_tmp-a(i,j)*b(j)
11        continue
    else if (sum_tmp.ne.0.) then
      ii=i
    endif
    b(i)=sum_tmp
12    continue
  do 14 i=n,1,-1
    sum_tmp=b(i)
    do 13 j=i+1,n
      sum_tmp=sum_tmp-a(i,j)*b(j)
13      continue
    b(i)=sum_tmp/a(i,i)
14    continue
end subroutine lubksb


subroutine ludcmp(a,n,np,indx,d)
  implicit none
  integer n,np,indx(n) !indx is an output vector that records the row permutation effected by the
                       !partial pivoting, see Numerical Recipes in FORTRAN77, chapter 3.2 for more details
  integer, parameter :: NMAX=500 !maximum size of the array a
  double precision d,a(np,np)
  double precision, parameter :: TINY=1.0d-20 !small number
  integer i,imax,j,k
  double precision aamax,dum,sum_tmp,vv(NMAX)
  d=1.0
  do 12 i=1,n
    aamax=0.0
    do 11 j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
    if (aamax.eq.0.) stop 'singular matrix in ludcmp, terminated'
    vv(i)=1./aamax
12    continue
  do 19 j=1,n
    do 14 i=1,j-1
      sum_tmp=a(i,j)
      do 13 k=1,i-1
        sum_tmp=sum_tmp-a(i,k)*a(k,j)
13        continue
      a(i,j)=sum_tmp
14      continue
    aamax=0.
    do 16 i=j,n
      sum_tmp=a(i,j)
      do 15 k=1,j-1
        sum_tmp=sum_tmp-a(i,k)*a(k,j)
15        continue
      a(i,j)=sum_tmp
      dum=vv(i)*abs(sum_tmp)
      if (dum.ge.aamax) then
        imax=i
        aamax=dum
      endif
16      continue
    if (j.ne.imax)then
      do 17 k=1,n
        dum=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=dum
17        continue
      d=-d
      vv(imax)=vv(j)
    endif
    indx(j)=imax
    if(a(j,j).eq.0.)a(j,j)=TINY
    if(j.ne.n)then
      dum=1./a(j,j)
      do 18 i=j+1,n
        a(i,j)=a(i,j)*dum
18        continue
    endif
19    continue
end subroutine ludcmp
