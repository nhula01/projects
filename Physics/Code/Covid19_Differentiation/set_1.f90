!define the function at a=1
subroutine function(x,f) !done
  implicit none
  double precision, INTENT(IN) :: x
  double precision, INTENT(OUT) :: f
  f=cos(x)*sinh(x)
end subroutine function

!define the exact derivative
subroutine exact_df_dx(h,N_max) !a
    implicit none
    double precision :: x_min=0, x_max !define the range of x
    double precision :: pi, tmp1, tmp2, df_dx, x_new, f, h
    integer :: N_max,n 
    pi = 4.D0*DATAN(1.D0) !define pi=3.14
    x_max = pi/2
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_exact.dat',unit=11) !file to record convergence
        do n=2,N_max
            x_new=x_min+(x_max-x_min)*float(n-1)/float(N_max-1) ! find new value of x
            df_dx=-sin(x_new)*sinh(x_new)+cos(x_new)*cosh(x_new)
            write(11,*) x_new, df_dx  !record the value x and df_dx
        enddo
    close(unit=11)
end subroutine exact_df_dx

subroutine forward(h,N_max) ! calculate forward done !a
  implicit none
  double precision :: h_pi 
  double precision :: x_min=0, x_max, x_new !define the range of x
  double precision :: tmp1, tmp2, df_dx, f, h
  integer :: n,N_max !N_max is the number of points for calculating derivative
  character(len=100) :: name_backward = '/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_forward.dat'
  h_pi = 4.D0*DATAN(1.D0)/2
  x_max = h_pi
  open(file=name_backward,unit=11) !file to record convergence
      do n=2,N_max ! 2 to prevent extrapolate
          x_new=x_min+(x_max-x_min)*float(n-1)/float(N_max-1) ! find new value of x
          call function(x_new+h, tmp1)          !function at x_new
          call function(x_new, tmp2)        !function at x_new-h
          df_dx=(tmp1-tmp2)/h                 !backward difference
          write(11,*) x_new, df_dx            !record the value x and df_dx
      enddo
  close(unit=11)
end subroutine forward

subroutine backward(h,N_max) ! calculate backward !a
    implicit none
    double precision :: h_pi 
    double precision :: x_min=0, x_max, x_new !define the range of x
    double precision :: tmp1, tmp2, df_dx, f, h
    integer :: n,N_max !N_max is the number of points for calculating derivative
    character(len=100) :: name_backward = '/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_backward.dat'
    h_pi = 4.D0*DATAN(1.D0)/2
    x_max = h_pi
    open(file=name_backward,unit=11) !file to record convergence
        do n=2,N_max ! 2 to prevent extrapolate
            x_new=x_min+(x_max-x_min)*float(n-1)/float(N_max-1) ! find new value of x
            call function(x_new, tmp1)          !function at x_new
            call function(x_new-h, tmp2)        !function at x_new-h
            df_dx=(tmp1-tmp2)/h                 !backward difference
            write(11,*) x_new, df_dx            !record the value x and df_dx
        enddo
    close(unit=11)
end subroutine backward

subroutine central(h,N_max) ! calculate central !a
    implicit none
    double precision :: h_pi 
    double precision :: x_min=0, x_max, x_new !define the range of x
    double precision :: tmp1, tmp2, df_dx, f, h
    integer :: n,N_max !N_max is the number of points for calculating derivative
    character(len=100) :: name_central = '/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_central.dat'
    h_pi = 4.D0*DATAN(1.D0)/2
    x_max = h_pi
    open(file=name_central,unit=11) !file to record convergence
        do n=2,N_max-1 !prevents extrapolation
            x_new=x_min+(x_max-x_min)*float(n-1)/float(N_max-1) ! find new value of x
            call function(x_new+h, tmp1)     !function at x_point+dx
            call function(x_new-h, tmp2)        !function at x_new
            df_dx=(tmp1-tmp2)/(2*h)      !forward difference
            write(11,*) x_new, df_dx  !record the value x and df_dx
        enddo
    close(unit=11)
end subroutine central

subroutine thirdorder(h,N_max) ! calculate third order !a
    implicit none
    double precision :: h_pi 
    double precision :: x_min=0, x_max, x_new !define the range of x
    double precision :: tmp1, tmp2, df_dx, f, h
    integer :: n,N_max !N_max is the number of points for calculating derivative
    character(len=100) :: name_thirdorder = '/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_thirdorder.dat'
    h_pi = 4.D0*DATAN(1.D0)/2
    x_max = h_pi
    open(file=name_thirdorder,unit=11) !file to record convergence
        do n=2,N_max
            x_new=x_min+(x_max-x_min)*float(n-1)/float(N_max-1) ! find new value of x
            call function(x_new+h, tmp1)     
            df_dx=6*tmp1
            call function(x_new-h, tmp2)       
            df_dx=df_dx-2*tmp2              
            call function(x_new+2*h, tmp2)
            df_dx=df_dx-tmp2 
            call function(x_new, tmp2)
            df_dx=(df_dx-3*tmp2)/(6*h)
            write(11,*) x_new, df_dx 
        enddo
    close(unit=11)
end subroutine thirdorder

subroutine fourthorder(h,N_max) ! calculate fourthorder !a
    implicit none
    double precision :: h_pi 
    double precision :: x_min=0, x_max, x_new !define the range of x
    double precision :: tmp1, tmp2, df_dx, f, h
    integer :: n,N_max !N_max is the number of points for calculating derivative
    character(len=100) :: name_fourthorder = '/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_fourthorder.dat'
    h_pi = 4.D0*DATAN(1.D0)/2
    x_max = h_pi
    open(file=name_fourthorder,unit=11) !file to record convergence
        do n=2,N_max
            x_new=x_min+(x_max-x_min)*float(n-1)/float(N_max-1) ! find new value of x
            call function(x_new-2*h, tmp1)     !function at x_point+dx
            df_dx=tmp1
            call function(x_new-h, tmp2)        !function at x_new
            df_dx=df_dx-8*tmp2              !forward difference
            call function(x_new+h, tmp2)
            df_dx=df_dx+8*tmp2 
            call function(x_new+2*h, tmp2)
            df_dx=(df_dx-tmp2)/(12*h)
            write(11,*) x_new, df_dx  !record the value x and df_dx
        enddo
    close(unit=11)
end subroutine fourthorder

subroutine df_dx_and_d2f_dx2(x,df_dx,d2f_dx2,a)
  implicit none
  double precision, INTENT(IN) :: x,a ! take input and the value of a
  double precision, INTENT(OUT) :: df_dx, d2f_dx2 !output first derivative and second derivative
  df_dx=-sin(x)*sinh(a*x)+cos(x)*cosh(a*x)*a
  d2f_dx2 = (a**2-1)*cos(x)*sinh(a*x) -2*sin(x)*cosh(a*x)*a
end subroutine df_dx_and_d2f_dx2

subroutine root_nr(a,x_iteration_old) !b find a root for any a input
    implicit none
    double precision, parameter :: tol=1.0d-6 !maximum relative error
    double precision, parameter :: x_initial=1.0 !initial guess
    double precision :: x_iteration,relative_error
    double precision, INTENT(IN) :: a
    double precision, INTENT(OUT) :: x_iteration_old
    double precision :: df_dx, d2f_dx2, h_pi
    double precision :: x1=0,x2
    integer n,k

    !define the range
    h_pi = 4.D0*DATAN(1.D0)/2 !define pi=3.14
    x2 = h_pi

    call df_dx_and_d2f_dx2(x_initial,df_dx,d2f_dx2,a) !call to find 1st der (tmp1) and 2nd (dfdx1) given a
    relative_error=1.0 !set inital error to any value higher than tol
    x_iteration_old=x_initial-df_dx/d2f_dx2 !zeroth iteration

    do while (relative_error>tol) !iterations untill relative error is less than tol
        n=n+1
        call df_dx_and_d2f_dx2(x_iteration_old,df_dx,d2f_dx2,a)
        x_iteration=x_iteration_old-df_dx/d2f_dx2
        relative_error=abs((x_iteration-x_iteration_old)/x_iteration)
        x_iteration_old=x_iteration
    enddo
end subroutine root_nr

subroutine bisection !c
    implicit none
    integer, parameter :: N_max = 4990
    integer, parameter :: n_akima = 10000
    integer n, index_1, index_2, temp_index
    double precision :: x_data(N_max), y_data(N_max)
    double precision :: x_akima(n_akima), y_akima(n_akima)
    double precision :: x_left ,pi, dx, sum=0, x_right
    double precision :: x_left_n, x_right_n, x_iteration, x_iteration_old
    double precision :: tol = 1.0d-3, relative_error, tmp1, tmp2, x_middle_n
    character(len=100) :: name_bisection = '/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_central.dat'
    !read the data /Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_central_ordered.dat'

    open(file=name_bisection,unit=11) ! read the data generated from central subroutine
        do n=1,N_max !but not all
            read(11,*) x_data(n), y_data(n)
        enddo
    close(unit=11)

    !create n akima values for new interpolation
    pi = 4.D0*DATAN(1.D0)
    dx = (pi/2)/n_akima
    do n=1,n_akima
        sum = sum + dx
        x_akima(n)= sum
    enddo

    !calculate the interpolation 
    call akima(3,N_max,x_data,y_data,n_akima,x_akima,y_akima)

    open(file='/Users/phihung/Documents/PHY/first/homework/hw_03/df_dx_akima.dat',unit=11)
        do n=1,n_akima
            write(11,*) x_akima(n), y_akima(n)
        enddo
    close(unit=11)

    !doing bisection technique
    n=0 !zeroth iteration
    index_1 = 1 !first index
    index_2 = n_akima !last index
    x_left_n=x_akima(index_1)    !initial point on the left
    x_right_n=x_akima(index_2)  !initial point on the right

    temp_index=INT((index_1+index_2)/2.0) !finding the middle index since they are uniformly generated
    x_middle_n = x_akima(temp_index)! find the value at that point

    relative_error=1.0 !set inital error to any value higher than tol
    x_iteration_old=x_middle_n !zeroth iteration
    
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_03/bisection.dat',unit=11)
    do while (relative_error>tol) !iterations untill relative error is less than tol
        n=n+1 ! update the iteration values
        !Update the boundary
        tmp1 = y_akima(index_1)
        tmp2 = y_akima(temp_index)
        if(tmp1*tmp2<=0.0)then 
            index_2=temp_index !<--- the root is in left subinterval
        else
            index_1=temp_index  !<--- the root is in right subinterval
        endif
        !update data
        temp_index=INT((index_1+index_2)/2.0) !finding the middle index
        x_middle_n = x_akima(temp_index) !find the new middle value
        x_iteration=x_middle_n 
        relative_error=abs((x_iteration-x_iteration_old)/x_iteration) !find the differences btw two values
        x_iteration_old=x_iteration
        write(11,*) n, x_iteration
    enddo
    close(unit=11)
end subroutine bisection

subroutine rootnr_a !d
    implicit none
    double precision :: a_initial, a_final, dx, x_convergence
    integer :: N = 2000, n1 !generate 2000 a data points

    a_initial = .1
    a_final = 5
    dx = (a_final-a_initial)/N
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_03/afx.dat',unit=11)
    do n1=1, N
        call root_nr(a_initial,x_convergence)
        write(11,*) a_initial, x_convergence
        a_initial = a_initial + dx
    enddo
    close(unit=11)
end subroutine rootnr_a


program hw3_1
    implicit none
    double precision :: h=.5, iteration_out, a=1.0
    integer :: N_max = 5000
    call forward(h,N_max)
    call backward(h,N_max)
    call central(h,N_max)
    call thirdorder(h,N_max)
    call fourthorder(h,N_max)
    call exact_df_dx(h,N_max)
    call root_nr(a,iteration_out) !b
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_03/a_1.dat',unit=11)
            write(11,*) iteration_out
    close(unit=11)
    call bisection
    call rootnr_a
end program hw3_1


Subroutine akima(NP,ND,XD,YD,NI,XI,YI)
   ! Univariate Interpolation (Improved Akima Method)

   ! Hiroshi Akima
   ! U.S. Department of Commerce, NTIA/ITS
   ! Version of 89/07/04

   ! This subroutine performs univariate interpolation.  It is based
   ! on the improved A method developed by Hiroshi Akima, 'A method
   ! of univariate interpolation that has the accuracy of a third-
   ! degree polynomial,' ACM TOMS, vol. xx, pp. xxx-xxx, 19xx.  (The
   ! equation numbers referred to in the comments below are those in
   ! the paper.)

   ! In this method, the interpolating function is a piecewise
   ! function composed of a set of polynomials applicable to
   ! successive intervals of the given data points.  This method
   ! uses third-degree polynomials as the default, but the user has
   ! an option to use higher-degree polynomial to reduce undulations
   ! in resulting curves.

   ! This method has the accuracy of a third-degree polynomial if
   ! the degree of the polynomials for the interpolating function is
   ! set to three.

   ! The input arguments are
   !   NP = degree of the polynomials for the interpolating
   !        function,
   !   ND = number of input data points
   !        (must be equal to 2 or greater),
   !   XD = array of dimension ND, containing the abscissas of
   !        the input data points
   !        (must be in a monotonic increasing order),
   !   YD = array of dimension ND, containing the ordinates of
   !        the input data points,
   !   NI = number of points for which interpolation is desired
   !        (must be equal to 1 or greater),
   !   XI = array of dimension NI, containing the abscissas of
   !        the desired points.

   ! The output argument is
   !   YI = array of dimension NI, where the ordinates of the
   !        desired points are to be stored.

   ! If an integer value smaller than 3 is given to the NP argument,
   ! this subroutine assumes NP = 3.

   ! The XI array elements need not be monotonic, but this
   ! subroutine interpolates faster if the XI array elements are
   ! given in a monotonic order.

   ! If the XI array element is less than XD(1) or greater than
   ! XD(ND), this subroutine linearly interpolates the YI value.

   ! Specification statement

      Implicit None

!      DOUBLE PRECISION, Dimension(:) :: XD,YD,XI,YI
      INTEGER, PARAMETER :: N_data_max=4990,N_spline_max=10000
      DOUBLE PRECISION XD(N_data_max),YD(N_data_max),XI(N_spline_max),YI(N_spline_max)

      DOUBLE PRECISION :: RENPM1,RENNM2,XII,X0,X1,X2,X3
      DOUBLE PRECISION :: Y0,Y1,Y2,Y3,DLT,A1,SMPEF,SMWTF,SMPEI,SMWTI
      DOUBLE PRECISION :: PE,SX,SY,SXX,SXY,DNM,B0,B1,DY0,DY1,DY2,DY3
      DOUBLE PRECISION :: VOL,EPSLN,WT,YP,YP0,YP1,DX,DY,A0,A2,A3
      DOUBLE PRECISION :: T0,T1,AA0,AA1,XX,U,UC,V,A12,A13

      INTEGER :: ND,NI,ID,NP0,NP,NPM1,II,IINTPV,IINT,IDMN
      INTEGER :: IDMX,IDMD,IEPT,ID0,IPE,ID1,ID2,ID3

   ! Error check
   if(ND>N_data_max)then
     write(*,*) 'check number of input data points and N_data_max'
     stop
   endif

   if(NI>N_spline_max)then
     write(*,*) 'check number of output data points and N_spline_max'
     stop
   endif


   10 IF (ND.LE.1)   GO TO 90
      IF (NI.LE.0)   GO TO 91
      DO 11  ID=2,ND
        IF (XD(ID).LE.XD(ID-1))     GO TO 92
   11 CONTINUE
   ! Branches off special cases.
      IF (ND.LE.4)   GO TO 50
   ! General case  --  Five data points of more
   ! Calculates some local variables.
   20 NP0=MAX(3,NP)
      NPM1=NP0-1
      RENPM1=NPM1
      RENNM2=NP0*(NP0-2)
   ! Main calculation for the general case
   ! First (outermost) DO-loop with respect to the desired points
   30 DO 39  II=1,NI
        IF (II.EQ.1)      IINTPV=-1
        XII=XI(II)
   ! Locates the interval that includes the desired point by binary
   ! search.
        IF (XII.LE.XD(1))  THEN
          IINT=0
        ELSE IF (XII.LT.XD(ND))  THEN
          IDMN=1
          IDMX=ND
          IDMD=(IDMN+IDMX)/2
   31     IF (XII.GE.XD(IDMD))  THEN
            IDMN=IDMD
          ELSE
            IDMX=IDMD
          END IF
          IDMD=(IDMN+IDMX)/2
          IF (IDMD.GT.IDMN)    GO TO 31
          IINT=IDMD
        ELSE
          IINT=ND
        END IF
   ! End of locating the interval of interest
   ! Interpolation or extrapolation in one of the three subcases
        IF (IINT.LE.0)  THEN
   ! Subcase 1  --  Linear extrapolation when the abscissa of the
   !                desired point is equal to that of the first data
   !                point or less.
   ! Estimates the first derivative when the interval is not the
   ! same as the one for the previous desired point.  --
   ! cf. Equation (8)
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
            X0=XD(1)
            X1=XD(2)-X0
            X2=XD(3)-X0
            X3=XD(4)-X0
            Y0=YD(1)
            Y1=YD(2)-Y0
            Y2=YD(3)-Y0
            Y3=YD(4)-Y0
            DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
            A1=(((X2*X3)**2)*(X3-X2)*Y1 &
&              +((X3*X1)**2)*(X1-X3)*Y2 &
&              +((X1*X2)**2)*(X2-X1)*Y3)/DLT
          END IF
   ! Evaluates the YI value.
          YI(II)=Y0+A1*(XII-X0)
   ! End of Subcase 1
        ELSE IF (IINT.GE.ND)  THEN
   ! Subcase 2  --  Linear extrapolation when the abscissa of the
   !                desired point is equal to that of the last data
   !                point or greater.
   ! Estimates the first derivative when the interval is not the
   ! same as the one for the previous desired point.  --
   ! cf. Equation (8)
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
            X0=XD(ND)
            X1=XD(ND-1)-X0
            X2=XD(ND-2)-X0
            X3=XD(ND-3)-X0
            Y0=YD(ND)
            Y1=YD(ND-1)-Y0
            Y2=YD(ND-2)-Y0
            Y3=YD(ND-3)-Y0
            DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
            A1=(((X2*X3)**2)*(X3-X2)*Y1 &
&              +((X3*X1)**2)*(X1-X3)*Y2 &
&              +((X1*X2)**2)*(X2-X1)*Y3)/DLT
          END IF
   ! Evaluates the YI value.
          YI(II)=Y0+A1*(XII-X0)
   ! End of Subcase 2
        ELSE
   ! Subcase 3  --  Interpolation when the abscissa of the desired
   !                point is  between those of the first and last
   !                data points.
   ! Calculates the coefficients of the third-degree polynomial (for
   ! NP.LE.3) or the factors for the higher-degree polynomials (for
   ! NP.GT.3), when the interval is not the same as the one for the
   ! previous desired point.
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
   ! The second DO-loop with respect to the two endpoints of the
   ! interval
            DO 37  IEPT=1,2
   ! Calculates the estimate of the first derivative at an endpoint.
   ! Initial setting for calculation
              ID0=IINT+IEPT-1
              X0=XD(ID0)
              Y0=YD(ID0)
              SMPEF=0.0d0
              SMWTF=0.0d0
              SMPEI=0.0d0
              SMWTI=0.0d0
   ! The third (innermost) DO-loop with respect to the four primary
   ! estimate of the first derivative
              DO 36  IPE=1,4
   ! Selects point numbers of four consecutive data points for
   ! calculating the primary estimate of the first derivative.
                IF (IPE.EQ.1)  THEN
                  ID1=ID0-3
                  ID2=ID0-2
                  ID3=ID0-1
                ELSE IF (IPE.EQ.2)  THEN
                  ID1=ID0+1
                ELSE IF (IPE.EQ.3)  THEN
                  ID2=ID0+2
                ELSE
                  ID3=ID0+3
                END IF
   ! Checks if any point number falls outside the legitimate range
   ! (between 1 and ND).  Skips calculation of the primary estimate
   ! if any does.
                IF (ID1.LT.1.OR.ID2.LT.1.OR.ID3.LT.1.OR. &
&                   ID1.GT.ND.OR.ID2.GT.ND.OR.ID3.GT.ND) &
&                    GO TO 36
   ! Calculates the primary estimate of the first derivative  --
   ! cf. Equation (8)
                X1=XD(ID1)-X0
                X2=XD(ID2)-X0
                X3=XD(ID3)-X0
                Y1=YD(ID1)-Y0
                Y2=YD(ID2)-Y0
                Y3=YD(ID3)-Y0
                DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
                PE=(((X2*X3)**2)*(X3-X2)*Y1 &
&                  +((X3*X1)**2)*(X1-X3)*Y2 &
&                  +((X1*X2)**2)*(X2-X1)*Y3)/DLT
   ! Calculates the volatility factor, VOL, and distance factor,
   ! SXX, for the primary estimate.  --  cf. Equations (9) and (11)
                SX=X1+X2+X3
                SY=Y1+Y2+Y3
                SXX=X1*X1+X2*X2+X3*X3
                SXY=X1*Y1+X2*Y2+X3*Y3
                DNM=4.0d0*SXX-SX*SX
                B0=(SXX*SY-SX*SXY)/DNM
                B1=(4.0d0*SXY-SX*SY)/DNM
                DY0=-B0
                DY1=Y1-(B0+B1*X1)
                DY2=Y2-(B0+B1*X2)
                DY3=Y3-(B0+B1*X3)
                VOL=DY0*DY0+DY1*DY1+DY2*DY2+DY3*DY3
   ! Calculates the EPSLN value, which is used to decide whether or
   ! not the volatility factor, VOL, is essentially zero.
                EPSLN=(YD(ID0)**2+YD(ID1)**2 &
&                     +YD(ID2)**2+YD(ID3)**2)*1.0D-12
   ! Accumulates the weighted primary estimates.  --
   ! cf. Equations (13) and (14)
                IF (VOL.GT.EPSLN)  THEN
   ! - For finite weight.
                  WT=1.0d0/(VOL*SXX)
                  SMPEF=SMPEF+PE*WT
                  SMWTF=SMWTF+WT
                ELSE
   ! - For infinite weight.
                  SMPEI=SMPEI+PE
                  SMWTI=SMWTI+1.0d0
                END IF
   36         CONTINUE
   ! End of the third DO-loop
   ! Calculates the final estimate of the first derivative.  --
   ! cf. Equation (14)
              IF (SMWTI.LT.0.5d0)  THEN
   ! - When no infinite weights exist.
                YP=SMPEF/SMWTF
              ELSE
   ! - When infinite weights exist.
                YP=SMPEI/SMWTI
              END IF
              IF (IEPT.EQ.1)  THEN
                YP0=YP
              ELSE
                YP1=YP
              END IF
   ! End of the calculation of the estimate of the first derivative
   ! at an endpoint
   37       CONTINUE
   ! End of the second DO-loop
            IF (NP0.LE.3)  THEN
   ! Calculates the coefficients of the third-degree polynomial
   ! (when NP.LE.3).  --  cf. Equation (4)
              DX=XD(IINT+1)-XD(IINT)
              DY=YD(IINT+1)-YD(IINT)
              A0=YD(IINT)
              A1=YP0
              YP1=YP1-YP0
              YP0=YP0-DY/DX
              A2=-(3.0d0*YP0+YP1)/DX
              A3= (2.0d0*YP0+YP1)/(DX*DX)
            ELSE
   ! Calculates the factors for the higher-degree polynomials
   ! (when NP.GT.3).  --  cf. Equation (20)
              DX=XD(IINT+1)-XD(IINT)
              DY=YD(IINT+1)-YD(IINT)
              T0=YP0*DX-DY
              T1=YP1*DX-DY
              AA0= (T0+RENPM1*T1)/RENNM2
              AA1=-(RENPM1*T0+T1)/RENNM2
            END IF
          END IF
   ! End of the calculation of the coefficients of the third-degree
   ! polynomial (when NP.LE.3) or the factors for the higher-degree
   ! polynomials (when NP.GT.3), when the interval is not the same
   ! as the one for the previous desired point.
   ! Evaluates the YI value.
          IF (NP0.LE.3)  THEN
   ! - With a third-degree polynomial.  --  cf. Equation (3)
            XX=XII-XD(IINT)
            YI(II)=A0+XX*(A1+XX*(A2+XX*A3))
          ELSE
   ! - With a higher-degree polynomial.  --  cf. Equation (19)
            U=(XII-XD(IINT))/DX
            UC=1.0d0-U
            V=AA0*((U**NP0)-U)+AA1*((UC**NP0)-UC)
            YI(II)=YD(IINT)+DY*U+V
          END IF
   ! End of Subcase 3
        END IF
   39 CONTINUE
   ! End of the first DO-loop
   ! End of general case
      RETURN
   ! Special cases  --  Four data points or less
   ! Preliminary processing for the special cases
   50 X0=XD(1)
      Y0=YD(1)
      X1=XD(2)-X0
      Y1=YD(2)-Y0
      IF (ND.EQ.2)   GO TO 60
      X2=XD(3)-X0
      Y2=YD(3)-Y0
      IF (ND.EQ.3)   GO TO 70
      X3=XD(4)-X0
      Y3=YD(4)-Y0
      GO TO 80
   ! Special Case 1  --  Two data points
   ! (Linear interpolation and extrapolation)
   60 A1=Y1/X1
      DO 61  II=1,NI
        YI(II)=Y0+A1*(XI(II)-X0)
   61 CONTINUE
   ! End of Special Case 1
      RETURN
   ! Special Case 2  --  Three data points
   ! (Quadratic interpolation and linear extrapolation)
   70 DLT=X1*X2*(X2-X1)
      A1=(X2*X2*Y1-X1*X1*Y2)/DLT
      A2=(X1*Y2-X2*Y1)/DLT
      A12=2.0d0*A2*X2+A1
      DO 71  II=1,NI
        XX=XI(II)-X0
        IF (XX.LE.0.0d0)  THEN
          YI(II)=Y0+A1*XX
        ELSE IF (XX.LT.X2) THEN
          YI(II)=Y0+XX*(A1+XX*A2)
        ELSE
          YI(II)=Y0+Y2+A12*(XX-X2)
        END IF
   71 CONTINUE
   ! End of Special Case 2
      RETURN
   ! Special Case 3  --  Four data points
   ! (Cubic interpolation and linear extrapolation)
   80 DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
      A1=(((X2*X3)**2)*(X3-X2)*Y1 &
&        +((X3*X1)**2)*(X1-X3)*Y2 &
&        +((X1*X2)**2)*(X2-X1)*Y3)/DLT
      A2=(X2*X3*(X2*X2-X3*X3)*Y1 &
&        +X3*X1*(X3*X3-X1*X1)*Y2 &
&        +X1*X2*(X1*X1-X2*X2)*Y3)/DLT
      A3=(X2*X3*(X3-X2)*Y1 &
&        +X3*X1*(X1-X3)*Y2 &
&        +X1*X2*(X2-X1)*Y3)/DLT
      A13=(3.0d0*A3*X3+2.0d0*A2)*X3+A1
      DO 81  II=1,NI
        XX=XI(II)-X0
        IF (XX.LE.0.0d0)  THEN
          YI(II)=Y0+A1*XX
        ELSE IF (XX.LT.X3) THEN
          YI(II)=Y0+XX*(A1+XX*(A2+XX*A3))
        ELSE
          YI(II)=Y0+Y3+A13*(XX-X3)
        END IF
   81 CONTINUE
   ! End of Special Case 3
      RETURN
   ! Error exit
   90 WRITE (*,99090) ND
      GO TO 99
   91 WRITE (*,99091) NI
      GO TO 99
   92 WRITE (*,99092) ID,XD(ID-1),XD(ID)
   99 WRITE (*,99099)
      RETURN
   ! Format statements for error messages
   99090 FORMAT (1X/ ' ***   Insufficient data points.', &
&                7X,'ND =',I3)
   99091 FORMAT (1X/ ' ***   No desired points.', &
&                7X,'NI =',I3)
   99092 FORMAT (1X/ ' ***   Two data points identical or out of ', &
&                'sequence.'/ &
&                7X,'ID, XD(ID-1), XD(ID) =',I5,2F10.3)
   99099 FORMAT (' Error detected in the akima subroutine'/)

End Subroutine Akima
