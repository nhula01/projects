subroutine average
    implicit none
    integer :: N_days = 7
    integer, parameter :: n_data=504, n_akima =20000

    integer :: n, isTrue ,n1,nn, n_update
    double precision :: x_data(n_data), y_data(n_data), df_dx, d2f_dx2, sum, y_avg(n_akima), x_avg(n_akima)
    character(len=100) :: name_1 = 'df_dx_7', name_2='d2f_dx2_7', name_3='avg7'
    character(len=100) :: file_name, file_name_1, file_name_2
    ! open the file
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw.dat',unit=10) !make sure data file is in the same directory as the source and executable
        !read data
        do n=1,n_data !but not all
            read(10,*) x_data(n), y_data(n)
        enddo

        !average the value
        n1=0
        open(file='/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_average_7.dat', unit=11)
            do n=1,n_data-N_days !for no extrapolation
                ! take the sum
                sum = 0 
                do nn=0, N_days-1
                    sum = sum + y_data(n+nn)
                enddo
                !write the value to the set
                sum = sum/N_days
                y_avg(n) = sum 
                write(11,*) x_data(n), y_avg(n)
            enddo
            
            !calculate the last avg
            do n=n_data-N_days+1, n_data
                sum = sum + y_data(n)
            enddo

            do n=n_data-N_days+1, n_data
                write(11,*) x_data(n), sum/N_days
            enddo

        close(unit=11)

        !open our average file and read the value
        open(file='/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_average_7.dat', unit=11)
            do n=1,n_data !but not all
                read(11,*) x_avg(n), y_avg(n)
            enddo
        close(unit=11)

        call akima_generation(n_data, x_avg, y_avg, name_3) ! generate the new akima 
        file_name = '/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_akima' // name_3
        !find first derivative
        open(file=file_name, unit=11)
            do n=1,n_akima !read the new akima value
                read(11,*) x_avg(n), y_avg(n)
            enddo
        close(unit=11)

        file_name_1 = '/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_akima' // name_1
        open(file=file_name_1,unit=11)
            do n=3,n_akima-3
                df_dx=(y_avg(n-2)-8*y_avg(n-1)+8*y_avg(n+1)-y_avg(n+2))/(12) !12* the time steps
                write(11,*) x_avg(n), df_dx  !record the value x and df_dx
            enddo
        close(unit=11)

        file_name_2 = '/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_akima' // name_2
        open(file=file_name_2,unit=11)
            do n=3,n_akima-3
                d2f_dx2=(-y_avg(n-2)+16.0*y_avg(n+1)-30.0*y_avg(n)+16.0*y_avg(n-1)-y_avg(n+2))/(12.0)
                write(11,*) x_avg(n), d2f_dx2  !record the value x and d2f_dx2 !time steps squared
            enddo
        close(unit=11)
    close(unit=10)
end subroutine average

subroutine akima_generation(n_data, x_data, y_data, name)
    implicit none
    integer, INTENT(IN) :: n_data !number of data points
    double precision, intent(in) :: x_data(n_data), y_data(n_data)
    integer, parameter :: n_akima=20000 !number of interpolation points
    double precision :: x_akima(n_akima),y_akima(n_akima), df_dx, ddf_ddx
    double precision :: x_new, h
    integer n,nn
    character(len=100), intent(in) :: name
    character(len=100) :: file_name='/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_akima', end_line = '.dat'

    ! generate points for x_akima
    do n=1,n_akima
        x_akima(n)=x_data(1)+(x_data(n_data)-x_data(1))*float(n-1)/float(n_akima-1) !points at which to calculate interpolation
    enddo

    call akima(3,n_data,x_data,y_data,n_akima,x_akima,y_akima) !cubic spline
    ! generate more data points
    file_name = '/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw_akima' // name 
    open(file=file_name,unit=11)
        do n=1,n_akima
            write(11,*) x_akima(n),y_akima(n) !writing generated data into a new file
        enddo
    close(unit=11)
end subroutine akima_generation


program hw3_3
    implicit none
    integer, parameter :: n_data=504 !number of data points
    double precision :: x_data(n_data), y_data(n_data)
    integer n
    character(len=100) :: name='main'
    !open the file
    !generate more data with akima
    open(file='/Users/phihung/Documents/PHY/first/homework/hw_03/CovidRaw.dat',unit=10) !make sure data file is in the same directory as the source and executable
        do n=1,n_data
            read(10,*) x_data(n),y_data(n) !reading data file
        enddo
    close(unit=10)
    call akima_generation(n_data, x_data, y_data, name)

    !calculate the average
    call average

end program hw3_3


!making extra akima points 
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
      INTEGER, PARAMETER :: N_data_max=504,N_spline_max=20000
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
