! output two files of position (t,x,y) and (t,vx,vy) and final time with proper flag
subroutine poncho_stat(dt, v0, angle, file_name, flag, t1)
    implicit none
        double precision, parameter :: pi=3.141592653589793
        double precision, parameter :: g=9.81 !free fall acceleration
        double precision, parameter :: m=.45,D=.02 !mass, air drag coefficient
        double precision, parameter :: x0=0.0,y0=100.0 !initial condition x(0), y(0)
        double precision, intent(in):: angle 
        character(len=100), intent(in) :: file_name
        character(len=100) :: file_name_2, file_name_3 ! store position, and velocity
        double precision, intent(in) :: v0 !initial speed
        double precision, intent(in) :: dt
        double precision, intent(out) :: t1
        integer, intent(in) :: flag
        double precision v0x,v0y !initial velocity projections
        double precision, parameter :: t_min=0.0d0,t_max=100.0d0 !interval of integration of ODE
        integer, parameter :: N_steps=1000 !number of steps total time of flight
        double precision t_ideal, d_vx, d_vy,t
        double precision x,vx,y,vy, alpha
        double precision tmp1(4)
        double precision k1(4),k2(4),k3(4),k4(4),tmp(4)
    !find initial velocity
        t=0
        alpha=angle*pi/180.0 !launch angle in rad ->55 degrees
        v0x=v0*cos(alpha)
        v0y=v0*sin(alpha)
        file_name_2 = '/Users/phihung/NumMethod/first/homework/hw_05/xy_RK4_' // file_name
        file_name_3 = '/Users/phihung/NumMethod/first/homework/hw_05/vxvy_RK4_' // file_name
    !<--- RK4 method
    !x_dot = v_x
    !v_x = -(D/m)*(v_x)*v
    !y_dot = v_y
    !v_y = -g-(D/m)*(v_y)*v
        t=t_min !initialize the time step derivative
        !initial condition for for deq
        x=x0 
        vx=v0x
        y=y0
        vy=v0y
        if (flag.eq.1)then
                open(file=file_name_2,unit=11)
            !record position with flag=2
            else if (flag.eq.2)then
                open(file=file_name_3,unit=12)
            !record both if flag=3
            else if (flag.eq.3)then
                open(file=file_name_2,unit=11)
                open(file=file_name_3,unit=12)
            endif
        do while(y>=0.0)!leveled projectile motion
            !first RK step 
            !updating the condition
            tmp(1)=x
            tmp(2)=vx
            tmp(3)=y
            tmp(4)=vy
            !for 4 deq
            k1(1)=tmp(2)
            k1(2)=(-D/m)*tmp(2)*sqrt(tmp(2)**2+tmp(4)**2)!d_vx(tmp(2),tmp(4))!(-D/m)*tmp(2)*sqrt(tmp(2)**2+tmp(4)**2) 
            k1(3)=tmp(4)
            k1(4)=-g+(-D/m)*tmp(4)*sqrt(tmp(2)**2+tmp(4)**2)!d_vy(tmp(2),tmp(4))
            !second RK step
            tmp(1)=x+0.5*k1(1)*dt
            tmp(2)=vx+0.5*k1(2)*dt
            tmp(3)=y+0.5*k1(3)*dt
            tmp(4)=vy+0.5*k1(4)*dt
            k2(1)=tmp(2)
            k2(2)=(-D/m)*tmp(2)*sqrt(tmp(2)**2+tmp(4)**2)
            k2(3)=tmp(4)
            k2(4)=-g+(-D/m)*tmp(4)*sqrt(tmp(2)**2+tmp(4)**2)
            !third RK step
            tmp(1)=x+0.5*k2(1)*dt
            tmp(2)=vx+0.5*k2(2)*dt
            tmp(3)=y+0.5*k2(3)*dt
            tmp(4)=vy+0.5*k2(4)*dt
            k3(1)=tmp(2)
            k3(2)=(-D/m)*tmp(2)*sqrt(tmp(2)**2+tmp(4)**2)
            k3(3)=tmp(4)
            k3(4)=-g+(-D/m)*tmp(4)*sqrt(tmp(2)**2+tmp(4)**2)
            !fourth RK step
            tmp(1)=x+k3(1)*dt
            tmp(2)=vx+k3(2)*dt
            tmp(3)=y+k3(3)*dt
            tmp(4)=vy+k3(4)*dt
            k4(1)=tmp(2)
            k4(2)=(-D/m)*tmp(2)*sqrt(tmp(2)**2+tmp(4)**2)!d_vx(tmp(2),tmp(4))
            k4(3)=tmp(4)
            k4(4)=-g+(-D/m)*tmp(4)*sqrt(tmp(2)**2+tmp(4)**2)!d_vy(tmp(2),tmp(4))
            !calculate the postition and velocities
            x=x+dt*(k1(1)+2.0*(k2(1)+k3(1))+k4(1))/6.0
            vx=vx+dt*(k1(2)+2.0*(k2(2)+k3(2))+k4(2))/6.0
            y=y+dt*(k1(3)+2.0*(k2(3)+k3(3))+k4(3))/6.0
            vy=vy+dt*(k1(4)+2.0*(k2(4)+k3(4))+k4(4))/6.0
            t=t+dt
            !record position with flag=1
            if (flag.eq.1)then
                write(11,*) t,x,y
            !record position with flag=2
            else if (flag.eq.2)then
                write(12,*) t,vx,vy
            !record both if flag=3
            else if (flag.eq.3)then
                write(11,*) t,x,y
                write(12,*) t,vx,vy
            endif
            !update time
        enddo
        if (flag.eq.1)then
            close(unit=11)
        !record position with flag=2
        else if (flag.eq.2)then
            close(unit=12)
        !record both if flag=3
        else if (flag.eq.3)then
            close(unit=11)
            close(unit=12)
        endif
        t1=t
end subroutine

program hw5_3
    implicit none
    double precision :: t1 = 3, t2 = 1, dt, v0=0
    double precision :: v1=-5, v2=-1, v3=5, v4 =20
    double precision :: t
    double precision :: angle3=30.0, angle4=160, angle1=60, angle2=90.0
    integer :: i, flag=3
    character(len=100) :: str
    double precision :: v_min = -10, v_max = 10, dv,v
    ! investigating convergence of a good timestep
        ! find the timesteps
        ! investigate convergence
        dt = (t1-t2)/4
        do i=1, 4
            write(str, '(I5.5)') i ! convert integer to string
            call poncho_stat(t2, v0, angle2, str, flag, t)
            t2 = t2 + dt
        enddo
        str = ".01"
        t2=.1
        call poncho_stat(t2, v0, angle2, str, flag, t)
    ! a good time step is dt=.01
    ! investigate the time convergence of velocity
        flag = 2
        t2 = .001
        str = '-5_90'
        call poncho_stat(t2, v1, angle2, str, flag, t) !v=-5m/s, angle at 90
        str = '-1_90'
        call poncho_stat(t2, v2, angle2, str, flag, t) !v=-1m/s, angle at 90
        str = '5_90'
        call poncho_stat(t2, v3, angle2, str, flag, t) !v=5m/s, angle at 90
        str = '5_60'
        call poncho_stat(t2, v3, angle1, str, flag, t) !v=5m/s, angle at 60
        str = '5_30'
        call poncho_stat(t2, v3, angle3, str, flag, t) !v=5m/s, angle at 30
        str = '20_160'
        call poncho_stat(t2, v4, angle4, str, flag, t) !v=20m/s, angle at 160

    ! find the time it takes to hit the ground
        flag=4
        dv=(v_max-v_min)/100
        t2=.001
        v = v_min
        open(file="/Users/phihung/NumMethod/first/homework/hw_05/time_velocity",unit=11)
            do i=1, 100
                call poncho_stat(t2,v,angle2,str,flag,t)
                write(11,*) t, v
                v = v + dv
            enddo
        close(unit=11)
end program hw5_3