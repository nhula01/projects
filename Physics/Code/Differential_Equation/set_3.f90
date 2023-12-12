subroutine poncho_stat(dt, v0, angle, file_name)
    implicit none
        double precision, parameter :: pi=3.141592653589793
        double precision, parameter :: g=9.81 !free fall acceleration
        double precision, parameter :: m=.45,D=.02 !mass, air drag coefficient
        double precision, parameter :: x0=0.0,y0=100.0 !initial condition x(0), y(0)
        double precision, intent(in):: angle 
        character(len=100), intent(in) :: file_name
        character(len=100) :: file_name_2, file_name_3
        double precision, intent(in) :: v0 !initial speed
        double precision, intent(in) :: dt
        double precision v0x,v0y !initial velocity projections
        double precision, parameter :: t_min=0.0d0,t_max=100.0d0 !interval of integration of ODE
        integer, parameter :: N_steps=1000 !number of steps total time of flight
        double precision t_ideal, d_vx, d_vy
        double precision t,x,vx,y,vy, alpha
        double precision tmp1(4)
        double precision k1(4),k2(4),k3(4),k4(4),tmp(4)
    !find initial velocity
        alpha=angle*pi/180.0
        v0x=v0*cos(alpha)
        v0y=v0*sin(alpha)
        file_name_2 = '/Users/phihung/NumMethod/first/homework/hw_05/xy_projectile_motion_RK4_1' // file_name
        file_name_3 = '/Users/phihung/NumMethod/first/homework/hw_05/vxvy_projectile_motion_RK4_1' // file_name
        !t_ideal=2.0*v0y/g
    !dt=t_ideal/N_steps
    !<--- RK4 method for four differential equations
    !x_dot = v_x
    !v_x = -(D/m)*(v_x)*v
    !y_dot = v_y
    !v_y = -g-(D/m)*(v_y)*v
    open(file=file_name_2,unit=11)
    open(file=file_name_3,unit=12)
        t=t_min
        x=x0
        vx=v0x
        y=y0
        vy=v0y
        do while(y>=0.0)!leveled projectile motion
            !first RK step for 
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
            k2(2)=(-D/m)*tmp(2)*sqrt(tmp(2)**2+tmp(4)**2)!d_vx(tmp(2),tmp(4))
            k2(3)=tmp(4)
            k2(4)=-g+(-D/m)*tmp(4)*sqrt(tmp(2)**2+tmp(4)**2)!d_vy(tmp(2),tmp(4))

            !third RK step
            tmp(1)=x+0.5*k2(1)*dt
            tmp(2)=vx+0.5*k2(2)*dt
            tmp(3)=y+0.5*k2(3)*dt
            tmp(4)=vy+0.5*k2(4)*dt
            k3(1)=tmp(2)
            k3(2)=(-D/m)*tmp(2)*sqrt(tmp(2)**2+tmp(4)**2)!d_vx(tmp(2),tmp(4))
            k3(3)=tmp(4)
            k3(4)=-g+(-D/m)*tmp(4)*sqrt(tmp(2)**2+tmp(4)**2)!d_vy(tmp(2),tmp(4))

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

            !update time
            t=t+dt
            write(11,*) t,x,y
            write(12,*) t,vx,vy
        enddo
    close(unit=11)
    close(unit=12)
end subroutine
subroutine time_vel(dt, v0, angle,t)
    implicit none
    double precision, parameter :: pi=3.141592653589793
    double precision, parameter :: g=9.81 !free fall acceleration
    double precision, parameter :: m=.45,D=.02 !mass, air drag coefficient
    double precision, parameter :: x0=0.0,y0=100.0 !initial condition x(0), y(0)
    double precision, intent(in):: angle 
    double precision, intent(in) :: v0 !initial speed
    double precision, intent(in) :: dt
    double precision, intent(out) :: t
    double precision v0x,v0y !initial velocity projections
    double precision, parameter :: t_min=0.0d0,t_max=100.0d0 !interval of integration of ODE
    double precision t_ideal, d_vx, d_vy
    double precision x,vx,y,vy, alpha
    double precision tmp1(4)
    double precision k1(4),k2(4),k3(4),k4(4),tmp(4)

    alpha=angle*pi/180.0 !launch angle in rad ->55 degrees
    v0x=v0*cos(alpha)
    v0y=v0*sin(alpha)
    !t_ideal=2.0*v0y/g
    !dt=t_ideal/N_steps
    !<--- RK4 method
    t=t_min
    x=x0
    vx=v0x
    y=y0
    vy=v0y
    do while(y>=0.0)!leveled projectile motion
        !first RK step
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
        k2(2)=(-D/m)*tmp(2)*sqrt(tmp(2)**2+tmp(4)**2)!d_vx(tmp(2),tmp(4))
        k2(3)=tmp(4)
        k2(4)=-g+(-D/m)*tmp(4)*sqrt(tmp(2)**2+tmp(4)**2)!d_vy(tmp(2),tmp(4))

        !third RK step
        tmp(1)=x+0.5*k2(1)*dt
        tmp(2)=vx+0.5*k2(2)*dt
        tmp(3)=y+0.5*k2(3)*dt
        tmp(4)=vy+0.5*k2(4)*dt
        k3(1)=tmp(2)
        k3(2)=(-D/m)*tmp(2)*sqrt(tmp(2)**2+tmp(4)**2)!d_vx(tmp(2),tmp(4))
        k3(3)=tmp(4)
        k3(4)=-g+(-D/m)*tmp(4)*sqrt(tmp(2)**2+tmp(4)**2)!d_vy(tmp(2),tmp(4))

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
        !update time
        t=t+dt
    enddo
end subroutine

program hw5_3
    implicit none
    character(len=100) :: file_name_2="1d-2", file_name_3=".5", file_name_4="1d-3"
    double precision :: v=0.0, angle=270.0, dt=1.0d-2,t, dx = .1
    integer :: i 
    call poncho_stat(dt,v,angle,file_name_2)
    open(file="/Users/phihung/NumMethod/first/homework/hw_05/time_velocity.dat",unit=10)
        do i=1, 1000
            v = v + dx * i
            call time_vel(dt,v,angle,t)
            write(10,*) v, t
        enddo
    close(unit=10)
end program hw5_3