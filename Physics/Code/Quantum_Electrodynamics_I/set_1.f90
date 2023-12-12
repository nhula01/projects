subroutine speed(dist,value)
    implicit none
    !~~~===== fundamental constants =====~~~!
    double precision, parameter :: pi=3.141592653589793,c=299792458.0
    double precision, parameter :: sqrt2=1.414213562373095,sqrt3=1.7320508075688772
    double precision, parameter :: mu0=4.0D-7*pi,eps0=1.0/(c*c*mu0)
    double precision, parameter :: h=1.054571628D-34
    double complex, parameter :: Im=(0.0,1.0)
    double precision, parameter :: ev_to_radsec=2.0*pi*2.4180d14
    !~~~============== grid =============~~~!
    double precision, parameter :: dz=1.0d-9 !spatial resolution is 1 nm here (everything is in SI units)
    integer, parameter :: Nz=3000 !number of steps along z
    double precision, parameter :: z0=-dz*Nz*0.5,zM=dz*Nz*0.5-dz !grid is between z0 and zM
    double precision z(Nz) !coordinate z
    double precision, parameter :: dt=dz/(2.0*c) !time step
    !~~~======== time iterations ========~~~!
    integer, parameter :: Nt_1fs=int(1.0d-15/dt)  !number of time steps per 1 fs
    integer, parameter :: n_fs=35                 !<--- number of fs for free time dynamics when detecting
    integer, parameter :: Nt=n_fs*Nt_1fs          !total number of time iterations
    double precision detection(Nt) !observable we are detecting
    !~~~======== full EM fields =========~~~!
    double precision Ex(Nz),Hy(Nz-1)
    double precision, allocatable :: Ext(:,:),Hyt(:,:) !electric & magnetic fields
    double precision, parameter :: eps_delectric=1.0 !permittivity of medium (1 for vacuum)
    double precision, parameter :: dt_eps0=dt/(eps0*eps_delectric),dt_mu0=dt/mu0 !constants for Maxwell's eq.
    !~~~======== spatial position of the source =========~~~!
    integer, parameter :: k_source=100 !placing source 100 steps to the right of z0
    !~~~======== laser pulse parameters =========~~~!
    double precision, parameter :: tau=1.0d-15 !pulse duration
    double precision, parameter :: E0=1.0d0    !pulse amplitude
    !double precision, parameter :: omega=ev_to_radsec*3.0d0 !carrier frequency (number is in eV)
    double precision, parameter, dimension (4) :: aBH=(/0.353222222d0,-0.488d0,0.145d0,-0.010222222d0/)
    double precision pulse
    integer :: n_save=100,flag=0
    !~~~= other parameters and variables =~~~!
    integer :: k,n,i,nt_step,nt_tmp, nt_psi, t_temp
    double precision t,tmp, det1, det2,det_temp,n_temp,test_n,value
    integer :: dist
    integer :: detector=500
    
    detector = detector-dist
    !<--- grid
    do k=1,Nz
        z(k)=z0+dz*(k-1)
    enddo

    !count times to save
    nt_tmp=0
    nt_step=0
    t=0.0
    do while(nt_step<=Nt)
        t=t+dt
        nt_step=nt_step+1
        if(mod(nt_step,n_save)==0)then
            nt_tmp=nt_tmp+1
        endif
    enddo

    !allocate the arrays
    nt_psi=nt_tmp
    allocate(Ext(Nz,nt_psi))
    allocate(Hyt(Nz-1,nt_psi)) !i used Nz
    Ex=0.0
    Hy=0.0

    t_temp=0
    flag=0
    det1=0
    det2=0
    do n=1,Nt
        !---- pulse ----!
        t=dt*dble(n)
        if(t<=tau)then
            tmp= &
                        aBH(1)+ &
                aBH(2)*cos(2.0*pi*t/tau)+ &
                aBH(3)*cos(2.0*pi*2.0*t/tau)+ &
                aBH(4)*cos(2.0*pi*3.0*t/tau)
        else
            tmp=0.0
        endif
        pulse=tmp*E0!*sin(omega*t)

        !~~~~~~ Hy ~~~~~~~~!
        do k=1,Nz-1
            Hy(k)=Hy(k)+dt_mu0*(Ex(k)-Ex(k+1))/dz
        enddo
        !~~~~~ Ex ~~~~~~~~~!
        do k=2,Nz-1
            if(k==k_source)then
            Ex(k)=Ex(k)+dt_eps0*(Hy(k-1)-Hy(k))/dz+pulse
            else
            Ex(k)=Ex(k)+dt_eps0*(Hy(k-1)-Hy(k))/dz
            endif
        enddo

        if(mod(n,n_save)==0)then
            t_temp = t_temp+1
            Hyt(:,t_temp)=Hy(:)
            Ext(:,t_temp)=Ex(:)
        endif
    !--------------------------------------------------------------------------!
    !~~~~~~~~~~~~~~~~~~~~~~~~         detection          ~~~~~~~~~~~~~~~~~~~~~~!
    !--------------------------------------------------------------------------!
        detection(n)=Ex(Nz-detector)!detecting Ex 500 nm to the left of zM
        if((detection(n)>.9964).and.(flag==0))then
            det1=n*dt*1.0d15 !determine the time
            write(*,*) det1, detection(n)
            flag=1
            n_temp=n*dt*1.0d15
        endif
        test_n=(n*dt*1.0d15)-n_temp
        if((detection(n)<(-.9964)).and.(flag==1).and.(test_n>2.0))then
            det2=n*dt*1.0d15
            write(*,*) det2, detection(n)
            flag=2
        endif
    enddo !Nt

    value=2*detector*1.0d-9/(det2-det1)*1.0d15
    !open(file='/Users/phihung/NumMethod/first/homework/hw_09/speed_right_10.dat',unit=32)
    !    write(32,*) (det2-det1),2*detector*1.0d-9,2*detector*1.0d-9/(det2-det1)*1.0d15
    !close(unit=32)

    !open(file='/Users/phihung/NumMethod/first/homework/hw_09/electric_t.dat',unit=32)
    !do i=1,Nz
    !    write(32,*) (Ext(i,n),n=1,nt_psi) !here time in fs
    !enddo
    !close(unit=32)

    !open(file='/Users/phihung/NumMethod/first/homework/hw_09/magnetic_t.dat',unit=32)
    !do i=1,Nz-1
    !    write(32,*) (Hyt(i,n),n=1,nt_psi) !here time in fs
    !enddo
    !close(unit=32)

    !open(file='/Users/phihung/NumMethod/first/homework/hw_09/detection.dat',unit=32)
    !do n=1,Nt
    !    if(mod(n,n_save)==0)then
    !        write(32,*) dt*n*1.0d15,detection(n) !here time in fs
    !    endif
    !enddo
    !close(unit=32)
    !open(file='/Users/phihung/NumMethod/first/homework/hw_09/test1.dat',unit=32)
    !do n=1,Nt
    !    write(32,*) dt*n*1.0d15,detection(n) !here time in fs
    !enddo
    !close(unit=32)
end subroutine speed

program set1
    implicit none
    integer::i_start, i_final=50
    double precision :: value
    open(file='/Users/phihung/NumMethod/first/homework/hw_09/speed_distance.dat',unit=32)
        do i_start=-50,i_final, 10
            call speed(i_start,value)
            write(32,*) i_start,value
        enddo
    close(unit=32)
end program set1