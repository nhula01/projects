subroutine problem5(delta1,E_max)
    implicit none
    !<----<----<----<----<----<----<----<----!
    !~~~===== fundamental constants =====~~~!
    !<----<----<----<----<----<----<----<----!
        double precision, parameter :: pi=3.141592653589793,c=299792458.0
        double precision, parameter :: sqrt2=1.414213562373095,sqrt3=1.7320508075688772
        double precision, parameter :: mu0=4.0D-7*pi,eps0=1.0/(c*c*mu0)
        double precision, parameter :: h=1.054571628D-34 !reduced Plank constant hbar
        double complex, parameter :: Im=(0.0,1.0)
        double precision, parameter :: ev_to_radsec=2.0*pi*2.4180d14 !eV to rad per second
        double precision, parameter :: Debye_to_Cm=3.33564d-30

    !<----<----<----<----<----<----<----<----
    !~~~======== time iterations ========~~~!
    !<----<----<----<----<----<----<----<----
        double precision, parameter :: dt=0.1d-15 !time step <-- 0.1 fs
        integer, parameter :: Nt_1fs=int(1.0d-15/dt)  !number of time steps per 1 fs
        integer, parameter :: n_fs=60                 !<--- number of fs for time dynamics
        integer, parameter :: Nt=n_fs*Nt_1fs          !total number of time iterations

    !~~~======== laser pulse parameters =========~~~!
        double precision, parameter :: tau=50.0d-15 !pulse duration, 50 fs
        double precision E0                         !pulse amplitude - we are scanning this parameter
        integer, parameter :: N_E0=1000              ! number of steps for E0
        double precision, parameter :: E0_min=0.5d8,E0_max=1.0d10 !min and max values for E0
        double precision detection(N_E0,3) !observable as a function of E0
        double precision temp
        double precision, parameter :: omega=ev_to_radsec*2.0d0 !carrier frequency (number is in eV)
                                                                !here we making it resonant with omega_eg below
        double precision pulse,tmp,BH_envelope

    !<----<----<----<----<----<----<----<------------

    !~~~=========== two-level emitters ==========~~~!
        !<----<----<----<----<----<----<----<------------
        double precision, intent(in) :: delta1
        double precision, intent(out) :: E_max

        double precision :: delta !coupling between emitters
        double precision, parameter :: dp=Debye_to_Cm*10.0          !transition dipole, 10 Debye - typical value
        double precision, parameter :: omega_eg=ev_to_radsec*2.0    !transition frequency, 2 eV
        double complex Campl(4)                        !quantum amplitudes
        double complex Campl_tmp(4)                    !for RK4
        double complex kk1(4),kk2(4),kk3(4),kk4(4)  !for RK4
        double complex Hamiltonian(4,4) !this is Hamiltonian/hbar
        integer n,m
        double precision t,t1,t2,t3,t4
        double precision :: cpu1,cpu2, pop_max

    delta=ev_to_radsec*delta1
    E_max = 0
    pop_max = 0
    do m=1,N_E0 ! loop through a range of electric field
        temp=log10(E0_min)+(log10(E0_max)-log10(E0_min))*float(m-1)/float(N_E0-1)
        E0=10.0d0**temp

        Campl=(0.0d0,0.0d0)
        Campl(1)=(1.0d0,0.0d0) !both emitters are in the ground state
        Hamiltonian=(0.0d0,0.0d0) !H
        Hamiltonian(2,2)=omega_eg
        Hamiltonian(2,3)=-delta
        Hamiltonian(3,2)=-delta
        Hamiltonian(3,3)=omega_eg
        Hamiltonian(4,4)=2.0*omega_eg
        t=0.0
        
        ! loop through an electric field
        do n=1,Nt
        !~==========~~~~~~~~~~ first RK step ~==========~~~~~~~~~~
            !---- pulse at t----!
            t1=t
            pulse=E0*sin(omega*t1)*BH_envelope(tau,t1)

            Hamiltonian(1,2)=-pulse*dp/h
            Hamiltonian(1,3)=Hamiltonian(1,2)
            Hamiltonian(2,1)=Hamiltonian(1,2)
            Hamiltonian(2,4)=Hamiltonian(1,2)
            Hamiltonian(3,1)=Hamiltonian(1,2)
            Hamiltonian(3,4)=Hamiltonian(1,2)
            Hamiltonian(4,2)=Hamiltonian(1,2)
            Hamiltonian(4,3)=Hamiltonian(1,2)

            Campl_tmp=Campl

            kk1=-Im*matmul(Hamiltonian,Campl_tmp)*dt

            !~==========~~~~~~~~~~ second RK step ~==========~~~~~~~~~~
            !---- pulse at t+dt/2----!
            t2=t+dt/2.0
            pulse=E0*sin(omega*t2)*BH_envelope(tau,t2)

            Hamiltonian(1,2)=-pulse*dp/h
            Hamiltonian(1,3)=Hamiltonian(1,2)
            Hamiltonian(2,1)=Hamiltonian(1,2)
            Hamiltonian(2,4)=Hamiltonian(1,2)
            Hamiltonian(3,1)=Hamiltonian(1,2)
            Hamiltonian(3,4)=Hamiltonian(1,2)
            Hamiltonian(4,2)=Hamiltonian(1,2)
            Hamiltonian(4,3)=Hamiltonian(1,2)

            Campl_tmp=Campl+kk1/2.0

            kk2=-Im*matmul(Hamiltonian,Campl_tmp)*dt

            !~==========~~~~~~~~~~ third RK step ~==========~~~~~~~~~~
            !---- pulse at t+dt/2----! is the same as in previous part so we won't need to calculate it again here
            Campl_tmp=Campl+kk2/2.0

            kk3=-Im*matmul(Hamiltonian,Campl_tmp)*dt

            !~==========~~~~~~~~~~ fourth RK step ~==========~~~~~~~~~~
            !---- pulse at t+dt----!
            t4=t+dt
            pulse=E0*sin(omega*t4)*BH_envelope(tau,t4)

            Hamiltonian(1,2)=-pulse*dp/h
            Hamiltonian(1,3)=Hamiltonian(1,2)
            Hamiltonian(2,1)=Hamiltonian(1,2)
            Hamiltonian(2,4)=Hamiltonian(1,2)
            Hamiltonian(3,1)=Hamiltonian(1,2)
            Hamiltonian(3,4)=Hamiltonian(1,2)
            Hamiltonian(4,2)=Hamiltonian(1,2)
            Hamiltonian(4,3)=Hamiltonian(1,2)

            Campl_tmp=Campl+kk3

            kk4=-Im*matmul(Hamiltonian,Campl_tmp)*dt

            !~==========~~~~~~~~~~ updates of C ~==========~~~~~~~~~~
            Campl=Campl+(kk1+2.0*(kk2+kk3)+kk4)/6.0

            t=t+dt
        enddo !Nt

        if (pop_max <= abs(Campl(4))**2) then !(if we dont the maxium populaton)
            E_max = E0 ! replace the new ma
            pop_max = abs(Campl(4))**2
        else if (pop_max >= abs(Campl(4))**2) then
            exit
        endif

    enddo !N_E0
end subroutine problem5

double precision function BH_envelope(tau,t)
  implicit none
  double precision, parameter :: pi=3.141592653589793
  double precision, parameter, dimension (4) :: aBH=(/0.353222222d0,-0.488d0,0.145d0,-0.010222222d0/) !BH envelope
  double precision tau,t,tmp

  if(t<=tau)then
    tmp= &
                  aBH(1)+ &
         aBH(2)*cos(2.0*pi*t/tau)+ &
         aBH(3)*cos(2.0*pi*2.0*t/tau)+ &
         aBH(4)*cos(2.0*pi*3.0*t/tau)
   else
    tmp=0.0
  endif

  BH_envelope=tmp
end function BH_envelope

program make
    implicit none
    double precision :: delta_max, delta_min, ddelta, e0
    integer :: i
    delta_min = 0
    delta_max= 5
    ddelta = (delta_max-delta_min)/200
    open(file="/Users/phihung/NumMethod/first/homework/hw_05/delta_E01.dat",unit=10)
        do i=1, 200
            call problem5(delta_min, e0)
            write(10,*) delta_min, e0
            delta_min = delta_min + ddelta
        enddo
    close(unit=10)
end program make