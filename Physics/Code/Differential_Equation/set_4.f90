subroutine problem4(tau,a_e_pop)
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
        integer, parameter :: n_fs=170 !20->150       !<--- number of fs for time dynamics 60 fs
        integer, parameter :: Nt=n_fs*Nt_1fs          !total number of time iterations over all thee time
        integer, parameter :: n_skip=20 !detect every n_skip steps to save some memory
        double precision detection(Nt/n_skip) !observable we are detecting

    !~~~======== laser pulse parameters =========~~~!
        double precision, intent(in) :: tau!=50.0d-15 !pulse duration, 50 fs
        double precision, parameter :: E0=5.0d8     !pulse amplitude 0.6 nm/V - very strong field
        double precision, parameter :: omega=ev_to_radsec*2.0d0 !carrier frequency (number is in eV)
                                                                !here we making it resonant with omega_eg below
        double precision pulse,tmp,BH_envelope

    !<----<----<----<----<----<----<----<------------
    !~~~=========== two-level emitter ===========~~~!
    !<----<----<----<----<----<----<----<------------

        double precision, parameter :: dp=Debye_to_Cm*10.0          !transition dipole, 10 Debye - typical value
        double precision, parameter :: omega_eg=ev_to_radsec*2.0    !transition frequency, 2 eV
        double complex a_g,a_e                        !quantum amplitudes
        double precision, intent(out) :: a_e_pop  
        double complex a_g_tmp,a_e_tmp                !for RK4
        double complex kk1_g,kk1_e,kk2_g,kk2_e,kk3_g,kk3_e,kk4_g,kk4_e  !for RK4
        integer n
        double precision t,t1,t2,t3,t4
        double precision cpu1,cpu2

    a_g=(1.0,0.0) !amplitude of the ground state at t=0 [set to 1] (real,complex)
    a_e=(0.0,0.0) !amplitude of the excited state at t=0 [set to 0]

    t=0.0 ! start at 20fs and end at 150fs

    ! find the value of constants ag, ae over time
    do n=1,Nt
    !~==========~~~~~~~~~~ first RK step ~==========~~~~~~~~~~
        !---- pulse at t----!
        t1=t
        ! epsilon naught
        pulse=E0*sin(omega*t1)*BH_envelope(tau,t1) !what the pulse is like from electricfield

        ! we have 2 deq
        a_g_tmp=a_g
        a_e_tmp=a_e ! these two lines are pretty useless
        
        !first rk
        kk1_g=Im*dp*pulse*a_e_tmp/h
        kk1_e=-Im*omega_eg*a_e_tmp+Im*dp*pulse*a_g_tmp/h
        kk1_g=dt*kk1_g
        kk1_e=dt*kk1_e

        !~==========~~~~~~~~~~ second RK step ~==========~~~~~~~~~~
        !---- pulse at t+dt/2----!
        t2=t+dt/2.0
        pulse=E0*sin(omega*t2)*BH_envelope(tau,t2)
        a_g_tmp=a_g+kk1_g/2.0
        a_e_tmp=a_e+kk1_e/2.0
        kk2_g=Im*dp*pulse*a_e_tmp/h
        kk2_e=-Im*omega_eg*a_e_tmp+Im*dp*pulse*a_g_tmp/h
        kk2_g=dt*kk2_g ! rescalee the value
        kk2_e=dt*kk2_e

        !~==========~~~~~~~~~~ third RK step ~==========~~~~~~~~~~
        !---- pulse at t+dt/2----! is the same as in previous part so we won't need to calculate it again here
        a_g_tmp=a_g+kk2_g/2.0
        a_e_tmp=a_e+kk2_e/2.0
        kk3_g=Im*dp*pulse*a_e_tmp/h
        kk3_e=-Im*omega_eg*a_e_tmp+Im*dp*pulse*a_g_tmp/h
        kk3_g=dt*kk3_g
        kk3_e=dt*kk3_e

        !~==========~~~~~~~~~~ fourth RK step ~==========~~~~~~~~~~
        !---- pulse at t+dt----!
        t4=t+dt
        pulse=E0*sin(omega*t4)*BH_envelope(tau,t4)

        a_g_tmp=a_g+kk3_g
        a_e_tmp=a_e+kk3_e

        kk4_g=Im*dp*pulse*a_e_tmp/h
        kk4_e=-Im*omega_eg*a_e_tmp+Im*dp*pulse*a_g_tmp/h

        kk4_g=dt*kk4_g
        kk4_e=dt*kk4_e

        !~==========~~~~~~~~~~ updates of a ~==========~~~~~~~~~~
        a_g=a_g+(kk1_g+2.0*(kk2_g+kk3_g)+kk4_g)/6.0
        a_e=a_e+(kk1_e+2.0*(kk2_e+kk3_e)+kk4_e)/6.0

        !--------------------------------------------------------------------------!
        !~~~~~~~~~~~~~~~~~~~~~~~~         detection          ~~~~~~~~~~~~~~~~~~~~~~!
        !--------------------------------------------------------------------------!
        if(mod(n,n_skip)==0)then
            detection(n/n_skip)=abs(a_e)**2 !population of the excited state vs time ! detect the value 
        endif

        t=t+dt ! update time 
    enddo !Nt
    a_e_pop = abs(a_e)**2
    
    open(file='/Users/phihung/NumMethod/first/homework/hw_05/2level_system_5.dat',unit=32)
    do n=1,Nt/n_skip
        write(32,*) dt*n*n_skip*1.0d15,detection(n) ! times 1d15 to show unit in fs
    enddo
    close(unit=32)

end subroutine problem4

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
endfunction BH_envelope


program hw5_4
    implicit none
    ! tau!=50.0d-15
    double precision :: tau_min = 20.0d-15, tau_max = 150.0d-15, dtau
    double precision ::  a_e_pop
    integer :: i
    ! adjust electric field and compare
    dtau = (tau_max - tau_min) / 200
    open(file="/Users/phihung/NumMethod/first/homework/hw_05/tau_e_pop_5", unit=10)
        do i=1, 200
            call problem4(tau_min,a_e_pop)
            write(10,*) tau_min,a_e_pop
            tau_min = tau_min + dtau
        enddo
    close(unit=10)
end program 