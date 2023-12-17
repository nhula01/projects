program main
    implicit none
    character(len=100) :: name1='0.dat', name2='diode1.dat', name3='doublewell01.dat'
    character(len=100) :: name4='10.dat', name5='50.dat', name6='alter.dat'
    integer :: option=9
    double complex :: initial_cond = (1.0d0,1.0d0)
    call scattering(name1, option)
    !call tunnel(initial_cond, name6)
end program main

subroutine scattering(name,potentialflag)
    implicit none
    character(len=100) :: name, file_name
    integer :: potentialflag !1=constant, 2=alternating
    double precision, parameter :: pi=3.141592653589793
    double complex, parameter :: Im=(0.0d0,1.0d0) !imaginary unit
    double precision, parameter :: hbar=1.0 !atomic units
    double precision, parameter :: mass=5.0 !1 in atomic units means it is electron
    double precision, parameter :: dt=1.0d-3 !time step for split-operator
    double precision, parameter :: t_final=20.0 !total propagation time
    integer, parameter :: nt_final=int(t_final/dt)
    double precision, parameter :: a=-80.0,b=80.0 !coordinate interval
    integer, parameter :: N=2048 !total number of points - for fast results make it 2**n
    double precision, parameter :: dx=(b-a)/(N+1), dp=2*pi/(dx*N) !space/momentum
    double precision x,t
    double precision wsave(4*N+15) !working array for FFT
    double precision tmp,tmp1,tmp2,potential,p, Apotential, potential2, oscillator, twopotential
    double complex psi(N),exp_T(N),exp_V,initial_wavepacket
    integer i,nt
    integer, parameter :: n_save=100 !save wavefunction every n_save steps
    integer, parameter :: i_min=1,i_max=N !interval for saving psi
    integer nt_tmp,nt_psi !number of time points for the wavefunction plot
    double precision, allocatable :: psi_save(:,:),psi_save_momentum(:,:) !i_min:i_max,nt_psi
    double complex, allocatable :: psi_save_parts(:,:),psi_save_momentum_parts(:,:) !i_min:i_max,nt_psi
    double precision av_x(nt_final),delta_x(nt_final),delta_p(nt_final),av_p(nt_final)
    double precision :: Ns=45.2548339959
    !<---initialize fftpack -- call it only once
    call zffti(N,wsave)
    !initial wavefunction
    do i=1,N
        x=a+dx*float(i)
        psi(i)=initial_wavepacket(x)
    enddo

    !<-- normalization
    tmp=0.0
    do i=1,N
        tmp=tmp+abs(psi(i))**2
    enddo
    tmp=tmp*dx
    psi=psi/sqrt(tmp)
    
    !initial state is ready!

    !<-- precalculate exponents with T for split-operator propagation
    !exp(-i*T*dt/(2*hbar)) * exp(-i*V*dt/hbar) * exp(-i*T*dt/(2*hbar))
    do i=1,N
        if(i<=(N/2))then !<-- momentum in Fourier space
            p=hbar*2.0*pi*(i-1)/(dx*N) !i*dp positive momentum
            else
            p=-hbar*2.0*pi*(N+1-i)/(dx*N) ! negative momentum
        endif
        exp_T(i)=exp(-Im*p**2*dt/(4.0*mass*hbar))
        !x=a+dx*float(i)
        !exp_V(i)=exp(-Im*potential(x)*dt/hbar)
    enddo

    !calculate number of time steps needed for wavefunction plot
    nt_tmp=0
    nt=0
    do while(nt<=nt_final)
        t=t+dt
        nt=nt+1
        if(mod(nt,n_save)==0)then
            nt_tmp=nt_tmp+1
        endif
    enddo

    !allocate the arrays
    nt_psi=nt_tmp
    allocate(psi_save(i_min:i_max,nt_psi))
    allocate(psi_save_parts(i_min:i_max,nt_psi))
    allocate(psi_save_momentum(i_min:i_max,nt_psi))
    allocate(psi_save_momentum_parts(i_min:i_max,nt_psi))
    psi_save=(0.0d0,0.0d0)
    psi_save_parts=(0.0d0,0.0d0)
    psi_save_momentum=(0.0d0,0.0d0)
    psi_save_momentum_parts=(0.0d0,0.0d0)

    t=0.0 !begin time propagation
    nt=0
    nt_tmp=0
    do while(nt<=nt_final)
        t=t+dt
        nt=nt+1

        !--- kinetic part ---!
        call zfftf(N,psi,wsave) !forward Fourier transform
        psi=exp_T*psi/float(N) !1/N is needed due to forward and backward fft
        call zfftb(N,psi,wsave) !backward Fourier transform

        !--- potential part ---!
        do i=1,N
            x=a+dx*float(i)
            if(potentialflag==2)then !ooption2=alterating potetial
                exp_V=exp(-Im*Apotential(x,t)*dt/hbar)
            elseif(potentialflag==3)then
                exp_V=exp(-Im*potential2(x)*dt/hbar) !double well potential
            elseif(potentialflag==4)then
                exp_V=exp(-Im*oscillator(t)*dt/hbar) !simple well potential
            elseif(potentialflag==5)then
                exp_V=exp(-Im*twopotential(x)*dt/hbar) !double well potential
            else
                exp_V=exp(-Im*potential(x)*dt/hbar) !option1 =  stationary coordinate
            endif
            psi(i)=exp_V*psi(i)
        enddo

        !--- kinetic part ---!
        call zfftf(N,psi,wsave) !forward Fourier transform
        psi=exp_T*psi/float(N) !1/N is needed due to fourth and back fft
        !~~~~~~~~~~~~Momentum Space~~~~~~~~~~~~~~~!
        !save momentum probability density space 
        if(mod(nt,n_save)==0)then !<-- save psi every n_save steps
            nt_tmp=nt_tmp+1
            do i=i_min,i_max
                psi_save_momentum(i,nt_tmp)=abs(psi(i)*Ns)**2!*sqrt(CMPLX(N,0.0d0)))**2
                psi_save_momentum_parts(i,nt_tmp)=psi(i)*Ns!*sqrt(CMPLX(N,0.0d0))
            enddo
        endif
        !--- calculate average p and deviation of p ---!
        tmp1=0.0
        tmp2=0.0
        do i=1,N
            if(i<=(N/2))then !<-- momentum in Fourier space
                p=hbar*2.0*pi*(i-1)/(dx*N) !positive
            else
                p=-hbar*2.0*pi*(N+1-i)/(dx*N) !negative
            endif
            tmp1=tmp1+p*abs(psi(i)*Ns)**2
            tmp2=tmp2+(p**2)*abs(psi(i)*Ns)**2
        enddo
        av_p(nt)=tmp1*dp !average p
        tmp2=tmp2*dp !average p**2
        delta_p(nt)=sqrt(abs(tmp2)-av_p(nt)**2)

        !psi=psi/Ns !1/sqrtN is needed due to fourth and back fft
        call zfftb(N,psi,wsave) !backward Fourier transform
        !~~~~~~~~~~~Coordinate Space~~~~~~~~~~~!
        !save the coordinate space
        if(mod(nt,n_save)==0)then !<-- save psi every n_save steps
            !nt_tmp=nt_tmp
            do i=i_min,i_max
                psi_save(i,nt_tmp)=abs(psi(i))**2
                psi_save_parts(i,nt_tmp)=psi(i)
            enddo
        endif

        !--- calculate average x and deviation of x ---!
        tmp1=0.0
        tmp2=0.0
        do i=1,N
            x=a+dx*float(i)
            tmp1=tmp1+x*abs(psi(i))**2
            tmp2=tmp2+(x**2)*abs(psi(i))**2
        enddo
        av_x(nt)=tmp1*dx !average x
        tmp2=tmp2*dx !average x**2
        delta_x(nt)=sqrt(abs(tmp2)-av_x(nt)**2) !uncertainty

        if(mod(nt,5000)==0)then
            write(*,*) nt !<--- print nt every 10000 steps
            !<-- test normalization
            tmp=0.0
            do i=1,N
                tmp=tmp+abs(psi(i))**2
            enddo
            tmp=tmp*dx
            write(*,*) 'norm',tmp
        endif
    enddo


    !~~~~~~~record average value
    file_name = '/Users/phihung/NumMethod/first/final/x_vs_t' // name
    open(file=file_name,unit=11)
        do nt=2,nt_final
            t=dt*nt
            write(11,*) t,av_x(nt),delta_x(nt)
        enddo
    close(unit=11)

    file_name = '/Users/phihung/NumMethod/first/final/p_vs_t' // name
    open(file=file_name,unit=11)
        do nt=2,nt_final
            t=dt*nt
            write(11,*) t,av_p(nt),delta_p(nt)
        enddo
    close(unit=11)

    !~~~~~~~~~~~~~~~density 
    file_name = '/Users/phihung/NumMethod/first/final/density_psi_xt' // name
    open(file=file_name,unit=11)
        do i=i_min,i_max
            write(11,*) (psi_save(i,nt),nt=1,nt_psi)
        enddo
    close(unit=11)
    file_name = '/Users/phihung/NumMethod/first/final/density_psi_pt' // name
    open(file=file_name,unit=11)
        do i=i_min,i_max
            write(11,*) (psi_save_momentum(i,nt),nt=1,nt_psi)
        enddo
    close(unit=11)

    !~~~~~~~~~~~~~real and im part in space coordinate
    file_name = '/Users/phihung/NumMethod/first/final/im_psi_xt' // name
    open(file=file_name,unit=12)
    file_name = '/Users/phihung/NumMethod/first/final/real_psi_xt' // name
    open(file=file_name,unit=11)
        do i=i_min,i_max
            write(11,*) (REAL(psi_save_parts(i,nt)),nt=1,nt_psi)
            write(12,*) (AIMAG(psi_save_parts(i,nt)),nt=1,nt_psi)
        enddo
    close(unit=11)
    close(unit=12)

    !~~~~~~~~~~~~~~~~~~~~real and im part in momentum coordinate
    file_name = '/Users/phihung/NumMethod/first/final/im_psi_pt' // name
    open(file=file_name,unit=12)
    file_name = '/Users/phihung/NumMethod/first/final/real_psi_pt' // name
    open(file=file_name,unit=11)
        do i=i_min,i_max
            write(11,*) (REAL(psi_save_momentum_parts(i,nt)),nt=1,nt_psi)
            write(12,*) (AIMAG(psi_save_momentum_parts(i,nt)),nt=1,nt_psi)
        enddo
    close(unit=11)
    close(unit=12)
end subroutine scattering

subroutine tunnel(initial_cond, name)
    implicit none
    double precision, parameter :: pi=3.141592653589793
    double complex, parameter :: Im=(0.0d0,1.0d0) !imaginary unit
    double precision, parameter :: hbar=1.0 !atomic units
    character(len=100) :: name, file_name
    double precision :: mass !1 in atomic units means it is electron
    double precision, parameter :: a=-5.0,b=5.0 !barrier interval
    integer, parameter :: N=500 !total number of points 
    double precision, parameter :: dx=(b-a)/(N-1) !N-1 segments in the barrier
    double precision :: x,t
    double precision :: potential2, E,dm_dx,tmp1,tmp2
    double complex :: psi(N),k1(2),k2(2),k3(2),k4(2),tmp(2)
    double complex :: dpsi,psi_0,initial_cond
    integer :: i,nt
    !find psi2
    E=1
    psi=(0d0,0d0)
    psi(1)=initial_cond
    file_name = '/Users/phihung/NumMethod/first/final/psi' // name
    open(file=file_name,unit=11)
        x=a
        psi_0=psi(1)
        dpsi=(0d0,0d0)
        do i=2,N
            !calculate the change of mass
            tmp1=mass(x+dx)                !function at x_point+dx
            tmp2=mass(x-dx)                !function at x_point-dx
            dm_dx=(tmp1-tmp2)/(2.0*dx)        !central difference
            !first RK step
            tmp(1)=psi_0
            tmp(2)=dpsi
            k1(1)=tmp(2)
            if(dm_dx==0.0)then
                k1(2)=2*mass(x)/hbar**2*tmp(1)*(potential2(x)-E)
            else
                k1(2)=-(mass(x)/dm_dx)*tmp(2)+2*mass(x)/hbar**2*tmp(1)*(potential2(x)-E)
            endif
            !second RK step
            tmp(1)=psi_0+0.5*k1(1)*dx
            tmp(2)=dpsi+0.5*k1(2)*dx
            k2(1)=tmp(2)
            if(dm_dx==0.0)then
                k2(2)=2*mass(x)/hbar**2*tmp(1)*(potential2(x)-E)
            else
                k2(2)=-(mass(x)/dm_dx)*tmp(2)+2*mass(x)/hbar**2*tmp(1)*(potential2(x)-E)
            endif
            !third RK step
            tmp(1)=psi_0+0.5*k2(1)*dx
            tmp(2)=dpsi+0.5*k2(2)*dx
            k3(1)=tmp(2)
            if(dm_dx==0.0)then
                k3(2)=2*mass(x)/hbar**2*tmp(1)*(potential2(x)-E)
            else
                k3(2)=-(mass(x)/dm_dx)*tmp(2)+2*mass(x)/hbar**2*tmp(1)*(potential2(x)-E)
            endif
            !fourth RK step
            tmp(1)=psi_0+k3(1)*dx
            tmp(2)=dpsi+k3(2)*dx
            k4(1)=tmp(2)
            if(dm_dx==0.0)then
                k4(2)=2*mass(x)/hbar**2*tmp(1)*(potential2(x)-E)
            else
                k4(2)=-(mass(x)/dm_dx)*tmp(2)+2*mass(x)/hbar**2*tmp(1)*(potential2(x)-E)
            endif
            psi_0=psi_0+dx*(k1(1)+2.0*(k2(1)+k3(1))+k4(1))/6.0
            dpsi=dpsi+dx*(k1(2)+2.0*(k2(2)+k3(2))+k4(2))/6.0
            x=a+dx*float(i-1)
            psi(i)=psi_0
            write(11,*) x,AIMAG(psi_0)
        enddo
    close(unit=11)

end subroutine tunnel


!temp1  = (hbar)**2/(2*dx**2*mass(a-dx))
!temp2 = (hbar)**2/(2*dx**2*mass(a))
!Kl = sqrt(2*mass(a-dx)*(E-potential2(a-dx)))/hbar
!psi(2) = (2*Im*temp1*sin(Kl*dx)-(temp2*exp(Im*Kl*dx)+E-potential2(x)-(temp1+temp2))*psi(1))/temp2
!do i=2, N-2
!    x=a+dx*float(i-1)
!    temp1 = (hbar)**2/(2*dx**2*mass(x))
!    temp2 = (hbar)**2/(2*dx**2*mass(x-dx))
!    psi(i+1) = (temp2*psi(i-1)+(E-potential2(x)-(temp1+temp2))*psi(i))/temp1
!enddo
!find psiN
!temp1 = (hbar)**2/(2*dx**2*mass(b-dx))
!temp2 = (hbar)**2/(2*dx**2*mass(b))
!KR = sqrt(2*mass(b+dx)*(E-potential2(b+dx)))/hbar
!psi(N)=-temp1*psi(N-1)/(E-potential2(b)-temp1-temp2+temp2*exp(Im*KR*dx))

double precision function mass(x)
    implicit none
    double precision :: x
    mass = 1d-2
end function mass

double precision function twopotential(x)
    implicit none
    double precision, parameter:: distance=20
    double precision, parameter :: width1=5.0d0 !width of well 1
    double precision, parameter :: width2=5.0d0 !width of well 2
    double precision a_V1,b_V1,x !first well
    double precision a_V2,b_V2  !second well
    double precision :: V01=5 !depth of the first well (all in atomic units)
    double precision :: V02=5 !depth of the second well (all in atomic units)
    double precision :: Vmiddle=0.0 !height of the barrier between wells
    a_V1=-distance/2.0-width1 
    b_V1=-distance/2.0 
    a_V2=distance/2.0
    b_V2=distance/2.0+width2
    if((x>=a_V1).and.(x<=b_V1))then
        twopotential=V01+1
    elseif((x>=a_V2).and.(x<=b_V2))then
        twopotential=V02+1
    elseif((x>b_V1).and.(x<a_V2))then !initially b_V2 ->error
        twopotential=0.0+1!(1/40)*x**2
    else
        twopotential=0.0+1
    endif
end function twopotential

double precision function potential2(x) ! potential for tuneeling
    implicit none
    double precision, parameter :: alpha=.005,r=2
    double precision, parameter :: V0=16
    double precision x,tmp,coshf
    potential2=-1/sqrt((x-10)**2/2000+alpha)-1/sqrt(x**2/2000+alpha)+V0 !shifr right 1 and energy in 
end function potential2

double precision function Apotential(x,t) ! the potential depends on time
    implicit none
    double precision, parameter :: pi=3.141592653589793
    double precision, parameter :: V0=20
    double precision, parameter :: omega=.3
    double precision x,tmp,coshf,t
    Apotential = V0 * abs(sin(2 * pi * omega * t) * exp(-0.2*(x-10.0)**2))
end function Apotential

double precision function oscillator(x) ! the potential depends on time
    implicit none
    double precision, parameter :: pi=3.141592653589793
    double precision, parameter :: V0=20
    double precision, parameter :: omega=.3
    double precision x,tmp,coshf,t
    oscillator =(1/40)*x**2-10
end function oscillator

double precision function potential(x) ! the potential doesnt depend on time
    implicit none
    double precision, parameter :: x0=10.0,width=7.5
    double precision, parameter :: V0=0
    double precision x,tmp,coshf
    if((x>=(x0-width/2.0)).and.(x<=(x0+width/2.0)))then
        potential=V0
        else
        potential=0.0
    endif
end function potential

double complex function initial_wavepacket(x)
    implicit none
    double precision, parameter :: pi=3.141592653589793
    double complex, parameter :: Im=(0.0d0,1.0d0)
    double precision, parameter :: hbar=1.0 !atomic units
    double precision, parameter :: beta=3.0 !see PHY314, Ch. 6, Gaussian wavepacket
    double precision, parameter :: p0=5.0 !initial momentum of the wavepacket
    double precision x,tmp
    tmp=sqrt(beta/(pi*hbar))
    initial_wavepacket=tmp*exp(Im*p0*x/hbar-(x*beta/hbar)**2)
end function initial_wavepacket

!Necessary subroutine
    SUBROUTINE ZFFTI (N,WSAVE)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       WSAVE(1)
    IF (N .EQ. 1) RETURN
    IW1 = N+N+1
    IW2 = IW1+N+N
    CALL CFFTI1 (N,WSAVE(IW1),WSAVE(IW2))
    RETURN
    END

    SUBROUTINE ZFFTF (N,C,WSAVE)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       C(1)       ,WSAVE(1)
    IF (N .EQ. 1) RETURN
    IW1 = N+N+1
    IW2 = IW1+N+N
    CALL CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
    RETURN
    END

    SUBROUTINE ZFFTB (N,C,WSAVE)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       C(1)       ,WSAVE(1)
    IF (N .EQ. 1) RETURN
    IW1 = N+N+1
    IW2 = IW1+N+N
    CALL CFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
    RETURN
    END

    SUBROUTINE CFFTI1 (N,WA,IFAC)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
    DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
    NL = N
    NF = 0
    J = 0
    101 J = J+1
    IF (J-4) 102,102,103
    102 NTRY = NTRYH(J)
    GO TO 104
    103 NTRY = NTRY+2
    104 NQ = NL/NTRY
    NR = NL-NTRY*NQ
    IF (NR) 101,105,101
    105 NF = NF+1
    IFAC(NF+2) = NTRY
    NL = NQ
    IF (NTRY .NE. 2) GO TO 107
    IF (NF .EQ. 1) GO TO 107
    DO 106 I=2,NF
    IB = NF-I+2
    IFAC(IB+2) = IFAC(IB+1)
    106 CONTINUE
    IFAC(3) = 2
    107 IF (NL .NE. 1) GO TO 104
    IFAC(1) = N
    IFAC(2) = NF
    TPI =  6.28318530717958647692D0
    ARGH = TPI/FLOAT(N)
    I = 2
    L1 = 1
    DO 110 K1=1,NF
    IP = IFAC(K1+2)
    LD = 0
    L2 = L1*IP
    IDO = N/L2
    IDOT = IDO+IDO+2
    IPM = IP-1
    DO 109 J=1,IPM
        I1 = I
        WA(I-1) = 1.0D0
        WA(I) = 0.0D0
        LD = LD+L1
        FI = 0.0D0
        ARGLD = FLOAT(LD)*ARGH
        DO 108 II=4,IDOT,2
            I = I+2
            FI = FI+1.D0
            ARG = FI*ARGLD
            WA(I-1) = COS(ARG)
            WA(I) = SIN(ARG)
    108       CONTINUE
        IF (IP .LE. 5) GO TO 109
        WA(I1-1) = WA(I-1)
        WA(I1) = WA(I)
    109    CONTINUE
    L1 = L2
    110 CONTINUE
    RETURN
    END

    SUBROUTINE CFFTF1 (N,C,CH,WA,IFAC)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
    NF = IFAC(2)
    NA = 0
    L1 = 1
    IW = 1
    DO 116 K1=1,NF
    IP = IFAC(K1+2)
    L2 = IP*L1
    IDO = N/L2
    IDOT = IDO+IDO
    IDL1 = IDOT*L1
    IF (IP .NE. 4) GO TO 103
    IX2 = IW+IDOT
    IX3 = IX2+IDOT
    IF (NA .NE. 0) GO TO 101
    CALL PASSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
    GO TO 102
    101    CALL PASSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
    102    NA = 1-NA
    GO TO 115
    103    IF (IP .NE. 2) GO TO 106
    IF (NA .NE. 0) GO TO 104
    CALL PASSF2 (IDOT,L1,C,CH,WA(IW))
    GO TO 105
    104    CALL PASSF2 (IDOT,L1,CH,C,WA(IW))
    105    NA = 1-NA
    GO TO 115
    106    IF (IP .NE. 3) GO TO 109
    IX2 = IW+IDOT
    IF (NA .NE. 0) GO TO 107
    CALL PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
    GO TO 108
    107    CALL PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
    108    NA = 1-NA
    GO TO 115
    109    IF (IP .NE. 5) GO TO 112
    IX2 = IW+IDOT
    IX3 = IX2+IDOT
    IX4 = IX3+IDOT
    IF (NA .NE. 0) GO TO 110
    CALL PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
    GO TO 111
    110    CALL PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
    111    NA = 1-NA
    GO TO 115
    112    IF (NA .NE. 0) GO TO 113
    CALL PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
    GO TO 114
    113    CALL PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
    114    IF (NAC .NE. 0) NA = 1-NA
    115    L1 = L2
    IW = IW+(IP-1)*IDOT
    116 CONTINUE
    IF (NA .EQ. 0) RETURN
    N2 = N+N
    DO 117 I=1,N2
    C(I) = CH(I)
    117 CONTINUE
    RETURN
    END

    SUBROUTINE CFFTB1 (N,C,CH,WA,IFAC)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
    NF = IFAC(2)
    NA = 0
    L1 = 1
    IW = 1
    DO 116 K1=1,NF
    IP = IFAC(K1+2)
    L2 = IP*L1
    IDO = N/L2
    IDOT = IDO+IDO
    IDL1 = IDOT*L1
    IF (IP .NE. 4) GO TO 103
    IX2 = IW+IDOT
    IX3 = IX2+IDOT
    IF (NA .NE. 0) GO TO 101
    CALL PASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
    GO TO 102
    101    CALL PASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
    102    NA = 1-NA
    GO TO 115
    103    IF (IP .NE. 2) GO TO 106
    IF (NA .NE. 0) GO TO 104
    CALL PASSB2 (IDOT,L1,C,CH,WA(IW))
    GO TO 105
    104    CALL PASSB2 (IDOT,L1,CH,C,WA(IW))
    105    NA = 1-NA
    GO TO 115
    106    IF (IP .NE. 3) GO TO 109
    IX2 = IW+IDOT
    IF (NA .NE. 0) GO TO 107
    CALL PASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
    GO TO 108
    107    CALL PASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
    108    NA = 1-NA
    GO TO 115
    109    IF (IP .NE. 5) GO TO 112
    IX2 = IW+IDOT
    IX3 = IX2+IDOT
    IX4 = IX3+IDOT
    IF (NA .NE. 0) GO TO 110
    CALL PASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
    GO TO 111
    110    CALL PASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
    111    NA = 1-NA
    GO TO 115
    112    IF (NA .NE. 0) GO TO 113
    CALL PASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
    GO TO 114
    113    CALL PASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
    114    IF (NAC .NE. 0) NA = 1-NA
    115    L1 = L2
    IW = IW+(IP-1)*IDOT
    116 CONTINUE
    IF (NA .EQ. 0) RETURN
    N2 = N+N
    DO 117 I=1,N2
    C(I) = CH(I)
    117 CONTINUE
    RETURN
    END

    SUBROUTINE PASSF (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          , &
                    C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP), &
                    CH2(IDL1,IP)
    IDOT = IDO/2
    NT = IP*IDL1
    IPP2 = IP+2
    IPPH = (IP+1)/2
    IDP = IP*IDO
    IF (IDO .LT. L1) GO TO 106
    DO 103 J=2,IPPH
    JC = IPP2-J
    DO 102 K=1,L1
        DO 101 I=1,IDO
            CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
            CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
    101       CONTINUE
    102    CONTINUE
    103 CONTINUE
    DO 105 K=1,L1
    DO 104 I=1,IDO
        CH(I,K,1) = CC(I,1,K)
    104    CONTINUE
    105 CONTINUE
    GO TO 112
    106 DO 109 J=2,IPPH
    JC = IPP2-J
    DO 108 I=1,IDO
        DO 107 K=1,L1
            CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
            CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
    107       CONTINUE
    108    CONTINUE
    109 CONTINUE
    DO 111 I=1,IDO
    DO 110 K=1,L1
        CH(I,K,1) = CC(I,1,K)
    110    CONTINUE
    111 CONTINUE
    112 IDL = 2-IDO
    INC = 0
    DO 116 L=2,IPPH
    LC = IPP2-L
    IDL = IDL+IDO
    DO 113 IK=1,IDL1
        C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
        C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
    113    CONTINUE
    IDLJ = IDL
    INC = INC+IDO
    DO 115 J=3,IPPH
        JC = IPP2-J
        IDLJ = IDLJ+INC
        IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
        WAR = WA(IDLJ-1)
        WAI = WA(IDLJ)
        DO 114 IK=1,IDL1
            C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
            C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
    114       CONTINUE
    115    CONTINUE
    116 CONTINUE
    DO 118 J=2,IPPH
    DO 117 IK=1,IDL1
        CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
    117    CONTINUE
    118 CONTINUE
    DO 120 J=2,IPPH
    JC = IPP2-J
    DO 119 IK=2,IDL1,2
        CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
        CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
        CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
        CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
    119    CONTINUE
    120 CONTINUE
    NAC = 1
    IF (IDO .EQ. 2) RETURN
    NAC = 0
    DO 121 IK=1,IDL1
    C2(IK,1) = CH2(IK,1)
    121 CONTINUE
    DO 123 J=2,IP
    DO 122 K=1,L1
        C1(1,K,J) = CH(1,K,J)
        C1(2,K,J) = CH(2,K,J)
    122    CONTINUE
    123 CONTINUE
    IF (IDOT .GT. L1) GO TO 127
    IDIJ = 0
    DO 126 J=2,IP
    IDIJ = IDIJ+2
    DO 125 I=4,IDO,2
        IDIJ = IDIJ+2
        DO 124 K=1,L1
            C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
            C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
    124       CONTINUE
    125    CONTINUE
    126 CONTINUE
    RETURN
    127 IDJ = 2-IDO
    DO 130 J=2,IP
    IDJ = IDJ+IDO
    DO 129 K=1,L1
        IDIJ = IDJ
        DO 128 I=4,IDO,2
            IDIJ = IDIJ+2
            C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
            C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
    128       CONTINUE
    129    CONTINUE
    130 CONTINUE
    RETURN
    END

    SUBROUTINE PASSF2 (IDO,L1,CC,CH,WA1)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           , &
                    WA1(1)
    IF (IDO .GT. 2) GO TO 102
    DO 101 K=1,L1
    CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
    CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
    CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
    CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
    101 CONTINUE
    RETURN
    102 DO 104 K=1,L1
    DO 103 I=2,IDO,2
        CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
        TR2 = CC(I-1,1,K)-CC(I-1,2,K)
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
        TI2 = CC(I,1,K)-CC(I,2,K)
        CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
        CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
    103    CONTINUE
    104 CONTINUE
    RETURN
    END

    SUBROUTINE PASSF3 (IDO,L1,CC,CH,WA1,WA2)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           , &
                    WA1(1)     ,WA2(1)
    !     *** TAUI IS -SQRT(3)/2 ***
    DATA TAUR,TAUI /-0.5D0,-0.86602540378443864676D0/
    IF (IDO .NE. 2) GO TO 102
    DO 101 K=1,L1
    TR2 = CC(1,2,K)+CC(1,3,K)
    CR2 = CC(1,1,K)+TAUR*TR2
    CH(1,K,1) = CC(1,1,K)+TR2
    TI2 = CC(2,2,K)+CC(2,3,K)
    CI2 = CC(2,1,K)+TAUR*TI2
    CH(2,K,1) = CC(2,1,K)+TI2
    CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
    CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
    CH(1,K,2) = CR2-CI3
    CH(1,K,3) = CR2+CI3
    CH(2,K,2) = CI2+CR3
    CH(2,K,3) = CI2-CR3
    101 CONTINUE
    RETURN
    102 DO 104 K=1,L1
    DO 103 I=2,IDO,2
        TR2 = CC(I-1,2,K)+CC(I-1,3,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,2,K)+CC(I,3,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
        CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
        CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
        CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
        CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
    103    CONTINUE
    104 CONTINUE
    RETURN
    END

    SUBROUTINE PASSF4 (IDO,L1,CC,CH,WA1,WA2,WA3)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           , &
                    WA1(1)     ,WA2(1)     ,WA3(1)
    IF (IDO .NE. 2) GO TO 102
    DO 101 K=1,L1
    TI1 = CC(2,1,K)-CC(2,3,K)
    TI2 = CC(2,1,K)+CC(2,3,K)
    TR4 = CC(2,2,K)-CC(2,4,K)
    TI3 = CC(2,2,K)+CC(2,4,K)
    TR1 = CC(1,1,K)-CC(1,3,K)
    TR2 = CC(1,1,K)+CC(1,3,K)
    TI4 = CC(1,4,K)-CC(1,2,K)
    TR3 = CC(1,2,K)+CC(1,4,K)
    CH(1,K,1) = TR2+TR3
    CH(1,K,3) = TR2-TR3
    CH(2,K,1) = TI2+TI3
    CH(2,K,3) = TI2-TI3
    CH(1,K,2) = TR1+TR4
    CH(1,K,4) = TR1-TR4
    CH(2,K,2) = TI1+TI4
    CH(2,K,4) = TI1-TI4
    101 CONTINUE
    RETURN
    102 DO 104 K=1,L1
    DO 103 I=2,IDO,2
        TI1 = CC(I,1,K)-CC(I,3,K)
        TI2 = CC(I,1,K)+CC(I,3,K)
        TI3 = CC(I,2,K)+CC(I,4,K)
        TR4 = CC(I,2,K)-CC(I,4,K)
        TR1 = CC(I-1,1,K)-CC(I-1,3,K)
        TR2 = CC(I-1,1,K)+CC(I-1,3,K)
        TI4 = CC(I-1,4,K)-CC(I-1,2,K)
        TR3 = CC(I-1,2,K)+CC(I-1,4,K)
        CH(I-1,K,1) = TR2+TR3
        CR3 = TR2-TR3
        CH(I,K,1) = TI2+TI3
        CI3 = TI2-TI3
        CR2 = TR1+TR4
        CR4 = TR1-TR4
        CI2 = TI1+TI4
        CI4 = TI1-TI4
        CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
        CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
        CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
        CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
        CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
        CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
    103    CONTINUE
    104 CONTINUE
    RETURN
    END

    SUBROUTINE PASSF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           , &
                    WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)
    !     *** TR11=COS(2*PI/5), TI11=-SIN(2*PI/5)
    !     *** TR12=-COS(4*PI/5), TI12=-SIN(4*PI/5)
    DATA TR11,TI11,TR12,TI12 /0.3090169943749474241D0, &
        -0.95105651629515357212D0, &
        -0.8090169943749474241D0, -0.58778525229247312917D0/
    IF (IDO .NE. 2) GO TO 102
    DO 101 K=1,L1
    TI5 = CC(2,2,K)-CC(2,5,K)
    TI2 = CC(2,2,K)+CC(2,5,K)
    TI4 = CC(2,3,K)-CC(2,4,K)
    TI3 = CC(2,3,K)+CC(2,4,K)
    TR5 = CC(1,2,K)-CC(1,5,K)
    TR2 = CC(1,2,K)+CC(1,5,K)
    TR4 = CC(1,3,K)-CC(1,4,K)
    TR3 = CC(1,3,K)+CC(1,4,K)
    CH(1,K,1) = CC(1,1,K)+TR2+TR3
    CH(2,K,1) = CC(2,1,K)+TI2+TI3
    CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
    CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
    CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
    CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
    CR5 = TI11*TR5+TI12*TR4
    CI5 = TI11*TI5+TI12*TI4
    CR4 = TI12*TR5-TI11*TR4
    CI4 = TI12*TI5-TI11*TI4
    CH(1,K,2) = CR2-CI5
    CH(1,K,5) = CR2+CI5
    CH(2,K,2) = CI2+CR5
    CH(2,K,3) = CI3+CR4
    CH(1,K,3) = CR3-CI4
    CH(1,K,4) = CR3+CI4
    CH(2,K,4) = CI3-CR4
    CH(2,K,5) = CI2-CR5
    101 CONTINUE
    RETURN
    102 DO 104 K=1,L1
    DO 103 I=2,IDO,2
        TI5 = CC(I,2,K)-CC(I,5,K)
        TI2 = CC(I,2,K)+CC(I,5,K)
        TI4 = CC(I,3,K)-CC(I,4,K)
        TI3 = CC(I,3,K)+CC(I,4,K)
        TR5 = CC(I-1,2,K)-CC(I-1,5,K)
        TR2 = CC(I-1,2,K)+CC(I-1,5,K)
        TR4 = CC(I-1,3,K)-CC(I-1,4,K)
        TR3 = CC(I-1,3,K)+CC(I-1,4,K)
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
        CH(I,K,1) = CC(I,1,K)+TI2+TI3
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
        CR5 = TI11*TR5+TI12*TR4
        CI5 = TI11*TI5+TI12*TI4
        CR4 = TI12*TR5-TI11*TR4
        CI4 = TI12*TI5-TI11*TI4
        DR3 = CR3-CI4
        DR4 = CR3+CI4
        DI3 = CI3+CR4
        DI4 = CI3-CR4
        DR5 = CR2+CI5
        DR2 = CR2-CI5
        DI5 = CI2-CR5
        DI2 = CI2+CR5
        CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
        CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
        CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
        CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
        CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
        CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
        CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
        CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
    103    CONTINUE
    104 CONTINUE
    RETURN
    END

    SUBROUTINE PASSB (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          , &
                    C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP), &
                    CH2(IDL1,IP)
    IDOT = IDO/2
    NT = IP*IDL1
    IPP2 = IP+2
    IPPH = (IP+1)/2
    IDP = IP*IDO
    IF (IDO .LT. L1) GO TO 106
    DO 103 J=2,IPPH
    JC = IPP2-J
    DO 102 K=1,L1
        DO 101 I=1,IDO
            CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
            CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
    101       CONTINUE
    102    CONTINUE
    103 CONTINUE
    DO 105 K=1,L1
    DO 104 I=1,IDO
        CH(I,K,1) = CC(I,1,K)
    104    CONTINUE
    105 CONTINUE
    GO TO 112
    106 DO 109 J=2,IPPH
    JC = IPP2-J
    DO 108 I=1,IDO
        DO 107 K=1,L1
            CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
            CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
    107       CONTINUE
    108    CONTINUE
    109 CONTINUE
    DO 111 I=1,IDO
    DO 110 K=1,L1
        CH(I,K,1) = CC(I,1,K)
    110    CONTINUE
    111 CONTINUE
    112 IDL = 2-IDO
    INC = 0
    DO 116 L=2,IPPH
    LC = IPP2-L
    IDL = IDL+IDO
    DO 113 IK=1,IDL1
        C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
        C2(IK,LC) = WA(IDL)*CH2(IK,IP)
    113    CONTINUE
    IDLJ = IDL
    INC = INC+IDO
    DO 115 J=3,IPPH
        JC = IPP2-J
        IDLJ = IDLJ+INC
        IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
        WAR = WA(IDLJ-1)
        WAI = WA(IDLJ)
        DO 114 IK=1,IDL1
            C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
            C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
    114       CONTINUE
    115    CONTINUE
    116 CONTINUE
    DO 118 J=2,IPPH
    DO 117 IK=1,IDL1
        CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
    117    CONTINUE
    118 CONTINUE
    DO 120 J=2,IPPH
    JC = IPP2-J
    DO 119 IK=2,IDL1,2
        CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
        CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
        CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
        CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
    119    CONTINUE
    120 CONTINUE
    NAC = 1
    IF (IDO .EQ. 2) RETURN
    NAC = 0
    DO 121 IK=1,IDL1
    C2(IK,1) = CH2(IK,1)
    121 CONTINUE
    DO 123 J=2,IP
    DO 122 K=1,L1
        C1(1,K,J) = CH(1,K,J)
        C1(2,K,J) = CH(2,K,J)
    122    CONTINUE
    123 CONTINUE
    IF (IDOT .GT. L1) GO TO 127
    IDIJ = 0
    DO 126 J=2,IP
    IDIJ = IDIJ+2
    DO 125 I=4,IDO,2
        IDIJ = IDIJ+2
        DO 124 K=1,L1
            C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
            C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
    124       CONTINUE
    125    CONTINUE
    126 CONTINUE
    RETURN
    127 IDJ = 2-IDO
    DO 130 J=2,IP
    IDJ = IDJ+IDO
    DO 129 K=1,L1
        IDIJ = IDJ
        DO 128 I=4,IDO,2
            IDIJ = IDIJ+2
            C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
            C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
    128       CONTINUE
    129    CONTINUE
    130 CONTINUE
    RETURN
    END

    SUBROUTINE PASSB2 (IDO,L1,CC,CH,WA1)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           , &
                    WA1(1)
    IF (IDO .GT. 2) GO TO 102
    DO 101 K=1,L1
    CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
    CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
    CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
    CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
    101 CONTINUE
    RETURN
    102 DO 104 K=1,L1
    DO 103 I=2,IDO,2
        CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
        TR2 = CC(I-1,1,K)-CC(I-1,2,K)
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
        TI2 = CC(I,1,K)-CC(I,2,K)
        CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
        CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
    103    CONTINUE
    104 CONTINUE
    RETURN
    END

    SUBROUTINE PASSB3 (IDO,L1,CC,CH,WA1,WA2)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           , &
                    WA1(1)     ,WA2(1)
    !     *** TAUI IS SQRT(3)/2 ***
    DATA TAUR,TAUI /-0.5D0,0.86602540378443864676D0/
    IF (IDO .NE. 2) GO TO 102
    DO 101 K=1,L1
    TR2 = CC(1,2,K)+CC(1,3,K)
    CR2 = CC(1,1,K)+TAUR*TR2
    CH(1,K,1) = CC(1,1,K)+TR2
    TI2 = CC(2,2,K)+CC(2,3,K)
    CI2 = CC(2,1,K)+TAUR*TI2
    CH(2,K,1) = CC(2,1,K)+TI2
    CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
    CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
    CH(1,K,2) = CR2-CI3
    CH(1,K,3) = CR2+CI3
    CH(2,K,2) = CI2+CR3
    CH(2,K,3) = CI2-CR3
    101 CONTINUE
    RETURN
    102 DO 104 K=1,L1
    DO 103 I=2,IDO,2
        TR2 = CC(I-1,2,K)+CC(I-1,3,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,2,K)+CC(I,3,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
        CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
        CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
        CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
        CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
    103    CONTINUE
    104 CONTINUE
    RETURN
    END

    SUBROUTINE PASSB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           , &
                    WA1(1)     ,WA2(1)     ,WA3(1)
    IF (IDO .NE. 2) GO TO 102
    DO 101 K=1,L1
    TI1 = CC(2,1,K)-CC(2,3,K)
    TI2 = CC(2,1,K)+CC(2,3,K)
    TR4 = CC(2,4,K)-CC(2,2,K)
    TI3 = CC(2,2,K)+CC(2,4,K)
    TR1 = CC(1,1,K)-CC(1,3,K)
    TR2 = CC(1,1,K)+CC(1,3,K)
    TI4 = CC(1,2,K)-CC(1,4,K)
    TR3 = CC(1,2,K)+CC(1,4,K)
    CH(1,K,1) = TR2+TR3
    CH(1,K,3) = TR2-TR3
    CH(2,K,1) = TI2+TI3
    CH(2,K,3) = TI2-TI3
    CH(1,K,2) = TR1+TR4
    CH(1,K,4) = TR1-TR4
    CH(2,K,2) = TI1+TI4
    CH(2,K,4) = TI1-TI4
    101 CONTINUE
    RETURN
    102 DO 104 K=1,L1
    DO 103 I=2,IDO,2
        TI1 = CC(I,1,K)-CC(I,3,K)
        TI2 = CC(I,1,K)+CC(I,3,K)
        TI3 = CC(I,2,K)+CC(I,4,K)
        TR4 = CC(I,4,K)-CC(I,2,K)
        TR1 = CC(I-1,1,K)-CC(I-1,3,K)
        TR2 = CC(I-1,1,K)+CC(I-1,3,K)
        TI4 = CC(I-1,2,K)-CC(I-1,4,K)
        TR3 = CC(I-1,2,K)+CC(I-1,4,K)
        CH(I-1,K,1) = TR2+TR3
        CR3 = TR2-TR3
        CH(I,K,1) = TI2+TI3
        CI3 = TI2-TI3
        CR2 = TR1+TR4
        CR4 = TR1-TR4
        CI2 = TI1+TI4
        CI4 = TI1-TI4
        CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
        CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
        CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
        CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
        CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
        CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
    103    CONTINUE
    104 CONTINUE
    RETURN
    END

    SUBROUTINE PASSB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           , &
                    WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)
    !     *** TR11=COS(2*PI/5), TI11=SIN(2*PI/5)
    !     *** TR12=COS(4*PI/5), TI12=SIN(4*PI/5)
    DATA TR11,TI11,TR12,TI12 /0.3090169943749474241D0, &
        0.95105651629515357212D0, &
        -0.8090169943749474241D0,0.58778525229247312917D0/
    IF (IDO .NE. 2) GO TO 102
    DO 101 K=1,L1
    TI5 = CC(2,2,K)-CC(2,5,K)
    TI2 = CC(2,2,K)+CC(2,5,K)
    TI4 = CC(2,3,K)-CC(2,4,K)
    TI3 = CC(2,3,K)+CC(2,4,K)
    TR5 = CC(1,2,K)-CC(1,5,K)
    TR2 = CC(1,2,K)+CC(1,5,K)
    TR4 = CC(1,3,K)-CC(1,4,K)
    TR3 = CC(1,3,K)+CC(1,4,K)
    CH(1,K,1) = CC(1,1,K)+TR2+TR3
    CH(2,K,1) = CC(2,1,K)+TI2+TI3
    CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
    CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
    CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
    CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
    CR5 = TI11*TR5+TI12*TR4
    CI5 = TI11*TI5+TI12*TI4
    CR4 = TI12*TR5-TI11*TR4
    CI4 = TI12*TI5-TI11*TI4
    CH(1,K,2) = CR2-CI5
    CH(1,K,5) = CR2+CI5
    CH(2,K,2) = CI2+CR5
    CH(2,K,3) = CI3+CR4
    CH(1,K,3) = CR3-CI4
    CH(1,K,4) = CR3+CI4
    CH(2,K,4) = CI3-CR4
    CH(2,K,5) = CI2-CR5
    101 CONTINUE
    RETURN
    102 DO 104 K=1,L1
    DO 103 I=2,IDO,2
        TI5 = CC(I,2,K)-CC(I,5,K)
        TI2 = CC(I,2,K)+CC(I,5,K)
        TI4 = CC(I,3,K)-CC(I,4,K)
        TI3 = CC(I,3,K)+CC(I,4,K)
        TR5 = CC(I-1,2,K)-CC(I-1,5,K)
        TR2 = CC(I-1,2,K)+CC(I-1,5,K)
        TR4 = CC(I-1,3,K)-CC(I-1,4,K)
        TR3 = CC(I-1,3,K)+CC(I-1,4,K)
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
        CH(I,K,1) = CC(I,1,K)+TI2+TI3
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
        CR5 = TI11*TR5+TI12*TR4
        CI5 = TI11*TI5+TI12*TI4
        CR4 = TI12*TR5-TI11*TR4
        CI4 = TI12*TI5-TI11*TI4
        DR3 = CR3-CI4
        DR4 = CR3+CI4
        DI3 = CI3+CR4
        DI4 = CI3-CR4
        DR5 = CR2+CI5
        DR2 = CR2-CI5
        DI5 = CI2-CR5
        DI2 = CI2+CR5
        CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
        CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
        CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
        CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
        CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
        CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
        CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
        CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
    103    CONTINUE
    104 CONTINUE
    RETURN
    END
