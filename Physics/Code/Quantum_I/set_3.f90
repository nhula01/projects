program main
    implicit none
    double precision, parameter :: pi=3.141592653589793
    double complex, parameter :: Im=(0.0d0,1.0d0) !imaginary unit
    double precision, parameter :: hbar=1.0 !atomic units

    double precision, parameter :: dt=5.0d-3 !time step for split-operator, need to be small but no super small
            !gotta experiment to find proper value
            !rule of thumb 1/dt (in a.u.) is energy - should be large enough to cover potential well
    double precision, parameter :: tol=1.0d-6 !relative error to stop iterations

    integer, parameter :: N=2048 !total number of points - for fast results make it 2**n
    double precision, parameter :: a=-50.0,b=50.0 !coordinate interval (must be greater than potental width!)
    double precision, parameter :: Bound1=-400.0,Bound2=400.0 !coordinate interval (must be greater than potental width!)
    double precision, parameter :: mass=1.0 !1 in atomic units means it is electron
    double precision, parameter :: dx=(Bound2-Bound1)/N
    double precision x
    character(len=100) :: name, file_name
    double complex psi(N) !wavefunction
    double complex psi_tmp(N),tmp_c !for reusing
    double precision tmp,exp_V(N),exp_T(N),Vpart,Tpart,En1,En2,eps,potential,p, potential1
    double precision wsave(4*N+15) !working array for FFT
    integer i,j
    character(len=100) :: name1="_2.dat", name2="_4.dat", name3="_1.dat", name4="_1imp2.dat"

    file_name = '/Users/phihung/NumMethod/first/homework/hw_07/wells1'//name4
    open(file=file_name,unit=10)
        do i=1,N-1
            x=Bound1+dx*float(i)
            write(10,*) x,potential(x)
        enddo
    close(unit=10)
    call k(name4)
end program main

subroutine k(name)
    implicit none
    double precision, parameter :: pi=3.141592653589793
    double complex, parameter :: Im=(0.0d0,1.0d0) !imaginary unit
    double precision, parameter :: hbar=1.0 !atomic units
    double precision, parameter :: dt=1.0d-2 !time step for split-operator, need to be small but not super small
            !gotta experiment to find proper value
            !rule of thumb 1/dt (in a.u.) is energy - should be large enough to cover potential well
    double precision, parameter :: tol=1.0d-6 !relative error to stop iterations
    integer, parameter :: N=2*8192 !total number of points - for the fastest results make it 2**n
    double precision, parameter :: a=-350.0,b=350.0 !coordinate interval
    double precision, parameter :: mass=1.0 !1 in atomic units means it is electron
    double precision, parameter :: dx=(b-a)/N
    double precision x,cc(N)

    double complex psi(N) !wavefunction
    character(len=100) :: name, file_name
    double complex psi_tmp(N),tmp_c !for reusing
    double precision tmp,exp_V(N),exp_T(N),Vpart,Tpart,En1,En2,eps,potential,p
    double precision wsave(4*N+15) !working array for FFT
    double complex psi0(N) !save the ground state wavefunction
    double precision E0 !save the ground state energy
    integer i,j

    !initial wavefunction, could be const but never 0
    psi=(1.0,0.0)
    !<-- normalize to 1
    tmp=0.0
    do i=1,N
    tmp=tmp+abs(psi(i))**2
    enddo
    tmp=tmp*dx
    psi=psi/sqrt(tmp) !now it is propely normalized

    !<-- precalculate exponents for split-operator propagation
    !exp(-i*T*dt/(2*hbar)) * exp(-i*V*dt/hbar) * exp(-i*T*dt/(2*hbar))
    do i=1,N
    x=a+dx*float(i)
    exp_V(i)=exp(-potential(x)*dt/hbar)

    if(i<=(N/2))then !<-- momentum in Fourier space
        p=hbar*2.0*pi*(i-1)/(dx*N)
        else
        p=-hbar*2.0*pi*(N+1-i)/(dx*N)
    endif
    exp_T(i)=exp(-p**2*dt/(4.0*mass*hbar))
    enddo

    !<---initialize fftpack -- call it only once
    call zffti(N,wsave)

    !<-- calculate energy first
    !--- potential part ---!
    tmp=0.0
    do i=1,N
    x=a+dx*float(i)
    tmp=tmp+(abs(psi(i))**2)*potential(x)
    enddo
    Vpart=tmp*dx

    !--- kinetic part ---!
    psi_tmp=psi
    call zfftf(N,psi_tmp,wsave) !forward Fourier transform

    do i=1,N
    if(i<=(N/2))then !<-- momentum in Fourier space
        p=hbar*2.0*pi*(i-1)/(dx*N)
        else
        p=-hbar*2.0*pi*(N+1-i)/(dx*N)
    endif
    psi_tmp(i)=psi_tmp(i)*(p**2/(2.0*mass))
    enddo
    !now go back to coordinate representation
    call zfftb(N,psi_tmp,wsave) !backward Fourier transform

    tmp_c=(0.0,0.0)
    do i=1,N
    tmp_c=tmp_c+conjg(psi(i))*psi_tmp(i)
    enddo

    Tpart=dx*dreal(tmp_c)/float(N) !1/N is needed due to forth and back fft

    En1=Vpart+Tpart
    write(*,*) En1

    write(*,*) '******************************'
    write(*,*) 'first we find the ground state'
    write(*,*) '******************************'

    eps=1.0
    j=0
    do while(eps>tol)
    j=j+1
    psi_tmp=psi
        !--- kinetic part ---!
        call zfftf(N,psi_tmp,wsave) !forward Fourier transform
        psi_tmp=exp_T*psi_tmp/float(N) !1/N is needed due to forth and back fft
        call zfftb(N,psi_tmp,wsave) !backward Fourier transform

        !--- potential part ---!
        psi_tmp=exp_V*psi_tmp

        !--- kinetic part ---!
        call zfftf(N,psi_tmp,wsave) !forward Fourier transform
        psi_tmp=exp_T*psi_tmp/float(N) !1/N is needed due to forth and back fft
        call zfftb(N,psi_tmp,wsave) !backward Fourier transform

        !<-- normalize to 1
        tmp=0.0
        do i=1,N
        tmp=tmp+abs(psi_tmp(i))**2
        enddo
        tmp=tmp*dx
        psi_tmp=psi_tmp/sqrt(tmp) !now it is propely normalized

        psi=psi_tmp !save both psi for energy calculation below

        !<-- calculate energy
        !--- potential part ---!
        tmp=0.0
        do i=1,N
        x=a+dx*float(i)
        tmp=tmp+(abs(psi(i))**2)*potential(x)
        enddo
        Vpart=tmp*dx
        !--- kinetic part ---!
        call zfftf(N,psi_tmp,wsave) !forward Fourier transform
        do i=1,N
        if(i<=(N/2))then !<-- momentum in Fourier space
            p=hbar*2.0*pi*(i-1)/(dx*N)
            else
            p=-hbar*2.0*pi*(N+1-i)/(dx*N)
        endif
        psi_tmp(i)=psi_tmp(i)*(p**2/(2.0*mass))
        enddo
        !now go back to coordinate representation
        call zfftb(N,psi_tmp,wsave) !backward Fourier transform
        tmp_c=(0.0,0.0)
        do i=1,N
        tmp_c=tmp_c+conjg(psi(i))*psi_tmp(i)
        enddo
        Tpart=dx*dreal(tmp_c)/float(N) !1/N is needed due to fourth and back fft

        En2=Tpart+Vpart
        eps=abs((En2-En1)/En1)
        En1=En2

        if(mod(j,100)==0)then !write energy every 100 time steps
            write(*,*) j,En1,eps
        endif
    enddo

    psi0=psi  !ground state wavefunction
    E0=En1    !ground state energy

    file_name ='/Users/phihung/NumMethod/first/homework/hw_07/converged_energy_E0'//name
    open(file=file_name,unit=11)
        write(11,*) j,E0,eps
    close(unit=11)
    file_name='/Users/phihung/NumMethod/first/homework/hw_07/psi0'//name
    open(file=file_name,unit=11)
        do i=1,N
        x=a+dx*float(i)
        write(11,*) x,dreal(psi0(i))
    enddo
    close(unit=11)








    write(*,*) '*****************************************************'
    write(*,*) 'the ground state is done, now the first excited state'
    write(*,*) '*****************************************************'

    !initial wavefunction, use Gaussian with sin to have an antisymmetric guess
    !you may use any antisymmetric guess with the node at 0
    !other guesses such as simple Gaussian may not converge
    do i=1,N
    x=a+dx*float(i)
    psi(i)=x !antisymmetric
    enddo

    !now we need to remove psi0 from the guess
    !the new guess should be normalized AND orthogonal to psi0
    tmp_c=(0.0,0.0) !calculate overlap <0|guess>
    do i=1,N
    tmp_c=tmp_c+psi(i)*conjg(psi0(i))
    enddo

    psi=psi-tmp_c*dx*psi0 !this guess is orthogonal to |0>
    !<-- now we normalize it to 1
    tmp=0.0
    do i=1,N
    tmp=tmp+abs(psi(i))**2
    enddo
    tmp=tmp*dx
    psi=psi/sqrt(tmp) !now it is propely normalized
    write(*,*) 'normalized initial guess'

    !<-- calculate energy first
    !--- potential part ---!
    tmp=0.0
    do i=1,N
    x=a+dx*float(i)
    tmp=tmp+(abs(psi(i))**2)*potential(x)
    enddo
    Vpart=tmp*dx

    !--- kinetic part ---!
    psi_tmp=psi
    call zfftf(N,psi_tmp,wsave) !forward Fourier transform

    do i=1,N
    if(i<=(N/2))then !<-- momentum in Fourier space
        p=hbar*2.0*pi*(i-1)/(dx*N)
        else
        p=-hbar*2.0*pi*(N+1-i)/(dx*N)
    endif
    psi_tmp(i)=psi_tmp(i)*(p**2/(2.0*mass))
    enddo
    !now go back to coordinate representation
    call zfftb(N,psi_tmp,wsave) !backward Fourier transform

    tmp_c=(0.0,0.0)
    do i=1,N
    tmp_c=tmp_c+conjg(psi(i))*psi_tmp(i)
    enddo

    Tpart=dx*dreal(tmp_c)/float(N) !1/N is needed due to forth and back fft

    En1=Vpart+Tpart
    write(*,*) 'first excited state energy, initial guess',En1
    write(*,*)
    write(*,*) '************************************'
    write(*,*) 'now we try to find the excited state'
    write(*,*) '************************************'
    write(*,*)
    eps=1.0
    j=0
    do while(eps>tol)
    j=j+1
    psi_tmp=psi
        !--- kinetic part ---!
        call zfftf(N,psi_tmp,wsave) !forward Fourier transform
        psi_tmp=exp_T*psi_tmp/float(N) !1/N is needed due to forth and back fft
        call zfftb(N,psi_tmp,wsave) !backward Fourier transform

        !--- potential part ---!
        psi_tmp=exp_V*psi_tmp

        !--- kinetic part ---!
        call zfftf(N,psi_tmp,wsave) !forward Fourier transform
        psi_tmp=exp_T*psi_tmp/float(N) !1/N is needed due to forth and back fft
        call zfftb(N,psi_tmp,wsave) !backward Fourier transform

        !remove the ground state & normalize <-------------------------------
        tmp_c=(0.0,0.0) !calculate overlap <0|guess>
        do i=1,N
        tmp_c=tmp_c+psi_tmp(i)*conjg(psi0(i))
        enddo
        psi_tmp=psi_tmp-tmp_c*dx*psi0 !this wf is orthogonal to |0>
        !<-- normalize to 1
        tmp=0.0
        do i=1,N
        tmp=tmp+abs(psi_tmp(i))**2
        enddo
        tmp=tmp*dx
        psi_tmp=psi_tmp/sqrt(tmp) !now it is propely normalized

        psi=psi_tmp !save both psi for energy calculation below

        !<-- calculate energy
        !--- potential part ---!
        tmp=0.0
        do i=1,N
        x=a+dx*float(i)
        tmp=tmp+(abs(psi(i))**2)*potential(x)
        enddo
        Vpart=tmp*dx
        !--- kinetic part ---!
        call zfftf(N,psi_tmp,wsave) !forward Fourier transform
        do i=1,N
        if(i<=(N/2))then !<-- momentum in Fourier space
            p=hbar*2.0*pi*(i-1)/(dx*N)
            else
            p=-hbar*2.0*pi*(N+1-i)/(dx*N)
        endif
        psi_tmp(i)=psi_tmp(i)*(p**2/(2.0*mass))
        enddo
        !now go back to coordinate representation
        call zfftb(N,psi_tmp,wsave) !backward Fourier transform
        tmp_c=(0.0,0.0)
        do i=1,N
        tmp_c=tmp_c+conjg(psi(i))*psi_tmp(i)
        enddo
        Tpart=dx*dreal(tmp_c)/float(N) !1/N is needed due to fourth and back fft

        En2=Tpart+Vpart
        eps=abs((En2-En1)/En1)
        En1=En2

        if(mod(j,100)==0)then !write energy every 100 time steps
            write(*,*) j,En1,eps
        endif
    enddo
    file_name="/Users/phihung/NumMethod/first/homework/hw_07/converged_energy_E1"//name
    open(file=file_name,unit=11)
        write(11,*) j,En1,eps
    close(unit=11)
    file_name='/Users/phihung/NumMethod/first/homework/hw_07/psi1'//name
    open(file=file_name,unit=11)
    do i=1,N
        x=a+dx*float(i)
        write(11,*) x,dreal(psi(i))
    enddo
    close(unit=11)
end subroutine k
double precision function potential(x)
  implicit none
  integer :: impurity
  double precision, parameter :: distance=10.0d0 !distance between wells
  double precision, parameter :: width=10.0d0 !width of each well
  double precision, parameter :: V0=0.5 !depth of each well
  double precision x, magnitude
  if((x>=(0.5*distance)).and.(x<=(0.5*distance+width)))then
      potential=0.0
    elseif((x>=(0.5*distance)+1.0*(width+distance)).and.(x<=(0.5*distance+width+1.0*(width+distance))))then
      potential=0.2
    elseif((x>=(0.5*distance)+2.0*(width+distance)).and.(x<=(0.5*distance+width+2.0*(width+distance))))then
      potential=0.0
    elseif((x>=(0.5*distance)+3.0*(width+distance)).and.(x<=(0.5*distance+width+3.0*(width+distance))))then
      potential=0.0
    elseif((x>=(0.5*distance)+4.0*(width+distance)).and.(x<=(0.5*distance+width+4.0*(width+distance))))then
      potential=0.0
    elseif((x>=(0.5*distance)-1.0*(width+distance)).and.(x<=(0.5*distance+width-1.0*(width+distance))))then
      potential=0.0
    elseif((x>=(0.5*distance)-2.0*(width+distance)).and.(x<=(0.5*distance+width-2.0*(width+distance))))then
      potential=0.0
    elseif((x>=(0.5*distance)-3.0*(width+distance)).and.(x<=(0.5*distance+width-3.0*(width+distance))))then
      potential=0.0
    elseif((x>=(0.5*distance)-4.0*(width+distance)).and.(x<=(0.5*distance+width-4.0*(width+distance))))then
      potential=0.0
    elseif((x>=(0.5*distance)-5.0*(width+distance)).and.(x<=(0.5*distance+width-5.0*(width+distance))))then
      potential=0.0
    else
      potential=V0
  endif
end function potential

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
