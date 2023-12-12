subroutine DVR(a,b, name)
    implicit none
    double precision, parameter :: pi=3.141592653589793
    double complex, parameter :: Im=(0.0d0,1.0d0) !imaginary unit
    double precision, parameter :: hbar=1.0, omg = 0.1 !atomic units
    integer, parameter :: N=1025 !total number of points
    double precision, intent(in) :: a,b!=-50.0,b=50.0 !coordinate interval (must be greater than potental width!)
    double precision, parameter :: mass=1.0 !1 in atomic units means it is electron
    double precision :: x, dx
    double precision H(N-1,N-1) !DVR representation of the Hamiltonian
    double precision T,V,potential, exact_eng
    double precision En(N-1)      !eigenenergies
    double precision psi(N-1,N-1) !eignefunctions
    double precision const,tmp,tmp1,tmp2
    character(len=100) :: name, file_name
    double precision fv1(N-1),fv2(N-1) !working arrays for rs subroutine
    integer i,j,ierr,matz

    !find dx
    dx=(b-a)/(N-1)
    !find constant before the T matrix
    const=(hbar*pi)**2/(4.0*mass*(b-a)**2)

    !Find hamiltonian
    H=0.0
    do i=1,N-1
        x=a+dx*float(i-1)
        do j=1,N-1
            if(i==j)then
                T=const*((2.0*(N**2)+1.0)/3.0-1.0/sin(pi*float(i)/float(N))**2)
                V=potential(x)
            else
                T=const*((-1)**(i-j))*(1.0/sin(pi*float(i-j)/(2.0*float(N)))**2- &
                                    1.0/sin(pi*float(i+j)/(2.0*float(N)))**2)
                V=0.0
            endif
            H(i,j)=T+V !Hamiltonian
        enddo
    enddo

    !H is ready
    !call subroutine RS to diagonalize H
    matz=1 !set this to 0 if ONLY eigenenergies needed
    call rs(N-1,N-1,H,En,matz,psi,fv1,fv2,ierr) !
    write(*,*) ierr !if ierr=0 => calculations have no errors, if not something went wrong

    !recored eigenenergies stored in En
    file_name='/Users/phihung/NumMethod/first/homework/hw_07/energies' // name
    open(file=file_name,unit=11)
        do i=1, 50 !print first 5 energies
            exact_eng = ((.5)+i-1)*omg*hbar
            write(11,*) i-1,En(i), exact_eng, (En(i)-exact_eng) !relative error
        enddo
    close(unit=11)

    !normalize all eignefunctions properly
    !note that for psi(i,j) i is coordinate and j is quantum number
    do i=1,N-1
        tmp=0.0
        do j=1,N-1
            tmp=tmp+abs(psi(j,i))**2
        enddo
        tmp=tmp*dx
        psi(:,i)=psi(:,i)/sqrt(tmp) !now it is properly normalized
    enddo

    !record the first three energies
    file_name='/Users/phihung/NumMethod/first/homework/hw_07/states123' // name
    open(file=file_name,unit=11)
        do i=1,N-1
            x=a+dx*float(i-1)
            write(11,*) x,psi(i,1),psi(i,2), psi(i,3) !first two states
        enddo
    close(unit=11)

end subroutine DVR

subroutine time_prop(a,b,name)
    implicit none
    double precision, parameter :: pi=3.141592653589793
    double complex, parameter :: Im=(0.0d0,1.0d0) !imaginary unit
    double precision, parameter :: hbar=1.0 !atomic units
    double precision :: dt=.5 !time step for split-operator, need to be small but not super small
            !gotta experiment to find proper value
            !rule of thumb 1/dt (in a.u.) is energy - should be large enough to cover potential well
    double precision, parameter :: tol=1.0d-4 !relative error to stop iterations
    integer, parameter :: N=2048 !total number of points - for the fastest results make it 2**n
    double precision, intent(in) :: a,b
    double precision, parameter :: mass=1.0 !1 in atomic units means it is electron
    double precision :: dx!, parameter :: dx=(b-a)/N
    double precision x
    character(len=100) :: name, file_name
    double complex psi(N) !wavefunction
    double complex psi_tmp(N),tmp_c !for reusing
    double precision tmp,exp_V(N),exp_T(N),Vpart,Tpart,En1,En2,eps,potential,p
    double precision wsave(4*N+15) !working array for FFT
    integer i,j

    dx = (b-a)/N
    !initial wavefunction, could be const but never 0
    psi=(1.0,0.0)
    !<-- normalize to 1
    tmp=0.0
    do i=1,N
        tmp=tmp+abs(psi(i))**2
    enddo
    tmp=tmp*dx
    psi=psi/sqrt(tmp) !now it is properly normalized

    !<-- precalculate the exponents for split-operator propagation
    !exp(-i*T*dt/(2*hbar)) * exp(-i*V*dt/hbar) * exp(-i*T*dt/(2*hbar))
    do i=1,N
        x=a+dx*float(i)
        exp_V(i)=exp(-potential(x)*dt/hbar)
        if(i<=(N/2))then !<-- momentum in Fourier space
            p=hbar*2.0*pi*(i-1)/(dx*N) ! use uncertainty principle
        else
            p=-hbar*2.0*pi*(N+1-i)/(dx*N)
        endif
        exp_T(i)=exp(-p**2*dt/(4.0*mass*hbar))
    enddo

    !<---initialize fftpack -- call it only once
    call zffti(N,wsave)

    !<-- calculate energy first
    !--- potential part ---!
    ! average value for the first guess
    tmp=0.0
    do i=1,N
        x=a+dx*float(i)
        tmp=tmp+(abs(psi(i))**2)*potential(x)
    enddo
    Vpart=tmp*dx

    !--- kinetic part ---!
    psi_tmp=psi
    call zfftf(N,psi_tmp,wsave) !forward Fourier transform to momentum (y-value)
    do i=1,N
        if(i<=(N/2))then !<-- momentum in Fourier space (x-value)
            p=hbar*2.0*pi*(i-1)/(dx*N)
            else
            p=-hbar*2.0*pi*(N+1-i)/(dx*N)
        endif
        psi_tmp(i)=psi_tmp(i)*(p**2/(2.0*mass)) 
    enddo
    !now go back to coordinate representation
    call zfftb(N,psi_tmp,wsave) !backward Fourier transform
    tmp_c=(0.0,0.0)
    do i=1,N !find the expected value in integral
        tmp_c=tmp_c+conjg(psi(i))*psi_tmp(i)
    enddo
    Tpart=dx*dreal(tmp_c)/float(N) !1/N is needed due to forth and back fft

    En1=Vpart+Tpart
    write(*,*) En1

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
    !    write(*,*) j,En1,eps
    enddo
    
    file_name = "/Users/phihung/NumMethod/first/homework/hw_07/converged_energy" // name
    open(file=file_name,unit=11)
        write(11,*) j,En1,eps
    close(unit=11)

    file_name = "/Users/phihung/NumMethod/first/homework/hw_07/ground_state" // name
    open(file=file_name,unit=11)
        do i=1,N
            x=a+dx*float(i)
            write(11,*) x,dreal(psi(i))
        enddo
    close(unit=11)
end subroutine time_prop

subroutine getDVR
    implicit none
    double precision :: a=-50, b=50, a3=-5, b3=5
    double precision :: a2=-10, b2=10, a4=-2, b4=2
    character(len=100) :: name = "50.dat", name2 = "10.dat", name3 = "5.dat", name4 = "2.dat"
    call DVR(a,b,name)
    call DVR(a2,b2,name2)
    call DVR(a3,b3,name3)
    call DVR(a4,b4,name4)
end subroutine getDVR
program set1
    implicit none
    double precision :: a=-50, b=50, a3=-10, b3=10
    double precision :: a2=-20, b2=20, a4=-5, b4=5
    character(len=100) :: name = "505.dat", name2 = "205.dat", name3 = "105.dat", name4 = "55.dat"
    !call time_prop(a,b,name)
    !call time_prop(a2,b2,name2)
    !call time_prop(a3,b3,name3)
    !call time_prop(a4,b4,name4)

    call getDVR
end program set1


double precision function potential(x)
    implicit none
    double precision :: x, mass=1
    double precision, parameter :: omg = 0.1! automic unit
    potential = (.5)*mass*((omg)**2)*(x**2)
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



subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
    integer n,nm,ierr,matz
    double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
    !
    !     this subroutine calls the recommended sequence of
    !     subroutines from the eigensystem subroutine package (eispack)
    !     to find the eigenvalues and eigenvectors (if desired)
    !     of a real symmetric matrix.
    !
    !     on input
    !
    !        nm  must be set to the row dimension of the two-dimensional
    !        array parameters as declared in the calling program
    !        dimension statement.
    !
    !        n  is the order of the matrix  a.
    !
    !        a  contains the real symmetric matrix.
    !
    !        matz  is an integer variable set equal to zero if
    !        only eigenvalues are desired.  otherwise it is set to
    !        any non-zero integer for both eigenvalues and eigenvectors.
    !
    !     on output
    !
    !        w  contains the eigenvalues in ascending order.
    !
    !        z  contains the eigenvectors if matz is not zero.
    !
    !        ierr  is an integer output variable set equal to an error
    !           completion code described in the documentation for tqlrat
    !           and tql2.  the normal completion code is zero.
    !
    !        fv1  and  fv2  are temporary storage arrays.
    !
    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory
    !
    !     this version dated august 1983.
    !
    !     ------------------------------------------------------------------
    !
    if (n .le. nm) go to 10
    ierr = 10 * n
    go to 50
    !
    10 if (matz .ne. 0) go to 20
    !     .......... find eigenvalues only ..........
    call  tred1(nm,n,a,w,fv1,fv2)
    !  tqlrat encounters catastrophic underflow on the Vax
    !     call  tqlrat(n,w,fv2,ierr)
    call  tql1(n,w,fv1,ierr)
    go to 50
    !    .......... find both eigenvalues and eigenvectors ..........
    20 call  tred2(nm,n,a,w,fv1,z)
    call  tql2(nm,n,w,fv1,z,ierr)
    50 return
    end

    double precision function pythag(a,b)
    double precision a,b
    !
    !     finds dsqrt(a**2+b**2) without overflow or destructive underflow
    !
    double precision p,r,s,t,u
    p = dmax1(dabs(a),dabs(b))
    if (p .eq. 0.0d0) go to 20
    r = (dmin1(dabs(a),dabs(b))/p)**2
    10 continue
    t = 4.0d0 + r
    if (t .eq. 4.0d0) go to 20
    s = r/t
    u = 1.0d0 + 2.0d0*s
    p = u*p
    r = (s/u)**2 * r
    go to 10
    20 pythag = p
    return
    end

    subroutine tql1(n,d,e,ierr)
    integer i,j,l,m,n,ii,l1,l2,mml,ierr
    double precision d(n),e(n)
    double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
    ierr = 0
    if (n .eq. 1) go to 1001
    !
    do 100 i = 2, n
    100 e(i-1) = e(i)
    !
    f = 0.0d0
    tst1 = 0.0d0
    e(n) = 0.0d0
    !
    do 290 l = 1, n
    j = 0
    h = dabs(d(l)) + dabs(e(l))
    if (tst1 .lt. h) tst1 = h
    !     .......... look for small sub-diagonal element ..........
    do 110 m = l, n
        tst2 = tst1 + dabs(e(m))
        if (tst2 .eq. tst1) go to 120
    !     .......... e(n) is always zero, so there is no exit
    !                through the bottom of the loop ..........
    110    continue
    !
    120    if (m .eq. l) go to 210
    130    if (j .eq. 30) go to 1000
    j = j + 1
    !     .......... form shift ..........
    l1 = l + 1
    l2 = l1 + 1
    g = d(l)
    p = (d(l1) - g) / (2.0d0 * e(l))
    r = pythag(p,1.0d0)
    d(l) = e(l) / (p + dsign(r,p))
    d(l1) = e(l) * (p + dsign(r,p))
    dl1 = d(l1)
    h = g - d(l)
    if (l2 .gt. n) go to 145
    !
    do 140 i = l2, n
    140    d(i) = d(i) - h
    !
    145    f = f + h
    !     .......... ql transformation ..........
    p = d(m)
    c = 1.0d0
    c2 = c
    el1 = e(l1)
    s = 0.0d0
    mml = m - l
    !    .......... for i=m-1 step -1 until l do -- ..........
    do 200 ii = 1, mml
        c3 = c2
        c2 = c
        s2 = s
        i = m - ii
        g = c * e(i)
        h = c * p
        r = pythag(p,e(i))
        e(i+1) = s * r
        s = e(i) / r
        c = p / r
        p = c * d(i) - s * g
        d(i+1) = h + s * (c * g + s * d(i))
    200    continue
    !
    p = -s * s2 * c3 * el1 * e(l) / dl1
    e(l) = s * p
    d(l) = c * p
    tst2 = tst1 + dabs(e(l))
    if (tst2 .gt. tst1) go to 130
    210    p = d(l) + f
    !     .......... order eigenvalues ..........
    if (l .eq. 1) go to 250
    !     .......... for i=l step -1 until 2 do -- ..........
    do 230 ii = 2, l
        i = l + 2 - ii
        if (p .ge. d(i-1)) go to 270
        d(i) = d(i-1)
    230    continue
    !
    250    i = 1
    270    d(i) = p
    290 continue
    !
    go to 1001
    !     .......... set error -- no convergence to an
    !                eigenvalue after 30 iterations ..........
    1000 ierr = l
    1001 return
    end

    subroutine tql2(nm,n,d,e,z,ierr)
    integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
    double precision d(n),e(n),z(nm,n)
    double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
    ierr = 0
    if (n .eq. 1) go to 1001
    do 100 i = 2, n
    100 e(i-1) = e(i)
    f = 0.0d0
    tst1 = 0.0d0
    e(n) = 0.0d0
    do 240 l = 1, n
    j = 0
    h = dabs(d(l)) + dabs(e(l))
    if (tst1 .lt. h) tst1 = h
    do 110 m = l, n
        tst2 = tst1 + dabs(e(m))
        if (tst2 .eq. tst1) go to 120
    110    continue
    120    if (m .eq. l) go to 220
    130    if (j .eq. 30) go to 1000
    j = j + 1
    l1 = l + 1
    l2 = l1 + 1
    g = d(l)
    p = (d(l1) - g) / (2.0d0 * e(l))
    r = pythag(p,1.0d0)
    d(l) = e(l) / (p + dsign(r,p))
    d(l1) = e(l) * (p + dsign(r,p))
    dl1 = d(l1)
    h = g - d(l)
    if (l2 .gt. n) go to 145
    do 140 i = l2, n
    140    d(i) = d(i) - h
    145    f = f + h
    p = d(m)
    c = 1.0d0
    c2 = c
    el1 = e(l1)
    s = 0.0d0
    mml = m - l
    do 200 ii = 1, mml
        c3 = c2
        c2 = c
        s2 = s
        i = m - ii
        g = c * e(i)
        h = c * p
        r = pythag(p,e(i))
        e(i+1) = s * r
        s = e(i) / r
        c = p / r
        p = c * d(i) - s * g
        d(i+1) = h + s * (c * g + s * d(i))
        do 180 k = 1, n
            h = z(k,i+1)
            z(k,i+1) = s * z(k,i) + c * h
            z(k,i) = c * z(k,i) - s * h
    180       continue
    200    continue
    p = -s * s2 * c3 * el1 * e(l) / dl1
    e(l) = s * p
    d(l) = c * p
    tst2 = tst1 + dabs(e(l))
    if (tst2 .gt. tst1) go to 130
    220    d(l) = d(l) + f
    240 continue
    do 300 ii = 2, n
    i = ii - 1
    k = i
    p = d(i)
    do 260 j = ii, n
        if (d(j) .ge. p) go to 260
        k = j
        p = d(j)
    260    continue
    if (k .eq. i) go to 300
    d(k) = d(i)
    d(i) = p
    do 280 j = 1, n
        p = z(j,i)
        z(j,i) = z(j,k)
        z(j,k) = p
    280    continue
    300 continue
    go to 1001
    1000 ierr = l
    1001 return
    end

    subroutine tred1(nm,n,a,d,e,e2)
    integer i,j,k,l,n,ii,nm,jp1
    double precision a(nm,n),d(n),e(n),e2(n)
    double precision f,g,h,scale
    do 100 i = 1, n
    d(i) = a(n,i)
    a(n,i) = a(i,i)
    100 continue
    do 300 ii = 1, n
    i = n + 1 - ii
    l = i - 1
    h = 0.0d0
    scale = 0.0d0
    if (l .lt. 1) go to 130
    do 120 k = 1, l
    120    scale = scale + dabs(d(k))
    if (scale .ne. 0.0d0) go to 140
    do 125 j = 1, l
        d(j) = a(l,j)
        a(l,j) = a(i,j)
        a(i,j) = 0.0d0
    125    continue
    130    e(i) = 0.0d0
    e2(i) = 0.0d0
    go to 300
    140    do 150 k = 1, l
        d(k) = d(k) / scale
        h = h + d(k) * d(k)
    150    continue
    e2(i) = scale * scale * h
    f = d(l)
    g = -dsign(dsqrt(h),f)
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g
    if (l .eq. 1) go to 285
    do 170 j = 1, l
    170    e(j) = 0.0d0
    do 240 j = 1, l
        f = d(j)
        g = e(j) + a(j,j) * f
        jp1 = j + 1
        if (l .lt. jp1) go to 220
        do 200 k = jp1, l
            g = g + a(k,j) * d(k)
            e(k) = e(k) + a(k,j) * f
    200       continue
    220       e(j) = g
    240    continue
    f = 0.0d0
    do 245 j = 1, l
        e(j) = e(j) / h
        f = f + e(j) * d(j)
    245    continue
    h = f / (h + h)
    do 250 j = 1, l
    250    e(j) = e(j) - h * d(j)
    do 280 j = 1, l
        f = d(j)
        g = e(j)
        do 260 k = j, l
    260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
    280    continue
    285    do 290 j = 1, l
        f = d(j)
        d(j) = a(l,j)
        a(l,j) = a(i,j)
        a(i,j) = f * scale
    290    continue
    300 continue
    return
    end

    subroutine tred2(nm,n,a,d,e,z)
    integer i,j,k,l,n,ii,nm,jp1
    double precision a(nm,n),d(n),e(n),z(nm,n)
    double precision f,g,h,hh,scale
    do 100 i = 1, n
    do 80 j = i, n
    80    z(j,i) = a(j,i)
    d(i) = a(n,i)
    100 continue
    if (n .eq. 1) go to 510
    do 300 ii = 2, n
    i = n + 2 - ii
    l = i - 1
    h = 0.0d0
    scale = 0.0d0
    if (l .lt. 2) go to 130
    do 120 k = 1, l
    120    scale = scale + dabs(d(k))
    if (scale .ne. 0.0d0) go to 140
    130    e(i) = d(l)
    do 135 j = 1, l
        d(j) = z(l,j)
        z(i,j) = 0.0d0
        z(j,i) = 0.0d0
    135    continue
    go to 290
    140    do 150 k = 1, l
        d(k) = d(k) / scale
        h = h + d(k) * d(k)
    150    continue
    f = d(l)
    g = -dsign(dsqrt(h),f)
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g
    do 170 j = 1, l
    170    e(j) = 0.0d0
    do 240 j = 1, l
        f = d(j)
        z(j,i) = f
        g = e(j) + z(j,j) * f
        jp1 = j + 1
        if (l .lt. jp1) go to 220
        do 200 k = jp1, l
            g = g + z(k,j) * d(k)
            e(k) = e(k) + z(k,j) * f
    200       continue
    220       e(j) = g
    240    continue
    f = 0.0d0
    do 245 j = 1, l
        e(j) = e(j) / h
        f = f + e(j) * d(j)
    245    continue
    hh = f / (h + h)
    do 250 j = 1, l
    250    e(j) = e(j) - hh * d(j)
    do 280 j = 1, l
        f = d(j)
        g = e(j)
        do 260 k = j, l
    260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
        d(j) = z(l,j)
        z(i,j) = 0.0d0
    280    continue
    290    d(i) = h
    300 continue
    do 500 i = 2, n
    l = i - 1
    z(n,l) = z(l,l)
    z(l,l) = 1.0d0
    h = d(i)
    if (h .eq. 0.0d0) go to 380
    do 330 k = 1, l
    330    d(k) = z(k,i) / h
    do 360 j = 1, l
        g = 0.0d0
        do 340 k = 1, l
    340       g = g + z(k,i) * z(k,j)
        do 360 k = 1, l
            z(k,j) = z(k,j) - g * d(k)
    360    continue
    380    do 400 k = 1, l
    400    z(k,i) = 0.0d0
    500 continue
    510 do 520 i = 1, n
    d(i) = z(n,i)
    z(n,i) = 0.0d0
    520 continue
    z(n,n) = 1.0d0
    e(1) = 0.0d0
    return
    end
