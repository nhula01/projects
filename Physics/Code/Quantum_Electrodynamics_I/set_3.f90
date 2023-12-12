subroutine thickness_func(thickness_mid, name)
    implicit none
    character(len=100), intent(in) :: name
    character(len=100) :: file_name
    !~~~===== fundamental constants =====~~~!
    double precision, parameter :: pi=3.141592653589793,c=299792458.0
    double precision, parameter :: sqrt2=1.414213562373095,sqrt3=1.7320508075688772
    double precision, parameter :: mu0=4.0D-7*pi,eps0=1.0/(c*c*mu0)
    double precision, parameter :: h=1.054571628D-34
    double complex, parameter :: Im=(0.0,1.0)
    double precision, parameter :: ev_to_radsec=2.0*pi*2.4180d14,Hz_to_ev=4.1356691d-15
    !~~~============== grid =============~~~!
    double precision, parameter :: dz=1.0d-9 !spatial resolution is 1 nm here (everything is in SI units)
    integer, parameter :: Nz=3000 !number of steps along z
    double precision, parameter :: z0=-dz*Nz*0.5,zM=dz*Nz*0.5-dz !grid is between z0 and zM
    double precision z(Nz) !coordinate z
    double precision, parameter :: dt=dz/(2.0*c) !time step
    !~~~======== time iterations ========~~~!
    integer, parameter :: Nt_1fs=int(1.0d-15/dt)  !number of time steps per 1 fs
    integer, parameter :: n_fs=1000                 !<--- number of fs for free time dynamics when detecting
    integer, parameter :: Nt=n_fs*Nt_1fs          !total number of time iterations
    !~~~======== full EM fields =========~~~!
    double precision Ex(Nz),Hy(Nz-1) !electric & magnetic fields
    double precision, parameter :: eps_dielectric=1.0 !permittivity of medium (1 for vacuum)
    double precision, parameter :: dt_eps0=dt/(eps0*eps_dielectric),dt_mu0=dt/mu0 !constants for Maxwell's eq.
    !~~~======== spatial position of the source =========~~~!
    integer, parameter :: k_source=50 !placing source 50 steps to the right of z0
    !~~~======== laser pulse parameters =========~~~!
    double precision, parameter :: tau=1.0d-15 !pulse duration
    double precision, parameter :: E0=1.0d0    !pulse amplitude
    double precision, parameter :: omega=ev_to_radsec*2.0d0 !carrier frequency (number is in eV)
    double precision, parameter, dimension (4) :: aBH=(/0.353222222d0,-0.488d0,0.145d0,-0.010222222d0/)
    double precision pulse
    !~~~= other parameters and variables =~~~!
    integer k,n
    double precision t,tmp
    !
    !~~~ absorbing boundaries - CPML ~~~!
    integer, parameter :: npml=19,m=3,ma=1 !npml corresponds to the total number of PM layers
    double precision, parameter :: alphaCPML=0.05,kappaCPML=5.0 !<--- optimal parameters per Taflove
    double precision sigmaCPML        !<--- this one is also set to its optimal value
    double precision psi_Exz_1(npml),psi_Exz_2(npml)
    double precision psi_Hyz_1(npml-1),psi_Hyz_2(npml-1)
    double precision be_z(npml),ce_z(npml),alphae_z(npml),sige_z(npml),kappae_z(npml)
    double precision bh_z(npml-1),ch_z(npml-1),alphah_z(npml-1),sigh_z(npml-1),kappah_z(npml-1)
    double precision den_ez(Nz),den_hz(Nz)
    integer kk
    !
    !~~~~======================== Drude model for Au ========================~~~~!
    double precision, parameter :: omegaD=ev_to_radsec*7.039005173444189 !plasma frequency
    double precision, parameter :: GammaD=ev_to_radsec*0.180924535083724 !damping
    double precision, parameter :: thickness=50.0345435d-9 !metal thickness
    double precision, parameter :: distance=800.0345098d-9 !distance between mirrors
    double precision,intent(in) :: thickness_mid!=20.0345435d-9 !metal thickness
    double precision :: thickness_middle!=20.0345435d-9 !metal thickness
    logical FB_metal(Nz) !logical array defining metal film on FDTD grid
    double precision PDx(Nz),tmpE !current density inside metal and tmp variable
    double precision, parameter :: eps_r=1.0 !<--- this and below parameters for the Drude model and ADE propagation
    double precision, parameter :: A1=(2.0-GammaD*dt)/(2.0+GammaD*dt),A2=eps0*omegaD*omegaD*dt/(2.0+GammaD*dt)
    double precision, parameter :: C1=(eps_r*eps0/dt-0.5*A2)/(eps_r*eps0/dt+0.5*A2)
    double precision, parameter :: C3=1.0/(eps_r*eps0/dt+0.5*A2)
    double precision, parameter :: C4=0.5*(A1+1.0)/(eps_r*eps0/dt+0.5*A2)

    !~~~~~~~ Fourier analysis !!!!!!!
    double precision, parameter :: omega_min=0.5,omega_max=4.0
    integer N_w
    integer, parameter :: Nt_skip=20
    integer n_new,nn
    integer, parameter :: Nt_new=int(Nt/Nt_skip)
    double precision wsave(4*Nt_new+15) !working array for FFT
    double complex detection_Ex(Nt_new),detection_Hy(Nt_new),Sz

    thickness_middle=thickness_mid*1.0d-9
    !<--- grid
    do k=1,Nz
        z(k)=z0+dz*(k-1)
    enddo

    !<--- setting geometry
    FB_metal=.false.
    do k=1,Nz
        if((z(k)>(distance/2.0)).and.(z(k)<(distance/2.0+thickness)))then !mirror on the right
            FB_metal(k)=.true.
            elseif((z(k)>(-distance/2.0-thickness)).and.(z(k)<(-distance/2.0)))then !mirror on the left
            FB_metal(k)=.true.
            elseif((z(k)>(-thickness_middle/2.0)).and.(z(k)<(thickness_middle/2.0)))then !mirror in the middle
            FB_metal(k)=.true.
        endif
    enddo

    !setting absorbing boundaries!
    sigmaCPML=0.8*(m+1)/(dz*(mu0/eps0*eps_dielectric)**0.5)
    do k=1,npml
        sige_z(k)=sigmaCPML*((npml-k)/(npml-1.0))**m
        alphae_z(k)=alphaCPML*((k-1)/(npml-1.0))**ma
        kappae_z(k)=1.0+(kappaCPML-1.0)*((npml-k)/(npml-1.0))**m
        be_z(k)=exp(-(sige_z(k)/kappae_z(k)+alphae_z(k))*dt/eps0)
        if ((sige_z(k)==0.0).and. &
            (alphae_z(k)==0.0).and. &
            (k==npml))then
            ce_z(k)=0.0
        else
            ce_z(k)=sige_z(k)*(be_z(k)-1.0)/(sige_z(k)+kappae_z(k)*alphae_z(k))/kappae_z(k)
        endif
    enddo

    do k=1,npml-1
        sigh_z(k)=sigmaCPML*((npml-k-0.5)/(npml-1.0))**m
        alphah_z(k)=alphaCPML*((k-0.5)/(npml-1.0))**ma
        kappah_z(k)=1.0+(kappaCPML-1.0)*((npml-k-0.5)/(npml-1.0))**m
        bh_z(k)=exp(-(sigh_z(k)/kappah_z(k)+alphah_z(k))*dt/eps0)
        ch_z(k)=sigh_z(k)*(bh_z(k)-1.0)/(sigh_z(k)+kappah_z(k)*alphah_z(k))/kappah_z(k)
    enddo

    kk=npml
    do k=1,Nz-1
        if(k<=npml)then
            den_ez(k)=1.0/(kappae_z(k)*dz)
        elseif(k>=(Nz+1-npml))then
            den_ez(k)=1.0/(kappae_z(kk)*dz)
            kk=kk-1
        else
            den_ez(k)=1.0/dz
        endif
    enddo

    kk=npml-1
    do k=1,Nz-1
        if(k<=(npml-1))then
            den_hz(k)=1.0/(kappah_z(k)*dz)
        elseif(k>=(Nz+1-npml))then
            den_hz(k)=1.0/(kappah_z(kk)*dz)
            kk=kk-1
        else
            den_hz(k)=1.0/dz
        endif
    enddo

    psi_Exz_1=0.0
    psi_Exz_2=0.0
    psi_Hyz_1=0.0
    psi_Hyz_2=0.0
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
            !~~~ end of CPML setup ~~~!
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    Ex=0.0
    Hy=0.0
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
        pulse=tmp*E0*sin(omega*t)

    !~~~~~~ Hy ~~~~~~~~!
    do k=1,Nz-1
        Hy(k)=Hy(k)+dt_mu0*(Ex(k)-Ex(k+1))*den_hz(k)
    enddo
    !  PML for the left side Hy
        do k=1,npml-1
        psi_Hyz_1(k)=bh_z(k)*psi_Hyz_1(k)+ch_z(k)*(Ex(k)-Ex(k+1))/dz
        Hy(k)=Hy(k)+dt_mu0*psi_Hyz_1(k)
    enddo
    !  PML for the right side Hy
        kk=npml-1
        do k=Nz+1-npml,Nz-1
        psi_Hyz_2(kk)=bh_z(kk)*psi_Hyz_2(kk)+ch_z(kk)*(Ex(k)-Ex(k+1))/dz
        Hy(k)=Hy(k)+dt_mu0*psi_Hyz_2(kk)
        kk=kk-1
        enddo

    !~~~~~ Ex ~~~~~~~~~!
    do k=2,Nz-1
        if(k==k_source)then
            Ex(k)=Ex(k)+dt_eps0*(Hy(k-1)-Hy(k))/dz+pulse
        elseif(FB_metal(k))then
            tmpE=C1*Ex(k)+C3*(Hy(k-1)-Hy(k))/dz-C4*PDx(k)
            PDx(k)=A1*PDx(k)+A2*(tmpE+Ex(k))
            Ex(k)=tmpE
        else
            Ex(k)=Ex(k)+dt_eps0*(Hy(k-1)-Hy(k))*den_ez(k)
        endif
    enddo
    !  PML for the left side Ex
        do k=2,npml
            psi_Exz_1(k)=be_z(k)*psi_Exz_1(k)+ce_z(k)*(Hy(k-1)-Hy(k))/dz
            Ex(k)=Ex(k)+dt_eps0*psi_Exz_1(k)
        enddo
    !  PML for right side Ex
        kk=npml
        do k=Nz+1-npml,Nz-1
            psi_Exz_2(kk)=be_z(kk)*psi_Exz_2(kk)+ce_z(kk)*(Hy(k-1)-Hy(k))/dz
            Ex(k)=Ex(k)+dt_eps0*psi_Exz_2(kk)
            kk=kk-1
        enddo
    !--------------------------------------------------------------------------!
    !~~~~~~~~~~~~~~~~~~~~~~~~         detection          ~~~~~~~~~~~~~~~~~~~~~~!
    !--------------------------------------------------------------------------!
        if(mod(n,Nt_skip)==0)then
            n_new=int(n/Nt_skip)
            detection_Ex(n_new)=Ex(Nz-npml-50) !detection occurs on the other side of the cavity
            detection_Hy(n_new)=(Hy(Nz-npml-50)-Hy(Nz-npml-50-1))/2.0 !i.e. we calculate transmission
        endif
    enddo !Nt

    call zffti(Nt_new,wsave)
    call zfftf(Nt_new,detection_Ex,wsave)
    detection_Ex=detection_Ex/sqrt(float(Nt_new))
    call zfftf(Nt_new,detection_Hy,wsave)
    detection_Hy=detection_Hy/sqrt(float(Nt_new))

    file_name = '/Users/phihung/NumMethod/first/homework/hw_09/' // name
    open(file=file_name,unit=32)
        do n=2,Nt_new/2
            tmp=dble(n-1)/(dt*Nt_skip*Nt_new)
            tmp=Hz_to_ev*tmp
            if((tmp>omega_min).and.(tmp<omega_max))then
            Sz=detection_Ex(n)*conjg(detection_Hy(n))
            write(32,*) tmp,abs(Sz)
            endif
        enddo
    close(unit=32)


end subroutine thickness_func

program main
    implicit none
    double precision :: thickness
    integer :: i
    character(len=100) :: str
    do i=10,70,10
        thickness = i 
        write(str, '(I5.5)') i ! convert to string
        call thickness_func(thickness,str)
    enddo 
end program main
!
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
