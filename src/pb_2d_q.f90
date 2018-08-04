    !****************************************************************************
    !
    !  PROGRAM: PB_2D_quenched
    !
    !  PURPOSE:  Entry point for the console application PB_2D_1
    ! scan partile charge and chain discretization
    ! non-restartable
    !****************************************************************************
    !
    module fields
    USE BOX2D
    IMPLICIT NONE
    REAL,DIMENSION(NT,0:N_comp-1),PUBLIC :: RHO,W,WP ! rhocp ! rho(:,0): A-chains,rho(:,1) B-chains rhocp(:,:) AB copolymer
    REAL,DIMENSION(NT),PUBLIC ::  PP! rhocp ! rho(:,0): A-chains,rho(:,1) B-chains rhocp(:,:) AB copolymer
    end module
    PROGRAM  PB_2D_1
    USE BOX2D
    use fields
    USE MYFFT2D
    USE DENSITY
    USE INOUT2d
    USE PBSOLVER2D
    USE RAN_NUMBERS
    USE TIME
    USE INIT_DENS
    IMPLICIT NONE
    INTEGER  ::  I,J,IMAX,PBIT                ! IMAX maximum number of SCF iterations in the current run
    real :: RHOIPB,RHOIMB !bulk equilibrium concentrations of the ions
    REAL :: CHIABN
    REAL ::  M_A, M_B,M_AB,F_AB      !  F_AB is the composition of the diblock
    REAL :: DELTA_WA,DELTA_WB,DELTA_WP,ERRW,ERRP,RHONORM,ERRP2,ERRW02,ERRW12
    REAL ::CHIAMN,CHIAPN, CHIBMN, CHIBPN
    LOGICAL :: SAVE_SNAP_RHO, SAVE_SNAP_PSI
    NAMELIST /START/ IMAX,DELTA_WA,DELTA_WB,DELTA_WP,DELTA_PSI,LX,LY,LZ,SAVE_SNAP_RHO,SAVE_SNAP_PSI
    NAMELIST /STARTPOLYMERBLEND/   RG,RGA,RGB,RGAB,FH,F_AB,CHIABN,RHO0 !RHO0 is the inverce occupied volume of a chain of of a size Rg


    NAMELIST /IONS/  LBX,BSP,FQ,M_PLUS, M_MINUS,RBplus,RBMnus,PBIT !FQ in the nymber of ions per A chain; LBX is the 4*pi*Bjerrum_length in vacuum in nano meters at 300K , PBIT -P-B iterations
    !BSP/(R_ion*Epsilon) is the Born energy of an ion
    NAMELIST /EPSILONS/ EPS_AI,EPS_BI

    REAL(8) :: DELTAPHI_MAX,FE,ESE
    REAL(8) :: TSPL,TSMIN
    REAL(8) :: ENTL, TS(0:1)

    REAL(8) ::L_D,Epsilion_Sol ! RESXX(24),


    !integer::ix,iy,iz,ir
    logical:: EXIST, EXISTD
    CHARACTER(*), PARAMETER :: SFN="scenariop.lst" ,SAFN="scenarioa.lst"
    CHARACTER(*), PARAMETER :: STARTFN="start.lst",POLYFN="startapolymer.lst",IONFN="ions.lst",epFN="epsilons.lst"

    TYPE(RAN_SAVE):: XSEED

    CALL SYSTEM_CLOCK (COUNT=J)
    CALL RAN_INIT(J,XSEED)


    CALL INIT_FFT()
    print*, "Stage 1"
    IF(INIT_READ_POYMER()) THEN
        PRINT*, "init configuration is red"
    ELSE
        PRINT*, "init CONFIGURATION CREATED"
        STOP
    END IF
    ! Body of test_diffusion
    print *, 'Begin computation of the Free Energy'


    CALL MAKE_IQ2(IQ2)
    call SET_CHARGE()
    RHOIPB=FH(0)*M_PLUS*LVI*(PI*4./3.)*RBplus**3
    RHOIMB=FH(0)*M_MINUS*LVI**(PI*4./3.)*RBplus**3
    L_D=1/SQRT((M_PLUS+M_MINUS)*LVI*EPSCI(1.,0.))
    Epsilion_Sol=LBX/EPSCI(1.,0.)
    Print '("Ions Debye-Huckel length is:",G11.3," Epsilon=",G11.3," fi=",G11.3)',  Real(L_D,4),REAL( Epsilion_Sol,4),REAL(M_IONS*LVI,4)
    print*, "Stage 3 ", "norm=",LVN*sum(rho)

    CALL SET_CPU_TIME()

    CALL SET_DS(64,RG,RG,0.)
    !W=0._8 !init field
    CALL SET_DELTA_PSI(1.,DELTA_PSI)

    INQUIRE (FILE='fields.dat', EXIST=EXISTD)
    INQUIRE (FILE='densxy.dat', EXIST=exist)
    READwrho: if(EXISTD) then
        OPEN (unit=15, FILE='fields.dat', FORM='FORMATTED', STATUS='OLD')
        READ(15,*) W,PSI,WP,PP
        CLOSE(15)
        CALL HOMOPOLYMER_DENSITY(W,RHO,QH,NS)
        print*, 'fields.dat has been read '
        CALL SET_DELTA_PSI(.1)
        CALL ITERATE_PB(RHO,1000,DELTAPHI_MAX,ESE,TSPL,TSMIN)
        print*, "dPhi=", real(DELTAPHI_MAX,4),"dPSI=",real(DELTA_PSI,4)
    ELSE IF (EXIST) then
        OPEN (unit=15, FILE='densxy.dat', FORM='FORMATTED', STATUS='OLD')
        READ(15,*) RHO
        CLOSE(15)
        PSI=0._8
        CALL SET_DELTA_PSI(1.)
        CALL ITERATE_PB(RHO,1)
        CALL SET_DELTA_PSI(.1)
        CALL ITERATE_PB(RHO,5000,DELTAPHI_MAX,ESE,TSPL,TSMIN)
        print*, "dPhi=", real(DELTAPHI_MAX,4),"dPSI=",real(DELTA_PSI,4)
        print*, "dPhi=", real(DELTAPHI_MAX,4),"dPSI=",real(DELTA_PSI),"ese=",real(ESE)
        print*, 'densxyz.dat has been read '
        CALL MAKE_WAB(RHO,PSI,E2,W,WP,PP)
    ELSE
        call MAKE_INIT(RHO,Real(LY*0.3989422804014327))
        OPEN (unit=15, FILE='dens.dat', FORM='FORMATTED', STATUS='UNKNOWN')
        WRITE(15,*) RHO
        CLOSE(15)
        PSI=0._8
        CALL SET_DELTA_PSI(1.)
        CALL ITERATE_PB(RHO,1)
        CALL SET_DELTA_PSI(.1)
        CALL ITERATE_PB(RHO,200000,DELTAPHI_MAX,ESE,TSPL,TSMIN)
        print*, "dPhi=", real(DELTAPHI_MAX,4),"dPSI=",real(DELTA_PSI,4)
        !pause
        CALL MAKE_WAB(RHO,PSI,E2,W,WP,PP)
    END IF READwrho


    OPEN (unit=15, FILE='fields.dat', FORM='FORMATTED', STATUS='UNKNOWN')
    WRITE(15,*) W,PSI
    CLOSE(15)
    print*, 'fields.dat has been written '

    DO I=1,IMAX
        DO J=1,100
            CALL HOMOPOLYMER_DENSITY(W,RHO,QH,NS)
            CALL ITERATE_PB(RHO,100)
            call EVOLVE_WAB(RHO,W,WP,PP,psi,E2,DELTA_WA,DELTA_WB,DELTA_WP,ERRW,ERRP,RHONORM,errp2,errw02,errw12)

        END DO
        CALL ITERATE_PB(RHO,1000,DELTAPHI_MAX,ESE,TSPL,TSMIN)

        OPEN (unit=14, FILE='dens.dat', FORM='FORMATTED', STATUS='UNKNOWN')
        WRITE(14,*) RHO
        CLOSE(14)
        OPEN (unit=15, FILE='fields.dat', FORM='FORMATTED', STATUS='UNKNOWN')
        WRITE(15,*) W,PSI,WP,PP
        CLOSE(15)
        OPEN (unit=16, FILE='ions.dat', FORM='FORMATTED', STATUS='UNKNOWN')
        WRITE(16,*) RHOI
        CLOSE(16)
        OPEN (unit=17, FILE='qdens.dat', FORM='FORMATTED', STATUS='UNKNOWN')
        WRITE(17,*) RHOQ
        CLOSE(17)
        print*, 'fields.dat has been written '

        CALL GET_F_HOMO(FE,TS,ENTL)

        Print'("ESE=",1PG24.15E3," TS+",1PG24.15E3," TS-",1PG24.15E3)',ESE,TSPL,TSMIN
        print'("Fe=",1PG24.15E3," TS",2(1PG24.15E3)," ENTL",1PG24.15E3)',FE,TS,ENTL
        Print*, "Fe_Total=",FE+ESE+TSPL+TSMIN
        print '("dPsi=",G9.3," Ers=",3G12.4," ErrP=",2G12.4)',DELTAPHI_MAX,errw,errw02,errw12,errp,errp2
        Print '("SumRho",2G12.4," SumWa",2G12.4," SumWp",2G12.4)', nti*sum(rho(:,0)),nti*sum(rho(:,1)), nti*sum(w(:,0)),nti*sum(w(:,1)), sum(wp(:,0))/NT, sum(wp(:,1))/NT
        print*,"steptime=",RTIME(),"cputime=",UTIME()/nut


        !IF(ABS(DW)<1.E-4.and.abs(DELTAPHI_MAX)<1.E-6) EXIT
    END DO


    L_D=1/SQRT((M_PLUS+M_MINUS)*LVI*LBX*EPSCI(1.,0.))
    Epsilion_Sol=LBX/(EPSCI(1.,0.))
    Print '("Ions Debye-Huckel length is:",G11.3," Epsilon=",G11.3," fi=",G11.3)', &
        Real(L_D,4),REAL( Epsilion_Sol,4),REAL(FQ*RHO0*M_A*LVI,4)
    !CALL MAKE_IQ2(IQ2)

    STOP  "calculations ended"
    CONTAINS
    !=============================================================================================

    LOGICAL FUNCTION INiT_READ_POYMER() RESULT(STARTED)
    LOGICAL :: EXIST !, FINISHED   !INTEGER :: IVAR

    INQUIRE (FILE=STARTFN,EXIST=STARTED)
    if(STARTED)then
        OPEN(UNIT=10,FILE=STARTFN,FORM ='FORMATTED', ACTION='READ',STATUS='UNKNOWN',DELIM='APOSTROPHE')
        READ(10,START)
        CLOSE(10)
        PRINT START
    ELSE
        IMAX=1000
        DELTA_WA=.1;DELTA_WB=.1;DELTA_WP=.1;
        DELTA_PSI=.1
        LX=173.20508075688772; LY=100.; LZ=100.
        SAVE_SNAP_RHO=.false. ;  SAVE_SNAP_PSI=.false.
    END IF
    call INIT_BOX()
    OPEN(UNIT=10,FILE=STARTFN,FORM ='FORMATTED', ACTION='WRITE',STATUS='UNKNOWN',DELIM='APOSTROPHE')
    WRITE(10,START)
    CLOSE(10)
    CALL SET_DELTA_PSI(1.,DELTA_PSI)


    INQUIRE (FILE=POLYFN,EXIST=exist)
    if(exist)then
        OPEN(UNIT=10,FILE=POLYFN,FORM ='FORMATTED', ACTION='READ',STATUS='UNKNOWN',DELIM='APOSTROPHE')
        READ(10,STARTPOLYMERBLEND)
        CLOSE(10)
        PRINT STARTPOLYMERBLEND
    ELSE
        RG=5.;   RGA=5. ;   RGB=5.
        RGAB=0.
        FH(0)=.5 ;   FH(1)=.5
        CHIABN=25
        RHO0=1000.*LVI
    END IF

    RG=MAX(RGA,RGB,RGAB)
    IF(RGAB>0.) THEN
        M_AB=(1.-FH(0)-FH(1))**RHO0*LX*LY*LZ*(RG/RGAB)**2
    ELSE
        M_AB=0.
        FH(1)=1.-FH(0)
    END IF

    M_A=FH(0)*RHO0*VOL*(RG/RGA)**2
    M_B=FH(1)*RHO0*VOL*(RG/RGA)**2

    OPEN(UNIT=10,FILE=POLYFN,FORM ='FORMATTED', ACTION='WRITE',STATUS='UNKNOWN',DELIM='APOSTROPHE')
    WRITE(10,STARTPOLYMERBLEND)
    CLOSE(10)

    INQUIRE (FILE=IONFN,EXIST=exist)
    if(exist)then
        OPEN(UNIT=10,FILE=IONFN,FORM ='FORMATTED', ACTION='READ',STATUS='UNKNOWN',DELIM='APOSTROPHE')
        READ(10,IONS)
        CLOSE(10)
        PRINT IONS
    ELSE
        LBX=700.! LBX=210000//T(Kelvins)
        BSP=28. !BSP=8400/T(Kelvins)
        FQ=10;
        M_MINUS=0.
        RBplus=0.3;RBMnus=0.3
        PBIT=1000
    END IF

    M_PLUS=FQ*(M_A+M_AB*F_AB*(RGAB/RGA)**2)+M_MINUS
    M_IONS=M_PLUS+M_MINUS

    IRP=BSP/(LBX*RBplus)
    IRM=BSP/(LBX*RBMnus)
    CHIAMN=bsp*EPS_AI/(rho0*RBMnus)
    CHIAPN=bsp*EPS_AI/(rho0*RBplus)
    CHIBMN=bsp*EPS_BI/(rho0*RBMnus)
    CHIBPN=bsp*EPS_BI/(rho0*RBplus)
    OPEN(UNIT=10,FILE=IONFN,FORM ='FORMATTED', ACTION='WRITE',STATUS='UNKNOWN',DELIM='APOSTROPHE')
    WRITE(10,IONS)
    CLOSE(10)

    INQUIRE (FILE=epFN,EXIST=exist)
    if(exist)then
        OPEN(UNIT=10,FILE=epFN,FORM ='FORMATTED', ACTION='READ',STATUS='UNKNOWN',DELIM='APOSTROPHE')
        READ(10,EPSILONS)
        CLOSE(10)
        PRINT EPSILONS
    ELSE
        EPS_AI=0.1
        EPS_BI=0.5
        !EPS_IMI=0.5
        !EPS_IPI=0.5
    END IF

    OPEN(UNIT=10,FILE=epFN,FORM ='FORMATTED', ACTION='WRITE',STATUS='UNKNOWN',DELIM='APOSTROPHE')
    WRITE(10,EPSILONS)
    CLOSE(10)
    STARTED=STARTED.AND.EXIST
    RETURN
3   STOP 'STARTED FILE ERROR '
    END FUNCTION INIT_READ_POYMER

    SUBROUTINE SET_CHARGE()
    FQA=FQ*(RGA/RG)**2
    M_PLUS=FQ*(M_A*(RGA/RG)**2+M_AB*F_AB*(RGAB/RG)**2)+M_MINUS
    M_IONS=M_PLUS+M_MINUS
    PQ_norm=FQA*RHO0

    END  SUBROUTINE SET_CHARGE
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE MAKE_WAB(RHO,PSI,E2,W,WP,PP)
    REAL, DIMENSION(NT,0:1) :: RHO,W,WP
    REAL(8), DIMENSION(NT) :: PSI,E2
    REAL, DIMENSION(NT) :: PP
    REAL :: WES,WESA,WESB,XP
    INTENT(IN) :: PSI,E2
    INTENT(INOUT) :: RHO
    INTENT(OUT):: W,WP
    INTEGER :: I
    REAL(8):: XW0,XW1, XPX,RHONORM
    xw0=0.; xw1=0.; XPX=0. ;RHONORM=0.
    !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I, WES,WESA,WESB,XP)
    !$OMP  DO   REDUCTION(+:RHONORM)
    DO I = 1, NT
        RHONORM =RHONORM+(RHO(I,1)+RHO(I,0))
    END DO
    !$OMP END  DO
    !$OMP SINGLE
    RHONORM=NT/RHONORM
    !$OMP END SINGLE
    !$OMP BARRIER
    !$OMP DO
    DO I = 1, NT
        RHO(I,0)=RHONORM*RHO(I,0)
        RHO(I,1)=RHONORM*RHO(I,1)
    END DO
    !$OMP END  DO
    !$OMP  DO REDUCTION(+:XW0,XW1,XPX)
    DO I = 1, NT
        XP=RHO(I,1)+RHO(I,0)-1.
        PP(i)=XP
        WES=-.5*LBX*E2(I)/(EPS_I(I)**2) ;  WESA=-PQ_NORM*PSI(I)+WES*EPS_AI ;  WESB=WES*EPS_BI
        XPX=XPX+XP
        WP(I,0) =RHO(I,1)*CHIABN +WESA+RHOI(0,i)*CHIAPN+RHOI(1,i)*CHIAMN !WPA
        W(I,0) =XP + WP(I,0)
        XW0=XW0+W(I,0)
        WP(I,1) =CHIABN*RHO(I,0)+WESB+RHOI(0,i)*CHIBPN+RHOI(1,i)*CHIBPN !WPB
        W(I,1) =XP+WP(I,1)
        XW1=XW1+W(I,1)
    END DO
    !$OMP END DO
    !$OMP SINGLE
    XW0=XW0/NT; XW1=XW1/NT; XPX=XPX/NT
    !$OMP END SINGLE
    !$OMP BARRIER
    !$OMP  DO
    DO I = 1, NT
        W(I,0)=W(I,0)-xw0
        W(I,1)=W(I,1)-xw1
        PP(I)=PP(I)-XPX
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    END SUBROUTINE MAKE_WAB

    SUBROUTINE EVOLVE_WAB(RHO,W,WP,PP,psi,E2,EPSWA,EPSWB,EPSP,ERRW,ERRP,RHONORM,errp2,errw02,errw12)
    REAL, DIMENSION(NT,0:1) :: RHO,W,WP
    REAL, DIMENSION(NT) :: PP
    REAL(8), DIMENSION(NT) :: PSI,E2
    REAL :: EPSWA,EPSWB,EPSP,ERRW,ERRP,RHONORM,errp2,errw02,errw12
    INTENT(IN) :: RHO,PSI,E2,EPSWA,EPSWB,EPSP
    INTENT(INOUT) :: W,WP,PP
    INTEGER :: I
    REAL(8)::  XW0,XW1,XP
    REAL(8)::  XN0,XN1,XPX,WESA,WESB,WES
    REAL(8)::  RHA,RHB
    INTENT(OUT) ::ERRW,ERRP,errp2,errw02,errw12,RHONORM

    XN0=0.;XN1=0.;XPX=0.;errp2=0.;errw02=0.;errw12=0.
    ERRW=0.; ERRP=0.;RHONORM=0.
    !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I,XW0,XW1,RHA,RHB,XP,WES,WESA,WESB)
    !$OMP  DO   REDUCTION(+:RHONORM)
    DO I = 1, NT
        RHONORM =RHONORM+(RHO(I,1)+RHO(I,0))
    END DO
    !$OMP END  DO
    !$OMP SINGLE
    RHONORM=NT/RHONORM
    !$OMP END SINGLE
    !!$OMP FLUSH(RHONORM)
    !$OMP  DO  REDUCTION(MAX:ERRW,ERRP),reduction(+:xn0,xn1,XPX,errp2,errw02,errw12)
    DO I = 1, NT
        RHA=RHO(I,0);RHB=RHO(I,1)
        XP=RHONORM*(RHA+RHB)-1.
        WES=-.5*LBX*E2(I)/(EPS_I(I)**2) ;  WESA=-PQ_NORM*PSI(I)+WES*EPS_AI ;  WESB=WES*EPS_BI
        ERRP=MAX(ERRP,ABS(XP))
        ERRP2=ERRP2+XP*XP
        PP(I)=PP(I)+EPSP*XP
        XW0 =RHB*CHIABN+WESA-WP(I,0)+RHOI(0,i)*CHIAPN+RHOI(1,i)*CHIAMN !WPA
        errw02=errw02+XW0**2
        XW1 =CHIABN*rha+WESB-WP(I,1)+RHOI(0,i)*CHIBPN+RHOI(1,i)*CHIBPN !WPB
        ERRW12=ERRW12+XW1**2
        ERRW=MAX(ERRW,ABS(XW0),ABS(XW1))
        WP(I,0)=WP(I,0)+EPSWA*XW0
        WP(I,1)=WP(I,1)+EPSWB*XW1
        XN0=XN0+WP(I,0)
        XN1=XN1+WP(I,1)
        XPX=XPX+PP(I)
    END DO
    !$OMP END  DO
    !$OMP SINGLE
    XN0=XN0/NT; XN1=XN1/NT; XPX=XPX/NT
    !$OMP END SINGLE
    !!$OMP FLUSH(XN0,XN1,XPX)
    !$OMP  DO
    DO I = 1, NT
        PP(I)=PP(I)-XPX
        W(I,0)=WP(I,0)+PP(I)-XN0
        W(I,1)=WP(I,1)+PP(I)-XN1
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    ERRP2=SQRT(ERRP2/NT)
    ERRW02=SQRT(ERRW02/NT)
    ERRW12=SQRT(ERRW12/NT)
    END SUBROUTINE evolve_WAB
    !-----------------------------------------------------------------------------------------------------

    SUBROUTINE SET_DELTA_PSI(X,BASE)
    REAL,INTENT(IN) :: X,BASE
    OPTIONAL :: BASE
    REAL, SAVE :: XBASE=1.
    IF(PRESENT(BASE))  XBASE=BASE
    DELTA_PSI=X*XBASE
    LAMBDA_PSI=1.-DELTA_PSI
    Print*, "delta psi set at",DELTA_PSI, LAMBDA_PSI
    END SUBROUTINE SET_DELTA_PSI

    SUBROUTINE GET_F_HOMO(F,ts,entl)
    REAL(8),INTENT(OUT) :: F,TS(0:1),entl
    REAL(8) :: TS0,TS1
    INTEGER::I
    TS0=0.; ts1=0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
    !$OMP DO   REDUCTION(+:TS0)
    DO I=1,NT
        TS0=TS0+RHO(I,0)*W(I,0)
    ENDDO
    !$OMP END  DO NOWAIT
    !$OMP DO   REDUCTION(+:TS1)
    DO I=1,NT
        TS1=TS0+RHO(I,1)*W(I,1)
    ENDDO
    !$OMP END  DO
    !$OMP SINGLE
    TS(0)=TS0*NTI+log(QH(0))
    TS(1)=TS1*NTI+log(QH(1))
    ts0=0
    !$OMP END SINGLE
    !$OMP BARRIER
    !$OMP DO REDUCTION(+:TS0)
    DO I=1,NT
        TS0=RHO(I,0)*(CHIABN*RHO(I,1)+RHOI(0,i)*CHIAPN+RHOI(1,i)*CHIAMN)+RHO(I,1)*(RHOI(0,i)*CHIBPN+RHOI(1,i)*CHIBPN)
    ENDDO
    !$OMP END  DO
    !$OMP END  PARALLEL
    ENTL=RHO0*VOL*TS0*NTI
    F=-TS(0)*M_A-TS(1)*M_B+ENTL
    RETURN
    END SUBROUTINE GET_F_HOMO
    END PROGRAM PB_2D_1

