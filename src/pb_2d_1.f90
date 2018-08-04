    !****************************************************************************
    !
    !  PROGRAM: PB_2D_1
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
    IMPLICIT NONE
    INTEGER  ::  I,J,K,IMAX,PBIT ,ls               ! IMAX maximum number of SCF iterations in the current run
    REAL :: F1,F2,RHO1,RHO2,RHON
    real :: RHOIPB,RHOIMB !bulk equilibrium concentrations of the ions
    REAL :: CHIABN, CHIAIM,CHIAIP,CHIBIM,CHIBIP
    REAL ::  M_POLY            ! NUMBER OF POLYMER CHAINS IN THE PERIODICAL CELL
    REAL :: DELTA_W,LAMBDA_W
    LOGICAL :: SAVE_SNAP_RHO, SAVE_SNAP_PSI, SAVE_SNAP_P
    NAMELIST /STARTPOLYMERBLEND/   RG,RGA,RGB,RGAB,M_POLY,CHIABN,FQ, M_PLUS, M_MINUS,PKA_P,RHO0

    NAMELIST /START/ IMAX,DELTA_W,DELTA_PSI,LX,LY,LZ,PBIT, & !, PBIT -P-B iterations
    SAVE_SNAP_RHO,SAVE_SNAP_PSI

    NAMELIST /IONS/  RBplus,RBMnus,LBX !LBX is the 4*pi*Bjerrum_length in vacuum in nano meters at 300K 

    NAMELIST /EPSILONS/ EPS_AI,EPS_BI,EPS_IMI,EPS_IPI

    REAL(8) :: WX,DW,DELTAPHI_MAX,FE,TS,TSP,TSPLUS,TSMINUS,EW,EPLUS,EMINUS,ESE,ENTLP,QTEST,FE0,Eads,ALPHAPS,ALPHBETCS
    REAL(8) :: EAPLUS,EAMINUS

    REAL(8) ::L_D,Epsilion_Sol ! RESXX(24),
    real(8) :: psir(0:nx/2-1),rhor(0:nx/2-1),rhoir(0:nx/2-1)
    real(8):: alphapbulk,wbulk,alphaFEbulk,LZX
    !integer::ix,iy,iz,ir

    !CHARACTER(LEN=16):: HEADXX(24)=["Ns","M_poly","Qp","M_+", "M_-","MI+","MI-","TSp","TS+","TS-","TS","Entlp","EsE","FE","Zx","Qz","Q+","Q-","Ea+","Ea-","F_DS_P","FDS_C","L_D","Epsilion_Sol"]

    CHARACTER(*), PARAMETER :: SFN="scenariop.lst" ,SAFN="scenarioa.lst"
    CHARACTER(*), PARAMETER :: STARTFN="start.lst",POLYFN="startapolymer.lst",epFN="epsilons.lst"

    INTEGER :: NDS0N,NDSI,NDSN !,IXQMAX,IXQI,IXQM ! initial, start, and
    REAL(8) :: MPMAX,MPI,MPM ! initial, start, and
    REAL(8) ::  ALPHA0,ALPHAM,ALPHAS
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
    
    L_D=1/SQRT((M_PLUS+M_MINUS)*LVI*EPSCI(1.,0.,RHOIPB,RHOIMB))
    Epsilion_Sol=1./EPSCI(1.,0.,RHOIPB,RHOIMB)
    Print '("Ions Debye-Huckel length is:",G11.3," Epsilon=",G11.3," fi=",G11.3)', &
        Real(L_D,4),REAL( Epsilion_Sol,4),REAL(FQ*RHO0*M_POLY*LVI,4)





    print*, "Stage 3 ", "norm=",LVN*sum(rho)

    CALL SET_CPU_TIME()


    !CALL INIT_DAT("results.tsv",HEAD=HEADXX,N=29)


    !DO WHILE(Q_PARTICLE<=IXQMAX)

    CALL SET_DS(NDS0N,RG,RG,0.)
    W=0._8 !init field
    PSI=0._8
    ZX=1.

    IF(M_POLY==0.)THEN
        RHO=LVI
        QZ=1.
    ELSE
        CALL HOMOPOLYMER_DENSITY(W,RHO,QZ,NS)
    END IF

    DO WHILE(M_POLY<=MPMAX)
        alphaPbulk=ALPHA0
        CALL SET_CHARGE()


        CALL SET_DELTA_PSI(.01)

        CALL ITERATE_PB(RHO,1000,QTEST,DELTAPHI_MAX,ESE)
        print*, "dPhi=", real(DELTAPHI_MAX,4),"dPSI=",real(DELTA_PSI,4),"Qt=",Real(QTEST,4),"Zx=",Real(ZX,4)
        DO WHILE(alphaPbulk<=ALPHAM)
            CALL SET_CHARGE()
            DO LS=1,NDSN
                call SET_DS((NDSI*2**LS)/2,RG,RG,0.)
                CALL SET_DELTA_PSI(.1)
                DO K=1,IMAX
                    DO J=1,100


                        IF(M_POLY==0.)THEN
                            RHO=LVI
                            QZ=1.
                            DW=0.
                        ELSE
                            DW=0.
                            !$OMP PARALLEL  DO DEFAULT(SHARED) PRIVATE(I,WX) REDUCTION(MAX:DW)
                            DO I=1,NT
                                WX=(RHO(I,1)-FQ*ALPHAP(I)*PSI(I)-0.5*LBX*EPS_AI*EPS_BI*FQ*RHO0*E2(I)/(EPS_I(I)**2)+(XLNX(ALPHAP(I))+XLNX(1.-ALPHAP(I))+ALPHAP(I)*LZX)*FQ)-WBULK
                                DW=MAX(DW,ABS(WX-W(I,1)))
                                W(I,1)=LAMBDA_W*W(I,1)+DELTA_W*WX
                            ENDDO
                            !$OMP END PARALLEL DO
                            CALL HOMOPOLYMER_DENSITY(W,RHO,QZ,NS)
                        END IF

                        CALL ITERATE_PB(RHO,max(1,(PBIT-1)/LS),QTEST)
                    END DO
                    CALL ITERATE_PB(RHO,max(1,(PBIT-1)/LS),QTEST,DELTAPHI_MAX,ESE)


                    CALL GET_TS_PB()
                    TSP=(1.+EW)
                    TSPLUS=MIPLUS*(1.+EPLUS)
                    TSMINUS=MIMINUS*(1.+EMINUS)
                    TS=TSP+TSPLUS+TSMINUS-ALPHAPS-ALPHBETCS
                    FE=(ENTLP+EsE+EAPLUS+EAMINUS-TS)-FE0
                    print*, "EsEtest-> rhofi=",ese,"E^2=",ese2




                    PRINT '("dPsi=",G9.3," dW=",G9.3,"DPsi=",G9.3," dRho=",G9.3," nmm=",G10.3," Zx=",G10.3," Qt=",G10.3,"spu:",2(I5,"c"))',&
                        DELTAPHI_MAX,dW,f1-f2,1._8-RHON,ZX,QTEST,Nint(rtime()), Nint(utime()/nut)
                    PRINT '("Qz=",G9.3," Mp=",G9.3," Q+=",G9.3," M+=",G9.3," Q-=",G9.3E3,"M-=",G10.3E3,"E+=",G10.3," E-=",G10.3," Ep=",G10.3)', &
                        QZ,Qplus,MIPLUS/M_PLUS,Qminus,MIMINUS/M_MINUS,EPLUS,EMINUS,EW
                    PRINT '(" Qp=",G10.3," Sp=",G10.3," S+=",G9.2," S-=",G9.2," S=",G10.3,"DP=",G11.3," DC=",G11.3)', Qpx,TSP,TSPLUS,TSMINUS,TS,-ALPHAPS,-ALPHBETCS
                    PRINT '("EsE=",G10.3,"Entlp=",G10.3," Ea+=",G10.3," Ea-=",G10.3," F=",G13.6)', ESE, ENTLP,EAPLUS,EAMINUS,FE

                    IF(ABS(DW)<1.E-4.and.abs(DELTAPHI_MAX)<1.E-6) EXIT
                END DO


                !RESXX=[REAL(NS,8),Qpx,M_PLUS,M_MINUS,MIPLUS,MIMINUS,TSP,TSPLUS,TSMINUS,TS,ENTLP,ESE,FE,ZX,QZ,Qplus,Qminus,EAPLUS,EAMINUS,ALPHAPS,ALPHBETCS,L_D,Epsilion_Sol]

                !CALL APPEND_DAT("results.tsv",dat=RESXX,N=26)
            END DO



            IF(SAVE_SNAP_RHO)THEN
                CALL DUMP_DENS(RHO,"Mp"//NUMCHAR(INT(M_POLY))//"Al"//NUMCHAR(INT(1000000*alphaPbulk))//"rho.dat")
            END IF





            IF(save_SNAP_psi)THEN
                CALL DUMP_DENS(PSI,"Mp"//NUMCHAR(INT(M_POLY))//"Al"//NUMCHAR(INT(1000000*alphaPbulk))//"psi.dat")
                CALL DUMP_IDENS(RHOI,"Mp"//NUMCHAR(INT(M_POLY))//"Al"//NUMCHAR(INT(1000000*alphaPbulk))//"plus.dat")
                
                !$OMP PARALLEL  DO DEFAULT(SHARED) PRIVATE(I)
                DO I=1,nt
                    PSI(I)=-FQ*ALPHAP(I)*RHO(I,1)
                ENDDO
                !$OMP END PARALLEL DO

                CALL DUMP_DENS(PSI,"Mp"//NUMCHAR(INT(M_POLY))//"Al"//NUMCHAR(INT(1000000*alphaPbulk))//"Qp.dat")

            END IF
            alphaPbulk=alphaPbulk+ ALPHAS
        ENDDO
        IF(M_POLY==0.AND.MPI==0 ) THEN
            M_POLY=MPM
        ELSE
            M_POLY=MPM*M_POLY+MPI
        END IF

        print*, "M_POLY=",M_POLY

        L_D=1/SQRT((M_PLUS+M_MINUS)*LVI*LBX*EPSCI(1.,0.,RHOIPB,RHOIMB))
        Epsilion_Sol=1./(EPSCI(1.,0.,RHOIPB,RHOIMB))
        Print '("Ions Debye-Huckel length is:",G11.3," Epsilon=",G11.3," fi=",G11.3)', &
            Real(L_D,4),REAL( Epsilion_Sol,4),REAL(FQ*RHO0*M_POLY*LVI,4)
        !CALL MAKE_IQ2(IQ2)
    END DO
    STOP  "calculations ended"
    CONTAINS
    !=============================================================================================
    SUBROUTINE TWUP()
    IF(DELTAPHI_MAX<1. ) THEN
        CALL SET_DELTA_PSI(3.)
    ELSEIF(DELTAPHI_MAX<10. ) THEN
        CALL SET_DELTA_PSI(2.5)
    ELSEIF(DELTAPHI_MAX<20. ) THEN
        CALL SET_DELTA_PSI(2.0)
    ELSEIF(DELTAPHI_MAX<50. ) THEN
        CALL SET_DELTA_PSI(1.5)
    ELSEIF(DELTAPHI_MAX<70. ) THEN
        CALL SET_DELTA_PSI(1.)
    ELSEIF(DELTAPHI_MAX<100. ) THEN
        CALL SET_DELTA_PSI(.5)
    ELSE
        CALL SET_DELTA_PSI(.2)
    END IF
    END SUBROUTINE TWUP
    !=============================================================================================
    LOGICAL FUNCTION INiT_READ_POYMER() RESULT(STARTED)

    LOGICAL :: EXIST !, FINISHED
    !INTEGER :: IVAR



    INQUIRE (FILE=STARTFN,EXIST=STARTED)
    if(STARTED)then
        OPEN(UNIT=10,FILE=STARTFN,FORM ='FORMATTED', ACTION='READ',STATUS='UNKNOWN',DELIM='APOSTROPHE')
        READ(10,START)
        CLOSE(10)
        PRINT START
    ELSEIF(started) then
        STOP "started, but start.lst does not exist"
    ELSE
        IMAX=100
        DELTA_W=.1
        DELTA_PSI=.1
        LBX=704.
       LX=100,LY,LZ,
        SAVE_SNAP_RHO=.false.
        SAVE_SNAP_PSI=.false.
    END IF
    call INIT_BOX()
    OPEN(UNIT=10,FILE=STARTFN,FORM ='FORMATTED', ACTION='WRITE',STATUS='UNKNOWN',DELIM='APOSTROPHE')
    WRITE(10,START)
    CLOSE(10)
    CALL SET_DELTA_PSI(1.,DELTA_PSI)
    LAMBDA_W  =1.-DELTA_W

    INQUIRE (FILE=POLYFN,EXIST=exist)
    if(exist)then
        OPEN(UNIT=10,FILE=POLYFN,FORM ='FORMATTED', ACTION='READ',STATUS='UNKNOWN',DELIM='APOSTROPHE')
        READ(10,STARTPOLYMERBLEND)
        CLOSE(10)
        PRINT STARTPOLYMERBLEND
    ELSEIF(started) then
        STOP "started, but startpolymer.lst does not exist"
    ELSE
        RG=5.
        M_POLY=100
        FQ=10
        M_MINUS=M_POLY*FQ/10
        M_PLUS=M_MINUS+M_POLY*FQ
        PKA_P=5.
       
    END IF
    STARTED=STARTED.AND.EXIST
    !ZPA=10.**(-PKA_P)
   


    PRINT*," FE0=",FE0,M_POLY*FQ*alphaFEbulk,M_PLUS,M_MINUS

    
    M_PLUS=M_POLY*FQ+M_MINUS
    M_IONS=M_PLUS+M_MINUS
    alphaPbulk=ALPHA0


    OPEN(UNIT=10,FILE=POLYFN,FORM ='FORMATTED', ACTION='WRITE',STATUS='UNKNOWN',DELIM='APOSTROPHE')
    WRITE(10,STARTPOLYMERBLEND)
    CLOSE(10)
    RETURN
3   STOP 'STARTED FILE ERROR '
    END FUNCTION INIT_READ_POYMER

    SUBROUTINE SET_CHARGE()
    !real(8) :: alphaPbulk
    alphaPbulk=min(.9999,max(.0001,alphaPbulk))
    
    ZPA=(1._8-ALPHAPBULK)/ALPHAPBULK
    LZX=log(ZPA)
    FQ=FQ/alphaPbulk
    PRINT*,"AlphapBulk=",Real(ALPHAPBULK,4)," FQ=",Real(FQ,4)," pKa_P=",Real(-LOG10(ZPA),4)


    !M_PLUS=M_POLY*FQ+M_MINUS
    M_PLUS=M_POLY*FQ+M_MINUS
    M_IONS=M_PLUS+M_MINUS

    alphaFEbulk=(XLNX(alphaPbulk)+XLNX(1.-alphaPbulk)+alphaPbulk*LZX)

    
    FE0=-M_POLY*(1.+M_POLY*LVI)- M_PLUS-M_MINUS
    !PRINT*,"QZ0=",QZ0I," FE0=",FE0
    PRINT*," FE0=",FE0
    !QZ0I=1._8/QZ0I
    !    PQN=MPARTICLES*Q_PARTICLE

    WBULK=M_POLY*LVI+FQ*ALPHAFEBULK
  
    Print*," wbulk=", wbulk
    ZPLUS=M_PLUS
    ZMINUS=M_MINUS
    ZX=1.
    END  SUBROUTINE SET_CHARGE
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE EVOLVE_WAB(RHO,W,WP,PP,EPSWA,EPSWB,EPSP,ERRW,ERRP,RHONORM,errp2,errw02,errw12)
    REAL(8), DIMENSION(NT,0:1) :: RHO
    REAL(8), DIMENSION(NT,0:1) :: W,WP
    REAL(8), DIMENSION(NT) :: PP
    REAL(8) :: EPSWA,EPSWB,EPSP,ERRW,ERRP
    INTENT(IN) :: RHO,EPSWA,EPSWB,EPSP
    INTENT(INOUT) :: W,WP,PP
    INTEGER :: I
    REAL(8)::  XW0,XW1,XP,RHONORM,LRHO,LFQ
    REAL(8)::  XN0,XN1,XPX
    REAL(8)::  errp2,errw02,errw12,rha,mu
    INTENT(OUT) ::ERRW,ERRP,errp2,errw02,errw12

    XN0=0.;XN1=0.;XPX=0.;errp2=0.;errw02=0.;errw12=0.
    ERRW=0.; ERRP=0.;RHONORM=0.
    !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I,XW0,XW1,lrho,rha,XP,mu)
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
        rha=max(RHO(I,0),real(epsilon(1._4),8))
        XP=RHONORM*(RHO(I,1)+RHO(I,0))-1.
        
        ERRP=MAX(ERRP,ABS(XP))
        ERRP2=ERRP2+XP*XP
        PP(I)=PP(I)+EPSP*XP
        XW0 =RHO(I,1)*CHIABN-WP(I,0)
        errw02=errw02+XW0**2
        XW1 =CHIABN*rha-WP(I,1)
        errw12=errw12+XW1**2
        ERRW=MAX(ERRW,ABS(XW0),ABS(XW1))
        WP(I,0)=WP(I,0)+EPSWA*XW0
        WP(I,1)=WP(I,1)+EPSWB*XW1
        xn0=xn0+WP(I,0)
        xn1=xn1+WP(I,1)
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
    errp2=sqrt(errp2/nt)
    errw02=sqrt(errw02/nt)
    errw12=sqrt(errw12/nt)
    END SUBROUTINE evolve_WAB
    !-----------------------------------------------------------------------------------------------------

    SUBROUTINE SET_DELTA_PSI(X,BASE)
    REAL :: X
    REAL(8) :: BASE,XBASE
    INTENT(IN)::X,BASE
    OPTIONAL :: BASE
    SAVE :: XBASE
    IF(PRESENT(BASE)) XBASE=BASE
    DELTA_PSI=X*XBASE
    LAMBDA_PSI=1.-DELTA_PSI
    !Print*, "delta psi set at",DELTA_PSI, LAMBDA_PSI
    END SUBROUTINE SET_DELTA_PSI

    SUBROUTINE GET_TS_PB()
    REAL(8) ::RHP,PS !,QQ
    REAL(8) :: WI,ROIP,ROIM
    INTEGER::I
    EPLUS=0.; EMINUS=0.
    RHON=0._8 ; EW=0.
    F2=HUGE(1.);F1=-F2
    RHO1=F1; RHO2=F2;  ENTLP=0.
    EAPLUS=0.;EAMINUS=0.
    ALPHAPS=0.; ALPHBETCS=0.
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,PS,RHP,ROIP,ROIM)
    !$OMP DO REDUCTION(MAX:F1,RHO1) REDUCTION(MIN:F2,RHO2) REDUCTION(+:RHON,EW,ENTLP)
    DO I=1,NT
        RHP=RHO(I,1)
        PS=PSI(I)
        RHON=RHON+RHP
        ENTLP=ENTLP+RHP*RHP
        RHO1=MAX(RHO1,RHP) ; RHO2=MIN(RHO2,RHP)
        F1=MAX(F1,PS) ; F2=MIN(F2,PS)
        EW=EW+RHP*W(I,1)

    ENDDO
    !$OMP END  DO
    !$OMP DO REDUCTION(+:EPLUS,EMINUS)
    DO I=1,NT
        ROIP=RHOI(0,I)
        ROIM=RHOI(1,I)
        EPLUS=EPLUS+XLNX(ROIP)
        EMINUS=EMINUS+XLNX(ROIM)
    ENDDO
    !$OMP END  DO
    !$OMP DO REDUCTION(+:ALPHAPS)
    DO I=1,NT
        ALPHAPS=ALPHAPS+(XLNX(ALPHAP(I))+XLNX(1.-ALPHAP(I))+ALPHAP(I)*LZX-alphaFEbulk)*RHO(I,1)
    ENDDO
    !$OMP END  DO

    !$OMP  DO  REDUCTION(+:ESE2)
    DO I=1,NT
        ESE2=ESE2+E2(I)/EPS_I(I)
    ENDDO
    !$OMP END DO
    !$OMP END  PARALLEL

    EW=EW*LVN
    EPLUS=-NTI*EPLUS/QPLUS
    EMINUS=-NTI*EMINUS/QMINUS
    ENTLP=0.5*LVN*ENTLP
    RHON=LVN*RHON

    EAPLUS=EAPLUS*NTI*MIPLUS/QPLUS
    EAMINUS=EAMINUS*NTI*MIMINUS/QMINUS
    ALPHAPS=ALPHAPS*LVN*FQ
    ALPHBETCS=ALPHBETCS*LVN
    ESE2=0.5*LVN*ESE2
    RETURN
    END SUBROUTINE GET_TS_PB

    

   

    END PROGRAM PB_2D_1

