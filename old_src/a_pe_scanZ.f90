    !  a_pe_scanZ.f90
    !
    !  
    !
    !
    !****************************************************************************
    !
    !  PROGRAM: PE_GCA_annealed
    !
    !  PURPOSE:  Entry point for the console application a_pe_scanZ.
    ! scan partile charge and chain discretization
    ! non-restartable
    !****************************************************************************
    !

    PROGRAM  A_PE_GCA
    USE BOX
    USE MYFFT
    USE DIFFUSION
    USE PARTITION
    USE INOUT
    USE PBSOLVER    
    USE RAN_NUMBERS
    USE TIME
    IMPLICIT NONE
    INTEGER  ::  I,J,K,IMAX,PBIT ,ls               ! IMAX maximum number of SCF iterations in the current run
    REAL :: F1,F2,RHO1,RHO2,RHON
    LOGICAL :: SAVE_SNAP_RHO, SAVE_SNAP_PSI, SAVE_SNAP_P
    NAMELIST /STARTPOLYMER/   RG,M_POLY,B2,M2, M_PLUS, M_MINUS,RHO0

    NAMELIST /START/ IMAX,DELTA_W,DELTA_PSI,LBX,LX,PBIT, & !LBX is the 4*pi*Bjerrum_length
    SAVE_SNAP_RHO,SAVE_SNAP_PSI                 

   

    NAMELIST /EPSILONS/ EPS_WI,EPS_CP,EPS_CC   

    REAL(8) :: WX,DW,DELTAPHI_MAX,FE,TS,TSP,TSPLUS,TSMINUS,EW,EPLUS,EMINUS,ESE,ENTLP,QTEST,FE0,ALPHAPS,ALPHBETCS
    REAL(8) :: EAPLUS,EAMINUS,XADS,bg

    REAL(8) :: RESXX(27),L_D,Epsilion_Sol
    real(8) :: psir(0:nx/2-1),rhor(0:nx/2-1),rhoir(0:nx/2-1)
    real(8):: alphapbulk,wbulk,alphaFEbulk,LZX
    !integer::ix,iy,iz,ir

    CHARACTER(LEN=16):: HEADXX(28)=["Ns","Qc","M_poly","Qp","B2","M_+", "M_-","Mpoly","MI+","MI-","TSp","TS+","TS-","TS","Entlp","EsE","FE","Zx","Qz","Q+","Q-","Ea+","Ea-","Xads","F_DS_P","FDS_C","L_D","Epsilion_Sol"]

    CHARACTER(*), PARAMETER :: APFN="aparticles.lst", PPFN="ppositions.start", SFN="scenariop.lst" 
    CHARACTER(*), PARAMETER :: STARTFN="start.lst",POLYFN="startapolymer.lst",epFN="epsilons.lst" 

    INTEGER :: NDS0N,NDSI,NDSN !,IXQMAX,IXQI,IXQM ! initial, start, and 
    REAL(8) :: MPMAX,MPI,MPM ! initial, start, and 
    TYPE(RAN_SAVE):: XSEED


    NAMELIST /SCENARIOP/ NDS0N,NDSI,NDSN, MPMAX,MPI,MPM

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


    !CALL INIT_PSM_DIFF(RG*SQRT(6.))
    !CALL MAKE_XYZQ2(IxQ2)
    CALL MAKE_IQ2(IQ2)
    call SET_CHARGE()
    !ZPOLY=EXP(B2*M_POLY*LVI+ (XLNX(alphaPbulk)+XLNX(1.-alphaPbulk)+alphaPbulk*LOG(ZPA/PLX))*M2)*M_POLY
    !Print*,"ZPOLY/M_POLY=",ZPOLY/M_POLY
    !ZPLUS=M_PLUS
    !ZMINUS=M_MINUS
    !ZX=1.
    L_D=1/SQRT((M_PLUS+M_MINUS)*LVI*EPSCI(M2*RHO0*M_POLY*LVI,0._8))
    Epsilion_Sol=LBX/(EPSCI(M2*RHO0*M_POLY*LVI,0._8))
    Print '("Ions Debye-Huckel length is:",G11.3," Epsilon=",G11.3," fi=",G11.3)', &
    Real(L_D,4),REAL( Epsilion_Sol,4),REAL(M2*RHO0*M_POLY*LVI,4)





    print*, "Stage 3 ", "norm=",LVN*sum(rho)

    CALL SET_CPU_TIME()


    CALL INIT_DAT("results.tsv",HEAD=HEADXX,N=29)


    !DO WHILE(Q_PARTICLE<=IXQMAX)
    DO WHILE(M_POLY<=MPMAX)
        CALL SET_DS(NDS0N,RG*SQRT(6.))
        W(:)=0._8
        PSI(:)=0._8
        ZX=1.

        IF(M_POLY==0.)THEN
            RHO=LVI
            QZ=1.
        ELSE
            CALL HOMOPOLYMER_DENSITY(W,RHO,QZ,NS)
        END IF


        CALL SET_DELTA_PSI(.01)

        CALL ITERATE_PB(1000,QTEST,DELTAPHI_MAX,ESE)
        print*, "dPhi=", real(DELTAPHI_MAX,4),"dPSI=",real(DELTA_PSI,4),"Qt=",Real(QTEST,4),"Zx=",Real(ZX,4)

        DO LS=1,NDSN
            call SET_DS((NDSI*2**LS)/2,RG*SQRT(6.))   
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
                        WX=(B2*MPOLY*RHO(I)-M2*ALPHAP(I)*PSI(I)-0.5*LBX*EPS_WI*EPS_CP*M2*RHO0*E2(I)/(EPS_I(I)**2)+&
                            (XLNX(ALPHAP(I))+XLNX(1.-ALPHAP(I))+ALPHAP(I)*LZX)*M2)-WBULK
                        DW=MAX(DW,ABS(WX-W(I)))
                        W(I)=LAMBDA_W*W(I)+DELTA_W*WX
                    ENDDO
                    !$OMP END PARALLEL DO
                    CALL HOMOPOLYMER_DENSITY(W,RHO,QZ,NS)
                END IF

                CALL ITERATE_PB(max(1,(PBIT-1)/LS),QTEST)
                END DO
                CALL ITERATE_PB(max(1,(PBIT-1)/LS),QTEST,DELTAPHI_MAX,ESE)


                CALL GET_TS_PB()
                TSP=MPOLY*(1.+EW)
                TSPLUS=MIPLUS*(1.+EPLUS)
                TSMINUS=MIMINUS*(1.+EMINUS)
                TS=TSP+TSPLUS+TSMINUS-ALPHAPS-ALPHBETCS
                FE=(ENTLP+EsE+EAPLUS+EAMINUS-TS)-FE0
                print*, "EsEtest->",ese,ese2

                BG=RHO(NX*(1+NY*(1+NZ)/2))


                PRINT '("dPsi=",G9.3," dW=",G9.3,"DPsi=",G9.3," dRho=",G9.3," nmm=",G10.3," Zx=",G10.3," Qt=",G10.3,"spu:",2(I5,"c"))',& 
                DELTAPHI_MAX,dW,f1-f2,MPOLY*(RHO1-RHO2),1._8-RHON,ZX,QTEST,Nint(rtime()), Nint(utime()/nut)
                PRINT '("Qz=",G9.3," Mp=",G9.3," Q+=",G9.3," M+=",G9.3," Q-=",G9.3E3,"M-=",G10.3E3,"E+=",G10.3," E-=",G10.3," Ep=",G10.3)', & 
                QZ,MPOLY/Max(M_POLY,Epsilon(1._8)),Qplus,MIPLUS/M_PLUS,Qminus,MIMINUS/M_MINUS,EPLUS,EMINUS,EW
                PRINT '("Sp=",G10.3," S+=",G9.2," S-=",G9.2," S=",G10.3,"DP=",G11.3," DC=",G11.3)',  & 
                TSP,TSPLUS,TSMINUS,TS,-ALPHAPS,-ALPHBETCS
                PRINT '("EsE=",G10.3,"Entlp=",G10.3," Ea+=",G10.3," Ea-=",G10.3," F=",G13.6)',  & 
                ESE, ENTLP,EAPLUS,EAMINUS,FE

                IF(ABS(DW)<1.E-4.and.abs(DELTAPHI_MAX)<1.E-6) EXIT
            END DO
            XADS=0._8
            !$OMP PARALLEL  DO DEFAULT(SHARED) PRIVATE(I) REDUCTION(+:XADS)
            DO I=1,NT
                XADS=XADS+IPHIP(I)*(RHO(I)-BG)
            ENDDO
            !$OMP END PARALLEL DO   
            XADS=LVN*M2*MPOLY*XADS
            PRINT '("Xads=",G26.12)', XADS

            RESXX=[REAL(NS,8),M_POLY,Qpx,B2, &
            M_PLUS,M_MINUS,MPOLY,MIPLUS,MIMINUS,TSP,TSPLUS,TSMINUS,TS,ENTLP,ESE,FE,ZX,QZ,Qplus,Qminus,EAPLUS,EAMINUS,XADS, &
            ALPHAPS,ALPHBETCS,L_D,Epsilion_Sol]

            CALL APPEND_DAT("results.tsv",dat=RESXX,N=28)
        END DO



        IF(SAVE_SNAP_RHO)THEN
            CALL DUMP_DENS(RHO,"Mp"//NUMCHAR(INT(M_POLY))//"rho.dat")
        END IF




       
        IF(save_SNAP_psi)THEN
            CALL DUMP_DENS(PSI,"Mp"//NUMCHAR(Int(M_POLY))//"psi.dat")
            CALL DUMP_DENS(RHOIPLUS,"Mp"//NUMCHAR(Int(M_POLY))//"plus.dat")
            CALL DUMP_DENS(RHOIMINUS,"Mp"//NUMCHAR(Int(M_POLY))//"minus.dat")
            !$OMP PARALLEL  DO DEFAULT(SHARED) PRIVATE(I)      
            DO I=1,nt
                RHOIMINUS(I)=-MPOLY*M2*ALPHAP(I)*RHO(I)
            ENDDO
            !$OMP END PARALLEL DO 

            CALL DUMP_DENS(RHOIMINUS,"Mp"//NUMCHAR(Int(M_POLY))//"Qp.dat")
            
        END IF

        !IF(Q_PARTICLE==0.and.IXQI==0 ) THEN
        !    Q_PARTICLE=1
        !ELSE
        !    Q_PARTICLE=IXQM*Q_PARTICLE+IXQI
        !END IF    
        !
        !print*, "Qp=",Q_PARTICLE

        IF(M_POLY==0.AND.MPI==0 ) THEN
            M_POLY=MPM
        ELSE
            M_POLY=MPM*M_POLY+MPI
        END IF    

        print*, "M_POLY=",M_POLY
        L_D=1/SQRT((M_PLUS+M_MINUS)*LVI*EPSCI(M2*RHO0*M_POLY*LVI,0._8))
        Epsilion_Sol=LBX/(EPSCI(M2*RHO0*M_POLY*LVI,0._8))
        Print '("Ions Debye-Huckel length is:",G11.3," Epsilon=",G11.3," fi=",G11.3)', &
        Real(L_D,4),REAL( Epsilion_Sol,4),REAL(M2*RHO0*M_POLY*LVI,4)

        !CALL MAKE_IQ2(IQ2)
        CALL SET_CHARGE()


    ENDDO


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
        LX=100.
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
        READ(10,STARTpolymer)
        CLOSE(10)
        PRINT STARTpolymer
    ELSEIF(started) then
        STOP "started, but startpolymer.lst does not exist"
    ELSE
        RG=5.
        M_POLY=100
        B2=1.
        M2=10
        M_MINUS=M_POLY*M2/10
        M_PLUS=M_MINUS+M_POLY*M2
        PKA_P=5.
        pH=6.
    END IF
    STARTED=STARTED.AND.EXIST
    ZPA=10.**(-PKA_P)
    PLX=(-1. + 55.5*10.**pH)/(2.76585231e8 + Sqrt(7.649939000732336e16 + 55.5*10.**pH))
    !alphaPbulk=PLX/(PLX+ZPA)
    !Print*,"alphaPbulk=",alphaPbulk
    M_PLUS=M_POLY*M2*alphaPbulk+M_MINUS
    M_IONS=M_PLUS+M_MINUS
    !MINDEX=(1./(1.+M2))
    !QZ0I=EXP(B2*M_POLY*LVI)
    !alphaFEbulk=(XLNX(alphaPbulk)+XLNX(1.-alphaPbulk)+alphaPbulk*LZX)
    !FE0=-M_POLY*(1.+.5*B2*M_POLY*LVI)+ M_POLY*m2*alphaFEbulk- M_PLUS-M_MINUS
    !PRINT*,"QZ0=",QZ0I," FE0=",FE0
    PRINT*," FE0=",FE0,M_POLY*m2*alphaFEbulk,M_PLUS,M_MINUS

    INQUIRE (FILE=SFN,EXIST=exist)
    if(exist)then
        OPEN(UNIT=10,FILE=SFN,FORM ='FORMATTED', ACTION='READ',STATUS='UNKNOWN',DELIM='APOSTROPHE')
        !READ(10,SCENARIO)
        READ(10,SCENARIOP)
        CLOSE(10)
        PRINT START
    ELSE
        NDS0N=16
        NDSI=32
        NDSN=3
        MPMAX=M_POLY
        MPI=1.
        MPM=1.
        OPEN(UNIT=10,FILE=SFN,FORM ='FORMATTED', ACTION='WRITE',STATUS='UNKNOWN',DELIM='APOSTROPHE')
        !WRITE(10,SCENARIO)
        WRITE(10,SCENARIOP)
    END IF 


    OPEN(UNIT=10,FILE=POLYFN,FORM ='FORMATTED', ACTION='WRITE',STATUS='UNKNOWN',DELIM='APOSTROPHE')
    WRITE(10,STARTPOLYMER)
    CLOSE(10) 
    RETURN 
3   STOP 'STARTED FILE ERROR '
    END FUNCTION INIT_READ_POYMER

    SUBROUTINE SET_CHARGE()
    !real(8) :: alphaPbulk
    alphaPbulk=PLX/(PLX+ZPA)
    LZX=log(ZPA/PLX)
    Print*,"alphaPbulk=",alphaPbulk

    !M_PLUS=M_POLY*M2+M_MINUS
    M_PLUS=M_POLY*M2*alphaPbulk+M_MINUS
    M_IONS=M_PLUS+M_MINUS
    
    alphaFEbulk=(XLNX(alphaPbulk)+XLNX(1.-alphaPbulk)+alphaPbulk*LZX)
    
    !FE0=-M_POLY*(1.+.5*B2*M_POLY*LVI)+  M_POLY*m2*alphaFEbulk- M_PLUS-M_MINUS
    FE0=-M_POLY*(1.+.5*B2*M_POLY*LVI)- M_PLUS-M_MINUS
    !PRINT*,"QZ0=",QZ0I," FE0=",FE0
    PRINT*," FE0=",FE0
    !QZ0I=1._8/QZ0I
    !    PQN=MPARTICLES*Q_PARTICLE
    
    WBULK=B2*M_POLY*LVI+M2*ALPHAFEBULK 
    ZPOLY=M_POLY
    Print*," wbulk=", wbulk
    ZPLUS=M_PLUS
    ZMINUS=M_MINUS
    ZX=1.
    END  SUBROUTINE SET_CHARGE
    !-----------------------------------------------------------------------------------------------------


    
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
    REAL(8) :: ROIP,ROIM
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
        RHP=RHO(I)
        PS=PSI(I)       
        RHON=RHON+RHP
        ENTLP=ENTLP+RHP*RHP
        RHO1=MAX(RHO1,RHP) ; RHO2=MIN(RHO2,RHP)
        F1=MAX(F1,PS) ; F2=MIN(F2,PS)
        EW=EW+RHP*W(I)
        ENDDO
    !$OMP END  DO
    !$OMP DO REDUCTION(+:EPLUS,EMINUS,EAPLUS,EAMINUS) 
    DO I=1,NT
        ROIP=RHOIPLUS(I)
        ROIM=RHOIMINUS(I)
        EPLUS=EPLUS+XLNX(ROIP)
        
        EMINUS=EMINUS+XLNX(ROIM)
        
    ENDDO
    !$OMP END  DO 
    !$OMP DO REDUCTION(+:ALPHAPS)
    DO I=1,NT
        ALPHAPS=ALPHAPS+(XLNX(ALPHAP(I))+XLNX(1.-ALPHAP(I))+ALPHAP(I)*LZX-alphaFEbulk)*RHO(I)
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
    ENTLP=0.5*LVN*B2*MPOLY*MPOLY*ENTLP
    
    RHON=LVN*RHON

    EAPLUS=EAPLUS*NTI*MIPLUS/QPLUS
    EAMINUS=EAMINUS*NTI*MIMINUS/QMINUS
    ALPHAPS=ALPHAPS*LVN*MPOLY*M2
    ALPHBETCS=ALPHBETCS*LVN
    ESE2=0.5*LVN*ESE2
    RETURN
    END SUBROUTINE GET_TS_PB

    SUBROUTINE PSI_R(PSI,PSIR)
    REAL(8):: PSI(0:NX-1,0:NY-1,0:NZ-1),PSIR(0:NX/2-1),PSI0
    INTEGER :: I,IX,IY,IZ,PSIC(0:NX/2-1)
    PSIR=0._8
    PSIC=0
    PSI0=PSI(NX/2-1,NY/2-1,NZ/2-1)
    !$OMP PARALLEL  DO DEFAULT(SHARED) PRIVATE(I,IX,IY,IZ) REDUCTION(+:PSIR,PSIC)      
    DO IZ=0,NZ-1
        DO IY=0,NY-1
            DO IX=0,NX-1
                I=FLOOR(DSQRT(REAL(RD2(IX,NX)+RD2(IY,NY)+RD2(IZ,NZ),8)))
                IF(I<=NX/2-1) THEN
                    PSIR(I)=PSIR(I)+(PSI(IX,IY,IZ)-PSI0)
                    PSIC(I)=PSIC(I)+1
                END IF
            ENDDO 
        ENDDO 
    ENDDO
    !$OMP END PARALLEL DO 
    !PRINT*, "PSI0=",PSI0
    !PRINT*,PSIR
    !PRINT*,PSIC
    DO I=1,NX/2-1
        PSIR(I)=PSIR(I)/PSIC(I)
    END DO
    END SUBROUTINE PSI_R

    SUBROUTINE PSI_RX(PSI,PSIR)
    REAL(8):: PSI(0:NX-1,0:NY-1,0:NZ-1),PSIR(0:NX/2-1)
    INTEGER :: I,IX,IY,IZ,PSIC(0:NX/2-1)
    PSIR=0._8
    PSIC=0
    !$OMP PARALLEL  DO DEFAULT(SHARED) PRIVATE(I,IX,IY,IZ) REDUCTION(+:PSIR,PSIC)      
    DO IZ=0,NZ-1
        DO IY=0,NY-1
            DO IX=0,NX-1
                I=FLOOR(DSQRT(REAL(RD2(IX,NX)+RD2(IY,NY)+RD2(IZ,NZ),8)))
                IF(I<=NX/2-1) THEN
                    PSIR(I)=PSIR(I)+(PSI(IX,IY,IZ))
                    PSIC(I)=PSIC(I)+1
                END IF
            ENDDO 
        ENDDO 
    ENDDO
    !$OMP END PARALLEL DO 
    !PRINT*, "PSI0=",PSI0
    !PRINT*,PSIR
    !PRINT*,PSIC
    DO I=1,NX/2-1
        PSIR(I)=PSIR(I)/PSIC(I)
    END DO
    END SUBROUTINE PSI_RX

    END PROGRAM A_PE_GCA

