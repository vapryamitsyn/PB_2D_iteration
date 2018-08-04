    MODULE PBSOLVER2D
    USE BOX2D
    USE MYFFT2D
    USE VECTORS3D
    USE INOUT2d
    IMPLICIT NONE
    REAL(8),PUBLIC:: FQ  !FQ - NUMBER OF CHARGEABLE MONOMERS PER "A" CHAIN
    REAL(8), PUBLIC ::QZ !,QZ0I !SINGLE CHAIN PARTITION FUNCTION
    REAL(8),PUBLIC:: QPLUS,QMINUS,M_PLUS,M_MINUS, M_IONS ! SINGLE ION PARTITION FUNCTIONS AND REFERENCE CONCENTRATIONS: M_IONS=M_PLU+M_MINUS
    REAL(8),PUBLIC ::    MIPLUS, MIMINUS ! INSTANT NUMBER OF POSITIVE AND NEGATIVE SALT IONS IN THE SYSTEM
    REAL(8),PUBLIC ::    ZPLUS, ZMINUS ! ACTIVITIES OF POLYMERS, AND IONS
    REAL(8),PUBLIC :: RBplus,RBMnus

    REAL,PUBLIC    ::  IQ2(NT2)
    REAL,DIMENSION(NT),PUBLIC :: PSI, RHOQ, E2, ALPHAP, EPS_I, RHOIPLUS, RHOIMINUS, EXPSI, GEGP, LEPSILON
    REAL,DIMENSION(0:1,NT),PUBLIC :: RHOI
    REAL,PUBLIC :: EPS_AI,EPS_BI,EPS_IMI,EPS_IPI !invere epsilons
    REAL(8),PUBLIC :: PKA_P,ZPA
    REAL(8),PUBLIC :: QC,QPX

    REAL(8),PUBLIC :: ZX
    REAL(8),PUBLIC :: EsE2

    REAL(8),PUBLIC :: NORMPLUS,NORMIMINUS

    REAL(8),PUBLIC ::LAMBDA_PSI,DELTA_PSI,LBX !,MINDEX  ! MINDEX=(1/(1+FQ))

    CONTAINS


    SUBROUTINE ITERATE_PB(RHO,N_IT,DQ,DELTAPHI_MAX,ESE)
    REAL,INTENT(IN)::RHO(NT,0:1) !densitie of the ionomer and the neutral polymer
    INTEGER::I,N_IT,J,K
    INTENT(IN) ::N_IT
    REAL(8),INTENT(OUT) ::DQ,DELTAPHI_MAX,ESE
    OPTIONAL:: DELTAPHI_MAX,ESE
    REAL(8) :: RH,QMIN,QPL,PSX,DLPT
    REAL(8) :: FV !,EPI !,PSX
    REAL(8) :: DUP,DDN,ZUP,ZDN,fip
    REAL(8), PARAMETER :: DZ=.01, MZ=1.-DZ,DZx=.001, MZx=1._8-DZx
    MAIN:  DO, K=1,N_IT
        QPL=0.; QMIN=0.
        !$OMP PARALLEL DO  DEFAULT (SHARED)  PRIVATE(I,FV,RH) REDUCTION(+:QMIN,QPL)
        DO I=1,NT
            EXPSI(i)=DEXP(REAL(PSI(I),8))
            FV=IPHIP(I)
            RH=FV*EXPSI(i)
            QMIN=QMIN+RH
            RHOI(1,I)=RH
            RH=FV/EXPSI(I)
            QPL=QPL+RH
            RHOI(0,I)=RH
        END DO
        !$OMP END PARALLEL DO
        QPLUS=QPL*NTI
        QMINUS=QMIN*NTI

        !electroneutrality
        DUP=0.;DDN=0.;zup=10. ;ZDN=.001
        zx_iter: DO J=1, 1000

            DQ=DQF(ZX)
            IF(DQ<-1.D-8)THEN
                ZDN=ZX ; DDN=DQ
            ELSEIF(DQ>1.D-8) THEN
                ZUP=ZX; DUP=DQ
            ELSE
                EXIT zx_iter
            ENDIF

            IF(ABS(ZUP-ZDN)<0.00000001_8) EXIT ZX_ITER

            IF((DUP==0.).OR.(DDN==0.))THEN
                IF(ABS(DQ)>10.)THEN
                    ZX=DEXP(MZX*LOG(ZX)+DZX*MAX(MIN(DQ,100.),-100.))
                ELSE
                    ZX=DEXP(MZ*LOG(ZX)+DZ*DQ)
                END IF

            ELSEIF(ABS(DQ)>1.)THEN
                ZX=0.5*(ZUP+ZDN)
            ELSE
                ZX=(DDN*ZUP-DUP*ZDN)/(DDN-DUP)
            END IF


        END DO zx_iter
        !end electroneutrality
        NORMIMINUS=LVI*MIMINUS/QMINUS
        NORMPLUS=LVI*MIPLUS/QPLUS
        FIP=FQ*RHO0
       
        !$OMP PARALLEL DEFAULT (SHARED)  PRIVATE(i,j,k)
        !$OMP  DO
        DO J=1,NT
            EPS_I(J)=EPSCI(RHO(J,0),RHO(J,1),RHOI(0,J), RHOI(1,J))
            LEPSILON(J)=-LOG(EPS_I(J))
        ENDDO
        !$OMP END DO NOWAIT
        !$OMP  DO
        DO I=1,NT
            RHOQ(I)=NORMPLUS*RHOI(0,I)-NORMIMINUS*RHOI(1,I)-FQ*ALPHAP(I)*RHO(I,0)
        ENDDO
        !$OMP END DO
        call GPDOTGPHI(LEPSILON,PSI,GEGP)
        !$OMP  DO
        DO I=1,NT
            EXPSI(I)=RHOQ(I)*EPS_I(I)+GEGP(I)
        ENDDO
        !$OMP END DO
        !$OMP END PARALLEL

        ! INVERSE LAPLACIAN
        CALL CONVOLUTE(EXPSI,IQ2,E2)
        ! END INVERSE LAPLACIAN

        IF (PRESENT(DELTAPHI_MAX).and.(K==N_IT)) THEN        !
            DELTAPHI_MAX=0.
            ESE=0.   ;EsE2=0.
            !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I,PSX)
            !$OMP  DO  REDUCTION(MAX:DELTAPHI_MAX)
            DO I=1,NT
                PSX=E2(I)-PSI(I)
                DELTAPHI_MAX=MAX(DELTAPHI_MAX,ABS(PSX))
                PSI(I)=PSI(I)+DELTA_PSI*PSX
            ENDDO
            !$OMP END  DO
            !$OMP  DO  REDUCTION(+:ESE)
            DO I=1,NT
                ESE=ESE+E2(I)*RHOQ(I)
            ENDDO
            !$OMP END DO
            call GRAD_PSI2(PSI,E2)
            !$OMP END PARALLEL
            ESE=0.5*LVN*ESE


            RETURN
        ELSE
            !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I)
            !$OMP  DO
            DO I=1,NT
                PSI(I)=PSI(I)+DELTA_PSI*(E2(I)-PSI(I))
            ENDDO
            !$OMP END  DO
            !$OMP END PARALLEL
        ENDIF
    END DO MAIN


    call GRAD_PSI2(PSI,E2)

    RETURN
    CONTAINS
    REAL(8) FUNCTION DQF(ZX) RESULT(DQ)
    REAL(8), INTENT(IN) ::ZX
    REAL(8) :: EXL
    INTEGER ::I,J
    DQ=0.; QPX=0.
    !$OMP PARALLEL DEFAULT (SHARED)  PRIVATE(i,j,k)
    !$OMP DO
    DO I=1,NT
        EXL=EXPSI(I)
        ALPHAP(I)=EXL/(EXL+ZPA*ZX)


    END DO
    !$OMP END DO
    !$OMP DO    REDUCTION(+:QPX)
    DO I=1,NT
        QPX=QPX+FQ*ALPHAP(I)*RHO(I,0)
    END DO
    !$OMP END DO NOWAIT

    !$OMP END PARALLEL
    QPX=QPX*LVN
    QC=LVN*DQ
    MIMINUS=ZMINUS*QMINUS*ZX
    MIPLUS=ZPLUS*QPLUS/ZX
    DQ=QC-QPX+MIPLUS-MIMINUS
    END FUNCTION DQF
    END SUBROUTINE ITERATE_PB
    !------------------------------------------------------------------------------------------------------------
    SUBROUTINE GRAD_PSI2(PSI,E2) !DIFFERENT FOR 2D AND 3D 
    REAL,INTENT(IN) :: PSI(0:NX-1,0:NY-1)
    REAL,INTENT(OUT) :: E2(0:NX-1,0:NY-1)
    REAL(8):: EX
    INTEGER :: X,Y,UP,DN
    !$OMP PARALLEL  DEFAULT(SHARED)  PRIVATE(X,Y,UP,DN,EX)
    !$OMP  DO
    DO Y=0,NY-1
        DO X=0,NX-1
            CALL UD(X,NX-1,DN,UP)
            EX=D2XI*(PSI(UP,Y)-PSI(DN,Y))
            E2(X,Y)=EX*EX
        ENDDO
    ENDDO
    !$OMP END DO
    !$OMP  DO
    DO Y=0,NY-1
        CALL UD(1,NY-1,DN,UP)
        DO X=0,NX-1
            EX=D2YI*(PSI(X,UP)-PSI(X,DN))
            E2(X,Y)=E2(X,Y)+EX*EX
        ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL

    !END GRAD
    END SUBROUTINE GRAD_PSI2
    !------------------------------------------------------------------------------------------------------------
    SUBROUTINE GPDOTGPHI(EPSI,PSI,GEGP) !DIFFERENT FOR 2D AND 3D 
    REAL,INTENT(IN) :: PSI(0:NX-1,0:NY-1)
    REAL,INTENT(IN) :: EPSI(0:NX-1,0:NY-1)
    REAL,INTENT(OUT) :: GEGP(0:NX-1,0:NY-1)
    REAL(8):: EX
    INTEGER :: X,Y,UP,DN
    !$OMP  DO  PRIVATE(X,Y,UP,DN,EX)
    DO Y=0,NY-1
        DO X=0,NX-1
            CALL UD(X,NX-1,DN,UP)
            EX=D2XI*(PSI(UP,Y)-PSI(DN,Y))
            GEGP(X,Y)=EX*D2XI*(EPSI(UP,Y)-EPSI(DN,Y))
        ENDDO
    ENDDO
    !$OMP END DO
    !$OMP  DO  PRIVATE(X,Y,UP,DN,EX)
    DO Y=0,NY-1
        CALL UD(1,NY-1,DN,UP)
        DO X=0,NX-1
            EX=D2YI*(PSI(X,UP)-PSI(X,DN))
            GEGP(X,Y)=GEGP(X,Y)+EX*D2YI*(EPSI(X,UP)-EPSI(X,DN))
        ENDDO
    ENDDO

    !$OMP END DO
    !END GRAD
    END SUBROUTINE GPDOTGPHI
    !------------------------------------------------------------------------------------------------------------
    ELEMENTAL SUBROUTINE UD(I,N,DN,UP)
    INTEGER, INTENT(IN)::I,N
    INTEGER, INTENT(OUT)::UP,DN
    IF(I==0) THEN
        DN=N ; UP=1
    ELSE IF(I==N) THEN
        DN=N-1;UP=0
    ELSE
        DN=I-1;UP=I+1
    END IF
    RETURN
    END SUBROUTINE UD
    !------------------------------------------------------------------------------------------------------------
    SUBROUTINE MAKE_IQ2(XQ2) !DIFFERENT FOR 2D AND 3D 
    REAL,INTENT(INOUT) :: XQ2(0:NX/2,0:NY-1)
    REAL(8) :: Q2,L1X,L1Y,y2
    INTEGER :: X,Y
    L1X=2*PI/LX
    L1Y=2*PI/LY
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(X,Y,Y2,Q2)
    DO Y=0,NY-1
        Y2=(L1Y*MIN(Y,NY-Y))**2
        DO X=0,NX/2
            Q2=(L1X*X)**2+Y2
            IF(Q2>0.0) THEN
                XQ2(X,Y)=1./Q2
            ELSE
                XQ2(X,Y)=0.
            END IF
        ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    END SUBROUTINE MAKE_IQ2

    ELEMENTAL REAL FUNCTION EPSCI(RHOA,RHOB,ROIP,ROIM)
    REAL,INTENT(IN) ::RHOA,RHOB,ROIP,ROIM
    EPSCI=EPS_AI*RHOA+EPS_BI*RHOB+EPS_IPI*ROIP+EPS_IMI*ROIM
    RETURN
    END FUNCTION EPSCI
    END MODULE PBSOLVER2D
