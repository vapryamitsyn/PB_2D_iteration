    MODULE PBSOLVER
    USE BOX
    USE MYFFT
    USE VECTORS3D
    USE INOUT
    USE DIFFUSION
    USE PARTITION
    IMPLICIT NONE
    !REAL,PUBLIC    ::  IXQ2(NT2,3)
    REAL,PUBLIC    ::  IQ2(NT2)
    REAL(8),PUBLIC :: PSI(NT),E2(NT)
    REAL(8),PUBLIC :: RHOIPLUS(NT),RHOIMINUS(NT),RHOQ(NT)
    REAL(8),PUBLIC :: ALPHAP(NT),ALPHAC(NT),BETAC(NT) !,BETAP(NT)
    REAL(8),PUBLIC :: EPS_I(NT)
    REAL(8),PUBLIC :: E(NT,3)
    REAL(8),PUBLIC :: EPS_WI,EPS_CP,EPS_CC
    REAL(8),PUBLIC :: PKA_P,PKA_C,PKB_C,APN,BPN,ZPA,ZCA,ZCB,PH,PLX
    REAL(8),PUBLIC :: QC,QPX
    REAL(8),PRIVATE :: EXPSI(NT)
    REAL(8),PUBLIC :: ZX
    REAL(8),PUBLIC :: EsE2

    REAL(8),PUBLIC :: NORMPLUS,NORMIMINUS

    REAL(8),PUBLIC ::LAMBDA_PSI,DELTA_PSI,LBX !,MINDEX  ! MINDEX=(1/(1+m2))

    CONTAINS


    SUBROUTINE ITERATE_PB(N_IT,DQ,DELTAPHI_MAX,ESE)
    INTEGER::I,N_IT,J,K
    INTENT(IN) ::N_IT
    REAL(8),INTENT(OUT) ::DQ,DELTAPHI_MAX,ESE
    OPTIONAL:: DELTAPHI_MAX,ESE
    REAL(8) :: RH,QMIN,QPL,PSX
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
            RHOIMINUS(I)=RH        
            RH=FV/EXPSI(I)        
            QPL=QPL+RH
            RHOIPLUS(I)=RH
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
        FIP=M2*RHO0*MPOLY
        !FIP=M2*RHO0*MPOLY*LVI
        !$OMP PARALLEL DEFAULT (SHARED)  PRIVATE(I,J)
        !$OMP  DO 
        DO J=1,NT
            EPS_I(J)=EPSCI(fip*RHO(J),(1.-IPHIP(J)))
        ENDDO
        !$OMP END DO NOWAIT
        !$OMP  DO 
        DO I=1,NT
            RHOQ(I)=NORMPLUS*RHOIPLUS(I)-NORMIMINUS*RHOIMINUS(I)-MPOLY*M2*ALPHAP(I)*RHO(I)
        ENDDO
        !$OMP END DO
        !$OMP END PARALLEL

        ! INVERSE LAPLACIAN
        call CONVOLUTE(RHOQ,IQ2,EXPSI)
        call PdivGRAD(E,EXPSI,EPS_I)
        CALL CONVOLUTE(EXPSI,IQ2,E2)
        ! END INVERSE LAPLACIAN

        IF (PRESENT(DELTAPHI_MAX).and.(K==N_IT)) THEN        !
            DELTAPHI_MAX=0.
            ESE=0.   ;EsE2=0.
            !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I,PSX)
            !$OMP  DO  REDUCTION(MAX:DELTAPHI_MAX) 
            DO I=1,NT
                PSX=E2(I)
                DELTAPHI_MAX=MAX(DELTAPHI_MAX,ABS(PSI(I)-PSX))
                PSI(I)=LAMBDA_PSI*PSI(I)+DELTA_PSI*PSX
            ENDDO
            !$OMP END  DO
            !$OMP  DO  REDUCTION(+:ESE)
            DO I=1,NT
                ESE=ESE+E2(I)*RHOQ(I)
            ENDDO
            !$OMP END DO
            !$OMP  DO  
            DO I=1,NT
                E2(I)=E(I,1)**2+E(I,2)**2+E(I,3)**2
            ENDDO
            !$OMP END DO    

            !$OMP END PARALLEL 
            ESE=0.5*LVN*ESE


            RETURN        
        ELSE
            !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I) 
            !$OMP  DO  
            DO I=1,NT
                PSI(I)=LAMBDA_PSI*PSI(I)+DELTA_PSI*E2(I)
            ENDDO
            !$OMP END  DO
            !$OMP END PARALLEL
        ENDIF
    END DO MAIN 

    !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I) 
    !$OMP  DO  
    DO I=1,NT
        E2(I)=E(I,1)**2+E(I,2)**2+E(I,3)**2
    ENDDO
    !$OMP END DO 
    !$OMP END PARALLEL
    RETURN
    CONTAINS
    REAL(8) FUNCTION DQF(ZX) RESULT(DQ)
    REAL(8), INTENT(IN) ::ZX
    REAL(8) :: EXL
    INTEGER ::I,J
    DQ=0.; QPX=0.
    !$OMP PARALLEL DEFAULT (SHARED)  PRIVATE(I,J)
    !$OMP DO   
    DO I=1,NT
        EXL=EXPSI(I)*PLX
        ALPHAP(I)=EXL/(EXL+ZPA*ZX)
        ALPHAC(I)=EXL/(EXL+ZCA*ZX)
        BETAC(I)=ZX/(EXL*ZCB+ZX)
    END DO
    !$OMP END DO
    !$OMP DO    REDUCTION(+:QPX)  
    DO I=1,NT
        QPX=QPX+M2*ALPHAP(I)*RHO(I)
    END DO
    !$OMP END DO NOWAIT

    !$OMP END PARALLEL
    QPX=QPX*LVN
    QC=LVN*DQ
    MPOLY=ZPOLY*QZ*ZX**QPX
    MIMINUS=ZMINUS*QMINUS*ZX
    MIPLUS=ZPLUS*QPLUS/ZX
    DQ=QC-MPOLY*QPX+MIPLUS-MIMINUS
    END FUNCTION DQF
    END SUBROUTINE ITERATE_PB

    SUBROUTINE PdivGRAD(E,PSI,epsi)
    REAL(8),INTENT(INOUT) :: PSI(0:NX-1,0:NY-1,0:NZ-1)
    REAL(8),INTENT(OUT) :: E(0:NX-1,0:NY-1,0:NZ-1,3)
    REAL(8),INTENT(IN) :: epsi(0:NX-1,0:NY-1,0:NZ-1)
    INTEGER :: X,Y,Z,UP,DN,X1,Y1,Z1,UP1,DN1,X2,Y2,Z2,UP2,DN2
    !$OMP PARALLEL  DEFAULT(SHARED)  PRIVATE(X,Y,Z,UP,DN,X1,Y1,Z1,UP1,DN1,X2,Y2,Z2,UP2,DN2)
    !$OMP  DO  
    DO Z=0,NZ-1
        DO Y=0,NY-1
            DO X=0,NX-1
                CALL UD(X,NX-1,DN,UP)
                E(X,Y,Z,1)=D2XI*epsi(X,Y,Z)*(PSI(UP,Y,Z)-PSI(DN,Y,Z))
            ENDDO
        ENDDO
    ENDDO 
    !$OMP END DO NOWAIT

    !$OMP  DO 
    DO Z1=0,NZ-1
        DO Y1=0,NY-1
            CALL UD(Y1,NY-1,DN1,UP1)
            DO X1=0,NX-1            
                E(X1,Y1,Z1,2)=D2YI*epsi(X1,Y1,Z1)*(PSI(X1,UP1,Z1)-PSI(X1,DN1,Z1))
            ENDDO
        ENDDO
    ENDDO 
    !$OMP END DO NOWAIT
    !$OMP  DO
    DO Z2=0,NZ-1
        CALL UD(Z2,NZ-1,DN2,UP2)
        DO Y2=0,NY-1
            DO X2=0,NX-1
                E(X2,Y2,Z2,3)=D2ZI*epsi(X2,Y2,Z2)*(PSI(X2,Y2,UP2)-PSI(X2,Y2,DN2))
            ENDDO
        ENDDO
    ENDDO 
    !$OMP END DO
    !end grad

    !$OMP  DO
    DO Z=0,NZ-1
        DO Y=0,NY-1
            DO X=0,NX-1
                CALL UD(X,NX-1,DN,UP)
                PSI(X,Y,Z)=-D2XI*(E(UP,Y,Z,1)-E(DN,Y,Z,1))
            ENDDO
        ENDDO
    ENDDO 
    !$OMP END DO 

    !$OMP  DO 
    DO Z1=0,NZ-1
        DO Y1=0,NY-1
            CALL UD(Y1,NY-1,DN1,UP1)
            DO X1=0,NX-1            
                PSI(X1,Y1,Z1)=PSI(X1,Y1,Z1)-D2YI*(E(X1,UP1,Z1,2)-E(X1,DN1,Z1,2))
            ENDDO
        ENDDO
    ENDDO 
    !$OMP END DO 

    !$OMP  DO  
    DO Z2=0,NZ-1
        CALL UD(Z2,NZ-1,DN2,UP2)
        DO Y2=0,NY-1
            DO X2=0,NX-1
                PSI(X2,Y2,Z2)=PSI(X2,Y2,Z2)-D2ZI*(E(X2,Y2,UP2,3)-E(X2,Y2,DN2,3))
            ENDDO
        ENDDO
    ENDDO 
    !$OMP END DO
    !$OMP END PARALLEL 
    RETURN
    CONTAINS 
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
    END SUBROUTINE PdivGRAD

    SUBROUTINE MAKE_IQ2(XQ2)
    REAL,INTENT(INOUT) :: XQ2(0:NX/2,0:NY-1,0:NZ-1)
    REAL(8) :: Q2,L1X,L1Y,L1Z,Y2,Z2
    INTEGER :: X,Y,Z
    L1X=2*PI/LX
    L1Y=2*PI/LY
    L1Z=2*PI/LZ
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(X,Y,Z,Y2,Z2,Q2)
    DO Z=0,NZ-1
        Z2=(L1Z*MIN(Z,NZ-Z))**2
        DO Y=0,NY-1
            Y2=(L1Y*MIN(Y,NY-Y))**2
            DO X=0,NX/2
                Q2=(L1X*X)**2+Y2+Z2
                IF(Q2>0.0) THEN
                    XQ2(X,Y,Z)=1./Q2
                ELSE
                    XQ2(X,Y,Z)=0.
                END IF

            ENDDO
        ENDDO
    ENDDO 
    !$OMP END PARALLEL DO
    END SUBROUTINE MAKE_IQ2

    ELEMENTAL REAL(8) FUNCTION EPSCI(FI_P,RHO_C)
    REAL(8),INTENT(IN) ::FI_P,RHO_C 
    EPSCI=LBX*EPS_WI*(1.+EPS_CP*FI_P+EPS_CC*RHO_C)    
    RETURN
    END FUNCTION EPSCI 



    END MODULE PBSOLVER
