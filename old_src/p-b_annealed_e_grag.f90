    MODULE PBSOLVER
    USE BOX
    USE MYFFT
    USE VECTORS3D
    USE INOUT
    USE PARTICLES
    USE DIFFUSION
    USE PARTITION
    IMPLICIT NONE
    !REAL,PUBLIC    ::  IXQ2(NT2,3)
    REAL,PUBLIC    ::  IQ2(NT2)
    REAL(8),PUBLIC :: PSI(NT),PSIX(NT),E2(NT)
    REAL(8),PUBLIC :: RHOIPLUS(NT),RHOIMINUS(NT),RHOQ(NT)
    REAL(8),PUBLIC :: ALPHAP(NT),ALPHAC(NT),BETAC(NT) !,BETAP(NT)
    REAL(8),PUBLIC :: EPS_I(NT)
    REAL(8),PUBLIC :: D(NT,3)
    REAL(8),PUBLIC :: EPS_WI,EPS_CP,EPS_CC
    REAL(8),PUBLIC :: PKA_P,PKA_C,PKB_C,APN,BPN,ZPA,ZCA,ZCB
    REAL(8),PUBLIC :: QC,QPX


    real(8),private :: EXPSI(NT)
    REAL(8),PUBLIC :: ZX

    REAL(8),PUBLIC :: NORMPLUS,NORMIMINUS

    REAL(8),PUBLIC ::LAMBDA_PSI,DELTA_PSI,LBX !,MINDEX  ! MINDEX=(1/(1+m2))

    CONTAINS


    SUBROUTINE ITERATE_PB(N_IT,DQ,DELTAPHI_MAX,ESE)
    INTEGER::I,N_IT,J,K
    INTENT(IN) ::N_IT
    REAL(8),INTENT(OUT) ::DQ,DELTAPHI_MAX,ESE
    OPTIONAL:: DELTAPHI_MAX,ESE
    REAL(8) :: RH,QMIN,QPL,PSX
    REAL(8) :: FV,EPI !,PSX
    REAL(8) :: DUP,DDN,ZUP,ZDN
    REAL(8), PARAMETER :: DZ=.01, MZ=1.-DZ,DZx=.000001, MZx=1._8-DZx
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
        DUP=0.;DDN=0.
        DO J=1, 1000
            DQ=DQF(ZX)
            IF(DQ<-1.D-8)THEN
                DDN=DQ ; ZDN=ZX
            ELSEIF(DQ>1.D-8) THEN 
                ZUP=ZX; DUP=DQ
            ELSE ; EXIT ;        ENDIF
                IF((DUP==0.).OR.(DDN==0.))THEN
                    IF(ABS(DQ)>10.)THEN
                        ZX=DEXP(MZX*LOG(ZX)+DZX*DQ)
                    ELSE
                        ZX=DEXP(MZ*LOG(ZX)+DZ*DQ)
                    END IF

                ELSEIF(ABS(DQ)>10.)THEN
                    ZX=0.5*(ZUP+ZDN)
                ELSE            
                    ZX=ZUP-DUP*(ZUP-ZDN)/(DUP-DDN)    
                END IF

            END DO
            !end electroneutrality
            NORMIMINUS=LVI*MIMINUS/QMINUS    
            NORMPLUS=LVI*MIPLUS/QPLUS

            !$OMP PARALLEL DEFAULT (SHARED)  PRIVATE(I)
            !$OMP  DO 
            DO I=1,NT
                EPS_I(I)=EPSCI(RHO(I),B2*MPOLY,(1.-IPHIP(I)))
            ENDDO
            !$OMP END DO NOWAIT
            !$OMP  DO 
            DO I=1,NT
                RHOQ(I)=NORMPLUS*RHOIPLUS(I)-NORMIMINUS*RHOIMINUS(I)-MPOLY*M2*ALPHAP(I)*RHO(I)+bpn*BETAC(I)*QP(I)-apn*ALPHAC(I)*QP(I)
            ENDDO
            !$OMP END DO
            !$OMP END PARALLEL

            ! INVERSE LAPLACIAN

            !CALL CONVOLUTE(EXPW,IQ2,PSIX)
            call CONVOLUTE(RHOQ,IQ2,PSIX)
            !$OMP PARALLEL  DEFAULT(SHARED) 
            call PGRAD(D,PSIX)


            !$OMP DO PRIVATE(I,epi) 
            DO I=1,NT
                EPI=-EPS_I(I) ! minus sign is to compensate minus in Laplacian=-q^2
                D(I,1)=EPI*D(I,1)
                D(I,2)=EPI*D(I,2)
                D(I,3)=EPI*D(I,3)
            ENDDO
            !$OMP END  DO
            CALL PDIVX(EXPSI,D(:,1))
            CALL PDIVY(EXPSI,D(:,2))
            CALL PDIVZ(EXPSI,D(:,3))
            !$OMP END PARALLEL 
            CALL CONVOLUTE(EXPSI,IQ2,PSIX)
            ! END INVERSE LAPLACIAN

            IF (PRESENT(DELTAPHI_MAX).and.(K==N_IT)) THEN        !
                DELTAPHI_MAX=0.
                ESE=0.   
                !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I,PSX)
                !$OMP  DO  REDUCTION(MAX:DELTAPHI_MAX) 
                DO I=1,NT
                    PSX=PSIX(I)
                    DELTAPHI_MAX=MAX(DELTAPHI_MAX,ABS(PSI(I)-PSX))
                    PSI(I)=LAMBDA_PSI*PSI(I)+DELTA_PSI*PSX
                ENDDO
                !$OMP END  DO
                !$OMP  DO  REDUCTION(+:ESE)
                DO I=1,NT
                    ESE=ESE+PSIX(I)*RHOQ(I)
                ENDDO
                !$OMP END DO
                !$OMP END PARALLEL 
                ESE=0.5*LVN*ESE
                RETURN        
            ELSE
                !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
                DO I=1,NT
                    PSI(I)=LAMBDA_PSI*PSI(I)+DELTA_PSI*PSIX(I)
                ENDDO
                !$OMP END PARALLEL DO
            ENDIF
        END DO MAIN 
        RETURN
    CONTAINS
    REAL(8) FUNCTION DQF(ZX) RESULT(DQ)
    REAL(8), INTENT(IN) ::ZX
    DQ=0.; QPX=0.
    !$OMP PARALLEL DEFAULT (SHARED)  PRIVATE(I,FV,RH)
    !$OMP DO   
    DO I=1,NT
        ALPHAP(I)=EXPSI(I)/(EXPSI(I)+ZPA*ZX)
        ALPHAC(I)=EXPSI(I)/(EXPSI(I)+ZCA*ZX)
        BETAC(I)=ZX/(EXPSI(I)*ZCB+ZX)
    END DO
    !$OMP END DO
    !$OMP DO    REDUCTION(+:QPX)  
    DO I=1,NT
        QPX=QPX+M2*ALPHAP(I)*RHO(I)
    END DO
    !$OMP END DO NOWAIT
    !$OMP DO    REDUCTION(+:DQ)  
    DO I=1,NT            
        DQ=DQ+(BPN*BETAC(I)-APN*ALPHAC(I))*QP(I)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    QPX=QPX*LVN
    QC=LVN*DQ
    MPOLY=ZPOLY*QZ*ZX**QPX
    MIMINUS=ZMINUS*QMINUS*ZX
    MIPLUS=ZPLUS*QPLUS/ZX
    DQ=QC-MPOLY*QPX+MIPLUS-MIMINUS
    END FUNCTION DQF
    END SUBROUTINE ITERATE_PB

    SUBROUTINE PdivGRAD(E,PSI)
    REAL(8),INTENT(INOUT) :: PSI(0:NX-1,0:NY-1,0:NZ-1)
    REAL(8),INTENT(OUT) :: E(0:NX-1,0:NY-1,0:NZ-1,3)
    INTEGER :: X,Y,Z,UP,DN,X1,Y1,Z1,UP1,DN1,X2,Y2,Z2,UP2,DN2
    !$OMP  DO  PRIVATE(X,Y,Z,UP,DN)
    DO Z=0,NZ-1
        DO Y=0,NY-1
            DO X=0,NX-1
                CALL UD(X,NX-1,DN,UP)
                E(X,Y,Z,1)=D2XI*(PSI(UP,Y,Z)-PSI(DN,Y,Z))
            ENDDO
        ENDDO
    ENDDO 
    !$OMP END DO NOWAIT

    !$OMP  DO PRIVATE(X1,Y1,Z1,UP1,DN1)
    DO Z1=0,NZ-1
        DO Y1=0,NY-1
            CALL UD(Y1,NY-1,DN1,UP1)
            DO X1=0,NX-1            
                E(X1,Y1,Z1,2)=D2YI*(PSI(X1,UP1,Z1)-PSI(X1,DN1,Z1))
            ENDDO
        ENDDO
    ENDDO 
    !$OMP END DO NOWAIT
    !$OMP  DO PRIVATE(X2,Y2,Z2,UP2,DN2)
    DO Z2=0,NZ-1
        CALL UD(Z2,NZ-1,DN2,UP2)
        DO Y2=0,NY-1
            DO X2=0,NX-1
                E(X2,Y2,Z2,3)=D2ZI*(PSI(X2,Y2,UP2)-PSI(X2,Y2,DN2))
            ENDDO
        ENDDO
    ENDDO 
    !$OMP END DO

    return
    END SUBROUTINE PdivGRAD

    SUBROUTINE PDIVX(E,PSI)
    REAL(8),INTENT(IN) :: PSI(0:NX-1,0:NY-1,0:NZ-1)
    REAL(8),INTENT(OUT) :: E(0:NX-1,0:NY-1,0:NZ-1)
    INTEGER :: X,Y,Z,UP,DN

    !$OMP  DO PRIVATE(X,Y,Z,UP,DN)
    DO Z=0,NZ-1
        DO Y=0,NY-1
            DO X=0,NX-1
                CALL UD(X,NX-1,DN,UP)
                E(X,Y,Z)=D2XI*(PSI(UP,Y,Z)-PSI(DN,Y,Z))
            ENDDO
        ENDDO
    ENDDO 
    !$OMP END DO 
    RETURN
    END SUBROUTINE PDIVX

    SUBROUTINE PDIVY(E,PSI)
    REAL(8),INTENT(IN) :: PSI(0:NX-1,0:NY-1,0:NZ-1)
    REAL(8),INTENT(INOUT) :: E(0:NX-1,0:NY-1,0:NZ-1)
    INTEGER :: X1,Y1,Z1,UP1,DN1
    !$OMP  DO PRIVATE(X1,Y1,Z1,UP1,DN1)
    DO Z1=0,NZ-1
        DO Y1=0,NY-1
            CALL UD(Y1,NY-1,DN1,UP1)
            DO X1=0,NX-1            
                E(X1,Y1,Z1)=E(X1,Y1,Z1)+D2YI*(PSI(X1,UP1,Z1)-PSI(X1,DN1,Z1))
            ENDDO
        ENDDO
    ENDDO 
    !$OMP END DO 
    RETURN
    END SUBROUTINE PDIVY


    SUBROUTINE PDIVZ(E,PSI)
    REAL(8),INTENT(IN) :: PSI(0:NX-1,0:NY-1,0:NZ-1)
    REAL(8),INTENT(INOUT) :: E(0:NX-1,0:NY-1,0:NZ-1)
    INTEGER :: X2,Y2,Z2,UP2,DN2
    !$OMP  DO  PRIVATE(X2,Y2,Z2,UP2,DN2)
    DO Z2=0,NZ-1
        CALL UD(Z2,NZ-1,DN2,UP2)
        DO Y2=0,NY-1
            DO X2=0,NX-1
                E(X2,Y2,Z2)=E(X2,Y2,Z2)+D2ZI*(PSI(X2,Y2,UP2)-PSI(X2,Y2,DN2))
            ENDDO
        ENDDO
    ENDDO 
    !$OMP END DO
    RETURN
    END SUBROUTINE PDIVZ

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





    SUBROUTINE MAKE_IQ2(XQ2)
    REAL,INTENT(INOUT) :: XQ2(0:NX/2,0:NY-1,0:NZ-1)
    REAL(8) :: Q2,L1X,L1Y,L1Z,Y2,Z2
    INTEGER :: X,Y,Z
    L1X=2*PI/LX
    L1Y=2*PI/LY
    L1Z=2*PI/LZ
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(X,Y,Z,Y2,Z2)
    DO Z=0,NZ-1
        Z2=(L1Z*MIN(Z,NZ-Z))**2
        DO Y=0,NY-14
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

    ELEMENTAL REAL(8) FUNCTION EPSCI(RHO_P,MP,RHO_C)
    REAL(8),INTENT(IN) ::MP,RHO_P,RHO_C 
    EPSCI=lbx*(EPS_WI+MP*EPS_CP*RHO_P+EPS_CC*RHO_C)    
    END FUNCTION EPSCI 



    END MODULE PBSOLVER
