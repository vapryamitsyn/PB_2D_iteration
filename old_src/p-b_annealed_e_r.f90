    MODULE PBSOLVER
    USE BOX
    USE MYFFT
    USE VECTORS3D
    USE INOUT
    USE PARTICLES
    USE DIFFUSION
    USE PARTITION
    IMPLICIT NONE
    REAL,PUBLIC    ::  IxQ2(NT2,3)
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
    REAL(8) :: FV,ERR,EPI !,PSX
    REAL(8) :: DUP,DDN,ZUP,ZDN
    REAL(8), PARAMETER :: DZ=.001, MZ=1.-DZ
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
        DO J=1, 100
            DQ=DQF(ZX)
            IF(DQ<-1.D-8)THEN
                DDN=DQ ; ZDN=ZX
            ELSEIF(DQ>1.D-8) THEN 
                ZUP=ZX; DUP=DQ
            ELSE ; EXIT ;        ENDIF
                IF((DUP==0.).OR.(DDN==0.))THEN             
                    ZX=DEXP(MZ*LOG(ZX)+DZ*DQ)
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
                EPS_I(I)=EPSCI(QP(I),MPOLY,RHO(I))
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
            call gCONVOLUTE(RHOQ,IxQ2,D)
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,epi) 
            DO I=1,NT
                EPI=EPS_I(I)
                D(I,1)=EPI*D(I,1)
                D(I,2)=EPI*D(I,2)
                D(I,3)=EPI*D(I,3)
            ENDDO
            !$OMP END PARALLEL DO
            call dCONVOLUTE(D,IxQ2,PSIX)
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





    !
    !SUBROUTINE MAKE_IQ2(Q2I,LAMBDA) !cubic cell only 
    !REAL,INTENT(INOUT) :: Q2I(0:NX/2,0:NY-1,0:NZ-1)
    !REAL(8),INTENT(IN):: LAMBDA
    !REAL(8) :: L2
    !INTEGER :: IX,JY,KZ
    !L2=0.25*LAMBDA*LX*LX/PI2
    !!$OMP PARALLEL DO DEFAULT(shared) private(ix,jy,kz)
    !DO KZ=0,NZ-1
    !    DO JY=0,NY-1
    !        DO IX=0,NX/2
    !            Q2I(IX,JY,KZ)=L2/AMAX0(1,RD2(KZ,NZ)+RD2(JY,NY)+IX**2)
    !        ENDDO
    !    ENDDO
    !ENDDO 
    !!$OMP END PARALLEL DO
    !Q2I(0,0,0)=0.
    !END SUBROUTINE MAKE_IQ2


    SUBROUTINE MAKE_XYZQ2(XQ2)
    REAL,INTENT(INOUT) :: XQ2(0:NX/2,0:NY-1,0:NZ-1,3)
    REAL(8) :: Q2,L1X,L1Y,L1Z,X2,Y2,Z2,VX,VY,VZ
    INTEGER :: IX,JY,KZ
    L1X=2*PI/LX
    L1Y=2*PI/LY
    L1Z=2*PI/LZ
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(IX,JY,KZ,X2,Y2,Z2,VY,VZ,Q2)
    DO KZ=0,NZ-1
        VZ=L1Z*ZGRAD(KZ,NZ/2)
        Z2=(L1Z*MIN(KZ,NZ-KZ))**2
        DO JY=0,NY-1
            VY=L1Y*ZGRAD(JY,NY/2)
            Y2=(L1Y*MIN(JY,NY-JY))**2
            DO IX=0,NX/2
                VX=L1X*ZGRAD(IX,NX/2)
                Q2=(L1X*IX)**2+Y2+Z2
                !IF(Q2>0.)THEN
                IF((IX/=NX/2).AND.(JY/=NY/2).AND.(KZ/=NZ/2).AND.(Q2>0.0))THEN
                    Q2=1./Q2
                    XQ2(IX,JY,KZ,1)=Q2*VX
                    XQ2(IX,JY,KZ,2)=Q2*VY
                    XQ2(IX,JY,KZ,3)=Q2*VZ
                ELSE
                    XQ2(IX,JY,KZ,1)=0.
                    XQ2(IX,JY,KZ,2)=0.
                    XQ2(IX,JY,KZ,3)=0.
                END IF

            ENDDO
        ENDDO
    ENDDO 
    !$OMP END PARALLEL DO
    CONTAINS 
    INTEGER ELEMENTAL FUNCTION ZGRAD(I,NH)
    INTEGER,INTENT(IN)::I,NH
    IF(I<NH) THEN
        ZGRAD=I
    ELSEIF(I==NH) THEN
        ZGRAD=0
    ELSE
        ZGRAD =I-2*NH
    END IF
    RETURN
    END FUNCTION ZGRAD
    END SUBROUTINE MAKE_XYZQ2



    SUBROUTINE MAKE_IQ2(XQ2)
    REAL,INTENT(INOUT) :: XQ2(0:NX/2,0:NY-1,0:NZ-1)
    REAL(8) :: Q2,L1X,L1Y,L1Z,X2,Y2,Z2
    INTEGER :: IX,JY,KZ
    L1X=2*PI/LX
    L1Y=2*PI/LY
    L1Z=2*PI/LZ
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(IX,JY,KZ,X2,Y2,Z2)
    DO KZ=0,NZ-1
        Z2=(L1Z*MIN(KZ,NZ-KZ))**2
        DO JY=0,NY-1
            Y2=(L1Y*MIN(JY,NY-JY))**2
            DO IX=0,NX/2
                Q2=(L1X*IX)**2+Y2+Z2
                IF(Q2>0.0) THEN
                    XQ2(IX,JY,KZ)=1./Q2
                ELSE
                    XQ2(IX,JY,KZ)=0.
                END IF

            ENDDO
        ENDDO
    ENDDO 
    !$OMP END PARALLEL DO
    END SUBROUTINE MAKE_IQ2

    ELEMENTAL REAL(8) FUNCTION EPSCI(RHO_P,MP,RHO_C)
    REAL(8),INTENT(IN) ::MP,RHO_P,RHO_C 
    EPSCI=EPS_WI+MP*EPS_CP*RHO_P+EPS_CC*RHO_C    
    END FUNCTION EPSCI 



    END MODULE PBSOLVER
