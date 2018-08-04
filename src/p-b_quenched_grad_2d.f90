    MODULE PBSOLVER2D
    USE BOX2D
    USE MYFFT2D
    USE VECTORS3D
    USE INOUT2d
    use LAPLACIAN
    IMPLICIT NONE
    REAL,PUBLIC:: FQ, FQA  !FQ - NUMBER OF CHARGEABLE MONOMERS PER a reference "A" CHAIN with RGA=RG
    REAL,PUBLIC:: PQ_norm  !PQ_norm - normalization parameter between rohoA and the polymer charge density 
    !REAL(8),PUBLIC:: QI(0),QI(1) ! SINGLE ION PARTITION FUNCTIONS 
    REAL,PUBLIC:: M_PLUS,M_MINUS, M_IONS !IONS CONCENTRATIONS: M_IONS=M_PLU+M_MINUS
    !REAL(8),PUBLIC ::    M_PLUS, M_MINUS ! INSTANT NUMBER OF POSITIVE AND NEGATIVE SALT IONS IN THE SYSTEM
    REAL(8),PUBLIC ::    ZPLUS, ZMINUS ! ACTIVITIES OF POLYMERS, AND IONS
    REAL(8),PUBLIC ::    QI(0:1) !partition fuction of ions
    REAL,PUBLIC :: RBplus,RBMnus, IRP, IRM

    REAL,PUBLIC    ::  IQ2(NT2)
    REAL,DIMENSION(NT),private :: qxin,psiout
    REAL(8),DIMENSION(NT),PUBLIC :: PSI, RHOQ, E2,  EPS_I, RHOIPLUS, RHOIMINUS, EXPSI, GEGP, LEPSILON
    REAL(8),DIMENSION(0:1,NT),PUBLIC :: RHOI
    REAL,PUBLIC :: EPS_AI,EPS_BI!,EPS_IMI,EPS_IPI !invere epsilons
    REAL(8),PUBLIC :: NORMPLUS,NORM_MINUS

    REAL,PUBLIC ::LAMBDA_PSI,DELTA_PSI,LBX,BSP !,MINDEX  ! MINDEX=(1/(1+FQ))

    CONTAINS


    SUBROUTINE ITERATE_PB(RHO,N_IT,DELTAPHI_MAX,ESE,TSPL,TSMIN)
    REAL,INTENT(IN)::RHO(NT,0:1) !density of the ionomer and the neutral polymer
    INTEGER::I,N_IT,J,K
    INTENT(IN) ::N_IT
    REAL(8),INTENT(OUT) ::DELTAPHI_MAX,ESE,TSPL,TSMIN
    OPTIONAL:: DELTAPHI_MAX,ESE,TSPL,TSMIN
    REAL(8) :: RH,QMIN,QPL,PSX,DLPT
    REAL(8) :: FV !,EPI !,PSX
    REAL(8) :: DUP,DDN,ZUP,ZDN,fip
    
    !$OMP PARALLEL DO DEFAULT (SHARED)  PRIVATE(I,FV)
    DO I=1,NT
        FV=EPSCI(RHO(I,0),RHO(I,1))
        EPS_I(I)=FV
        LEPSILON(I)=-LOG(FV)
     ENDDO
    !$OMP END PARALLEL DO
    MAIN:  DO, K=1,N_IT
        QPL=0.; QMIN=0.
        !$OMP PARALLEL DEFAULT (SHARED)  PRIVATE(i,j,k,FV,RH)
        !$OMP  DO  REDUCTION(+:QMIN,QPL)
        DO I=1,NT
            RH=DEXP(-IRP*EPS_I(I)-PSI(I))
            QPL=QPL+RH
            RHOI(0,I)=RH
            RH=DEXP(-IRM*EPS_I(I)+PSI(I))
            QMIN=QMIN+RH
            RHOI(1,I)=RH
        END DO
        !$OMP END  DO
        !$OMP SINGLE
        QI(0)=LVN*QPL
        QI(1)=LVN*QMIN
        NORMPLUS=NT*LVI*M_PLUS/QPL
        NORM_MINUS=NT*LVI*M_MINUS/QMIN
        !$OMP END SINGLE
        !$OMP  DO
        DO I=1,NT
            RHOI(0,I)=NORMPLUS*RHOI(0,I)
            RHOI(1,I)=NORM_MINUS*RHOI(1,I)
            RHOQ(I)=RHOI(0,I)-RHOI(1,I)-PQ_norm*RHO(I,0)
        ENDDO
        !$OMP END DO
        call GPDOTGPHI(LEPSILON,PSI,GEGP)
        !$OMP  DO
        DO I=1,NT
            qxin(I)=RHOQ(I)*EPS_I(I)+GEGP(I)
        ENDDO
        !$OMP END DO
        !$OMP END PARALLEL

        ! INVERSE LAPLACIAN
        CALL CONVOLUTE(qxin,IQ2,PSIOUT)
        ! END INVERSE LAPLACIAN

        IF (PRESENT(DELTAPHI_MAX).and.(K==N_IT)) THEN        !
            DELTAPHI_MAX=0.
            ESE=0.;TSPL=0.;TSMIN=0.   
            !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I,PSX)
            !$OMP  DO  REDUCTION(MAX:DELTAPHI_MAX)
            DO I=1,NT
                PSX=psiout(I)-PSI(I)
                DELTAPHI_MAX=MAX(DELTAPHI_MAX,ABS(PSX))
                PSI(I)=PSI(I)+DELTA_PSI*PSX
            ENDDO
            !$OMP END  DO
            !$OMP  DO  REDUCTION(+:ESE)
            DO I=1,NT
                ESE=ESE+psiout(I)*RHOQ(I)
            ENDDO
            !$OMP END DO
            !$OMP  DO  REDUCTION(+:TSPL,TSMIN)
            DO I=1,NT
                TSPL =TSPL +RHOI(0,I)*(IRM*EPS_I(I)+PSI(I))
                TSMIN=TSMIN+RHOI(1,I)*(IRM*EPS_I(I)-PSI(I))
            ENDDO
            !$OMP END DO
            call GRAD_PSI2(PSI,E2)
            !$OMP END PARALLEL
            ESE=0.5*LVN*ESE
            TSPL=LVN*TSPL+M_PLUS*dlog(QI(0))
            TSMIN=LVN*TSMIN+M_MINUS*dlog(QI(1))
            RETURN
            ELSE if (K==N_IT) THEN 
            !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I)
            !$OMP  DO
            DO I=1,NT
                PSI(I)=PSI(I)+DELTA_PSI*(PSIOUT(I)-PSI(I))
            ENDDO
            !$OMP END  DO
                call GRAD_PSI2(PSI,E2)
            !$OMP END PARALLEL
            
        ELSE
            !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(I)
            !$OMP  DO
            DO I=1,NT
                PSI(I)=PSI(I)+DELTA_PSI*(PSIOUT(I)-PSI(I))
            ENDDO
            !$OMP END  DO
            !$OMP END PARALLEL
        ENDIF
    END DO MAIN
    RETURN
    END SUBROUTINE ITERATE_PB
    
    

    ELEMENTAL REAL FUNCTION EPSCI(RHOA,RHOB)
    REAL,INTENT(IN) ::RHOA,RHOB
    
    EPSCI=LBX*(EPS_AI*RHOA+EPS_BI*RHOB)/(RHOA+RHOB)
    RETURN
    END FUNCTION EPSCI
    END MODULE PBSOLVER2D

    !ELEMENTAL REAL FUNCTION EPSCI(RHOA,RHOB,ROIP,ROIM)
    !REAL,INTENT(IN) ::RHOA,RHOB,ROIP,ROIM
    !EPSCI=EPS_AI*RHOA+EPS_BI*RHOB+EPS_IPI*ROIP+EPS_IMI*ROIM
    !RETURN
    !END FUNCTION EPSCI
    !END MODULE PBSOLVER2D
