    module density
    USE OMP_LIB
    USE BOX2D
    USE MYFFT2D
    !USE PARAM3d
    USE INOUT2d
    IMPLICIT NONE
    !REAL(8),PRIVATE,ALLOCATABLE,TARGET :: QD(:,:)
    REAL,ALLOCATABLE,TARGET :: Q(:,:) !PRIVATE?
    REAL,PUBLIC::    RG,RGA,RGB,RGAB !,KH   
    INTEGER,PUBLIC::    NS,NSA,NSB,NSAB !,KH   
    REAL,PRIVATE ::QQOUT(NT), QQpr(NT),EXP_W(NT,0:N_comp-1)
    REAL(8),PUBLIC::QH(0:N_COMP-1)
    REAL(8),PRIVATE ::RH(NT)
    REAL,PRIVATE :: EXPDELTA(NT2)
    REAL, PRIVATE :: DS,DSA
    

    CONTAINS
    !DIFFERENT FOR 2D AND 3D==================================================================================================================================
    SUBROUTINE MAKE_EXPDELTA(ED,Rg,DS)   
    REAL,INTENT(IN) ::Rg,DS
    REAL,INTENT(OUT) :: ED(0:NX/2,0:NY-1)
    INTEGER :: I,J
    REAL(8) ::A,X2,Y2
    A=-4*Rg**2*PI2*DS
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,X2,Y2)    
        DO J=0,NY-1
            Y2=(LYI*MIN(J,NY-J))**2
            DO I=0,NX/2
                X2=(LXI*I)**2
                ED(I,J)=DEXP(A*(X2+Y2))
            END DO
        END DO
    !$OMP END PARALLEL DO
    END SUBROUTINE MAKE_EXPDELTA
    !END of DIFFERENT FOR 2D AND 3D================================================================================================================================
    
    SUBROUTINE SET_DS(NSX,RGA,RGB,RGAB) !Setting global dS, NS and Exp[Delta dS ; Rg should be the Rg of the largest chain in the system (for polydisperse systems)
    INTEGER,INTENT(IN)::NSX
    REAL,INTENT(IN) ::RGA,RGB,RGAB
    NS=NSX !NS must be an even number!
    !PRINT*,"NS=",NS
    RG=MAX(RGA,RGB,RGAB)
    NSA=2*NINT(0.5*REAL(NSX)*(RGA/RG)**2)
    NSB=2*NINT(0.5*REAL(NSX)*(RGB/RG)**2)
    NSAB=2*NINT(0.5*REAL(NSX)*(RGAB/RG)**2)
    NS=max(NSA,NSB,NSAB) !NS must be an even number!
    !PRINT*,"NS=",NS
    DS=1./REAL(NS) ; DSA=-0.5_8*DS
    Print*,"NX=",NX,"NY=",NY,"NS=",NS
    CALL MAKE_EXPDELTA(EXPDELTA,Rg,DS)
    IF(ALLOCATED(Q)) DEALLOCATE(Q)
    ALLOCATE(Q(NT,0:MAX0(NSAB,NSA/2,NSB/2)))
    END SUBROUTINE SET_DS



    SUBROUTINE MAKE_EXP_W_AB(W,EXP_W) !FIX FOR ABC copolymer
    REAL,INTENT(IN) ::W(N_comp*NT)
    REAL,INTENT(OUT) ::EXP_W(N_comp*NT)
    INTEGER:: I
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
    DO I = 1,N_comp*NT
        EXP_W(I)=EXP(DSA*W(I))
    END DO
    !$OMP  END PARALLEL DO
    END SUBROUTINE MAKE_EXP_W_AB



    SUBROUTINE GET_AB_DENSITY(EXP_W,RHO,t_m,Z,F,NB,NS,dmp) !Both NS and NB should be even numbers
    LOGICAL,OPTIONAL,INTENT(IN) :: DMP
    INTEGER,INTENT(INOUT):: NS,NB
    INTEGER(1)::T_M(NS)
    REAL,DIMENSION(NT,0:N_comp-1) :: RHO,EXP_W
    INTENT(OUT)::Z,RHO
    REAL(8) ::F,Z,Zdag
    INTENT(IN) ::F,T_M,EXP_W


    !CALL MAKE_EXP_W_AB(W,EXP_W)
    CALL GET_COPOLYMER_PARTITION(EXP_W,T_M,Z,ZDAG,NS)
    If(Abs(z/zdag-1.)>.0001) print*, "Z=",z," Z^dag=",zdag, " err=",1.-Z/ZDAG
    Z=0.5*(z+zdag)
    call SIMPSON(rho(1:NT,0),Q(1:NT,0:NB),NB,F/Z)
    call SIMPSON(rho(1:NT,1),Q(1:NT,NB:NS),NS-NB,(1.-F)/Z)
    IF (PRESENT(dmp)) then
        if(dmp) call dump_big2d(q,nt*(ns+1),[nx,ny],'mdens.dat')
    end if
  
    CONTAINS

    !computing polymer density  density Simpson's rule
    SUBROUTINE SIMPSON(RHO,QQ,N,F) !HERE Q IS LOCAL
    INTEGER::N,I,K
    REAL(8)::  NORM,F
    REAL:: RHO(NT)
    REAL:: QQ(NT,0:N)
    INTENT(IN)::QQ,F,N
    INTENT(OUT)::RHO
    NORM=F/REAL(3*N,8)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
    DO I = 1,NT
        RH(I)=QQ(I,0)
    END DO
    !$OMP  END PARALLEL DO
    DO K=1,N-3,2
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
        !$OMP  DO
        DO I = 1,NT
            RH(I)=RH(I)+4*QQ(I,K)
        END DO
        !$OMP  end DO
        !$OMP  DO
        DO I = 1,NT
            RH(I)=RH(I)+2*QQ(I,K+1)
        END DO
        !$OMP end  DO
        !$OMP  END PARALLEL
    END DO
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
    !$OMP  DO
    DO I = 1,NT
        RH(I)=RH(I)+4*QQ(I,N-1)
    END DO
    !$OMP  end DO
    !$OMP  DO
    DO I = 1,NT
        RH(I)=NORM*(RH(I)+QQ(I,N)) !LAST STEP AND NORMALIZATION
    END DO
    !$OMP  end DO
    !$OMP WORKSHARE
    rho=rh
    !$OMP END WORKSHARE
    !$OMP  END  PARALLEL
    END SUBROUTINE SIMPSON
    end SUBROUTINE GET_AB_DENSITY

    SUBROUTINE GET_COPOLYMER_PARTITION(EXP_W,T_M,Z,ZDAG,NS) ! here local NS<=Ns global, subroutie returns Q(x,s)=q(x,s)q^dag(x,N-s), quality test: z=zdag
    INTEGER,INTENT(IN):: NS
    INTEGER(1)::T_M(NS)
    REAL(8),INTENT(OUT)::Z,ZDAG
    REAL,INTENT(IN)::EXP_W(NT,0:N_comp-1)
    REAL(8):: Qint
    INTEGER:: I,S,IS,IT,IT0

    !compute partition function q
    it=t_m(1)  !
    CALL CONVOLUTE(EXP_W(:,it),EXPDELTA,QQOUT) !!Pre-computing and convoluting Q(:,1)->qqout
    DO S = 1, Ns-1
        it0=it ; it=t_m(S+1)
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Qint)
        DO i = 1,NT
            Qint= EXP_W(i,it0)*QQOUT(i)  !!Computing Q(I,S)
            Q(I,S)=QINT                   !Storing Q(I,S)
            QQpr(i)=EXP_W(i,it)*Qint       !!Pre-computing Q(I,S+1)
        END DO
        !$OMP END PARALLEL DO
        CALL CONVOLUTE(QQpr,EXPDELTA,QQOUT) ! !!convoluting Q(:,S+1)
    END DO
    Z=0.
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Qint) REDUCTION(+:Z)
    DO i = 1,NT
        Qint=EXP_W(i,it)*QQOUT(i) !!finalizing and Computing Q(I,NS) and Z
        Q(i,NS) = Qint
        Z=Z+Qint
    END DO
    !$OMP END PARALLEL DO


    ! calculate Qdagger

    CALL CONVOLUTE(EXP_W(:,it),EXPDELTA,QQOUT) !!convoluting QQpr(:,1) it=t_m(NS)

    DO IS = 1,NS-1 !backward do loop
        S=NS-IS
        it0=it
        it=t_m(S)
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Qint)
        DO i = 1,NT
            Qint= EXP_W(i,it0)*QQOUT(i)
            Q(I,S)=Qint*Q(I,S) !!!!Computing Q(i,S)*Qdagger(i,NS-S) and storing it in Q(I,S)
            QQpr(i)=EXP_W(i,it)*Qint       !!Pre-computing Qdagger(:,S)
        END DO
        !$OMP END PARALLEL DO
        CALL CONVOLUTE(QQpr,EXPDELTA,QQOUT) ! calling fft
    END DO

    Zdag=0.
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Qint) REDUCTION(+:Zdag)
    DO i = 1,NT
        Qint=EXP_W(i,it)*QQOUT(i) !!finalizing and Computing Qdagger(NS) and Zdag
        Q(I,0) = Qint !! sending Qdagger(I,NS) to Q(I,0)
        Zdag=Zdag+Qint
    END DO
    !$OMP END PARALLEL DO
    z=z*nti
    Zdag=Zdag*nti
    END SUBROUTINE GET_COPOLYMER_PARTITION


    SUBROUTINE HOMOPOLYMER_DENSITY(W,RHO,Z,NS,dmp) !Both NS and NB should be even numbers
    LOGICAL,OPTIONAL,INTENT(IN) :: DMP
    INTEGER:: NS, HNS
    REAL,DIMENSION(NT,0:N_comp-1) :: RHO,W
    INTENT(OUT)::Z,RHO
    REAL(8) :: Z(0:N_COMP-1)
    INTENT(IN) :: W,NS
    CALL MAKE_EXP_W_AB(W,EXP_W)

    CALL GET_HOMO_DENSITY(EXP_W(:,0),RHO(:,0),Z(0),NSA)
    CALL GET_HOMO_DENSITY(EXP_W(:,1),RHO(:,1),Z(1),NSB)


    IF (PRESENT(DMP)) THEN
        if(dmp) call dump_big2d(q,nt*(ns+1),[nx,ny],'mdens.dat')
    END IF

    END SUBROUTINE HOMOPOLYMER_DENSITY

    SUBROUTINE GET_HOMO_DENSITY(EXP_W,RHO,Z,NS) !Both NS and NB should be even numbers
    INTEGER:: NS, HNS
    REAL,DIMENSION(NT) :: RHO,EXP_W
    INTENT(OUT)::Z,RHO
    REAL(8) ::Z
    INTENT(IN) :: EXP_W,NS
    HNS=NS/2
    !CALL MAKE_EXP_W_AB(W,EXP_W)
    Call GET_HOMOPOLYMER_PARTITION(EXP_W,Z,HNS)
    call SIMPSONH(rho(1:NT),Q(1:NT,0:HNS),HNS,Z)
    CONTAINS
    !computing polymer density  density Simpson's rule
    SUBROUTINE SIMPSONH(RHO,QQ,N,Z) !HERE Q IS LOCAL
    INTEGER::N,I,K,HN
    REAL(8)::  NORM,F,LAST,Z
    REAL:: RHO(NT)
    REAL:: QQ(NT,0:N)
    INTENT(IN)::QQ,N,Z
    INTENT(OUT)::RHO
    NORM=1._8/(Z*3.*N)
    HN=N/2
    IF(2*hn==n)then
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
        DO I = 1,NT
            RH(I)=2*QQ(I,0)
        END DO
        !$OMP  END PARALLEL DO
        DO K=2,N-2,2
            !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
            !$OMP  DO
            DO I = 1,NT
                RH(I)=RH(I)+8*QQ(I,K-1)
            END DO
            !$OMP  end DO
            !$OMP  DO
            DO I = 1,NT
                RH(I)=RH(I)+4*QQ(I,K)
            END DO
            !$OMP end  DO
            !$OMP  END PARALLEL
        END DO
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
        !$OMP  DO
        DO I = 1,NT
            RH(I)=RH(I)+8*QQ(I,N-1)
        END DO
        !$OMP  end DO
        !$OMP  DO
        DO I = 1,NT
            RH(I)=NORM*(RH(I)+2*QQ(I,N)) !LAST STEP AND NORMALIZATION
        END DO
        !$OMP  end DO
        !$OMP  END  PARALLEL
    ELSE
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
        DO I = 1,NT
            RH(I)=2*QQ(I,0)
        END DO
        !$OMP  END PARALLEL DO
        DO K=2,N-1,2
            !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
            !$OMP  DO
            DO I = 1,NT
                RH(I)=RH(I)+8*QQ(I,K-1)
            END DO
            !$OMP  end DO
            !$OMP  DO
            DO I = 1,NT
                RH(I)=RH(I)+4*QQ(I,K)
            END DO
            !$OMP end  DO
            !$OMP  END PARALLEL
        END DO
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I)
        !$OMP  DO
        DO I = 1,NT
            RH(I)=NORM*(RH(I)+4*QQ(I,N)) !LAST STEP AND NORMALIZATION
        END DO
        !$OMP  END DO
        !$OMP WORKSHARE
        RHO=RH
        !$OMP END WORKSHARE
        !$OMP  END  PARALLEL
    END IF
    END SUBROUTINE SIMPSONH
    end SUBROUTINE GET_HOMO_DENSITY


    SUBROUTINE GET_HOMOPOLYMER_PARTITION(EXP_W,Z,HNS) ! NS should be an even numbershere local NS<=Ns global, subroutie returns Q(x,s)=q(x,s)q^dag(x,Ns-s), 0<S<N/2
    REAL,INTENT(in) :: EXP_W(NT)
    INTEGER,INTENT(IN):: HNS
    INTEGER(1)::T
    REAL(8),INTENT(OUT)::Z
    REAL(8):: Qint
    INTEGER:: I,S,IS,NS
    NS=2*HNS

    !compute partition function q
    CALL CONVOLUTE(EXP_W,EXPDELTA,QQOUT) !!Pre-computing and convoluting Q(:,1)->qqout
    DO S = 1, HNS-1
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Qint)
        DO i = 1,NT
            Qint= EXP_W(i)*QQOUT(i)  !!Computing Q(I,S)
            Q(I,S)=QINT                   !Storing Q(I,S)
            QQpr(i)=EXP_W(i)*Qint       !!Pre-computing Q(I,S+1)
        END DO
        !$OMP END PARALLEL DO
        CALL CONVOLUTE(QQpr,EXPDELTA,QQOUT) ! !!convoluting Q(:,S+1)
    END DO
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Qint)
    DO i = 1,NT
        Qint= EXP_W(i)*QQOUT(i)  !!Computing Q(I,S)
        Q(I,HNS)=QINT*QINT                   !Storing q(I,HNS)*q(I,HNS)
        QQpr(i)=EXP_W(i)*Qint       !!Pre-computing Q(I,S+1)
    END DO
    !$OMP END PARALLEL DO
    CALL CONVOLUTE(QQpr,EXPDELTA,QQOUT) ! !!convoluting Q(:,S+1)
    DO S =HNS+1,NS-1
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Qint)
        DO i = 1,NT
            Qint= EXP_W(i)*QQOUT(i)  !!Computing Q(I,S)
            Q(I,NS-S)=QINT*Q(I,NS-S)        !Storing q(I,S)*q(I,NS-S) in Q(I,NS-S)
            QQpr(i)=EXP_W(i)*Qint       !!Pre-computing Q(I,S+1)
        END DO
        !$OMP END PARALLEL DO
        CALL CONVOLUTE(QQpr,EXPDELTA,QQOUT) ! !!convoluting Q(:,S+1)
    END DO
    Z=0.
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Qint) REDUCTION(+:Z)
    DO i = 1,NT
        Qint=EXP_W(i)*QQOUT(i) !!finalizing and Computing Q(I,NS) and Z
        Q(i,0) = Qint
        Z=Z+Qint
    END DO
    !$OMP END PARALLEL DO



    END SUBROUTINE GET_HOMOPOLYMER_PARTITION

    subroutine dump_big2d(q,n,d,fn)
    integer,intent(in) :: n,d(2)
    character(*),intent(in) ::fn
    real,intent(in) :: q(n)
    integer :: I
    OPEN(UNIT=299,FILE=fn,FORM ='FORMATTED',ACTION='WRITE',STATUS='REPLACE',ASYNCHRONOUS='YES',BUFFERED='YES')
    WRITE(299,'(3I22)')  d
    WRITE(299,'(G29.22)')  (q(i),i=1,n)
    CLOSE(299)
    end subroutine dump_big2d

    end module density
