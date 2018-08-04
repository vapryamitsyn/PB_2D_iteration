    MODULE LAPLACIAN
    USE  BOX2D
    IMPLICIT NONE
    REAL(8),PRIVATE, PARAMETER :: PI=3.1415926535897932384626433833
    PUBLIC:: GRAD_PSI2,GPDOTGPHI,MAKE_IQ2
    PRIVATE:: UD
    CONTAINS
    !------------------------------------------------------------------------------------------------------------
    SUBROUTINE GRAD_PSI2(PSI,E2) !DIFFERENT FOR 2D AND 3D , call from inside OMP region
    REAL(8),INTENT(IN) :: PSI(0:NX-1,0:NY-1)
    REAL(8),INTENT(OUT) :: E2(0:NX-1,0:NY-1)
    REAL(8):: EX
    INTEGER :: X,Y,UP,DN
    !$OMP  DO PRIVATE(X,Y,UP,DN,EX) ! SHARED(PSI,E2)
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
    END SUBROUTINE GRAD_PSI2
    !------------------------------------------------------------------------------------------------------------
    SUBROUTINE GPDOTGPHI(EPSI,PSI,GEGP) !DIFFERENT FOR 2D AND 3D, call from inside OMP region
    REAL(8),INTENT(IN) :: PSI(0:NX-1,0:NY-1)
    REAL(8),INTENT(IN) :: EPSI(0:NX-1,0:NY-1)
    REAL(8),INTENT(OUT) :: GEGP(0:NX-1,0:NY-1)
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
    REAL(8) :: Q2,L1X,L1Y,Y2
    INTEGER :: X,Y
    L1X=2*PI/LX
    L1Y=2*PI/LY
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(X,Y,Y2,Q2)
    DO Y=0,NY-1
        Y2=(L1Y*MIN(Y,NY-Y))**2
        DO X=0,NX/2
            Q2=(L1X*X)**2+Y2
            IF(Q2>0.0) THEN
                XQ2(X,Y)=REAL(1._8/Q2)
            ELSE
                XQ2(X,Y)=0.
            END IF
        ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    END SUBROUTINE MAKE_IQ2
    END MODULE LAPLACIAN