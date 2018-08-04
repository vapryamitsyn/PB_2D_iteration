    MODULE INIT_DENS
    USE BOX2D
    IMPLICIT NONE
    PUBLIC::MAKE_INIT
    CONTAINS
    SUBROUTINE MAKE_INIT(RHO,R)
    REAL, INTENT(OUT):: RHO(0:NX-1,0:NY-1,0:1)
    REAL, INTENT(IN):: R
    INTEGER :: I,J
    real :: X,Y
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,X,Y)
    DO I = 0,NX-1
        X=LX*MIN(I,ABS(NX-I))/REAL(NX)
        DO J=0,NY-1
            Y=LY*MIN(J,ABS(NY-J))/REAL(NY)
            IF(((X-0.5*LX)**2 + (Y-0.5*LY)**2)<=R*R.OR.(X**2 + Y**2)<=R*R) THEN
                RHO(I,J,0)=1;RHO(I,J,1)=0
            ELSE
                RHO(I,J,0)=0.;RHO(I,J,1)=1.
            END IF
        END DO
    END DO
    !$OMP END PARALLEL DO
    END SUBROUTINE MAKE_INIT

    END MODULE INIT_DENS