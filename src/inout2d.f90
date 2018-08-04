    MODULE INOUT2d
    USE BOX2d
    USE MYFFT2D
    IMPLICIT NONE
    CHARACTER(LEN=5), PUBLIC :: IN_EXT,OUTEXT
    CHARACTER(1),PARAMETER,PRIVATE :: TAB=CHAR(9)
    INTEGER,PUBLIC ::ITERATION
    PUBLIC :: SETEXT,DUMP_DENS,NUMCHAR,IEXT
    CONTAINS
    !===========================================================================================
    CHARACTER(LEN=5)FUNCTION IEXT(I)RESULT(IT)
    INTEGER,INTENT(IN)::I
    WRITE(IT,'(".",I4.4)') I
    END  FUNCTION IEXT
    !===========================================================================================
    CHARACTER(LEN=7)FUNCTION NUMCHAR(I)RESULT(IT)
    INTEGER,INTENT(IN)::I
    WRITE(IT,'(I7.7)') I
    END  FUNCTION NUMCHAR
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE SETEXT(I)
    INTEGER, INTENT(INOUT) :: I
    IN_EXT=IEXT(I)
    I=I+1
    OUTEXT=IEXT(I)
    END SUBROUTINE
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE DEL_FILE(FILE)
    CHARACTER(*), INTENT(IN) ::FILE
    LOGICAL :: EXIST
    INTEGER :: IVAR
    INQUIRE (FILE=FILE , ERR=1, IOSTAT=IVAR ,  EXIST=EXIST)
    IF(EXIST) THEN 
        OPEN(UNIT=1,FILE=FILE,ERR=1,IOSTAT=IVAR, ACTION='WRITE',STATUS='REPLACE')
        CLOSE(1,STATUS='DELETE')
    END IF
    RETURN
1   STOP 'error in del_file'
    END SUBROUTINE DEL_FILE
    !----------------------------------------------------------------------------------------------
    CHARACTER(LEN=128) FUNCTION PRINTFILE(N)RESULT(FL)
    INTEGER,INTENT(IN)::N
    INQUIRE(UNIT=N,NAME=FL)
    END FUNCTION PRINTFILE
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE DUMP_DENS(QZ,FN)
    REAL,INTENT(IN) :: QZ(NT)
    CHARACTER(*),INTENT(IN) ::FN 
    OPEN(UNIT=299,FILE=FN,FORM ='FORMATTED',ACTION='WRITE',STATUS='REPLACE')
    WRITE(299,'(" ", 2G15.8E3," ")')  QZ
    CLOSE(299)
    END SUBROUTINE DUMP_DENS
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE DUMP_IDENS(QZ,FN)
    REAL,INTENT(IN) :: QZ(0:1,NT)
    CHARACTER(*),INTENT(IN) ::FN 
    OPEN(UNIT=299,FILE=FN,FORM ='FORMATTED',ACTION='WRITE',STATUS='REPLACE')
    WRITE(299,'(" ", 2G15.8E3," ")')  QZ
    CLOSE(299)
    END SUBROUTINE DUMP_IDENS
    !-----------------------------------------------------------------------------------------------------
   SUBROUTINE DUMP_DENS8(QZ,FN)
    REAL(8),INTENT(IN) :: QZ(NT)
    CHARACTER(*),INTENT(IN) ::FN 
    OPEN(UNIT=299,FILE=FN,FORM ='FORMATTED',ACTION='WRITE',STATUS='REPLACE')
    WRITE(299,'(" ", 2G15.8E3," ")')  QZ
    CLOSE(299)
    END SUBROUTINE DUMP_DENS8
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE DUMP_IDENS8(QZ,FN)
    REAL(8),INTENT(IN) :: QZ(0:1,NT)
    CHARACTER(*),INTENT(IN) ::FN 
    OPEN(UNIT=299,FILE=FN,FORM ='FORMATTED',ACTION='WRITE',STATUS='REPLACE')
    WRITE(299,'(" ", 2G15.8E3," ")')  QZ
    CLOSE(299)
    END SUBROUTINE DUMP_IDENS8
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE INIT_dat(FN,HEAD,N)
    CHARACTER(*) ::FN
    CHARACTER(*):: HEAD(N)
    INTEGER::I,N
    INTENT(IN) :: FN,N,HEAD
    OPEN(UNIT=513,FILE=FN,FORM ='FORMATTED',ACTION='WRITE',STATUS='REPLACE')
    WRITE(513,'(A24,<n-1>(A2,A24))')  HEAD(1),(TAB,HEAD(I),I=2,n)
    CLOSE(513,STATUS='SAVE')
    END SUBROUTINE INIT_dat

    SUBROUTINE APPEND_DAT(FN,DAT,N)
    CHARACTER(*) ::FN
    REAL(8):: DAT(N)
    INTEGER::I,N
    INTENT(IN) :: FN,N ,DAT
    OPEN(UNIT=513,FILE=FN,FORM ='FORMATTED',ERR=1,ACTION='WRITE',STATUS='OLD',POSITION='APPEND')
    WRITE(513,'(G24.15E3,<N-1>(A2,G24.15E3))',ERR=1)  DAT(1),(TAB,DAT(I),I=2,N)
    CLOSE(513,STATUS='SAVE')
    ! PRINT*,FN
    RETURN
1   STOP 'ERR' 
    END SUBROUTINE APPEND_DAT




    end module inout2d