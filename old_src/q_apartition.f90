    MODULE PARTITION
    IMPLICIT NONE
    REAL(8), PUBLIC ::  M_POLY            ! NUMBER OF POLYMER CHAINS IN THE PERIODICAL CELL
    REAL(8),PUBLIC :: DELTA_W,LAMBDA_W           
    REAL(8),PUBLIC::    RG,B2 !,KH   !b2-colil-coil virial coefficient, m2 - number of charges per chain
    INTEGER,PUBLIC:: M2
    REAL(8), PUBLIC ::QZ !,QZ0I !SINGLE CHAIN PARTITION FUNCTION
    REAL(8),PUBLIC:: QPLUS,QMINUS,M_PLUS,M_MINUS, M_IONS ! SINGLE ION PARTITION FUNCTIONS AND REFERENCE CONCENTRATIONS: M_IONS=M_PLU+M_MINUS
    REAL(8),PUBLIC ::    MPOLY,MiPLUS, MiMINUS ! instant Number of positive and negative salt ions in the system 
    REAL(8),PUBLIC ::    ZPOLY,ZPLUS, ZMINUS ! activities of polymers, and ions 
    REAL(8):: PQN !PQN=MPARTICLES*Q_PARTICLE

    !CONTAINS


    END MODULE partition