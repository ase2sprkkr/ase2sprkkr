C====================================================================
C  DEBUG_DUMP module (fixed-form .f) - write each variable to
C  <name>.txt with deterministic, diff-friendly formatting.
C  Supports: INTEGER, INTEGER(KIND=8), REAL*8, LOGICAL, CHARACTER
C            Scalars and 1D/2D/3D arrays (where applicable)
C
C  Usage (fixed-form):
C    USE DEBUG_DUMP
C    CALL DUMP_VAR('energy', energy)
C    CALL DUMP_VAR('forces', forces)
C
C  Notes:
C    - Files are overwritten (STATUS='REPLACE').
C    - REAL*8 uses ES25.16E3 format.
C    - Arrays dumped in column-major order with indices increasing
C      fastest in the first dimension.
C====================================================================
      MODULE DEBUG_DUMP
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT64, REAL64
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: DUMP_VAR
C
      INTERFACE DUMP_VAR
        MODULE PROCEDURE DD_INT_SCAL
        MODULE PROCEDURE DD_INT64_SCAL
        MODULE PROCEDURE DD_REAL64_SCAL
        MODULE PROCEDURE DD_LOG_SCAL
        MODULE PROCEDURE DD_CHAR_SCAL
C       1D
        MODULE PROCEDURE DD_INT_1D
        MODULE PROCEDURE DD_INT64_1D
        MODULE PROCEDURE DD_REAL64_1D
        MODULE PROCEDURE DD_LOG_1D
        MODULE PROCEDURE DD_CHAR_1D
C       2D
        MODULE PROCEDURE DD_INT_2D
        MODULE PROCEDURE DD_INT64_2D
        MODULE PROCEDURE DD_REAL64_2D
        MODULE PROCEDURE DD_LOG_2D
C       3D
        MODULE PROCEDURE DD_INT_3D
        MODULE PROCEDURE DD_INT64_3D
        MODULE PROCEDURE DD_REAL64_3D
        MODULE PROCEDURE DD_LOG_3D
      END INTERFACE
C
      CONTAINS
C--------------------------------------------------------------------
C  Helpers
C--------------------------------------------------------------------
      SUBROUTINE DD_OPEN(U, BASENAME)
        INTEGER            U
        CHARACTER*(*)      BASENAME
        OPEN(NEWUNIT=U, FILE=TRIM(BASENAME)//'.txt',
     &       ACTION='WRITE', STATUS='REPLACE', FORM='FORMATTED')
      END SUBROUTINE DD_OPEN
C
      SUBROUTINE DD_WRITE_SHAPE(U, RANK1, D1, D2, D3)
        INTEGER U, RANK1, D1, D2, D3
        CHARACTER*64 SH
        IF (RANK1 .EQ. 0) THEN
          WRITE(U,'(A)') '# shape: ()'
        ELSE IF (RANK1 .EQ. 1) THEN
          WRITE(SH,'(I0)') D1
          WRITE(U,'(A,A,A)') '# shape: (', TRIM(SH), ')'
        ELSE IF (RANK1 .EQ. 2) THEN
          WRITE(U,'(A,I0,A,I0,A)') '# shape: (', D1, ',', D2, ')'
        ELSE
          WRITE(U,'(A,I0,A,I0,A,I0,A)') '# shape: (', D1, ',', D2,
     &                                   ',', D3, ')'
        END IF
      END SUBROUTINE DD_WRITE_SHAPE
C--------------------------------------------------------------------
C  Scalars
C--------------------------------------------------------------------
      SUBROUTINE DD_INT_SCAL(NAME, X)
        CHARACTER*(*) NAME
        INTEGER       X
        INTEGER U
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer'
        CALL DD_WRITE_SHAPE(U,0,0,0,0)
        WRITE(U,'(I0)') X
        CLOSE(U)
      END SUBROUTINE DD_INT_SCAL
C
      SUBROUTINE DD_INT64_SCAL(NAME, X)
        CHARACTER*(*)   NAME
        INTEGER(INT64)  X
        INTEGER U
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer(kind=8)'
        CALL DD_WRITE_SHAPE(U,0,0,0,0)
        WRITE(U,'(I0)') X
        CLOSE(U)
      END SUBROUTINE DD_INT64_SCAL
C
      SUBROUTINE DD_REAL64_SCAL(NAME, X)
        CHARACTER*(*) NAME
        REAL(REAL64)  X
        INTEGER U
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: real(kind=8)'
        CALL DD_WRITE_SHAPE(U,0,0,0,0)
        WRITE(U,'(ES25.16E3)') X
        CLOSE(U)
      END SUBROUTINE DD_REAL64_SCAL
C
      SUBROUTINE DD_LOG_SCAL(NAME, X)
        CHARACTER*(*) NAME
        LOGICAL       X
        INTEGER U
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: logical'
        CALL DD_WRITE_SHAPE(U,0,0,0,0)
        WRITE(U,'(L1)') X
        CLOSE(U)
      END SUBROUTINE DD_LOG_SCAL
C
      SUBROUTINE DD_CHAR_SCAL(NAME, X)
        CHARACTER*(*) NAME
        CHARACTER*(*) X
        INTEGER U
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: character'
        CALL DD_WRITE_SHAPE(U,0,0,0,0)
        WRITE(U,'(A)') TRIM(X)
        CLOSE(U)
      END SUBROUTINE DD_CHAR_SCAL
C--------------------------------------------------------------------
C  Integer arrays
C--------------------------------------------------------------------
      SUBROUTINE DD_INT_1D(NAME, X)
        CHARACTER*(*) NAME
        INTEGER       X(:)
        INTEGER U, I
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer'
        CALL DD_WRITE_SHAPE(U,1,SIZE(X),0,0)
        DO I = 1, SIZE(X)
          WRITE(U,'(I0)') X(I)
        END DO
        CLOSE(U)
      END SUBROUTINE DD_INT_1D
C
      SUBROUTINE DD_INT_2D(NAME, X)
        CHARACTER*(*) NAME
        INTEGER       X(:,:)
        INTEGER U, I, J
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer'
        CALL DD_WRITE_SHAPE(U,2,SIZE(X,1),SIZE(X,2),0)
        DO J = 1, SIZE(X,2)
          DO I = 1, SIZE(X,1)
            WRITE(U,'(I0)') X(I,J)
          END DO
        END DO
        CLOSE(U)
      END SUBROUTINE DD_INT_2D
C
      SUBROUTINE DD_INT_3D(NAME, X)
        CHARACTER*(*) NAME
        INTEGER       X(:,:,:)
        INTEGER U, I, J, K
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer'
        CALL DD_WRITE_SHAPE(U,3,SIZE(X,1),SIZE(X,2),SIZE(X,3))
        DO K = 1, SIZE(X,3)
          DO J = 1, SIZE(X,2)
            DO I = 1, SIZE(X,1)
              WRITE(U,'(I0)') X(I,J,K)
            END DO
          END DO
        END DO
        CLOSE(U)
      END SUBROUTINE DD_INT_3D
C--------------------------------------------------------------------
C  INT64 arrays
C--------------------------------------------------------------------
      SUBROUTINE DD_INT64_1D(NAME, X)
        CHARACTER*(*)   NAME
        INTEGER(INT64)  X(:)
        INTEGER U, I
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer(kind=8)'
        CALL DD_WRITE_SHAPE(U,1,SIZE(X),0,0)
        DO I = 1, SIZE(X)
          WRITE(U,'(I0)') X(I)
        END DO
        CLOSE(U)
      END SUBROUTINE DD_INT64_1D
C
      SUBROUTINE DD_INT64_2D(NAME, X)
        CHARACTER*(*)   NAME
        INTEGER(INT64)  X(:,:)
        INTEGER U, I, J
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer(kind=8)'
        CALL DD_WRITE_SHAPE(U,2,SIZE(X,1),SIZE(X,2),0)
        DO J = 1, SIZE(X,2)
          DO I = 1, SIZE(X,1)
            WRITE(U,'(I0)') X(I,J)
          END DO
        END DO
        CLOSE(U)
      END SUBROUTINE DD_INT64_2D
C
      SUBROUTINE DD_INT64_3D(NAME, X)
        CHARACTER*(*)   NAME
        INTEGER(INT64)  X(:,:,:)
        INTEGER U, I, J, K
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer(kind=8)'
        CALL DD_WRITE_SHAPE(U,3,SIZE(X,1),SIZE(X,2),SIZE(X,3))
        DO K = 1, SIZE(X,3)
          DO J = 1, SIZE(X,2)
            DO I = 1, SIZE(X,1)
              WRITE(U,'(I0)') X(I,J,K)
            END DO
          END DO
        END DO
        CLOSE(U)
      END SUBROUTINE DD_INT64_3D
C--------------------------------------------------------------------
C  REAL*8 arrays
C--------------------------------------------------------------------
      SUBROUTINE DD_REAL64_1D(NAME, X)
        CHARACTER*(*) NAME
        REAL(REAL64)  X(:)
        INTEGER U, I
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: real(kind=8)'
        CALL DD_WRITE_SHAPE(U,1,SIZE(X),0,0)
        DO I = 1, SIZE(X)
          WRITE(U,'(ES25.16E3)') X(I)
        END DO
        CLOSE(U)
      END SUBROUTINE DD_REAL64_1D
C
      SUBROUTINE DD_REAL64_2D(NAME, X)
        CHARACTER*(*) NAME
        REAL(REAL64)  X(:,:)
        INTEGER U, I, J
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: real(kind=8)'
        CALL DD_WRITE_SHAPE(U,2,SIZE(X,1),SIZE(X,2),0)
        DO J = 1, SIZE(X,2)
          DO I = 1, SIZE(X,1)
            WRITE(U,'(ES25.16E3)') X(I,J)
          END DO
        END DO
        CLOSE(U)
      END SUBROUTINE DD_REAL64_2D
C
      SUBROUTINE DD_REAL64_3D(NAME, X)
        CHARACTER*(*) NAME
        REAL(REAL64)  X(:,:,:)
        INTEGER U, I, J, K
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: real(kind=8)'
        CALL DD_WRITE_SHAPE(U,3,SIZE(X,1),SIZE(X,2),SIZE(X,3))
        DO K = 1, SIZE(X,3)
          DO J = 1, SIZE(X,2)
            DO I = 1, SIZE(X,1)
              WRITE(U,'(ES25.16E3)') X(I,J,K)
            END DO
          END DO
        END DO
        CLOSE(U)
      END SUBROUTINE DD_REAL64_3D
C--------------------------------------------------------------------
C  LOGICAL arrays
C--------------------------------------------------------------------
      SUBROUTINE DD_LOG_1D(NAME, X)
        CHARACTER*(*) NAME
        LOGICAL       X(:)
        INTEGER U, I
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: logical'
        CALL DD_WRITE_SHAPE(U,1,SIZE(X),0,0)
        DO I = 1, SIZE(X)
          WRITE(U,'(L1)') X(I)
        END DO
        CLOSE(U)
      END SUBROUTINE DD_LOG_1D
C
      SUBROUTINE DD_LOG_2D(NAME, X)
        CHARACTER*(*) NAME
        LOGICAL       X(:,:)
        INTEGER U, I, J
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: logical'
        CALL DD_WRITE_SHAPE(U,2,SIZE(X,1),SIZE(X,2),0)
        DO J = 1, SIZE(X,2)
          DO I = 1, SIZE(X,1)
            WRITE(U,'(L1)') X(I,J)
          END DO
        END DO
        CLOSE(U)
      END SUBROUTINE DD_LOG_2D
C
      SUBROUTINE DD_LOG_3D(NAME, X)
        CHARACTER*(*) NAME
        LOGICAL       X(:,:,:)
        INTEGER U, I, J, K
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: logical'
        CALL DD_WRITE_SHAPE(U,3,SIZE(X,1),SIZE(X,2),SIZE(X,3))
        DO K = 1, SIZE(X,3)
          DO J = 1, SIZE(X,2)
            DO I = 1, SIZE(X,1)
              WRITE(U,'(L1)') X(I,J,K)
            END DO
          END DO
        END DO
        CLOSE(U)
      END SUBROUTINE DD_LOG_3D
C--------------------------------------------------------------------
C  CHARACTER arrays (1D)
C--------------------------------------------------------------------
      SUBROUTINE DD_CHAR_1D(NAME, X)
        CHARACTER*(*) NAME
        CHARACTER*(*) X(:)
        INTEGER U, I
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: character'
        CALL DD_WRITE_SHAPE(U,1,SIZE(X),0,0)
        DO I = 1, SIZE(X)
          WRITE(U,'(A)') TRIM(X(I))
        END DO
        CLOSE(U)
      END SUBROUTINE DD_CHAR_1D
C--------------------------------------------------------------------
      END MODULE DEBUG_DUMP
C*==emptyspheres.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE EMPTYSPHERES(AVC,BVC,CVC,BAS,IS,ALAT,NATOM,NSORT,S,
     &                        NSORTES,TAUES,SES,MAX2SORT,ITOP,IQA,NG,
     &                        ISNEW,BASNEW,NATOMNEW,ROTATIONS,
     &                        TRANSLATIONS,N_SYMMETRY_OPS)
      USE DEBUG_DUMP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NPATB1
      PARAMETER (NPATB1=6000)
C
C Dummy arguments
C
      REAL*8 ALAT
      INTEGER MAX2SORT,NATOM,NATOMNEW,NG,NSORT,NSORTES
      INTEGER N_SYMMETRY_OPS
      REAL*8 AVC(3),BAS(3,MAX2SORT),BASNEW(3,MAX2SORT),BVC(3),
     &       CVC(3),S(MAX2SORT),
     &       SES(MAX2SORT),
     &       TAUES(3,MAX2SORT)
      INTEGER IQA(48),IS(MAX2SORT),ISNEW(MAX2SORT),ITOP(48),N_SYMMETRY_OPS
      DOUBLE PRECISION,  DIMENSION(3,3,N_SYMMETRY_OPS)
     &       :: ROTATIONS
      DOUBLE PRECISION,  DIMENSION(3,N_SYMMETRY_OPS)
     &       :: TRANSLATIONS
C
C Local variables
C
      REAL*8 ANG(3),GEN(3,3,48),S2NEW(:),SNEW(:),TAU(:,:),VC(3,48),VV(3)
      INTEGER I,IER,IG,IGA,IGC,ISB(:),ISKIP(:),I_ISB,I_ISKIP,I_NHSORTES,
     &        I_S2NEW,I_SNEW,I_TAU,J,MAXNPAT,MAXNPATG,NHSORTES(:)

      COMMON /IPRINT/ IPRINT
      INTEGER IPRINT
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE SNEW,S2NEW,TAU,NHSORTES,ISKIP,ISB
C
      MAXNPAT = MAX2SORT
      MAXNPATG = NPATB1
C
      DO I = 1,3
         VV(I) = -BAS(I,1)
      END DO
C
CCC      DO IGC = 1,NG
CCC         IG = ITOP(IGC)
CCC         IGA = IG
CCC         IF ( IGA.GT.32 ) IGA = IGA - 32
CCC         CALL EILANG(ANG,IGA)
CCC         CALL TURNM(ANG,GEN(1,1,IGC))
CCCC                              !DEF TURN MATRIX
CCC         IF ( IG.NE.IGA ) THEN
CCCC                              !ADD INVERSION IF NEED
CCC            DO I = 1,3
CCC               DO J = 1,3
CCC                  GEN(I,J,IGC) = -GEN(I,J,IGC)
CCC               END DO
CCC            END DO
CCC         END IF
CCC         write(*,*) "IGC",IGC,ITOP(IGC),IQA(IGC)
CCC         write(*,*) "VEC",VV,BAS(1:3,IQA(IGC))
CCC         write(*,*) "ROTMAT1",GEN(1:3,1,IGC)
CCC         write(*,*) "ROTMAT2",GEN(1:3,2,IGC)
CCC         write(*,*) "ROTMAT3",GEN(1:3,3,IGC)
CCC         CALL VEC(VC(1,IGC),BAS(1,IQA(IGC)),VV,GEN(1,1,IGC))
CCC         write(*,*) "RES",VC(1:3,IGC)
CCC         END DO
CCC
                 
         NG=N_SYMMETRY_OPS 

         CALL READ_IQA1("IQA1.txt", IQA, NG)
         GEN=ROTATIONS

         DO IGC=1,N_SYMMETRY_OPS
           CALL VEC(VC(1,IGC),BAS(1,IQA(IGC)),VV,GEN(1,1,IGC))
         END DO
C
C
C
      ALLOCATE (SNEW(MAXNPAT),S2NEW(MAXNPAT),TAU(3,MAXNPATG))
      ALLOCATE (NHSORTES(MAXNPAT),ISKIP(MAXNPAT),ISB(MAXNPATG))
C
      CALL EMPTY(IER,AVC,BVC,CVC,BAS,IS,ALAT,NATOM,NSORT,S,BASNEW,ISNEW,
     &           SNEW,S2NEW,ISB,TAU,MAXNPAT,MAXNPATG,NSORTES,TAUES,SES,
     &           NHSORTES,ISKIP,GEN,VC,NG,NATOMNEW)
C
C     nsortes=iter-1
      IF( IPRINT > 0 ) THEN
        WRITE (6,*) ' IRR. POSITIONS:'
        DO I = 1,NSORTES
           WRITE (6,'(3f21.15)') (TAUES(J,I),J=1,3)
        END DO
      END IF
C      stop
      DEALLOCATE (SNEW,S2NEW,TAU)
      DEALLOCATE (NHSORTES,ISKIP,ISB)

      END
C*==empty.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE EMPTY(IER,AVC,BVC,CVC,BAS,IS,ALAT,NATOM,NSORT,S,BASNEW,
     &                 ISNEW,SNEW,S2NEW,ISB,TAU,MAXNPAT,MAXNPATG,
     &                 NSORTES,TAUES,SES,NHSORTES,ISKIP,GM,VC,NG,
     &                 NATOMNEW)

      USE DEBUG_DUMP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C COMMON variables
C
      REAL*8 BOA,COA,RMAXES,RMINES,U(3,3),UM1(3,3)
      INTEGER IDATASTR,N1ES,N2ES,N3ES,NATES
      COMMON /EMPTYS/ NATES,N1ES,N2ES,N3ES,RMINES,RMAXES
      COMMON /TAU2  / BOA,COA
      COMMON /U     / U,UM1,IDATASTR
C
C Dummy arguments
C
      REAL*8 ALAT
      INTEGER IER,MAXNPAT,MAXNPATG,NATOM,NATOMNEW,NG,NSORT,NSORTES
      REAL*8 AVC(3),BAS(3,*),BASNEW(3,*),BVC(3),CVC(3),GM(3,3,48),S(*),
     &       S2NEW(*),SES(*),SNEW(*),TAU(3,*),TAUES(3,*),VC(3,48)
      INTEGER IS(*),ISB(*),ISKIP(*),ISNEW(*),NHSORTES(*)
C
C Local variables
C
      REAL*8 ANG(3),CS,CVOL,DEMIN,EPS_ANG,EPS_ANGMAX,PA,PB,PC,PI43,RAD,
     &       RADNEW,SM,SMN,SSM,SVOLA,SVOLE,VE(3,48),VSHIFT(3),VV(3)
      REAL*8 DPI,TRNT
      INTEGER I,IATOM,IG,IIN,IN,INVMAX,IODN,ISORT,ITAV,ITER,IUN,J,K,
     &        NATB,NATS,NSORTNEW

      INTEGER IPRINT
      COMMON /IPRINT/ IPRINT
C
C*** End of declarations rewritten by SPAG
C
C
C COMMON variables
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
C
C COMMON variables
C
C
C Dummy arguments
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
C ---- Dump each variable
      DO I = 1,3
         DO J = 1,3
            U(I,J) = 0
            UM1(I,J) = 0
            GM(I,J,1) = 0
         END DO
         U(I,I) = 1
         UM1(I,I) = 1
         GM(I,I,1) = 1
         VSHIFT(I) = 0
         VC(I,1) = 0
      END DO
C      it(1)=1
C     ng=1
      BOA = 1
      COA = 1
C---------------------------------------
      IUN = 6
      NATES = 0
      SSM = 0
      DO I = 1,NSORT
         SNEW(I) = S(I)/ALAT
         SSM = MAX(SSM,SNEW(I))
         S2NEW(I) = SNEW(I)**2
      END DO
      INVMAX = 0
      IF ( RMAXES.LT.0 ) INVMAX = 1
      IF ( RMINES.LT.0 ) INVMAX = -INVMAX
      RMINES = ABS(RMINES)
      RMAXES = ABS(RMAXES)
      SM = RMAXES/ALAT + SSM
      SMN = RMINES/ALAT
      DO IATOM = 1,NATOM
         ISNEW(IATOM) = IS(IATOM)
         DO I = 1,3
            BASNEW(I,IATOM) = BAS(I,IATOM)
         END DO
      END DO
      NATOMNEW = NATOM
      NSORTNEW = NSORT
      ITER = 0
      IF(IPRINT >0) THEN
        PRINT *,'BOA,COA:',BOA,COA,ALAT
        WRITE (IUN,*) '  Searching of Empty spheres  '
        WRITE (*,*) '  Searching of Empty spheres  '
        WRITE (IUN,'(a,3i3)') ' Mesh: ',N1ES,N2ES,N3ES
        WRITE (IUN,99001) 'Smin, Smax (a.u.): ',RMINES,RMAXES
      END IF
 100  CONTINUE
      ITER = ITER + 1
      IF(IPRINT >0) THEN
        WRITE (IUN,*) '  Iteration ',ITER
        WRITE (*,*) '  Iteration ',ITER
      END IF
C
      CALL GETEPOS(MAXNPATG,AVC,BVC,CVC,BASNEW,NATOMNEW,ISNEW,TAU,NATB,
     &             ISB,SNEW,S2NEW,N1ES,N2ES,N3ES,PA,PB,PC,RAD,IER,SM)
C
      IF(IPRINT >0) THEN
        WRITE (IUN,99001) ' Maximum possible sphere is ',RAD*ALAT
        WRITE (IUN,99001) ' Position in a,b,c:',PA,PB,PC
      END IF
      IF ( RAD.LT.SMN ) THEN

         IF(IPRINT >0) THEN
           WRITE (IUN,*) 'skipped'
           PRINT *,'RAD=',RAD*ALAT,' RMINES=',RMINES
         END IF
         GOTO 300
      END IF
C
      IF ( IER.EQ.0 ) THEN
         ISKIP(ITER) = 0
         IF ( RAD.GT.RMAXES/ALAT ) THEN
            IF ( INVMAX.EQ.0 ) THEN
              IF(IPRINT >0) WRITE (IUN,99001) 'Contracted to be ',RMAXES
               RAD = RMAXES/ALAT        ! rad is used later
            ELSE
               ISKIP(ITER) = INVMAX
               IF (IPRINT > 0 ) THEN
                 WRITE (IUN,99001) 'Sphere is too large',RAD*ALAT
                 WRITE (IUN,99001)
     &                      'Small auxiliary sphere is inserted instead'
               END IF
            END IF
         END IF
C
C multiplication :
         VE(1,1) = PA*AVC(1) + PB*BVC(1) + PC*CVC(1) - VSHIFT(1)
         VE(2,1) = PA*AVC(2) + PB*BVC(2) + PC*CVC(2) - VSHIFT(2)
         VE(3,1) = PA*AVC(3) + PB*BVC(3) + PC*CVC(3) - VSHIFT(3)
C
         ITAV = 0
 150     CONTINUE
         ITAV = ITAV + 1
         IN = 1
         EPS_ANGMAX = 0.D0
         CALL SHORTN(VE,VE)
C       print*,(ve(i,1),i=1,3)
         DO IG = 1,NG
            CALL VEC(VV,VC(1,IG),VE,GM(1,1,IG))
C         print*,vv
            DO IIN = 1,IN
               ANG(1) = VE(1,IIN) - VV(1)
               ANG(2) = VE(2,IIN) - VV(2)
               ANG(3) = VE(3,IIN) - VV(3)
               CALL SHORTN(ANG,ANG)
               EPS_ANG = ANG(1)**2 + ANG(2)**2 + ANG(3)**2
               IF ( EPS_ANG.LT.1.D-5 ) THEN
                  EPS_ANGMAX = MAX(EPS_ANG,EPS_ANGMAX)
C              print*,iin,eps_ang
                  GOTO 200
               END IF
            END DO
            IN = IN + 1
            VE(1,IN) = VV(1)
            VE(2,IN) = VV(2)
            VE(3,IN) = VV(3)
C         print*,'in=',in,vv
 200     END DO
         VV(1) = 0
         VV(2) = 0
         VV(3) = 0
         IODN = 0
         DEMIN = 1.D10
         DO IIN = 2,IN
            ANG(1) = VE(1,IIN) - VE(1,1)
            ANG(2) = VE(2,IIN) - VE(2,1)
            ANG(3) = VE(3,IIN) - VE(3,1)
            CALL SHORTN(ANG,ANG)
            EPS_ANG = ANG(1)**2 + ANG(2)**2 + ANG(3)**2
            EPS_ANG = SQRT(EPS_ANG)
            IF ( DEMIN.GT.EPS_ANG ) DEMIN = EPS_ANG
            IF ( EPS_ANG.LT.RAD*1.5D0 ) THEN
               VV(1) = VV(1) + ANG(1)
               VV(2) = VV(2) + ANG(2)
               VV(3) = VV(3) + ANG(3)
               IODN = IODN + 1
            END IF
         END DO                         ! iin
         IF ( IODN.NE.0 ) THEN
            IF (IPRINT > 0 ) THEN
              WRITE (IUN,'(a,i3,2f21.15)') 'a bit shifted position'
     &             ,IODN, REAL(DEMIN),REAL(RAD)
              WRITE (IUN,99001) 'OLD:',VE(1,1),VE(2,1)/BOA,VE(3,1)/COA
            END IF
            VE(1,1) = VE(1,1) + VV(1)/(IODN+1)
            VE(2,1) = VE(2,1) + VV(2)/(IODN+1)
            VE(3,1) = VE(3,1) + VV(3)/(IODN+1)
            IF (IPRINT > 0 ) WRITE (IUN,99001)
     &              'new:',VE(1,1),VE(2,1)/BOA,VE(3,1)/COA
            CALL REGET_RE(TAU,ISB,SNEW,NATB,VSHIFT,VE(1,1),RADNEW)
            IF (IPRINT > 0 ) WRITE (IUN,99001)
     &              'rad,radtst',RAD*ALAT,RADNEW*ALAT
            RAD = RADNEW
            IF ( ITAV.LE.5 ) GOTO 150   ! to avoide infinite cycle
            IF (IPRINT > 0 )
     >          WRITE (IUN,*) 'can not find averaged position'
            STOP
         ELSE IF ( DEMIN.LT.2.D0*RAD ) THEN
                                        ! overlapping spheres
            IF (IPRINT > 0 )
     &          WRITE (IUN,*) 'radius is decreased from',RAD,' to',
     &                    DEMIN/2.D0
            RAD = DEMIN/2.D0            ! decrease radius to touching
            IF ( RAD.LT.SMN ) THEN
               IF (IPRINT > 0 ) THEN
                 WRITE (IUN,*) 'sphere is less then rmines'
                 WRITE (IUN,*) 'rad=',RAD*ALAT,' rmines=',RMINES
               END IF
               ISKIP(ITER) = -1         ! will be skipped later
            END IF
         END IF
         IF (IPRINT > 0 ) THEN
           WRITE (IUN,'(a,i4,a)') 'it will be ',IN,' ES per unit cell:'
         END IF
         IF ( NATOMNEW+IN.GT.MAXNPAT ) THEN
            PRINT *,'EMPTY: too many atoms',NATOMNEW + IN,' maxnpat=',
     &            MAXNPAT
            PRINT *,'try to increase array sizes in empty.f'
            STOP
         END IF
         DO IATOM = NATOMNEW + 1,NATOMNEW + IN
            ISNEW(IATOM) = NSORTNEW + 1
            DO I = 1,3
               BASNEW(I,IATOM) = VE(I,IATOM-NATOMNEW) + VSHIFT(I)
            END DO
            IF (IPRINT > 0 )
     &          WRITE (IUN,'(2x,3f15.8)') (BASNEW(I,IATOM),I=1,3)
         END DO
         write(*,*) "MUMU",VE,UM1
C
         DO K = 1,3
            VV(K) = 0
            DO J = 1,3
               VV(K) = VV(K) + VE(J,1)*UM1(J,K)
            END DO
            TAUES(K,ITER) = VV(K)
         END DO
C      write(iun,'(i3,3f15.8)')iter,vv
         NSORTNEW = NSORTNEW + 1
         NATOMNEW = NATOMNEW + IN
         NHSORTES(ITER) = IN
         IF ( ISKIP(ITER).EQ.0 ) THEN
            SNEW(NSORTNEW) = MIN(RAD,RMAXES/ALAT)
         ELSE
            SNEW(NSORTNEW) = RMINES/ALAT
         END IF
         S2NEW(NSORTNEW) = SNEW(NSORTNEW)**2
         SES(ITER) = SNEW(NSORTNEW)*ALAT
         IF (IPRINT > 0 )
     &      PRINT *,IN,' Empty spheres added, S=',SNEW(NSORTNEW)*ALAT
         NATES = NATES + IN
         GOTO 100
      END IF
C
 300  CONTINUE
      NSORTES = ITER - 1
      IF ( INVMAX.GE.0 ) THEN
C check if some spheres should be skipped due to small radius
         DO ISORT = 1,NSORTES
            IF ( ISKIP(ISORT).LT.0 ) INVMAX = ISKIP(ISORT)
         END DO
      END IF
      IF ( INVMAX.LT.0 ) THEN
C skip auxiliary spheres
C$$$        write(iun,*)'nsortes,natomnew',nsortes,natomnew
C$$$        do isort=1,nsortes
C$$$          write(iun,'(i3,3f21.15,2i4,f21.15)')isort,(taues(k,isort),k=1,3)
C$$$     $         ,nhsortes(isort),iskip(isort),ses(isort)
C$$$        enddo
C$$$        do iatom=1,natomnew
C$$$          write(iun,'(i3,3f21.15)')iatom,(basnew(k,iatom),k=1,3)
C$$$        enddo
C$$$        write(iun,*)'rmines',rmines
         IATOM = 1
         ISORT = 1
         DO WHILE ( ISORT.LE.NSORTES )
            IF (IPRINT > 0 )
     &         WRITE (IUN,*) 'isort,nsortes,s',ISORT,NSORTES,SES(ISORT)
            NATS = NHSORTES(ISORT)
            IF ( ISKIP(ISORT).LT.0 ) THEN
               IF (IPRINT > 0 )
     &            WRITE (IUN,99001) 'sphere is skipped',
     &                           (TAUES(K,ISORT),K=1,3)
               DO I = ISORT + 1,NSORTES
                  NHSORTES(I-1) = NHSORTES(I)
                  SES(I-1) = SES(I)
                  ISKIP(I-1) = ISKIP(I)
                  DO K = 1,3
                     TAUES(K,I-1) = TAUES(K,I)
                  END DO
               END DO
               DO I = IATOM + NATS,NATOMNEW
                  DO K = 1,3
                     BASNEW(K,I-NATS) = BASNEW(K,I)
                  END DO
               END DO
               NSORTES = NSORTES - 1
               NATOMNEW = NATOMNEW - NATS
            ELSE
               IATOM = IATOM + NATS
               ISORT = ISORT + 1
            END IF
         END DO
C$$$        write(iun,*)'nsortes,natomnew',nsortes,natomnew
C$$$        do isort=1,nsortes
C$$$          write(iun,'(i3,3f21.15,i4,f21.15)')isort,(taues(k,isort),k=1,3)
C$$$     $         ,nhsortes(isort),ses(isort)
C$$$        enddo
C$$$        do iatom=1,natomnew
C$$$          write(iun,'(i3,3f21.15)')iatom,(basnew(k,iatom),k=1,3)
C$$$        enddo
      END IF
C$$$      if(nsortes.gt.0)then
C$$$        open(24,status='scratch',form='unformatted')
C$$$        do iatom=natom+1,natomnew
C$$$          write(24)(basnew(i,iatom),i=1,3)
C$$$c$$$        write(iun,'(2x,3f15.8,i4)')(basnew(i,iatom),i=1,3),isnew(iatom)
C$$$        enddo
C$$$      endif
C MT -> ASA spheres
      PI43 = 4.D0*DPI()/3.D0
      SVOLA = 0.D0
      DO IATOM = 1,NATOM
         SVOLA = SVOLA + PI43*S(IS(IATOM))**3
      END DO
      SVOLE = 0.D0
      DO ISORT = 1,NSORTES
         SVOLE = SVOLE + NHSORTES(ISORT)*PI43*SES(ISORT)**3
C        print*,'SES::',ses(isort),isort
      END DO
      CVOL = TRNT(AVC,BVC,CVC)*ALAT**3
      CS = CVOL/(SVOLA+SVOLE)
      SVOLA = SVOLA*CS
      SVOLE = SVOLE*CS
      CS = CS**(1.D0/3.D0)
      DO ISORT = 1,NSORT
         S(ISORT) = CS*S(ISORT)
      END DO
      DO ISORT = 1,NSORTES
         SES(ISORT) = CS*SES(ISORT)
      END DO
      IF (IPRINT > 0 ) THEN
        WRITE (IUN,99001) 'Vat, Ves, V',SVOLA,SVOLE,SVOLA + SVOLE
        WRITE (IUN,99001) 'Vat/V, Ves/V',SVOLA/CVOL,SVOLE/CVOL
      END IF
99001 FORMAT (1x,a,3F21.15)
      END
C*==getepos.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      SUBROUTINE GETEPOS(NPAT,A,B,C,BAS,NATOM,IS,TAU,NATB,ISB,S,S2,N1,
     &                   N2,N3,PAE,PBE,PCE,RAD,IER,SM)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IER,N1,N2,N3,NATB,NATOM,NPAT
      REAL*8 PAE,PBE,PCE,RAD,SM
      REAL*8 A(3),B(3),BAS(3,*),C(3),S(*),S2(*),TAU(3,NPAT)
      INTEGER IS(NPAT),ISB(NPAT)
C
C Local variables
C
      REAL*8 ANA,ANB,ANC,DD,DDMAX,DIAG,DIST2,DMAX,DMIN,DX,DY,DZ,PA,PB,
     &       PC,X,Y,Z
      INTEGER I,IA,IATOM,IB,IC,J,K
C
C*** End of declarations rewritten by SPAG
C
C
      DIAG = MAX((A(1)+B(1)+C(1))**2+(A(2)+B(2)+C(2))
     &       **2+(A(3)+B(3)+C(3))**2,(A(1)-B(1)+C(1))
     &       **2+(A(2)-B(2)+C(2))**2+(A(3)-B(3)+C(3))**2,
     &       (A(1)+B(1)-C(1))**2+(A(2)+B(2)-C(2))**2+(A(3)+B(3)-C(3))
     &       **2,(A(1)-B(1)-C(1))**2+(A(2)-B(2)-C(2))
     &       **2+(A(3)-B(3)-C(3))**2)
C
      IER = 0
      DIAG = SQRT(DIAG)/2
      DDMAX = (SM+DIAG)**2
      TAU=0.0D0
C
      NATB = 0
      DO I = -3,3
         DO J = -3,3
            DO K = -3,3
               DX = I*A(1) + J*B(1) + K*C(1)
               DY = I*A(2) + J*B(2) + K*C(2)
               DZ = I*A(3) + J*B(3) + K*C(3)
               DO IATOM = 1,NATOM
                  X = DX + BAS(1,IATOM)
                  Y = DY + BAS(2,IATOM)
                  Z = DZ + BAS(3,IATOM)
                  DD = X**2 + Y**2 + Z**2
                  IF ( DD.LT.DDMAX ) THEN
                     NATB = NATB + 1
                     TAU(1,NATB) = X
                     TAU(2,NATB) = Y
                     TAU(3,NATB) = Z
                     IF ( NATB.GT.NPAT ) THEN
                        IER = 1
                        PRINT *,'IER=',IER
                        PRINT *,
     &                        'GETEPOS: too many atoms are generated ',
     &                        NATB,' npat=',NPAT
                        PRINT *,'try to increase npatb1 in empty.f'
                        RETURN
                     END IF
                     ISB(NATB) = IS(IATOM)
                  END IF
               END DO                   ! iatom
            END DO                      ! k
         END DO                         ! j
      END DO                            ! i
      do i=1,NATB
      write(999,*) TAU(1,i),TAU(2,i),TAU(3,i)
      end do
C
      PAE = 0
      PBE = 0
      PCE = 0
C
      ANA = 0.5D0/N1
      ANB = 0.5D0/N2
      ANC = 0.5D0/N3
      DMAX = 0
C      npnt=0
      DO IA = -N1,N1 - 1,2
         PA = ANA*IA
         DO IB = -N2,N2 - 1,2
            PB = ANB*IB
            DO IC = -N3,N3 - 1,2
C            npnt=npnt+1
               PC = ANC*IC
C$$$            lg=ia.eq.-24.and.ib.eq.0.and.ic.eq.-72
               X = PA*A(1) + PB*B(1) + PC*C(1)
               Y = PA*A(2) + PB*B(2) + PC*C(2)
               Z = PA*A(3) + PB*B(3) + PC*C(3)
C           if(lg)print*,'XYZ:',x,y,z
               DMIN = 1.D10
               DO I = 1,NATB
                  DIST2 = (X-TAU(1,I))**2 + (Y-TAU(2,I))
     &                    **2 + (Z-TAU(3,I))**2
C             if(lg.and.isb(i).eq.8)print*,'TTT:',tau(1,i),tau(2,i)
C    $             ,tau(3,i)
C$$$              print *,'ia,ib,ic,i',ia,ib,ic,i,natb
C$$$              print *,'dist2,isb(i),sq',dist2,isb(i),s2(isb(i))
C$$$                PRINT*,I,DIST2,ISB(I),S2(ISB(I))
C$$$                print*, x,y,z,ia,ib,ic
C$$$                print*,tau(1,i),tau(2,i),tau(3,i)
                  IF ( DIST2.LT.S2(ISB(I)) ) GOTO 20
                  DD = SQRT(DIST2) - S(ISB(I))
                  DMIN = MIN(DD,DMIN)
               END DO                   ! i
               IF ( DMIN.GT.DMAX ) THEN
C            print*,"HUHU",PA,PB,PC,DMIN
                  DMAX = DMIN
                  PAE = PA
                  PBE = PB
                  PCE = PC
C$$$              iae=ia
C$$$              ibe=ib
C$$$              ice=ic
               END IF
 20         END DO                      ! ic
         END DO                         ! ib
      END DO                            ! ia
C      print*,iter,'MAX is:',DMAX,iae,ibe,ice,pae,pbe,pce
C      print*,iter,'npnt=',npnt
      RAD = DMAX
      END
C*==reget_re.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      SUBROUTINE REGET_RE(TAU,ISB,S,NATB,VSHIFT,VS,RAD)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NATB
      REAL*8 RAD
      INTEGER ISB(*)
      REAL*8 S(*),TAU(3,*),VS(*),VSHIFT(*)
C
C Local variables
C
      REAL*8 D,DMIN,V(3)
      INTEGER I,K
C
C*** End of declarations rewritten by SPAG
C
      DO K = 1,3
         V(K) = VS(K) + VSHIFT(K)
      END DO
      DMIN = 1.D10
      DO I = 1,NATB
         D = (V(1)-TAU(1,I))**2 + (V(2)-TAU(2,I))**2 + (V(3)-TAU(3,I))
     &       **2
         D = SQRT(D) - S(ISB(I))
         IF ( D.LT.0.D0 ) THEN
            PRINT *,'REGET_RE: error'
            PRINT *,I,D + S(ISB(I)),ISB(I),S(ISB(I))
            PRINT *,(V(K),K=1,3)
            PRINT *,TAU(1,I),TAU(2,I),TAU(3,I)
            RETURN
         END IF
         DMIN = MIN(D,DMIN)
      END DO                            ! natb
      RAD = DMIN
      END


C====================================================================
C  DEBUG DUMP (Fortran 77-style, no modules)
C  Write each variable to its own <name>.txt for diffing across runs.
C  - Fixed form (.f)
C  - No modules, no generics. Call the routine that matches TYPE/RANK.
C  - For arrays, pass the dimensions explicitly.
C
C  Provided routines (scalars):
C    CALL DUMP_INT_SCAL(NAME, X)
C    CALL DUMP_INT8_SCAL(NAME, X)           ! INTEGER*8
C    CALL DUMP_REAL_SCAL(NAME, X)           ! default REAL (e.g., REAL*4)
C    CALL DUMP_REAL8_SCAL(NAME, X)          ! REAL*8 (double)
C    CALL DUMP_LOG_SCAL(NAME, X)
C    CALL DUMP_CHAR_SCAL(NAME, X)
C
C  Arrays (1D/2D/3D). You must pass lengths (N), (N1,N2), (N1,N2,N3):
C    CALL DUMP_INT_1D(NAME, X, N)
C    CALL DUMP_INT_2D(NAME, X, N1, N2)
C    CALL DUMP_INT_3D(NAME, X, N1, N2, N3)
C
C    CALL DUMP_INT8_1D/2D/3D(...)
C    CALL DUMP_REAL_1D/2D/3D(...)
C    CALL DUMP_REAL8_1D/2D/3D(...)
C    CALL DUMP_LOG_1D/2D/3D(...)
C
C  Character arrays:
C    CALL DUMP_CHAR_1D(NAME, X, N)          ! X is CHARACTER*(*) X(N)
C
C  Notes:
C    - Files are overwritten if compiler supports STATUS='REPLACE'.
C      If not, STATUS='UNKNOWN' also works for most compilers.
C    - Real values use ES25.16E3 format for determinism.
C    - Arrays are written column-major: i changes fastest.
C====================================================================
      SUBROUTINE DD_OPEN(U, BASENAME)
        INTEGER            U
        CHARACTER*(*)      BASENAME
        OPEN(NEWUNIT=U, FILE=TRIM(BASENAME)//'.txt',
     &       ACTION='WRITE', STATUS='REPLACE', FORM='FORMATTED')
      RETURN
      END

      SUBROUTINE DD_WRITE_SHAPE0(U)
        INTEGER U
        WRITE(U,'(A)') '# shape: ()'
      RETURN
      END

      SUBROUTINE DD_WRITE_SHAPE1(U, N)
        INTEGER U, N
        WRITE(U,'(A,I0,A)') '# shape: (', N, ')'
      RETURN
      END

      SUBROUTINE DD_WRITE_SHAPE2(U, N1, N2)
        INTEGER U, N1, N2
        WRITE(U,'(A,I0,A,I0,A)') '# shape: (', N1, ',', N2, ')'
      RETURN
      END

      SUBROUTINE DD_WRITE_SHAPE3(U, N1, N2, N3)
        INTEGER U, N1, N2, N3
        WRITE(U,'(A,I0,A,I0,A,I0,A)') '# shape: (', N1, ',', N2,
     &                                  ',', N3, ')'
      RETURN
      END
C--------------------------------------------------------------------
C  Scalars
C--------------------------------------------------------------------
      SUBROUTINE DUMP_INT_SCAL(NAME, X)
        CHARACTER*(*) NAME
        INTEGER       X
        INTEGER U
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer'
        CALL DD_WRITE_SHAPE0(U)
        WRITE(U,'(I0)') X
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_INT8_SCAL(NAME, X)
        CHARACTER*(*) NAME
        INTEGER*8     X
        INTEGER U
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer*8'
        CALL DD_WRITE_SHAPE0(U)
        WRITE(U,'(I0)') X
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_REAL_SCAL(NAME, X)
        CHARACTER*(*) NAME
        REAL          X
        INTEGER U
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: real'
        CALL DD_WRITE_SHAPE0(U)
        WRITE(U,'(ES25.16E3)') X
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_REAL8_SCAL(NAME, X)
        CHARACTER*(*) NAME
        REAL*8        X
        INTEGER U
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: real*8'
        CALL DD_WRITE_SHAPE0(U)
        WRITE(U,'(ES25.16E3)') X
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_LOG_SCAL(NAME, X)
        CHARACTER*(*) NAME
        LOGICAL       X
        INTEGER U
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: logical'
        CALL DD_WRITE_SHAPE0(U)
        WRITE(U,'(L1)') X
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_CHAR_SCAL(NAME, X)
        CHARACTER*(*) NAME
        CHARACTER*(*) X
        INTEGER U
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: character'
        CALL DD_WRITE_SHAPE0(U)
        WRITE(U,'(A)') TRIM(X)
        CLOSE(U)
      RETURN
      END
C--------------------------------------------------------------------
C  Integer arrays
C--------------------------------------------------------------------
      SUBROUTINE DUMP_INT_1D(NAME, X, N)
        CHARACTER*(*) NAME
        INTEGER       X(N)
        INTEGER U, I
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer'
        CALL DD_WRITE_SHAPE1(U,N)
        DO 10 I = 1, N
          WRITE(U,'(I0)') X(I)
   10   CONTINUE
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_INT_2D(NAME, X, N1, N2)
        CHARACTER*(*) NAME
        INTEGER       X(N1,N2)
        INTEGER U, I, J
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer'
        CALL DD_WRITE_SHAPE2(U,N1,N2)
        DO 30 J = 1, N2
          DO 20 I = 1, N1
            WRITE(U,'(I0)') X(I,J)
   20     CONTINUE
   30   CONTINUE
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_INT_3D(NAME, X, N1, N2, N3)
        CHARACTER*(*) NAME
        INTEGER       X(N1,N2,N3)
        INTEGER U, I, J, K
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer'
        CALL DD_WRITE_SHAPE3(U,N1,N2,N3)
        DO 60 K = 1, N3
          DO 50 J = 1, N2
            DO 40 I = 1, N1
              WRITE(U,'(I0)') X(I,J,K)
   40       CONTINUE
   50     CONTINUE
   60   CONTINUE
        CLOSE(U)
      RETURN
      END
C--------------------------------------------------------------------
C  INTEGER*8 arrays
C--------------------------------------------------------------------
      SUBROUTINE DUMP_INT8_1D(NAME, X, N)
        CHARACTER*(*) NAME
        INTEGER*8     X(N)
        INTEGER U, I
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer*8'
        CALL DD_WRITE_SHAPE1(U,N)
        DO 70 I = 1, N
          WRITE(U,'(I0)') X(I)
   70   CONTINUE
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_INT8_2D(NAME, X, N1, N2)
        CHARACTER*(*) NAME
        INTEGER*8     X(N1,N2)
        INTEGER U, I, J
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer*8'
        CALL DD_WRITE_SHAPE2(U,N1,N2)
        DO 90 J = 1, N2
          DO 80 I = 1, N1
            WRITE(U,'(I0)') X(I,J)
   80     CONTINUE
   90   CONTINUE
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_INT8_3D(NAME, X, N1, N2, N3)
        CHARACTER*(*) NAME
        INTEGER*8     X(N1,N2,N3)
        INTEGER U, I, J, K
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: integer*8'
        CALL DD_WRITE_SHAPE3(U,N1,N2,N3)
        DO 120 K = 1, N3
          DO 110 J = 1, N2
            DO 100 I = 1, N1
              WRITE(U,'(I0)') X(I,J,K)
  100       CONTINUE
  110     CONTINUE
  120   CONTINUE
        CLOSE(U)
      RETURN
      END
C--------------------------------------------------------------------
C  REAL (default) arrays
C--------------------------------------------------------------------
      SUBROUTINE DUMP_REAL_1D(NAME, X, N)
        CHARACTER*(*) NAME
        REAL          X(N)
        INTEGER U, I
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: real'
        CALL DD_WRITE_SHAPE1(U,N)
        DO 130 I = 1, N
          WRITE(U,'(ES25.16E3)') X(I)
  130   CONTINUE
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_REAL_2D(NAME, X, N1, N2)
        CHARACTER*(*) NAME
        REAL          X(N1,N2)
        INTEGER U, I, J
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: real'
        CALL DD_WRITE_SHAPE2(U,N1,N2)
        DO 150 J = 1, N2
          DO 140 I = 1, N1
            WRITE(U,'(ES25.16E3)') X(I,J)
  140     CONTINUE
  150   CONTINUE
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_REAL_3D(NAME, X, N1, N2, N3)
        CHARACTER*(*) NAME
        REAL          X(N1,N2,N3)
        INTEGER U, I, J, K
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: real'
        CALL DD_WRITE_SHAPE3(U,N1,N2,N3)
        DO 180 K = 1, N3
          DO 170 J = 1, N2
            DO 160 I = 1, N1
              WRITE(U,'(ES25.16E3)') X(I,J,K)
  160       CONTINUE
  170     CONTINUE
  180   CONTINUE
        CLOSE(U)
      RETURN
      END
C--------------------------------------------------------------------
C  REAL*8 arrays
C--------------------------------------------------------------------
      SUBROUTINE DUMP_REAL8_1D(NAME, X, N)
        CHARACTER*(*) NAME
        REAL*8        X(N)
        INTEGER U, I
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: real*8'
        CALL DD_WRITE_SHAPE1(U,N)
        DO 190 I = 1, N
          WRITE(U,'(ES25.16E3)') X(I)
  190   CONTINUE
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_REAL8_2D(NAME, X, N1, N2)
        CHARACTER*(*) NAME
        REAL*8        X(N1,N2)
        INTEGER U, I, J
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: real*8'
        CALL DD_WRITE_SHAPE2(U,N1,N2)
        DO 210 J = 1, N2
          DO 200 I = 1, N1
            WRITE(U,'(ES25.16E3)') X(I,J)
  200     CONTINUE
  210   CONTINUE
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_REAL8_3D(NAME, X, N1, N2, N3)
        CHARACTER*(*) NAME
        REAL*8        X(N1,N2,N3)
        INTEGER U, I, J, K
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: real*8'
        CALL DD_WRITE_SHAPE3(U,N1,N2,N3)
        DO 240 K = 1, N3
          DO 230 J = 1, N2
            DO 220 I = 1, N1
              WRITE(U,'(ES25.16E3)') X(I,J,K)
  220       CONTINUE
  230     CONTINUE
  240   CONTINUE
        CLOSE(U)
      RETURN
      END
C--------------------------------------------------------------------
C  LOGICAL arrays
C--------------------------------------------------------------------
      SUBROUTINE DUMP_LOG_1D(NAME, X, N)
        CHARACTER*(*) NAME
        LOGICAL       X(N)
        INTEGER U, I
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: logical'
        CALL DD_WRITE_SHAPE1(U,N)
        DO 250 I = 1, N
          WRITE(U,'(L1)') X(I)
  250   CONTINUE
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_LOG_2D(NAME, X, N1, N2)
        CHARACTER*(*) NAME
        LOGICAL       X(N1,N2)
        INTEGER U, I, J
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: logical'
        CALL DD_WRITE_SHAPE2(U,N1,N2)
        DO 270 J = 1, N2
          DO 260 I = 1, N1
            WRITE(U,'(L1)') X(I,J)
  260     CONTINUE
  270   CONTINUE
        CLOSE(U)
      RETURN
      END

      SUBROUTINE DUMP_LOG_3D(NAME, X, N1, N2, N3)
        CHARACTER*(*) NAME
        LOGICAL       X(N1,N2,N3)
        INTEGER U, I, J, K
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: logical'
        CALL DD_WRITE_SHAPE3(U,N1,N2,N3)
        DO 300 K = 1, N3
          DO 290 J = 1, N2
            DO 280 I = 1, N1
              WRITE(U,'(L1)') X(I,J,K)
  280       CONTINUE
  290     CONTINUE
  300   CONTINUE
        CLOSE(U)
      RETURN
      END
C--------------------------------------------------------------------
C  CHARACTER arrays (1D list of strings)
C--------------------------------------------------------------------
      SUBROUTINE DUMP_CHAR_1D(NAME, X, N)
        CHARACTER*(*) NAME
        INTEGER       N
        CHARACTER*(*) X(N)
        INTEGER U, I
        CALL DD_OPEN(U, NAME)
        WRITE(U,'(A)') '# type: character'
        CALL DD_WRITE_SHAPE1(U,N)
        DO 310 I = 1, N
          WRITE(U,'(A)') TRIM(X(I))
  310   CONTINUE
        CLOSE(U)
      RETURN
      END
C====================================================================
C  END OF FILE
C====================================================================

      SUBROUTINE READ_IQA1(FNAME, IQA, NG)
C-----------------------------------------------------------------------
C  Reads IQA1.txt:
C    Line 1 : NG  (number of symmetry operations)
C    Next NG items: 1-based indices (one per line is fine; list-directed
C                   read also accepts multiple per line)
C
C  Inputs:
C    FNAME  - file name (CHARACTER*(*))
C    MAXNG  - max length of IQA array
C  Outputs:
C    IQA    - integer array of length NG with 1-based indices
C    NG     - number of operations read
C    IERR   - 0 OK; 1 open fail; 2 header read fail; 3 NG>MAXNG; 4 data read fail
C-----------------------------------------------------------------------
      CHARACTER*(*) FNAME
      INTEGER IQA(*), NG, MAXNG, IERR,NGINP
      INTEGER U, IOS, I

      IERR = 0
      U = 41

      OPEN(U, FILE=FNAME, STATUS='OLD', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
         IERR = 1
         STOP "CAN NOT READ IQA1.txt"
         RETURN
      ENDIF

C     Read number of operations
      READ(U, *, IOSTAT=IOS) NGINP
      IF (IOS .NE. 0) THEN
         IERR = 2
         STOP "CAN NOT READ IQA1.txt"
         GOTO 900
      ENDIF

      IF (NG .NE. NGINP) THEN
         IERR = 3
         STOP "CAN NOT READ IQA1.txt"
         GOTO 900
      ENDIF

C     Read NG integers (list-directed: spans multiple lines if needed)
      READ(U, *, IOSTAT=IOS) (IQA(I), I = 1, NG)
      IF (IOS .NE. 0) THEN
         IERR = 4
         GOTO 900
      ENDIF

  900 CONTINUE
      CLOSE(U)
      RETURN
      END

