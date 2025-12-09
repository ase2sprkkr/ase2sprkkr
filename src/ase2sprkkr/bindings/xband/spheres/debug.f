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

