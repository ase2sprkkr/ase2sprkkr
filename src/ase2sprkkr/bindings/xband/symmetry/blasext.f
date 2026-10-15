C*==lengths.f    processed by SPAG 6.05Rc at 10:56 on  6 Mar 2001
      INTEGER FUNCTION LENGTHS(A)
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) A
C
C Local variables
C
      INTEGER J,L
      INTEGER LEN
C
C*** End of declarations rewritten by SPAG
C
C actual length of the string
      L = LEN(A)
      J = 0
      DO WHILE ( J.LT.L .AND. A(J+1:).NE.' ' )
         J = J + 1
      END DO
      LENGTHS = J
      END
C*==dpi.f    processed by SPAG 6.05Rc at 10:56 on  6 Mar 2001
C
      DOUBLE PRECISION FUNCTION DPI()
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      LOGICAL FIRST
      DOUBLE PRECISION PI
      SAVE PI
C
C*** End of declarations rewritten by SPAG
C
      DATA FIRST/.TRUE./
C
      IF ( FIRST ) THEN
         FIRST = .FALSE.
         PI = ATAN(1.D0)*4.D0
C$$$        print *,pi
      END IF
      DPI = PI
      END
C*==dinv33.f    processed by SPAG 6.05Rc at 10:56 on  6 Mar 2001
      SUBROUTINE DINV33(MATRIX,IOPT,INVRSE,DET)
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION DET
      INTEGER IOPT
      DOUBLE PRECISION INVRSE(3,3),MATRIX(3,3)
C
C Local variables
C
      DOUBLE PRECISION DATAN
      DOUBLE PRECISION DDOT
      INTEGER I,J
      DOUBLE PRECISION XX
C
C*** End of declarations rewritten by SPAG
C
C- INVERTS 3X3 MATRIX
C ----------------------------------------------------------------
CI INPUTS
CI   INVERSE: INPUT MATRIX
CI   IOPT:  IF 0, USUAL INVERSE
CI             1, TRANSPOSE OF INVERSE
CI             2, 2*PI*INVERSE
CI             3, 2*PI*TRANSPOSE OF INVERSE
CO OUTPUTS
CO   INVERSE, AS MODIFIED ACCORDING TO IOPT
CO   DET:      DETERMINANT, OR DET/2*PI (SIGN OK ??)
CR REMARKS
CR   TO GENERATE RECIPROCAL LATTICE VECTORS, CALL DINV33(PLAT,3,RLAT)
C ----------------------------------------------------------------
C     IMPLICIT NONE
      CALL CROSS(INVRSE,MATRIX(1,2),MATRIX(1,3))
      CALL CROSS(INVRSE(1,2),MATRIX(1,3),MATRIX)
      CALL CROSS(INVRSE(1,3),MATRIX,MATRIX(1,2))
      DET = DDOT(3,MATRIX,1,INVRSE,1)
      IF ( DET.EQ.0.D0 ) CALL ENDJOB(10,'INV33: VANISHING DETERMINANT')
      IF ( IOPT.GE.2 ) DET = DET/(8*DATAN(1.D0))
      IF ( MOD(IOPT,2).EQ.0 ) THEN
         DO I = 1,3
            DO J = I + 1,3
               XX = INVRSE(I,J)
               INVRSE(I,J) = INVRSE(J,I)
               INVRSE(J,I) = XX
            END DO
         END DO
      END IF
      CALL DSCAL(9,1.D0/DET,INVRSE,1)
      END
C*==cross.f    processed by SPAG 6.05Rc at 10:56 on  6 Mar 2001
      SUBROUTINE CROSS(A,B,C)
C
C  cross product (ax,ay,az)=(bx,by,bz)*(cx,cy,cz)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION A(3),B(3),C(3)
C
C*** End of declarations rewritten by SPAG
C
C
      A(1) = B(2)*C(3) - B(3)*C(2)
      A(2) = B(3)*C(1) - B(1)*C(3)
      A(3) = B(1)*C(2) - B(2)*C(1)
C
      END
C*==dscal.f    processed by SPAG 6.05Rc at 10:56 on  6 Mar 2001
      SUBROUTINE DSCAL(N,DA,DX,INCX)
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION DA
      INTEGER INCX,N
      DOUBLE PRECISION DX(*)
C
C Local variables
C
      INTEGER I,M,MP1,NINCX
C
C*** End of declarations rewritten by SPAG
C
C
C     scales a vector by a constant.
C     uses unrolled loops for increment equal to one.
C     jack dongarra, linpack, 3/11/78.
C     modified 3/93 to return if incx .le. 0.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
C
      IF ( N.LE.0 .OR. INCX.LE.0 ) RETURN
      IF ( INCX.EQ.1 ) THEN
C
C        code for increment equal to 1
C
C
C        clean-up loop
C
         M = MOD(N,5)
         IF ( M.NE.0 ) THEN
            DO I = 1,M
               DX(I) = DA*DX(I)
            END DO
            IF ( N.LT.5 ) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
         END DO
      ELSE
C
C        code for increment not equal to 1
C
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DX(I) = DA*DX(I)
         END DO
         RETURN
      END IF
      END
C*==ddot.f    processed by SPAG 6.05Rc at 10:56 on  6 Mar 2001
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER INCX,INCY,N
      DOUBLE PRECISION DX(*),DY(*)
C
C Local variables
C
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
C
C*** End of declarations rewritten by SPAG
C
C
C     forms the dot product of two vectors.
C     uses unrolled loops for increments equal to one.
C     jack dongarra, linpack, 3/11/78.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
C
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF ( N.LE.0 ) RETURN
      IF ( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
C
C        code for both increments equal to 1
C
C
C        clean-up loop
C
         M = MOD(N,5)
         IF ( M.NE.0 ) THEN
            DO I = 1,M
               DTEMP = DTEMP + DX(I)*DY(I)
            END DO
            IF ( N.LT.5 ) THEN
               DDOT = DTEMP
               GOTO 99999
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)
     &              *DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
         END DO
         DDOT = DTEMP
      ELSE
C
C        code for unequal increments or equal increments
C          not equal to 1
C
         IX = 1
         IY = 1
         IF ( INCX.LT.0 ) IX = (-N+1)*INCX + 1
         IF ( INCY.LT.0 ) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
         DDOT = DTEMP
         RETURN
      END IF
99999 CONTINUE
      END
C*==ilun.f    processed by SPAG 6.05Rc at 10:56 on  6 Mar 2001
      FUNCTION ILUN(I)
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I
      INTEGER ILUN

      COMMON /output/ output
      integer output
C
C*** End of declarations rewritten by SPAG
C
      ILUN = output
      END
C*==endjob.f    processed by SPAG 6.05Rc at 10:56 on  6 Mar 2001
      SUBROUTINE ENDJOB(J,TEXT)
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER J
      CHARACTER*(*) TEXT
C
C*** End of declarations rewritten by SPAG
C
C-----------------
      WRITE (*,*) ' === STOP === ',J,' ',TEXT
      STOP
      END
      SUBROUTINE WKINIT(NSIZE)
      IMPLICIT  INTEGER(O)
c      include 'COMPUTER.H'
      double precision
     $ qdummy8
      COMMON /q_LMTO_/ qdummy8,ofree,LIMIT
      data  L_WORD/4/, L_WORD_I/1/, L_WORD_R/1/, L_WORD_RR/2/,
     $     L_WORD_C/2/, L_WORD_CC/4/

C ----- DEFINE STORAGE SIZE ------
C  START OF FIRST ARRAY AND MAX NUMBER TO BE DEFINED:
      LIMIT=NSIZE
      OFREE=5
      RETURN
C ------ SUBROUTINES TO DEFINE ARRAYS OF VARIOUS TYPES -----
      ENTRY DEFI(ONAME,LENG)
         LENGTH=LENG*L_WORD_I
         i8=0
         GOTO 10
      ENTRY DEFR(ONAME,LENG)
         LENGTH=LENG*L_WORD_R
         i8=0
         GOTO 10
      ENTRY DEFC(ONAME,LENG)
         LENGTH=LENG*L_WORD_C
         i8=0
         GOTO 10
      ENTRY DEFRR(ONAME,LENG)
         LENGTH=LENG*L_WORD_RR
         i8=1
         GOTO 10
      ENTRY DEFDR(ONAME,LENG)
         LENGTH=LENG*L_WORD_RR
         i8=1
         GOTO 10
      ENTRY DEFCC(ONAME,LENG)
         LENGTH=LENG*L_WORD_CC
         i8=1
  10  IF(LENGTH.LT.0) call endjob(10,'**** LENGTH OF ARRAY NEGATIVE')
      IF(LENGTH.EQ.0) LENGTH=1
      if(i8.eq.1.and.mod(ofree,2).eq.0)ofree=ofree+1
      oname=ofree
      ofree=ofree+length
c$$$      if(mod(ofree,2).eq.0)ofree=ofree+1
      if(ofree.ge.LIMIT)then
        print 1,'def0',ofree,LIMIT
        call endjob(5,' NO MEMORY')
      endif
1     format(/45('*')/'  In ',a,' memory pool exhausted',/
     >     '  N___=',I7,' NMAX_= ',i7,/45('*'))
      return
c
      ENTRY RLSE(ONAME)
      IF(ONAME.GT.LIMIT)
     .         call endjob(10,'**** RESET POINTER GT LIMIT')
      IF(ONAME.LT.3) call endjob(10,'**** RESET POINTER LT 3')
      OFREE=ONAME
      RETURN
      END
 
