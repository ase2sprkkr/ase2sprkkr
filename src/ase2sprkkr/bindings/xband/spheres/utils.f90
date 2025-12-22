!*==dneville.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
FUNCTION DNEVILLE(X, Y, T, X1, N)
    !
    !     Neville -Algorithmus
    !
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER N
    DOUBLE PRECISION X1
    DOUBLE PRECISION DNEVILLE
    DOUBLE PRECISION T(N), X(N), Y(N)
    !
    ! Local variables
    !
    INTEGER I, J
    !
    !*** End of declarations rewritten by SPAG
    !
    DO I = 1, N
        T(I) = Y(I)
    END DO
    !
    ! It is not clear for me, what does it mean t(0)
    !
    DO I = 1, N - 1
        DO J = N, I + 1, -1  ! <=== My changes!!!
            T(J) = T(J) + (X1 - X(J))*(T(J) - T(J - 1))/(X(J) - X(J - I))
        END DO                 ! j
    END DO                    ! i
    DNEVILLE = T(N)
END
!*==dnevmod.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
FUNCTION DNEVMOD(X, Y, X1, N, M)
    !
    !     Neville -Algorithmus (modified)
    !
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! PARAMETER definitions
    !
    INTEGER MMAX
    PARAMETER(MMAX=20)
    !
    ! Dummy arguments
    !
    INTEGER M, N
    DOUBLE PRECISION X1
    DOUBLE PRECISION DNEVMOD
    DOUBLE PRECISION X(N), Y(N)
    !
    ! Local variables
    !
    DOUBLE PRECISION DNEVILLE
    INTEGER J, JJ, MM
    DOUBLE PRECISION T(MMAX)
    !
    !*** End of declarations rewritten by SPAG
    !
    IF (M .GT. MMAX) STOP ' DNEVMOD: m>mmax'
    CALL BISEC(J, N, X1, X)
    MM = MIN(M, 2*J, 2*(N - J))
    MM = MAX(MM, 8)
    JJ = MIN(MAX(J - MM/2, 1), N - MM + 1)
    !     print*,'JJ=',jj,x(jj),' N=',N
    DNEVMOD = DNEVILLE(X(JJ), Y(JJ), T, X1, MM)
END
!*==bisec.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE BISEC(K, N, R, RMESH)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER K, N
    DOUBLE PRECISION R
    DOUBLE PRECISION RMESH(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION HR
    INTEGER IB, IE, IH
    !
    !*** End of declarations rewritten by SPAG
    !
    K = 0
    IF (R .LT. RMESH(1) .OR. R .GT. RMESH(N)) RETURN
    !
    IB = 1
    IE = N
    !
100 CONTINUE
    IH = IB + (IE - IB)/2
    !
    HR = RMESH(IH)
    IF (IE - IB .EQ. 1) THEN
        K = IB
        RETURN
    END IF
    IF (R .GT. HR) THEN
        IB = IH
        GOTO 100
    ELSE IF (R .LT. HR) THEN
        IE = IH
        GOTO 100
    END IF
    K = IH
END
!*==qd.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE QD(H, Y, Z, NDIM)
    !
    !       THE SUBROUTINE CALCULATES NDIM INTEGRALS
    !       BETWEEN X1 AND XN OF FUNCTION Z ( N.LE.NDIM
    !        XN=X1+(N-1)*H ) FOR THE SPECIAL CASE:
    !       NDIM=2*ND+1 Y(1)=0 . Z(2*I+1) IS CALCULATED
    !       BY THE SIMPSON'S METHOD AND Z(2*I+2)=Z(2*I+1)+
    !       DELT(2*I+2) . TO INCREASE A ACCURACY IT IS
    !       USED DELT(2*I+2)=DELT4(2*I+2)*(Z(2*I+3)-Z(2*I+1))/
    !       (DELT4(2*I+2)+DELT4(2*I+3)) WHERE DELT4(2*I+2)
    !       IS CALCULATED BY USING A CUBIC INTERPOLATION
    !       BETWEEN Y(2*I),Y(2*I+1),Y(2*I+2),Y(2*I+3).
    !
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION H
    INTEGER NDIM
    DOUBLE PRECISION Y(*), Z(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION C24, C3, DELT, DELT3, DELT4, DELT5, SUM2, Y1, Y2, Y3, Y4, Y5
    INTEGER I
    !
    !*** End of declarations rewritten by SPAG
    !
    IF ((NDIM/2 + NDIM/2) .EQ. NDIM) &
      &      STOP 'THE SUBPROGRAM QD WORKS ONLY IF NDIM=2*ND+1'
    !
    C3 = H/3.D0
    C24 = H/24.D0
    Z(1) = 0.D0
    IF (NDIM .NE. 1) THEN
        Y2 = Y(1)
        Y3 = Y(2)
        Y4 = Y(3)
        Z(2) = C24*(10*Y2 + 16*Y3 - 2*Y4)
        SUM2 = Y2 + 4*Y3 + Y4
        Z(3) = C3*SUM2
        IF (NDIM .NE. 3) THEN
            Y5 = Y(4)
            !
            ! the main loop of a integration
            !
            DO I = 4, NDIM - 1, 2
                Y1 = Y3
                Y2 = Y4
                Y3 = Y5
                Y4 = Y(I + 1)
                DELT3 = Y2 + 4*Y3 + Y4
                SUM2 = SUM2 + DELT3
                Z(I + 1) = C3*SUM2
                DELT4 = -Y1 + 13*Y2 + 13*Y3 - Y4
                IF (I .LT. NDIM - 1) THEN
                    Y5 = Y(I + 2)
                    DELT5 = -Y1 + 12*Y2 + 26*Y3 + 12*Y4 - Y5
                ELSE
                    DELT5 = DELT3*8
                END IF
                IF (DELT5 .GT. 1.D-20) THEN
                    DELT = 8*DELT4*DELT3/DELT5
                    Z(I) = Z(I - 1) + C24*DELT
                ELSE
                    Z(I) = Z(I - 1) + C24*DELT4
                END IF
            END DO
        END IF
    END IF
END
!*==intlog.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE INTLOG(R1, DX, Y, Y0, S, YINT, RI)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION DX, R1, S, Y0, YINT
    DOUBLE PRECISION RI(*), Y(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION A, ALF, BET, DEL, F0, F1, F2, FM1, FM2, FSR1, GAM, P, PNMA, &
      &                 PS, RALF, SS, SUM, SUM1, SUM2, X0, XNMA, XS
    DOUBLE PRECISION FI
    INTEGER I, I0, I1, I2, IM1, IMA, N, NMA
    !
    !*** End of declarations rewritten by SPAG
    !
    FI(P) = F0*P + .5D0*ALF*P**2 + BET/3.D0*P**3 + GAM/4.D0*P**4 + &
      &        DEL/5.D0*P**5
    !  PERFORMS THE INTEGRAL OF Y IN THE INTERVAL (0,S)
    !  USING SIMPSON'S RULE
    ! LOGARITHMIC NET RI(I)=R1*EXP((I-1)*DX)
    !  IF (Y0.EQ.0.) Y(R)=A*R**ALF IS ASSUMED AND CALCULATED FROM
    !  THE FIRST 2 NET POINTS
    !     WRITE(6,*) S,R1,(Y(I),I=1,10)
    XS = LOG(S/R1)
    SS = 1.D0 + XS/DX + .00001D0
    I2 = 1 + SS
    I1 = I2 - 1
    I0 = I1 - 1
    IM1 = I1 - 2
    N = I1/4
    NMA = 4*N + 1
    IMA = NMA - 3
    SUM = (Y(1)*R1 - Y(NMA)*RI(NMA))*7.D0
    !      print*,'START'
    DO I = 2, IMA, 4
        SUM = SUM + 32.D0*Y(I)*RI(I) + 12.D0*Y(I + 1)*RI(I + 1) &
          &         + 32.D0*Y(I + 2)*RI(I + 2) + 14.D0*Y(I + 3)*RI(I + 3)
    END DO
    SUM = SUM*DX*2.D0/45.D0
    !      print*,'CONT'
    IF (Y0 .EQ. 0.) THEN
        SUM1 = 0.D0
        IF (Y(1) .NE. 0.) THEN
            ALF = LOG(Y(2)/Y(1))/DX
            RALF = R1**ALF
            A = Y(1)/RALF
            SUM1 = A/(ALF + 1.D0)*RALF*R1
        END IF
    ELSE
        FSR1 = 1.D0/DX*(-11.D0*Y(1) + 18.D0*Y(2) - 9.D0*Y(3) + 2.D0*Y(4)) &
          &          /6.D0/RI(1)
        SUM1 = RI(1)*(Y0 + Y(1))*.5D0 - RI(1)**2/12.D0*FSR1
    END IF
    !  19 FORMAT(' SUB. INTLOG: Y(1),Y(2),A,ALF,SUM1 =',1P5E21.10)
    !     PRINT 19,Y(1),Y(2),A,ALF,SUM1
    !      print*,'C33',I2,Y(I2),RI(I2)
    F2 = Y(I2)*RI(I2)
    !      print*,'C--'
    F1 = Y(I1)*RI(I1)
    F0 = Y(I0)*RI(I0)
    !      print*,'C--'
    !
    FM1 = Y(IM1)*RI(IM1)
    FM2 = Y(IM1 - 1)*RI(IM1 - 1)
    ALF = (2.D0*FM2 - 16.D0*FM1 + 16.D0*F1 - 2.D0*F2)/24.D0
    BET = (-FM2 + 16.D0*FM1 - 30.D0*F0 + 16.D0*F1 - F2)/24.D0
    GAM = (-2.D0*FM2 + 4.D0*FM1 - 4.D0*F1 + 2.D0*F2)/24.D0
    DEL = (FM2 - 4.D0*FM1 + 6.D0*F0 - 4.D0*F1 + F2)/24.D0
    XNMA = (NMA - 1)*DX
    X0 = (I0 - 1)*DX
    PNMA = (XNMA - X0)/DX
    PS = (XS - X0)/DX
    SUM2 = (FI(PS) - FI(PNMA))*DX
    YINT = SUM + SUM1 + SUM2
    !      print*,'C44'
    !     PRINT 101,XS,I2,NMA
    !  101 FORMAT(' XS,I2,NMA : ',F21.15,2I5)
    !  102 FORMAT(' RADII :',//,50(/5E21.10))
    !  103 FORMAT(' INTEGRAND :',//,50(/5E21.10))
    !     PRINT 102,RI
    !     PRINT 103,Y
    !     PRINT 104,F0,ALF,BET,GAM
    !  104 FORMAT(' F0,ALF,BET,GAM,DEL : ',5E21.10)
    !  105 FORMAT(' PNMA,PS ;',2E21.10)
    !     PRINT 105,PNMA,PS
    !     WRITE(6,100) SUM,SUM1,SUM2,YINT
    !  100 FORMAT('SUM,SUM1,SUM2,YINT:',4D18.10)
END
!*==au2a.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
DOUBLE PRECISION FUNCTION AU2A()
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
    DOUBLE PRECISION A0
    PARAMETER(A0=0.529177D0)
!
!*** End of declarations rewritten by SPAG
!
    AU2A = A0
!
END
!*==incopy.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE INCOPY(N, X, INCX, Y, INCY)
    !- Copies an integer vector, x, to an integer vector, y.
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   n     :lenght of x
    !i   x     :vector to copy
    !i   incx  :incrementation for x
    !i   incy  :incrementation for y
    !o Outputs:
    !o   y     :result vector is stored in y
    !r Remarks:
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER INCX, INCY, N
    INTEGER X(*), Y(*)
    !
    ! Local variables
    !
    INTEGER I, IX, IY
    !
    !*** End of declarations rewritten by SPAG
    !
    ! Passed parameters:
    ! Local parameters:
    !
    IF (N .LE. 0) RETURN
    IF (INCX .NE. 1 .OR. INCY .NE. 1) THEN
        ! ----  code for unequal increments or equal increments not equal to 1
        IX = 1
        IY = 1
        IF (INCX .LT. 0) IX = (-N + 1)*INCX + 1
        IF (INCY .LT. 0) IY = (-N + 1)*INCY + 1
        DO I = 1, N
            Y(IY) = X(IX)
            IX = IX + INCX
            IY = IY + INCY
        END DO
    ELSE
        ! ----  code for both increments equal to 1
        DO I = 1, N
            Y(I) = X(I)
        END DO
    END IF
END
!*==errmsg.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
! ------ ERRMSG ------
SUBROUTINE ERRMSG(TEXT, I)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER I
    CHARACTER*(*) TEXT
    !
    ! Local variables
    !
    CHARACTER*12 FATAL, WARN
    !
    !*** End of declarations rewritten by SPAG
    !
    ! i<0 info
    ! i=0 warning
    ! i>0 fatal error
    DATA FATAL/'Fatal error:'/, WARN/'Warning:'/
    !
    IF (I .LE. 0) THEN
        PRINT 99001, WARN, TEXT
        !        if(iun.ne.6)write(iun,1000)warn,text
    ELSE
        PRINT 99001, FATAL, TEXT
        !        if(iun.ne.6)write(iun,1000)fatal,text
        CALL ENDJOB(-11, TEXT)
    END IF
99001 FORMAT(/2x, a12, 1x, a)
END
!*==infomsg.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
! ------ INFOMSG ------
SUBROUTINE INFOMSG(TEXT)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    CHARACTER*(*) TEXT
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    PRINT 99001, TEXT
    !      if(iun.ne.6)write(iun,1000)text
99001 FORMAT(2x, a)
END
!*==warning.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
! ------ WARNING ------
SUBROUTINE WARNING(A)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    CHARACTER*(*) A
    !
    !*** End of declarations rewritten by SPAG
    !
    PRINT '(1x,a)', A
END
!*==ddet33.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
DOUBLE PRECISION FUNCTION DDET33(A)
!- Determinant of 3x3 matrix
! ----------------------------------------------------------------------
!i Inputs:
!i   a :input matrix
!o Outputs:
!o   ddet33 :determinant
! ----------------------------------------------------------------------
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
    DOUBLE PRECISION A(*)
!
! Local variables
!
    DOUBLE PRECISION V(3)
!
!*** End of declarations rewritten by SPAG
!
! Passed parameters:
! Local parameters:
!
    CALL CROSS(V, A(4), A(7))
    DDET33 = A(1)*V(1) + A(2)*V(2) + A(3)*V(3)
END
!*==trnt.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
DOUBLE PRECISION FUNCTION TRNT(A, B, C)
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
    DOUBLE PRECISION A(3), B(3), C(3)
!
!*** End of declarations rewritten by SPAG
!
    TRNT = A(1)*B(2)*C(3) + A(2)*B(3)*C(1) + A(3)*B(1)*C(2) - A(3) &
          &       *B(2)*C(1) - A(2)*B(1)*C(3) - A(1)*B(3)*C(2)
END
!*==eilang.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE EILANG(ANG, IG)
    !
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! PARAMETER definitions
    !
    DOUBLE PRECISION P, Z, H, H3, T, T2, T4, T5
    PARAMETER(P=3.141592653589793D0, Z=0.D0, H=P/2, H3=3*H, T=P/3, T2=T*2, &
      &           T4=T*4, T5=T*5)
    !
    ! Dummy arguments
    !
    INTEGER IG
    DOUBLE PRECISION ANG(3)
    !
    ! Local variables
    !
    DOUBLE PRECISION C(3, 32)
    !
    !*** End of declarations rewritten by SPAG
    !
    DATA C/Z, Z, Z, Z, P, Z, Z, P, P, P, Z, Z, P, H, H, Z, H, H3, P, H, H3, Z, H, H, H, H, Z, H, &
      &     H, P, H3, H, P, H3, H, Z, H3, P, Z, H, Z, Z, H3, Z, Z, H, P, Z, Z, H, P, P, H, Z, Z, H, &
      &     Z, P, H, P, H3, H, H3, H3, H, H, H, H, H, H, H, H3, T, Z, Z, T2, Z, Z, T4, Z, Z, T5, Z, &
      &     Z, T4, P, Z, T5, P, Z, T, P, Z, T2, P, Z/
    ANG(1) = C(1, IG)
    ANG(2) = C(2, IG)
    ANG(3) = C(3, IG)
END
!*==turnm.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE TURNM(OM, U)
    !      INCLUDE 'PREC.H'
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION OM(3), U(3, 3)
    !
    ! Local variables
    !
    DOUBLE PRECISION COSA, COSB, COSG, SINA, SINB, SING
    !
    !*** End of declarations rewritten by SPAG
    !
    SINA = SIN(OM(1))
    SINB = SIN(OM(2))
    SING = SIN(OM(3))
    COSA = COS(OM(1))
    COSB = COS(OM(2))
    COSG = COS(OM(3))
    U(1, 1) = COSA*COSB*COSG - SINA*SING
    U(1, 2) = -COSA*COSB*SING - SINA*COSG
    U(1, 3) = COSA*SINB
    U(2, 1) = SINA*COSB*COSG + COSA*SING
    U(2, 2) = -SINA*COSB*SING + COSA*COSG
    U(2, 3) = SINA*SINB
    U(3, 1) = -SINB*COSG
    U(3, 2) = SINB*SING
    U(3, 3) = COSB
END
!*==dnrm23.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
DOUBLE PRECISION FUNCTION DNRM23(R)
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
    DOUBLE PRECISION R(3)
!
!*** End of declarations rewritten by SPAG
!
!-  'norm'-square
!
    DNRM23 = R(1)*R(1) + R(2)*R(2) + R(3)*R(3)
!
END
!*==lengths.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
INTEGER FUNCTION LENGTHS(A)
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
    CHARACTER*(*) A
!
! Local variables
!
    INTEGER J, L
!
!*** End of declarations rewritten by SPAG
!
! actual length of the string
    L = LEN(A)
    J = 0
    DO WHILE (J .LT. L .AND. A(J + 1:) .NE. ' ')
        J = J + 1
    END DO
    LENGTHS = J
END
!*==dpi.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
DOUBLE PRECISION FUNCTION DPI()
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
    LOGICAL FIRST
    DOUBLE PRECISION PI
    SAVE PI
!
!*** End of declarations rewritten by SPAG
!
    DATA FIRST/.TRUE./
!
    IF (FIRST) THEN
        FIRST = .FALSE.
        PI = ATAN(1.D0)*4.D0
        !$$$        print *,pi
    END IF
    DPI = PI
END
!*==dinv33.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE DINV33(MATRIX, IOPT, INVRSE, DET)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION DET
    INTEGER IOPT
    DOUBLE PRECISION INVRSE(3, 3), MATRIX(3, 3)
    !
    ! Local variables
    !
    DOUBLE PRECISION DDOT
    INTEGER I, J
    DOUBLE PRECISION XX
    !
    !*** End of declarations rewritten by SPAG
    !
    !- INVERTS 3X3 MATRIX
    ! ----------------------------------------------------------------
    !I INPUTS
    !I   INVERSE: INPUT MATRIX
    !I   IOPT:  IF 0, USUAL INVERSE
    !I             1, TRANSPOSE OF INVERSE
    !I             2, 2*PI*INVERSE
    !I             3, 2*PI*TRANSPOSE OF INVERSE
    !O OUTPUTS
    !O   INVERSE, AS MODIFIED ACCORDING TO IOPT
    !O   DET:      DETERMINANT, OR DET/2*PI (SIGN OK ??)
    !R REMARKS
    !R   TO GENERATE RECIPROCAL LATTICE VECTORS, CALL DINV33(PLAT,3,RLAT)
    ! ----------------------------------------------------------------
    !     IMPLICIT NONE
    CALL CROSS(INVRSE, MATRIX(1, 2), MATRIX(1, 3))
    CALL CROSS(INVRSE(1, 2), MATRIX(1, 3), MATRIX)
    CALL CROSS(INVRSE(1, 3), MATRIX, MATRIX(1, 2))
    DET = DDOT(3, MATRIX, 1, INVRSE, 1)
    IF (DET .EQ. 0.D0) CALL ENDJOB(10, 'INV33: VANISHING DETERMINANT')
    IF (IOPT .GE. 2) DET = DET/(8*DATAN(1.D0))
    IF (MOD(IOPT, 2) .EQ. 0) THEN
        DO I = 1, 3
            DO J = I + 1, 3
                XX = INVRSE(I, J)
                INVRSE(I, J) = INVRSE(J, I)
                INVRSE(J, I) = XX
            END DO
        END DO
    END IF
    CALL DSCAL(9, 1.D0/DET, INVRSE, 1)
END
!*==cross.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE CROSS(A, B, C)
    !
    !  cross product (ax,ay,az)=(bx,by,bz)*(cx,cy,cz)
    !
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION A(3), B(3), C(3)
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    A(1) = B(2)*C(3) - B(3)*C(2)
    A(2) = B(3)*C(1) - B(1)*C(3)
    A(3) = B(1)*C(2) - B(2)*C(1)
    !
END
!*==dscal.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE DSCAL(N, DA, DX, INCX)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION DA
    INTEGER INCX, N
    DOUBLE PRECISION DX(*)
    !
    ! Local variables
    !
    INTEGER I, M, MP1, NINCX
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    !     scales a vector by a constant.
    !     uses unrolled loops for increment equal to one.
    !     jack dongarra, linpack, 3/11/78.
    !     modified 3/93 to return if incx .le. 0.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    !
    IF (N .LE. 0 .OR. INCX .LE. 0) RETURN
    IF (INCX .EQ. 1) THEN
        !
        !        code for increment equal to 1
        !
        !
        !        clean-up loop
        !
        M = MOD(N, 5)
        IF (M .NE. 0) THEN
            DO I = 1, M
                DX(I) = DA*DX(I)
            END DO
            IF (N .LT. 5) RETURN
        END IF
        MP1 = M + 1
        DO I = MP1, N, 5
            DX(I) = DA*DX(I)
            DX(I + 1) = DA*DX(I + 1)
            DX(I + 2) = DA*DX(I + 2)
            DX(I + 3) = DA*DX(I + 3)
            DX(I + 4) = DA*DX(I + 4)
        END DO
    ELSE
        !
        !        code for increment not equal to 1
        !
        NINCX = N*INCX
        DO I = 1, NINCX, INCX
            DX(I) = DA*DX(I)
        END DO
        RETURN
    END IF
END
!*==ddot.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
DOUBLE PRECISION FUNCTION DDOT(N, DX, INCX, DY, INCY)
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
    INTEGER INCX, INCY, N
    DOUBLE PRECISION DX(*), DY(*)
!
! Local variables
!
    DOUBLE PRECISION DTEMP
    INTEGER I, IX, IY, M, MP1
!
!*** End of declarations rewritten by SPAG
!
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
    DDOT = 0.0D0
    DTEMP = 0.0D0
    IF (N .LE. 0) RETURN
    IF (INCX .EQ. 1 .AND. INCY .EQ. 1) THEN
        !
        !        code for both increments equal to 1
        !
        !
        !        clean-up loop
        !
        M = MOD(N, 5)
        IF (M .NE. 0) THEN
            DO I = 1, M
                DTEMP = DTEMP + DX(I)*DY(I)
            END DO
            IF (N .LT. 5) THEN
                DDOT = DTEMP
                GOTO 99999
            END IF
        END IF
        MP1 = M + 1
        DO I = MP1, N, 5
            DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) + DX(I + 2) &
                   &              *DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
        END DO
        DDOT = DTEMP
    ELSE
        !
        !        code for unequal increments or equal increments
        !          not equal to 1
        !
        IX = 1
        IY = 1
        IF (INCX .LT. 0) IX = (-N + 1)*INCX + 1
        IF (INCY .LT. 0) IY = (-N + 1)*INCY + 1
        DO I = 1, N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
        END DO
        DDOT = DTEMP
        RETURN
    END IF
99999 CONTINUE
END
!*==endjob.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE ENDJOB(J, TEXT)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER J
    CHARACTER*(*) TEXT
    !
    !*** End of declarations rewritten by SPAG
    !
    !-----------------
    WRITE (*, *) ' === STOP === ', J, ' ', TEXT
    STOP
END
!*==wkinit.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE WKINIT(NSIZE)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! COMMON variables
    !
    INTEGER LIMIT, OFREE
    DOUBLE PRECISION QDUMMY8
    COMMON/Q_LMTO_/QDUMMY8, OFREE, LIMIT
    !
    ! Dummy arguments
    !
    INTEGER LENG, NSIZE, ONAME
    !
    ! Local variables
    !
    INTEGER I8, LENGTH, L_WORD_C, L_WORD_CC, L_WORD_I, L_WORD_R, L_WORD_RR
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! COMMON variables
    !
    !
    ! Dummy arguments
    !
    !
    ! Local variables
    !
    !
    !*** End of declarations rewritten by SPAG
    !
    DATA L_WORD_I/1/, L_WORD_R/1/, L_WORD_RR/2/, L_WORD_C/2/, L_WORD_CC/4/
    !
    ! ----- DEFINE STORAGE SIZE ------
    !  START OF FIRST ARRAY AND MAX NUMBER TO BE DEFINED:
    LIMIT = NSIZE
    OFREE = 5
    RETURN
    ! ------ SUBROUTINES TO DEFINE ARRAYS OF VARIOUS TYPES -----
    ENTRY DEFI(ONAME, LENG)
    LENGTH = LENG*L_WORD_I
    I8 = 0
    GOTO 100
    ENTRY DEFR(ONAME, LENG)
    LENGTH = LENG*L_WORD_R
    I8 = 0
    GOTO 100
    ENTRY DEFC(ONAME, LENG)
    LENGTH = LENG*L_WORD_C
    I8 = 0
    GOTO 100
    ENTRY DEFRR(ONAME, LENG)
    LENGTH = LENG*L_WORD_RR
    I8 = 1
    GOTO 100
    ENTRY DEFDR(ONAME, LENG)
    LENGTH = LENG*L_WORD_RR
    I8 = 1
    GOTO 100
    ENTRY DEFCC(ONAME, LENG)
    LENGTH = LENG*L_WORD_CC
    I8 = 1
100 CONTINUE
    IF (LENGTH .LT. 0) CALL ENDJOB(10, '**** LENGTH OF ARRAY NEGATIVE')
    IF (LENGTH .EQ. 0) LENGTH = 1
    IF (I8 .EQ. 1 .AND. MOD(OFREE, 2) .EQ. 0) OFREE = OFREE + 1
    ONAME = OFREE
    OFREE = OFREE + LENGTH
    IF (OFREE .GE. LIMIT) THEN
        PRINT 99001, 'def0', OFREE, LIMIT
        CALL ENDJOB(5, ' NO MEMORY')
    END IF
    RETURN
    !
    ENTRY RLSE(ONAME)
    IF (ONAME .GT. LIMIT) CALL ENDJOB(10, '**** RESET POINTER GT LIMIT' &
      &                                  )
    IF (ONAME .LT. 3) CALL ENDJOB(10, '**** RESET POINTER LT 3')
    OFREE = ONAME
99001 FORMAT(/45('*')/'  In ', a, ' memory pool exhausted', /'  N___=', I7, &
&        ' NMAX_= ', i7, /45('*'))
END
!*==rhfds.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
! ! 03-03-94 - PERLOV : NDP cann'n calculate (0.0)**(1./3.)
! !                         changed by (0.0 + 1.d-42)**(1./3.)
! **********************************************************************
!
! HARTREE FOCK DIRAC SLATER          J P DESCLAUX     CEA PARIS 1969
! MODIFIE JUILLET 1970
!
! **********************************************************************
SUBROUTINE RHFDS(IZ, IZCORE, RO, ROM, DR00)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! COMMON variables
    !
    CHARACTER*4 BAR(10), TITRE(30)
    DOUBLE PRECISION D(251), DC(251), DCOP, DEN(30), DEXE, DEXV, DFL(30), &
      &                 DGC(251, 30), DP(251), DPAS, DPC(251, 30), DQ(251), &
      &                 DQ1(30), DR(251), DV(251), DVF(251), DVN(251), Z
    INTEGER ICUT, ION, IPRAT, MAG(30), NCORB(30), NEL(30), NES, NITER, NK(30), &
      &        NMAX(30), NORB, NP, NQL(30), NQN(30), NSTOP, NUC, NWF(30)
    DOUBLE PRECISION TEST, TESTE, TESTV, TESTY, TETS
    COMMON DEN, DQ1, DFL, NQN, NQL, NK, NMAX, NEL, NORB, MAG
    COMMON/DEUX/DVN, DVF, D, DC, DGC, DPC
    COMMON/DIRA/DV, DR, DP, DQ, DPAS, Z, NSTOP, NES, TETS, NP, NUC
    COMMON/PS2/DEXV, DEXE, DCOP, TITRE, BAR, TEST, TESTE, TESTY, TESTV, &
      &                NITER, ION, ICUT, IPRAT
    COMMON/UKP/NCORB, NWF
    COMMON/IPRINT/IPRINT
    INTEGER IPRINT
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION DR00
    INTEGER IZ, IZCORE
    DOUBLE PRECISION RO(251), ROM(251)
    !
    ! Local variables
    !
    DOUBLE PRECISION DE, DENER, DVAL
    DOUBLE PRECISION EMAX, VAL, VMAX, Y, YMAX, YN
    INTEGER I, ICEL, IM, IMAX, ITER, J, K, KK, L, N, NPCH
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! COMMON variables
    !
    !
    ! Dummy arguments
    !
    !
    ! Local variables
    !
    !
    !*** End of declarations rewritten by SPAG
    !
    NSTOP = 1
100 CONTINUE
    IF (IZ .GE. 0) THEN
        CALL INSLD(NPCH, IZ, IZCORE)
    ELSE
        CALL INSLD(NPCH, IZ, IZCORE)
    END IF
    ITER = 1
    DO I = 1, NP
        DO J = 1, NORB
            DGC(I, J) = 0.D0
            DPC(I, J) = 0.D0
        END DO
    END DO
    !WRITE (8,99001) BAR
    !WRITE (8,99002)
    N = -(ION + 1)
200 CONTINUE
    DO I = 1, NP
        D(I) = 0.D0
    END DO
    TETS = TEST
    YMAX = 0.
    VMAX = 0.
    EMAX = 0.
    ! PE{EHiE uPABHEHiq diPAKA dlq KAvdOj OPbiTAli
    ! **********************************************************************
    DO J = 1, NORB
        DE = DEN(J)
250     CONTINUE
        CALL RESLD(NQN(J), NQL(J), NK(J), IMAX, DEN(J), DFL(J), DQ1(J), J)
        IF (NSTOP .EQ. 0) THEN
            VAL = ABS((DEN(J) - DE)/DE)
            IF (VAL .GT. EMAX) EMAX = VAL
            NMAX(J) = IMAX
            DO I = 1, NP
                VAL = DGC(I, J) - DP(I)
                IF (ABS(DP(I)) .GT. 1.D0) VAL = VAL/DP(I)
                IF (ABS(VAL) .GE. ABS(YMAX)) THEN
                    YMAX = VAL
                    Y = DP(I)
                    YN = DGC(I, J)
                END IF
                VAL = DPC(I, J) - DQ(I)
                IF (ABS(DQ(I)) .GT. 1.D0) VAL = VAL/DQ(I)
                IF (ABS(VAL) .GE. ABS(YMAX)) THEN
                    YMAX = VAL
                    Y = DQ(I)
                    YN = DPC(I, J)
                END IF
                DGC(I, J) = DP(I)
                DPC(I, J) = DQ(I)
                D(I) = D(I) + NEL(J)*(DP(I)*DP(I) + DQ(I)*DQ(I))
            END DO
        ELSE IF (NSTOP .NE. 362 .OR. ITER .GE. 10 .OR. TETS .GT. TEST) THEN
            !WRITE (8,99008) NSTOP,NQN(J),TITRE(J)
            GOTO 100
        ELSE
            TETS = TESTV
            GOTO 250
        END IF
    END DO
    CALL POTSL(DC, D, DP, DR, DPAS, DEXV, Z, NP, ION, ICUT)
    IF (NUC .GT. 0) THEN
        DO I = 1, NUC
            DC(I) = DC(I) + Z/DR(I) + Z*((DR(I)/DR(NUC))**2 - 3.D0) &
              &              /(DR(NUC) + DR(NUC))
        END DO
    END IF
    DO I = 1, NP
        DVAL = ABS(DC(I) - DV(I))
        IF ((DR(I)*DC(I)) .LE. N) DVAL = -DVAL/DC(I)
        IF (DVAL .GT. VMAX) THEN
            VMAX = DVAL
            J = I
        END IF
    END DO
    !WRITE (8,99003) ITER,VMAX,DR(J),DV(J),DC(J),EMAX,YMAX,YN,Y
    !$$$      call textout(20,20,0,7,40,buf)
    IF (TETS .GT. TEST .OR. EMAX .GT. TESTE .OR. VMAX .GT. TESTV .OR. &
      &     YMAX .GT. TESTY) THEN
        ITER = ITER + 1
        IF (ITER .LE. NITER) THEN
            ! pOTEHciAl dlq ClEdu`}Ej iTEPAcii
            ! **********************************************************************
            DVAL = 1.D0 - DCOP
            DO I = 1, NP
                DVN(I) = DV(I)
                DVF(I) = DC(I)
                DV(I) = DVAL*DV(I) + DCOP*DC(I)
            END DO
            GOTO 200
        ELSE
            !WRITE (8,99004) NITER
            NSTOP = 2
        END IF
    END IF
    IF (NPCH .NE. 0) THEN
        !WRITE (7,99009) BAR
        !WRITE (7,99010) (DR(KK),D(KK),KK=1,251)
    END IF
    IF (IPRINT > 0) WRITE (6, 99005)
!WRITE (8,99005)
! CPEdHiE zHA~EHiq R
! **********************************************************************
    DO I = 1, NP
        DVF(I) = DC(I)
        DQ(I) = 0.
    END DO
    DVAL = 0.
    DO I = 1, NORB
        IM = NMAX(I)
        DVAL = DVAL + NEL(I)*DEN(I)
        DO J = 1, IM
            DC(J) = DGC(J, I)*DGC(J, I) + DPC(J, I)*DPC(J, I)
        END DO
        L = 5
        IF (IABS(NK(I)) .EQ. 1) L = L - 1
        DO J = 1, L
            DP(J) = DFL(I) + DFL(I)
            IF (J .EQ. 1) THEN
                N = 4
            ELSE IF (J .EQ. 2) THEN
                N = 2
            ELSE IF (J .EQ. 3) THEN
                N = 1
            ELSE IF (J .EQ. 4) THEN
                N = -1
            ELSE
                N = -3
            END IF
            CALL SOMM(DR, DC, DQ, DPAS, DP(J), N, IM)
        END DO
        !         WRITE (8,99006) NQN(I),TITRE(I),DEN(I)*2,NEL(I),NCORB(I),
        !     &                   (DP(J),J=2,3)
        IF (IPRINT > 0) &
          &    WRITE (6, 99006) NQN(I), TITRE(I), DEN(I)*2, NEL(I), NCORB(I), &
          &                   (DP(J), J=2, 3)
    END DO
! pOlHAq |HEPgiq uCPEdHEHHAq pO CfEPE
! **********************************************************************
    DC(1) = 1.D0
    DO I = 1, NP
        DP(I) = D(I)/DR(I)
    END DO
    IF (NUC .GT. 0) THEN
        DO I = 1, NUC
            DP(I) = D(I)*(3.D0 - DR(I)*DR(I)/(DR(NUC)*DR(NUC))) &
              &              /(DR(NUC) + DR(NUC))
        END DO
        DC(1) = 4.D0
    END IF
    CALL SOMM(DR, DP, DQ, DPAS, DC(1), 0, NP)
    DO I = 1, NP
        DP(I) = D(I)*DVF(I)
        D(I) = D(I)*((D(I)*DR(I) + 1.D-42)**(1.D0/3.D0))
    END DO
    DC(2) = 3.D0
    DC(3) = 1.D0
    IF (NUC .NE. 0) DC(3) = 4.D0
    CALL SOMM(DR, DP, DQ, DPAS, DC(3), 0, NP)
    CALL SOMM(DR, D, DQ, DPAS, DC(2), -1, NP)
    DC(2) = -3.D0*DC(2)/(105.27578D0**(1.D0/3.D0))
    DC(1) = -Z*DC(1)
    DC(4) = DVAL - DC(3)
    DVAL = DVAL + (DC(1) - DC(3) + (DEXE - DEXV)*DC(2))/2.D0
    DC(3) = (DC(3) - DC(1) - DEXV*DC(2))/2.D0
    DC(2) = DC(2)*DEXE/2.D0
!     WRITE(8,1007) DVAL,DC(4),DC(3),DC(2),DC(1)
    IF (NORB .NE. 1) THEN
        !     WRITE(8,1001) BAR
        !     WRITE(8,1008)
        ! iHTEgPAly pEPEKPyTiq
        ! **********************************************************************
        DO I = 2, NORB
            K = I - 1
            DO J = 1, K
                IF (NQL(I) .EQ. NQL(J) .AND. NK(I) .EQ. NK(J)) THEN
                    IM = NMAX(J)
                    IF (NMAX(I) .LT. IM) IM = NMAX(I)
                    DO L = 1, IM
                        DQ(L) = DPC(L, I)*DPC(L, J)
                        DC(L) = DGC(L, I)*DGC(L, J)
                    END DO
                    DVAL = DFL(I) + DFL(J)
                    CALL SOMM(DR, DC, DQ, DPAS, DVAL, 0, IM)
                    !                  IF ( DVAL.GT.1.D-3 ) WRITE (8,99007) NQN(I),TITRE(I),
                    !     &                 NQN(J),TITRE(J),DVAL
                END IF
            END DO
        END DO
    END IF
!     WRITE(8,4003)
    DO J = 1, NORB
        IF (NWF(J) .NE. 0) THEN
            DENER = DEN(J)*2
            !            WRITE (41,99011) NQN(J),NQL(J),DENER,DR(1)
            !            WRITE (41,99012) (DGC(I,J),I=1,250)
        END IF
    END DO
!
    ICEL = 0
    DO J = 1, NORB
        IF (NCORB(J) .NE. 0) ICEL = ICEL + NEL(J)
    END DO
!
    DR00 = DR(1)
    DO I = 1, 251
        RO(I) = 0.D0
        ROM(I) = 0.D0
    END DO
    DO J = 1, NORB
        !
        DO I = 1, 251
            RO(I) = RO(I) + NEL(J)*(DGC(I, J)*DGC(I, J) + DPC(I, J)*DPC(I, J))
            ROM(I) = ROM(I) + MAG(J) &
                    &               *(DGC(I, J)*DGC(I, J) + DPC(I, J)*DPC(I, J))
        END DO
    END DO
!
99001 FORMAT(5X, 10A4/)
99002 FORMAT(' ITER', 4X, 'DVMAX', 10X, 'R', 14X, 'VN-1', 13X, 'VN', 10X, &
&        'DEMAX', 6X, 'DPMAX', 9X, 'PN-1', 13X, 'PN')
99003 FORMAT(I5, 1PE21.10, 3(E21.10), 2(1PE21.10), 2(E21.10))
99004 FORMAT(' Number of iterations is exceeded', i4)
99005 FORMAT('Level  Energy      Occup. Val.      R**2          R')
99006 FORMAT(I2, A2, 1X, E21.10, 2X, I2, 3X, I2, 2X, E21.10, 2X, E21.10)
99007 FORMAT(34X, I1, A2, I3, A2, F21.15)
99008 FORMAT('  NSTOP=', I4, '  dlq OPbiTAli   ', I3, A2)
99009 FORMAT(10A4)
99010 FORMAT(4D15.8)
99011 FORMAT(I3, 1X, I3, 1X, G14.7, 1X, G14.7)
99012 FORMAT(5(E21.10))
END
!*==insld.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE INSLD(NPCH, IZ, IZCORE)
    !
    ! ~TEHiE CTAPTOBOgO pOTEHciAlA
    ! pOClEdHqq KAPTA dOlvHA COdEPvATx 4 zBEzdO~Ki
    ! **********************************************************************
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! COMMON variables
    !
    CHARACTER*4 BAR(10), TITRE(30)
    DOUBLE PRECISION DCOP, DEN(30), DEXE, DEXV, DFL(30), DP(251), DPAS, &
      &                 DQ(251), DQ1(30), DR(251), DV(251), Z
    INTEGER ICUT, IEX, ION, IPRAT, IWAT, MAG(30), NCORB(30), NEL(30), NES, &
      &        NITER, NK(30), NMAX(30), NORB, NP, NQL(30), NQN(30), NSTOP, NUC, &
      &        NWF(30)
    DOUBLE PRECISION RWAT, TEST, TESTE, TESTV, TESTY, TETS
    COMMON DEN, DQ1, DFL, NQN, NQL, NK, NMAX, NEL, NORB, MAG
    COMMON/DIRA/DV, DR, DP, DQ, DPAS, Z, NSTOP, NES, TETS, NP, NUC
    COMMON/HLS/RWAT, IWAT, IEX
    COMMON/PS2/DEXV, DEXE, DCOP, TITRE, BAR, TEST, TESTE, TESTY, TESTV, &
      &                NITER, ION, ICUT, IPRAT
    COMMON/UKP/NCORB, NWF
    !
    ! Dummy arguments
    !
    INTEGER IZ, IZCORE, NPCH
    !
    ! Local variables
    !
    DOUBLE PRECISION D1, DR1, DVAL, DVC
    CHARACTER*4 ENDA, ITXCH(3, 2), TTIRE(9)
    DOUBLE PRECISION FPOT
    INTEGER I, ICH, IDEP, IMCH, ION1, IZ1, J, K, L, NKAI, NUC1
    DOUBLE PRECISION R, VAL
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! COMMON variables
    !
    !
    ! Dummy arguments
    !
    !
    ! Local variables
    !
    !
    !*** End of declarations rewritten by SPAG
    !
    DATA TTIRE/'S   ', 'P*  ', 'P   ', 'D*  ', 'D   ', 'F*  ', 'F   ', &
      &     'G*  ', 'G   '/, ENDA/'****'/
    DATA ITXCH/'BART', 'H-HE', 'DIN ', 'X-AL', 'PHA ', '    '/
    !      IF (NSTOP.EQ.0) GO TO 2
    IF (NSTOP .EQ. 1) THEN
        DPAS = 0.05D0
        DR1 = 1D-2
        !0.01D0
        NES = 15
        NITER = 50
        TESTE = 5.D-06
        TESTY = 1.D-05
        TESTV = 1.D-05
        TEST = 1.D-07
        NP = 251
        NSTOP = 30
        DEXV = 1.D0
        DEXE = 1.5D0
        DCOP = 0.3D0
        DVC = 137.0359895D0            !137.0373D0
        IZ1 = 0
        ION1 = 0
        NUC1 = -1
        ION = 0
        IDEP = 0
        ICUT = 0
        IPRAT = 1
        IEX = 0
        I = 0
        J = 0
        L = 0
        K = 0
        NUC = 0
        NPCH = 0
        IWAT = 0
        DVAL = 0.D0
        DO I = 1, 7
            BAR(I) = '----'
        END DO
        !      READ (5,1000) (BAR(I),I=1,7)
        !
        ! BAR TITRE COdEPvAT        24 CiMBOlA
        ! **********************************************************************
        IF (BAR(1) .EQ. ENDA) CALL EXIT_A
        I = 0
        CALL DEFORB(NQN, NK, NEL, NCORB, NORB, IZ, IZCORE, MAG)
        !
        ! IZ ATOMHyj HOMEP      ION=IZ-~iClO |lEKTPOHOB
        ! NORB ~iClO OPbiTAlEj    IDEP dOlvH byTx PABEH 0 ili 1
        ! IDEP=0 CTAPTOByj pOTHciAl =pOTEHciAl TOMACA-fEPMi
        ! IDEP=1 CTAPTOByj pOTEHciAl ~iTAETCq C KAPT
        ! ECli ICUT HOlx UL ON CORRIGE LE POTENTIEL EN -(ION+1)/R
        ! ECli IPRAT HOlx iCpOlxzyETCq pPOcEduPA pPATTA
        ! ECli IEX HOlx - ObMEH  BARTH-HEDIN'A
        ! I ~iClO TO~EK dlq iHTEgPiPOBAHiq (251 ECli I=0
        ! J ~iClO pOpyTOK pOdOgHATx |HEPgi`(15 ECli J=0
        ! K ~iClO iTEPAcij (50 ECli K=0
        ! L=0 CTAHdAPTHO BybPAHHAq TO~HOCTx
        ! KOHE~HyE PAzMEPy qdPA ECli NUC pOlOviTElEH
        ! **********************************************************************
        DO ICH = 8, 10
            IMCH = ICH - 7
            BAR(ICH) = ITXCH(IMCH, IEX + 1)
        END DO
        IF (NORB .LE. NSTOP) THEN
            IF (I .GT. 0) THEN
                I = 2*(I/2) + 1
                IF (I .LE. NP) THEN
                    NP = I
                ELSE
                    !                 WRITE (8,99011) I
                    GOTO 100
                END IF
            END IF
            IF (J .GT. 0) NES = J
            IF (K .GT. 0) NITER = K
            IF (IEX .EQ. 0) THEN
            END IF
            !      READ (5,1005) DEXV,DEXE
            !
            ! DEXV KO|fficiEHT pPi ObMEHHOM pOTEHciAlE      DEXV=1. dlq  SLATER
            ! DEXE KO|fficiEHT pPi ObMEHHOj |HEPgii
            ! DEXV dOlvEH byTx PABEH 2.*DEXE/3. dlq TEOPEMy BiPiAlA
            ! **********************************************************************
            IF (L .EQ. 0) THEN
            END IF
            !      READ (5,1002) DPAS,DR1,TEST,TESTE,TESTY,TESTV
            !
            ! DPAS |KCp.{Ag    DR1 OpPEdElqET pEPBu` TO~Ku  =DR1/IZ
            ! TEST TO~HOCTx pO |HEPgii dlq RESLD
            ! TESTE KPiTEPij CAMOCOglACOBAHiq OdHO|lEKTPOHHyX |HEPgij
            ! TESTE KPiTEPij CAMOCOglACOBAHiq BOlHOByX fuHKcij
            ! TESTE KPiTEPij CAMOCOglACOBAHiq pOTEHciAlA
            ! **********************************************************************
            IF (IPRAT .NE. 0) DCOP = 0.3D0
            !
            ! VI(N+1)=(1.-DCOP)*VI(N)+DCOP*VF(N)
            ! **********************************************************************
            Z = DFLOAT(IZ)
            IF (NUC .GT. 0) THEN
                !
                ! DVAL ATOMHAq MACCA ECli NUC pOlOviTElEH
                ! **********************************************************************
                DVAL = Z*(DVAL**(1.D0/3.D0))*2.2677D-05/EXP(4.D0*DPAS)
                IF (DVAL .LE. DR1) THEN
                    DR1 = DVAL
                    NUC = 5
                    GOTO 20
                ELSE
                    DVAL = DVAL*EXP(4.D0*DPAS)
                    DO I = 6, NP
                        D1 = DR1*EXP(DFLOAT(I - 1)*DPAS)
                        IF (D1 .GE. DVAL) GOTO 10
                    END DO
                    !                  WRITE (8,99018)
                    GOTO 100
                END IF
10              CONTINUE
                NUC = I
                DR1 = DR1*DVAL/D1
            END IF
20          CONTINUE
            !            WRITE (8,99002) IZ,ION,NITER,TESTE,TESTY,TESTV
            !            WRITE (8,99009) NP,DR1,IZ,DPAS
            !            WRITE (8,99010) TEST,NES
            !            IF ( IEX.EQ.0 ) WRITE (8,99001)
            !            IF ( IEX.NE.0 ) WRITE (8,99013) DEXV,DEXE
            K = 0
            DVAL = Z*Z/(DVC*DVC)
            !            IF ( NUC.GT.0 ) WRITE (8,99017)
            !
            ! dAHHyE pO OPbiTAlqM
            ! **********************************************************************
            DO I = 1, NORB
                !
                ! DEN |HEPgiq OPbiTAli B Ed.lAHdAu (<0)
                ! NQN glABHOE KBAHTOBOE ~iClO      NK KBAHTOBOE ~iClO KAppA
                ! NEL zAHqTOCTx OPbiTAli NCORB pPizHAK BAlEHTHOCTi(=1 dlq BAl.
                ! **********************************************************************
                K = K + NEL(I)
                DEN(I) = -Z*Z/(4.D0*NQN(I)*NQN(I))
                NQL(I) = IABS(NK(I))
                IF (NK(I) .LT. 0) NQL(I) = NQL(I) - 1
                IF (NUC .GT. 0) THEN
                    NKAI = IABS(NK(I))
                    DFL(I) = DFLOAT(NKAI)
                ELSE
                    DFL(I) = NK(I)*NK(I)
                    DFL(I) = SQRT(DFL(I) - DVAL)
                END IF
                L = 2*IABS(NK(I))
                IF (NQL(I) .LT. NQN(I) .AND. NEL(I) .LE. L .AND. NQN(I) &
                  &              .GT. 0 .AND. NQL(I) .LE. 4) THEN
                    J = NQL(I) + IABS(NK(I))
                    TITRE(I) = TTIRE(J)
                    !                  WRITE (8,99006) NQN(I),TITRE(I),NCORB(I),NEL(I),DEN(I)
                ELSE
                    !                  WRITE (8,99004) DEN(I),NQN(I),NQL(I),J,NEL(I)
                    GOTO 100
                END IF
            END DO
            IF (K .EQ. (IZ - ION)) THEN
                IF (IPRAT .NE. 0) THEN
                    !                  WRITE (8,99016) DCOP
                ELSE
                    !                  WRITE (8,99015)
                END IF
                IF (NUC .EQ. NUC1) THEN
                    IF (IZ .EQ. IZ1 .AND. ION .EQ. ION1) GOTO 50
                    IF (IZ .EQ. IZ1) GOTO 40
                END IF
                DR(1) = DR1/Z
                DO I = 2, NP
                    DR(I) = DR(1)*EXP(DFLOAT(I - 1)*DPAS)
                END DO
            ELSE
                !               WRITE (8,99005)
                GOTO 100
            END IF
            !
            ! CTAPTOByj pOTEHciAl
            ! **********************************************************************
40          CONTINUE
            VAL = -ION - 1
            IF (IDEP .EQ. 1) THEN
                !      READ(5,1004) (DV(I),I=1,NP)
                !
                ! ~TEHiE CTAPTOBOgO pOTEHciAlA   (EN U.A. ET NEGATIF) ECli IDEP=1
                ! **********************************************************************
                !               WRITE (8,99012) BAR,(DV(I),I=1,NP)
                DVAL = -Z/DV(1)
                IF (NUC .GT. 0) DVAL = 1.D0
                DO I = 1, NP
                    DV(I) = DV(I)*DVAL/DR(I)
                END DO
            ELSE IF (IDEP .EQ. 0) THEN
                IF (IZ .NE. IZ1 .OR. ION .LE. ION1 .OR. NUC .NE. NUC1) THEN
                    DO I = 1, NP
                        R = DR(I)
                        DV(I) = FPOT(R, Z, VAL)
                    END DO
                    IF (IWAT .NE. 0) THEN
                        DO I = 1, NP
                            IF (DR(I) .GT. RWAT) THEN
                                DV(I) = DV(I) + ION/DR(I)
                            ELSE
                                DV(I) = DV(I) + ION/RWAT
                            END IF
                        END DO
                    END IF
                    IF (NUC .GT. 0) THEN
                        DO I = 1, NUC
                            DV(I) = DV(I) + Z/DR(I) &
                              &                          + Z*((DR(I)/DR(NUC))**2 - 3.D0) &
                              &                          /(DR(NUC) + DR(NUC))
                        END DO
                    END IF
                    GOTO 50
                END IF
            ELSE
                !               WRITE (8,99014)
                GOTO 100
            END IF
            IF (ICUT .EQ. 0) THEN
                DO I = 1, NP
                    IF ((DR(I)*DV(I)) .GT. VAL) DV(I) = VAL/DR(I)
                END DO
            END IF
            VAL = Z + DV(1)*DR(1)
            IF (NUC .GT. 0) VAL = Z + DV(NUC)*DR(NUC)
            IF (ABS(VAL) .GE. 0.1D0) THEN
                !               WRITE (8,99007)
                GOTO 100
            END IF
        ELSE
            !            WRITE (8,99003) NORB
            GOTO 100
        END IF
50      CONTINUE
        IF (NORB .NE. 1) THEN
            DO I = 2, NORB
                K = I - 1
                DO J = 1, K
                    IF (NQN(I) .EQ. NQN(J) .AND. NK(I) .EQ. NK(J)) THEN
                        !                    WRITE (8,99008)
                        GOTO 100
                    END IF
                END DO
            END DO
        END IF
        IZ1 = IZ
        ION1 = ION
        NUC1 = NUC
        DO I = 1, NORB
            NMAX(I) = NP
            L = 1
            J = NQN(I) - NQL(I)
            IF ((J - 2*(J/2)) .EQ. 0) L = -L
            NKAI = IABS(NK(I))
            DQ1(I) = DFLOAT(L*NK(I)/NKAI)
            IF (NUC .NE. 0) THEN
                IF (NK(I) .LT. 0) DQ1(I) = DQ1(I)*(NK(I) - DFL(I))*DVC/Z
            END IF
        END DO
        RETURN
    END IF
100 CONTINUE
!      WRITE (8,99019)
!    READ (5,1000) BAR(1)
    IF (BAR(1) .EQ. ENDA) THEN
    END IF
    NSTOP = 1
99001 FORMAT(' Exchange: Barth-Hedin   ')
99002 FORMAT(' Atomic number  ', I3, '   Ionicity  ', I2/1X, &
&        'Maximal number of iterations ', &
&        I4/' Precision in energy ', E21.10/14x, 'wave function ', &
&        E21.10/14X, 'potential', E21.10/)
99003 FORMAT(' NORB=', I3, 'too large ************** ')
99004 FORMAT(' Input error       ', E21.10, I1, 2I2)
99005 FORMAT(' Erroneous number of electrons  **************')
99006 FORMAT(7X, I1, A2, 2I8, 1PE21.10)
99007 FORMAT(' Error in potential   ')
99008 FORMAT(' Bad configuration ')
99009 FORMAT(' Integration will be in ', i3, &
&        ' points'/' the first point is ', f7.4, '/', i2, '  step is ', &
&        f7.4,/)
!$$$ 2008 FORMAT ('  iHTEgPiPOBAHiE budET pPOBOdiTxCq B  ',I3,
!$$$     = ' TO4KAX ',/,'    pEPByj 4lEH PABEH    ',F7.4,'/'
!$$$     =,I2,' {Ag PABEH   ',F7.4,/)
99010 FORMAT(' In the program RESLD the relative precision in energy', &
&        ' is ', 1pe9.2, /' the number of attempts ', i3,/)
!/,'
!    2IDEP=',I3,2X,'ICUT=',I3,2X,'IPRAT=',I3,2X,'ION=',I3,2X,'IWAT=',
!    1I3,2X,'RWAT=',E10.3,/)
!$$$ 2009 FORMAT ('       pPi  PE{EHii uPABHEHiq diPAKA ',
!$$$     */,'    OTHOCiTElxHAq TO4HOCTx pO |HEPgii ',E21.10,
!$$$     */,'                4iClO pOpyTOK         ',I3,/)
99011 FORMAT('  NP=', I3, ' is too large *******')
99012 FORMAT(20X, 10A4, //, 5X, 'Starting potential * R'/, 5(2X, F21.15))
99013 FORMAT('  Exchange: Slater X-alpha, DEXV=', f21.15, '     DEXE=', &
&        f21.15,/)
99014 FORMAT(' Error in  IDEP *********************')
99015 FORMAT(' Pratt''s mixing is used  '/)
99016 FORMAT(' Potential is mixed with  ADMIX=', E21.10,/)
99017 FORMAT(10X, ' Finite nucleus'/)
99018 FORMAT(' Error in atomic weight  **************** ')
99019 FORMAT(' The next case ')
!      GO TO 1
END
!*==somm.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE SOMM(DR, DP, DQ, DPAS, DA, M, NP)
    !
    ! iHTEgPiPOBAHiE METOdOM CiMpCOHA  (DP+DQ)*DR**M OT 0 dO R=DR(NP)
    ! DPAS |KCpOHEHciAlxHyj {Ag dlq R->0   (DP+DQ)=CTE*R**DA
    ! **********************************************************************
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION DA, DPAS
    INTEGER M, NP
    DOUBLE PRECISION DP(*), DQ(*), DR(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION D1, DB, DC, DL
    INTEGER I, MM
    !
    !*** End of declarations rewritten by SPAG
    !
    MM = M + 1
    D1 = DA + DFLOAT(MM)
    DA = 0.D0
    DB = 0.D0
    DO I = 1, NP
        DL = DR(I)**MM
        IF (I .NE. 1 .AND. I .NE. NP) THEN
            DL = DL + DL
            IF ((I - 2*(I/2)) .EQ. 0) DL = DL + DL
        END IF
        DC = DP(I)*DL
        IF (DC .LT. 0) THEN
            DB = DB + DC
        ELSE IF (DC .NE. 0) THEN
            DA = DA + DC
        END IF
        DC = DQ(I)*DL
        IF (DC .LT. 0) THEN
            DB = DB + DC
        ELSE IF (DC .NE. 0) THEN
            DA = DA + DC
        END IF
    END DO
    DA = DPAS*(DA + DB)/3.D0
    DC = EXP(DPAS) - 1.D0
    DB = D1*(D1 + 1.D0)*DC*EXP((D1 - 1.D0)*DPAS)
    DB = DR(1)*(DR(2)**M)/DB
    DC = (DR(1)**MM)*(1.D0 + 1.D0/(DC*(D1 + 1.D0)))/D1
    DA = DA + DC*(DP(1) + DQ(1)) - DB*(DP(2) + DQ(2))
END
!*==fpot.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
FUNCTION FPOT(R, Z, WA)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION R, WA
    DOUBLE PRECISION Z
    DOUBLE PRECISION FPOT
    !
    ! Local variables
    !
    DOUBLE PRECISION S, WC, WD, WE
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! pOTEHciAl TOMACA-fEPMi B TO~KE  R  Z ATOMHyj HOMEP
    ! WA ~iClO |lEKTPOHOB-Z-1
    ! **********************************************************************
    S = Z
    WC = SQRT((R*(S + WA)**(1.D0/3.D0))/0.8853D0)
    WD = WC*(0.60112*WC + 1.81061D0) + 1.D0
    WE = WC*(WC*(WC*(WC*(0.04793D0*WC + 0.21465D0) + 0.77112D0) + 1.39515D0) &
      &     + 1.81061D0) + 1.D0
    WC = (Z + WA)*(WD/WE)**2 - WA
    FPOT = -WC/R
END
!*==potsl.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE POTSL(DV, D, DP, DR, DPAS, DEXV, Z, NP, ION, ICUT)
    !
    ! iHTEgPiPOBAHiE pOTEHciAlA pO 4 TO~KAM
    ! DV pOTEHciAl   D plOTHOCTx DP BLOC DE TRAVAIL  DR PAdiAlxHAq {KAlA
    ! DPAS |KCp.{Ag          DEXV MHOviTElx dlq ObMEHA
    ! Z ATOMHyj HOMEP     NP ~iClO TO~EK        ION=Z-~iClO |lEKTPOHOB
    ! SI ICUT EST NUL ON CORRIGE EVENTUELLEMENT LE POTENTIEL EN -(ION+1)/R
    ! **********************************************************************
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! COMMON variables
    !
    INTEGER IEX, IWAT
    DOUBLE PRECISION RWAT
    COMMON/HLS/RWAT, IWAT, IEX
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION DEXV, DPAS, Z
    INTEGER ICUT, ION, NP
    DOUBLE PRECISION D(*), DP(*), DR(*), DV(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION DARS, DAS, DBRS, DCF, DCP, DLO, DLO2, DNY, DR2, DS, DS1, &
      &                 DSF, DSF2, DSF3, DSP, DSP2, DSP3, DVXC
    INTEGER I, J, K
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! COMMON variables
    !
    !
    ! Dummy arguments
    !
    !
    ! Local variables
    !
    !
    !*** End of declarations rewritten by SPAG
    !
    DAS = DPAS/24.D0
    DO I = 1, NP
        DV(I) = D(I)*DR(I)
    END DO
    DLO = EXP(DPAS)
    DLO2 = DLO*DLO
    DP(2) = DR(1)*(D(2) - D(1)*DLO2)/(12.D0*(DLO - 1.D0))
    DP(1) = DV(1)/3.D0 - DP(2)/DLO2
    DP(2) = DV(2)/3.D0 - DP(2)*DLO2
    J = NP - 1
    DO I = 3, J
        DP(I) = DP(I - 1) + DAS*(13.D0*(DV(I) + DV(I - 1)) - (DV(I - 2) + DV(I + 1)))
    END DO
    DP(NP) = DP(J)
    DV(J) = DP(J)
    DV(NP) = DP(J)
    DO I = 3, J
        K = NP + 1 - I
        DV(K) = DV(K + 1) &
          &           /DLO + DAS*(13.D0*(DP(K + 1)/DLO + DP(K)) - (DP(K + 2)/DLO2 + &
          &           DP(K - 1)*DLO))
    END DO
    DV(1) = DV(3)/DLO2 + DPAS*(DP(1) + 4.D0*DP(2)/DLO + DP(3)/DLO2)/3.D0
    DLO = -DFLOAT(ION + 1)
    IF (IEX .EQ. 0) THEN
        DO I = 1, NP
            !
            !     pOTEHciAl pO BARTH-HEDIN
            !  *********************************************************************
            DR2 = DR(I)**2
            DS1 = (D(I)/(3.D0*DR2) + 1.D-42)**(1.D0/3.D0)
            !
            IF (ABS(DS1) .LT. 1.D-10) THEN
                DVXC = 0.D0
            ELSE
                DS = 1.D0/DS1
                DSF = DS/75.D0
                DSF2 = DSF*DSF
                DSF3 = DSF2*DSF
                DSP = DS/30.D0
                DSP2 = DSP*DSP
                DSP3 = DSP2*DSP
                DCF = (1.D0 + DSF3)*LOG(1.D0 + 1.D0/DSF) + 0.5D0*DSF - DSF2 - &
                  &               0.3333333333D0
                DCP = (1.D0 + DSP3)*LOG(1.D0 + 1.D0/DSP) + 0.5D0*DSP - DSP2 - &
                  &               0.3333333333D0
                DNY = 5.1297628D0*(0.0504D0*DCP - 0.0254D0*DCF)
                DARS = -1.22177412D0/DS + DNY
                DBRS = -0.0504D0*LOG(1.D0 + 30.D0/DS) - DNY
                DVXC = DARS + DBRS
            END IF
            DV(I) = DV(I) - (Z - 0.5D0*DR(I)*DVXC)
            !  *********************************************************************
            IF (ICUT .EQ. 0) THEN
                IF (DV(I) .GT. DLO) DV(I) = DLO
            END IF
            DV(I) = DV(I)/DR(I)
        END DO
        IF (IWAT .EQ. 0) RETURN
        DO I = 1, NP
            IF (DR(I) .GT. RWAT) THEN
                DV(I) = DV(I) + ION/DR(I)
            ELSE
                DV(I) = DV(I) + ION/RWAT
            END IF
        END DO
        GOTO 99999
    END IF
    !
    !     pOTEHciAl pO Cl|TTEPu
    !
    DO I = 1, NP
        DV(I) = DV(I) &
          &           - (Z + 3.D0*DEXV*((DR(I)*D(I)/105.27578D0 + 1.D-42)**(1.D0/ &
          &           3.D0)))
        IF (ICUT .EQ. 0) THEN
            IF (DV(I) .GT. DLO) DV(I) = DLO
        END IF
        DV(I) = DV(I)/DR(I)
    END DO
    IF (IWAT .EQ. 0) RETURN
    DO I = 1, NP
        IF (DR(I) .GT. RWAT) THEN
            DV(I) = DV(I) + ION/DR(I)
        ELSE
            DV(I) = DV(I) + ION/RWAT
        END IF
    END DO
    RETURN
99999 CONTINUE
END
!*==deforb.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE DEFORB(XN, NK, XZ, IVAL, NORB, IZN, IZC, MAG)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER IZC, IZN, NORB
    INTEGER IVAL(30), MAG(30), NK(30), XN(30), XZ(30)
    !
    ! Local variables
    !
    INTEGER I, IXT, IXV, IXZ1(105), IXZ2(105), IXZ3(105), IXZ4(105), &
      &        IXZ5(105), IZ, IZV, J, JF1, JFUL, JFULR(105), JREL(105), NMAG
    DOUBLE PRECISION XJ(35), XL(35), XMAG, XML, XXJ(31), XXL(31), XXN(31), XXZ(31), ZN
    !
    !*** End of declarations rewritten by SPAG
    !
    !=====================================================================
    !
    !     STORED DATA FOR FULL ORBITALS (AS THEY ARE NUMBERED IN 'ATOMGL'):
    !
    !     STORED DATA FOR RELATIVISTIC ORBITALS:
    !
    DATA XXN/1.0D0, 2.0D0, 2.0D0, 2.0D0, 3.0D0, 3.0D0, 3.0D0, 4.0D0, 3.0D0, &
      &     3.0D0, 4.0D0, 4.0D0, 5.0D0, 4.0D0, 4.0D0, 5.0D0, 5.0D0, 6.0D0, 4.0D0, &
      &     4.0D0, 5.0D0, 5.0D0, 6.0D0, 6.0D0, 7.0D0, 5.0D0, 5.0D0, 6.0D0, 6.0D0, &
      &     7.0D0, 7.0D0/
    DATA XXL/0.0D0, 0.0D0, 1.0D0, 1.0D0, 0.0D0, 1.0D0, 1.0D0, 0.0D0, 2.0D0, &
      &     2.0D0, 1.0D0, 1.0D0, 0.0D0, 2.0D0, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 3.0D0, &
      &     3.0D0, 2.0D0, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 3.0D0, 3.0D0, 2.0D0, 2.0D0, &
      &     1.0D0, 1.0D0/
    DATA XXJ/0.5D0, 0.5D0, 0.5D0, 1.5D0, 0.5D0, 0.5D0, 1.5D0, 0.5D0, 1.5D0, &
      &     2.5D0, 0.5D0, 1.5D0, 0.5D0, 1.5D0, 2.5D0, 0.5D0, 1.5D0, 0.5D0, 2.5D0, &
      &     3.5D0, 1.5D0, 2.5D0, 0.5D0, 1.5D0, 0.5D0, 2.5D0, 3.5D0, 1.5D0, 2.5D0, &
      &     0.5D0, 1.5D0/
    DATA XXZ/2.0D0, 2.0D0, 2.0D0, 4.0D0, 2.0D0, 2.0D0, 4.0D0, 2.0D0, 4.0D0, &
      &     6.0D0, 2.0D0, 4.0D0, 2.0D0, 4.0D0, 6.0D0, 2.0D0, 4.0D0, 2.0D0, 6.0D0, &
      &     8.0D0, 4.0D0, 6.0D0, 2.0D0, 4.0D0, 2.0D0, 6.0D0, 8.0D0, 4.0D0, 6.0D0, &
      &     2.0D0, 4.0D0/
    DATA JFULR/0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 8, 8, 7, 9, 9, 9, &
      &     9, 7, 10, 10, 11, 11, 11, 11, 12, 12, 13, 13, 13, 12, 12, 14, 12, 12, 12, 12, 15, &
      &     15, 16, 16, 16, 16, 17, 17, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, &
      &     19, 20, 20, 20, 20, 21, 21, 21, 21, 17, 17, 22, 22, 23, 23, 23, 23, 24, 24, 25, &
      &     25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27/
    DATA JREL/1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 9, 10, 10, 10, &
      &     10, 10, 10, 10, 11, 11, 12, 12, 12, 12, 13, 13, 14, 14, 14, 15, 15, 15, 15, 15, &
      &     15, 15, 16, 16, 17, 17, 17, 17, 18, 18, 21, 21, 19, 19, 19, 19, 20, 21, 20, 20, &
      &     20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 23, 23, 24, 24, 24, 24, &
      &     25, 25, 28, 28, 28, 28, 28, 26, 27, 28, 28, 27, 27, 27, 27, 27, 28, 28, 28/
    DATA IXZ1/1, 0, 1, 0, 1, 0, 1, 2, 3, 0, 1, 0, 1, 0, 1, 2, 3, 0, 1, 0, 1, 2, 3, 1, 1, 2, 3, 4, &
      &     1, 0, 1, 0, 1, 2, 3, 0, 1, 0, 1, 2, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 2, 3, 0, 1, 0, 0, 1, &
      &     3, 4, 5, 0, 1, 1, 3, 4, 5, 6, 7, 0, 1, 2, 3, 0, 1, 2, 3, 1, 1, 0, 1, 0, 1, 2, 3, 0, 1, 0, &
      &     0, 0, 2, 3, 4, 0, 1, 1, 2, 4, 5, 6, 7, 0, 1, 2, 3/
    DATA IXZ2/23*0, 4, 4*0, 4, 11*0, 4, 4, 0, 4, 4, 4, 4, 16*0, 1, 13*0, 6, 6, 16*0, 1, &
      &     1, 8*0/
    DATA IXZ3/23*0, 1, 4*0, 6, 12*0, 1, 0, 3, 4, 6, 6, 9*0, 1, 1, 19*0, 8, 8, 9*0, 1, 2, &
      &     1, 1, 1, 12*0/
    DATA IXZ4/77*0, 4, 4, 26*0/
    DATA IXZ5/77*0, 5, 6, 26*0/
    !=====================================================================
    ZN = IZN
    IZ = INT(ZN + 0.000001D0)
    IZV = IZN - IZC
    !------ J-REPRESENTATION ------
    JFUL = JFULR(IZ)
    !
    J = JREL(IZ)
    DO I = 1, JFUL
        XN(I) = NINT(XXN(I))
        XL(I) = XXL(I)
        XJ(I) = XXJ(I)
        XZ(I) = NINT(XXZ(I))
    END DO
    !
    XZ(JFUL + 1) = IXZ1(IZ)
    XZ(JFUL + 2) = IXZ2(IZ)
    XZ(JFUL + 3) = IXZ3(IZ)
    XZ(JFUL + 4) = IXZ4(IZ)
    XZ(JFUL + 5) = IXZ5(IZ)
    !
    !--------------------------------
    IF (JFUL .NE. J) THEN
        JF1 = JFUL + 1
        DO I = JF1, J
            XN(I) = NINT(XXN(I))
            XL(I) = XXL(I)
            XJ(I) = XXJ(I)
            !--------------------------------
        END DO
    END IF
    !=====================================================================
    IXV = 0
    IXT = 0
    DO I = J, 1, -1
        IXT = IXT + XZ(I)
        NK(I) = NINT(XL(I))
        MAG(I) = 0
        XMAG = 0
        IF (XJ(I) .GT. XL(I)) NK(I) = NINT(-XL(I) - 1)
        IF (IXV + XZ(I) .LE. IZV) THEN
            IF (I .GE. 31) THEN
                WRITE (6, *) '########### POTSL: I = ', I, ' > 30'
                STOP
            END IF
            IVAL(I) = 1
            IXV = IXV + XZ(I)
            IF (NK(I) .GT. 1) THEN
                XML = XZ(I) + XZ(I + 1)
                NMAG = NK(I)*2 + 1
                XMAG = NMAG - ABS(XML - NMAG)
                MAG(I) = NINT(XMAG)
                !            XMAG1=(nmag-1)*xmag/(2.*nmag)
                !            XMAG2=xmag-xmag1
                !            print*,xz(i),xz(i+1)
                !            print*, i,nk(i),' l=',nint(xl(i)),' n='
                !     $           ,xn(i),' val=',nint(XMAG),xmag1,xmag2
            END IF
        ELSE
            IVAL(I) = 0
        END IF
    END DO
    IF (IXV .NE. IZV) STOP ' error - I can not make initial data '
    NORB = J
END
!*==exit_a.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE EXIT_A
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    !*** End of declarations rewritten by SPAG
    !
END
!*==resld.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE RESLD(NQN, NQL, NK, IMAX, DE, DFL, DQ1, JC1)
    !
    ! PE{EHiE uPABHEHiq diPAKA
    ! NQN NglABHOE KBAHTOBOE ~iClO   NQL OPbiTAlxHOE KBAHTOBOE ~iClO
    ! NK KBAHTOBOE ~iClO KAppA    IMAX pOClEdHqq TO~KA B TAbulqcii BOl-
    ! HOBOj fuHKcii     DE |HEPgiq   DFL pOKAzATElx CTEpEHi B PAzlOvEHii
    ! BOlHOBOj fuHKcii       DQ1 OTHO{EHiE DP K DQ B HA~AlE KOOPdiHAT
    ! **********************************************************************
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! COMMON variables
    !
    DOUBLE PRECISION DB, DEP(5), DEQ(5), DK, DM, DP(251), DPAS, DPNO(4, 30), &
      &                 DQ(251), DQNO(4, 30), DR(251), DSAL, DV(251), DVC, Z
    INTEGER NES, NP, NSTOP, NUC
    DOUBLE PRECISION TEST
    COMMON/DIRA/DV, DR, DP, DQ, DPAS, Z, NSTOP, NES, TEST, NP, NUC
    COMMON/PS1/DEP, DEQ, DB, DVC, DSAL, DK, DM
    COMMON/TROIS/DPNO, DQNO
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION DE, DFL, DQ1
    INTEGER IMAX, JC1, NK, NQL, NQN
    !
    ! Local variables
    !
    DOUBLE PRECISION DBE, DD, DKOEF, DPM, DPQ, DQM, DSUM, DVAL
    DOUBLE PRECISION ELIM, VAL
    INTEGER I, IES, IMAT, IMM, J, JC, K, LLL, M, ND, NOEUD
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! COMMON variables
    !
    !
    ! Dummy arguments
    !
    !
    ! Local variables
    !
    !
    !*** End of declarations rewritten by SPAG
    !
    ! DV pOTEHciAl B AT.Ed.(<0)        DR PAdiAlxHAq CETKA
    ! DP bOlx{Aq KOMpOHEHTA  DQ MAlAq KOMpOHEHTA    DPAS |KCpOHEHc.{Ag
    ! Z ATOMHyj HOMEP     NSTOP KOHTPOlx ~iClEHHOgO iHTEgPiPOBAHiq
    ! NES MAKCiMAlxHOE ~iClO iTEPAcij pO OpPEdElEHi` |HEPgii
    ! TEST TO~HOCTx pOlu~EHiq |HEPgii          NP NMAKCiMAlxHOE ~iClO TO~EK
    ! KOHE~HyE PAzMEPy qdPA ECli    NUC HE   0
    ! **********************************************************************
    !
    ! DEP,DEQ pPOizBOdHyE DP i  DQ   DB=ENERGIE/DVC    DVC CKOPOCTx CBETA
    !  B AT.Ed.         DSAL=2.*DVC   DK KBAHTOBOE ~iClO KAppA
    ! DM=|KCpOHEHciAlxHyj {Ag/720., DKOEF=1./720.
    ! **********************************************************************
    DATA DKOEF/.1388888888888888D-2/
    JC = 1
    IF (JC1 .GT. 0) JC = JC1
    NSTOP = 0
    DVC = 137.0373D0
    DSAL = DVC + DVC
    IMM = 0
    IES = 0
    DK = NK
    LLL = (NQL*(NQL + 1))/2
    ND = 0
    NOEUD = NQN - NQL
    IF (LLL .NE. 0) THEN
        ELIM = DV(1) + LLL/(DR(1)*DR(1))
        DO I = 2, NP
            VAL = DV(I) + LLL/(DR(I)*DR(I))
            IF (VAL .LE. ELIM) ELIM = VAL
        END DO
        IF (ELIM .LT. 0) THEN
            IF (DE .LE. ELIM) DE = ELIM*0.5D0
        ELSE
            NSTOP = 17
            ! 2*V+L*(L+1)/R**2   BC`du pOlOviTElEH
            ! **********************************************************************
            RETURN
        END IF
    ELSE
        ELIM = -Z*Z/(1.5*NQN*NQN)
        IF (DE .LE. ELIM) DE = ELIM*0.5D0
    END IF
100 CONTINUE
    IF (IMM .NE. 1) THEN
        DO I = 7, NP, 2
            IMAT = NP + 1 - I
            IF ((DV(IMAT) + DFLOAT(LLL)/(DR(IMAT)*DR(IMAT)) - DE) .LE. 0.D0) &
              &           EXIT
        END DO
        IF (IMAT .LE. 5) THEN
            DE = DE*0.5D0
            IF (DE .LT. -TEST .AND. ND .LE. NOEUD) GOTO 100
            NSTOP = 28
            ! 2*V+L*(L+1)/R**2-2*E BC`du pOlOviTElEH
            ! **********************************************************************
            RETURN
        END IF
    END IF
    ! HA~AlxHyE zHA~EHiq dlq iHTEgPiPOBAHiq pO BHuTPEHHEj OblACTi
    ! **********************************************************************
    DB = DE/DVC
    CALL INOUH(DP, DQ, DR, DQ1, DFL, DV(1), Z, TEST, NUC, NSTOP, JC)
    IF (NSTOP .NE. 0) GOTO 99999
    !     NSTOP=45
    ! HET CXOdiMOCTi B HA~AlE KOOPdiHAT
    ! **********************************************************************
    ND = 1
    DO I = 1, 5
        DVAL = DR(I)**DFL
        IF (I .NE. 1) THEN
            IF (DP(I - 1) .NE. 0.D0) THEN
                IF ((DP(I)/DP(I - 1)) .LE. 0.D0) ND = ND + 1
            END IF
        END IF
        DP(I) = DP(I)*DVAL
        DQ(I) = DQ(I)*DVAL
        DEP(I) = DEP(I)*DVAL
        DEQ(I) = DEQ(I)*DVAL
    END DO
    K = -1 + 2*(NOEUD - 2*(NOEUD/2))
    IF ((DP(1)*DFLOAT(K)) .GT. 0.D0) THEN
        IF ((DFLOAT(K)*DFLOAT(NK)*DQ(1)) .GE. 0.D0) THEN
            DM = DPAS*DKOEF
            ! IiHTEgPiPOBAHiE pO BHuTPEHHEj OblACTi
            ! **********************************************************************
            DO I = 6, IMAT
                DP(I) = DP(I - 1)
                DQ(I) = DQ(I - 1)
                CALL INTH(DP(I), DQ(I), DV(I), DR(I))
                IF (DP(I - 1) .NE. 0.D0) THEN
                    IF ((DP(I)/DP(I - 1)) .LE. 0.D0) THEN
                        ND = ND + 1
                        IF (ND .GT. NOEUD) GOTO 200
                    END IF
                END IF
            END DO
            IF (ND .EQ. NOEUD) THEN
                ! HA~AlHyE zHA~EHiq dlq iHTgPiPOBAHiq pO BHE{HEj OblACTi
                ! **********************************************************************
                DQM = DQ(IMAT)
                DPM = DP(IMAT)
                IF (IMM .NE. 1) THEN
                    DO I = 1, NP, 2
                        IMAX = NP + 1 - I
                        IF (((DV(IMAX) - DE)*DR(IMAX)*DR(IMAX)) .LE. 300.D0) &
                          &                    EXIT
                    END DO
                END IF
                DD = SQRT(-DE*(2.D0 + DB/DVC))
                DPQ = -DD/(DSAL + DB)
                DM = -DM
                DO I = 1, 5
                    J = IMAX + 1 - I
                    DP(J) = EXP(-DD*DR(J))
                    DEP(I) = -DD*DP(J)*DR(J)
                    DQ(J) = DPQ*DP(J)
                    DEQ(I) = DPQ*DEP(I)
                END DO
                M = IMAX - 5
                ! iHTEgPiPOBAHiE pO BHE{HEj OblACTi
                !***********************************************************************
                DO I = IMAT, M
                    J = M + IMAT - I
                    DP(J) = DP(J + 1)
                    DQ(J) = DQ(J + 1)
                    CALL INTH(DP(J), DQ(J), DV(J), DR(J))
                END DO
                ! C{iBKA bOlx{Oj KOMpOHEHTy
                ! **********************************************************************
                DVAL = DPM/DP(IMAT)
                IF (DVAL .GT. 0.D0) THEN
                    DO I = IMAT, IMAX
                        DP(I) = DP(I)*DVAL
                        DQ(I) = DQ(I)*DVAL
                    END DO
                    ! By~iClEHiE HOPMy
                    ! **********************************************************************
                    DSUM = 3.D0*DR(1)*(DP(1)**2 + DQ(1)**2) &
                      &                   /(DPAS*(DFL + DFL + 1.D0))
                    DO I = 3, IMAX, 2
                        DSUM = DSUM + DR(I)*(DP(I)**2 + DQ(I)**2) &
                          &                      + 4.D0*DR(I - 1)*(DP(I - 1)**2 + DQ(I - 1)**2) &
                          &                      + DR(I - 2)*(DP(I - 2)**2 + DQ(I - 2)**2)
                    END DO
                    DSUM = DPAS*(DSUM + DR(IMAT)*(DQM*DQM - DQ(IMAT)*DQ(IMAT)) &
                      &                   )*0.3333333333333333D0
                    ! izMEHEHiE |HEPgii
                    ! **********************************************************************
                    DBE = DP(IMAT)*(DQM - DQ(IMAT))*DVC/DSUM
                    IMM = 0
                    VAL = ABS(DBE/DE)
                    IF (VAL .LE. TEST) THEN
                        IF (JC1 .LE. 0 .AND. IMAX .GT. (-JC1)) THEN
                            IMAX = -JC1
                            DSUM = 3.D0*DR(1)*(DP(1)**2 + DQ(1)**2) &
                              &                         /(DPAS*(DFL + DFL + 1.D0))
                            DO I = 3, IMAX, 2
                                DSUM = DSUM + DR(I)*(DP(I)**2 + DQ(I)**2) &
                                  &                            + 4.D0*DR(I - 1)*(DP(I - 1)**2 + DQ(I - 1)**2) &
                                  &                            + DR(I - 2)*(DP(I - 2)**2 + DQ(I - 2)**2)
                            END DO
                            DSUM = DPAS*(DSUM + DR(IMAT) &
                                        &                         *(DQM*DQM - DQ(IMAT)*DQ(IMAT))) &
                                  &                         *0.3333333333333333D0
                        END IF
                        DSUM = SQRT(DSUM)
                        DQ1 = DQ1/DSUM
                        DO I = 1, IMAX
                            DP(I) = DP(I)/DSUM
                            DQ(I) = DQ(I)/DSUM
                        END DO
                        DO I = 1, 4
                            DPNO(I, JC) = DPNO(I, JC)/DSUM
                            DQNO(I, JC) = DQNO(I, JC)/DSUM
                        END DO
                        IF (IMAX .NE. NP) THEN
                            J = IMAX + 1
                            DO I = J, NP
                                DP(I) = 0.D0
                                DQ(I) = 0.D0
                            END DO
                        END IF
                        NSTOP = 0
                        GOTO 99999
                    ELSE
102                     CONTINUE
                        DVAL = DE + DBE
                        IF (DVAL .LT. 0.D0) THEN
                            DE = DVAL
                            IF (VAL .LE. 0.1D0) IMM = 1
                            IES = IES + 1
                            IF (IES .LE. NES) GOTO 100
                            NSTOP = 362
                            ! pPEBy{EH ~iClO iTEPAcij
                            ! **********************************************************************
                            RETURN
                        ELSE
                            DBE = DBE*0.5D0
                            VAL = VAL*0.5D0
                            IF (VAL .GT. TEST) GOTO 102
                            NSTOP = 345
                            ! HulEBAq |HEPgiq
                            ! **********************************************************************
                            RETURN
                        END IF
                    END IF
                ELSE
                    NSTOP = 312
                    ! O{ibKA B zHAKE bOlx{Oj KOMpOHEHTy
                    ! **********************************************************************
                    RETURN
                END IF
            ELSE
                DE = 0.8D0*DE
                IF (DE .LT. -TEST) GOTO 100
                NSTOP = 206
                ! ~iClO uzlOB Cli{KOM MAlO
                ! **********************************************************************
                RETURN
            END IF
        END IF
    END IF
    NSTOP = 53
    ! O{ibKA B PAzlOvEHii B HEA~AlE KOOPdiHAT
    ! **********************************************************************
    RETURN
200 CONTINUE
    DE = 1.2D0*DE
    IF (DE .GT. ELIM) GOTO 100
    NSTOP = 210
    ! ~iClO uzlOB Cli{KOM BEliKO
    ! **********************************************************************
    RETURN
99999 CONTINUE
END
!*==inouh.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE INOUH(DP, DQ, DR, DQ1, DFL, DV, Z, TEST, NUC, NSTOP, JC)
    !
    ! HA~AlxHyE zHA~EHiq dlq iHTEgPiPOBAHiq pO BHuTPEHHEj OblACTi
    ! DP bOlx{Aq KOMpOHEHTA   DQ MAlAq KOMpOHEHTA     DR PAdiAlxHAq.CETKA
    ! DQ1 OTHO{EHiE DP K DQ B HA~AlE KOOPd.DFL pOKAzATElx CTEpEHi glABHOgO
    ! ~lEHA PAzlOvEHiq B HA~AlE DV pOTEHciAl B pEPBOj TO~KE
    ! Z ATOMHyj HOMEP      TEST TO~HOCTx
    ! KOHE~HyE PAzMEPy qdPA ECli    NUC HE 0
    ! NSTOP KOHTPOlx CXOdiMOCTi pPi PAzlOvEHii B Pqd
    ! **********************************************************************
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! COMMON variables
    !
    DOUBLE PRECISION DD, DEP(5), DEQ(5), DK, DM, DPNO(4, 30), DQNO(4, 30), &
      &                 DSAL, DVC
    COMMON/PS1/DEP, DEQ, DD, DVC, DSAL, DK, DM
    COMMON/TROIS/DPNO, DQNO
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION DFL, DQ1, DV, Z
    INTEGER JC, NSTOP, NUC
    DOUBLE PRECISION TEST
    DOUBLE PRECISION DP(*), DQ(*), DR(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION DBE, DEVA1, DEVA2, DEVA3, DPR, DQR, DSUM, DVAL
    INTEGER I, J, M
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! COMMON variables
    !
    !
    ! Dummy arguments
    !
    !
    ! Local variables
    !
    !
    !*** End of declarations rewritten by SPAG
    !
    ! DEP,DEQ pPOizBOdHyE DP i  DQ   DD=|HEPgiq/DVC    DVC CKOPOCTx CBETA
    !  B AT. Ed.        DSAL=2.*DVC   DK KBAHTOBOE ~iClO KAppA
    ! DM=|KCp.{Ag/720.
    ! **********************************************************************
    DO I = 1, 10
        DP(I) = 0.D0
        DQ(I) = 0.D0
    END DO
    IF (NUC .LE. 0) THEN
        DVAL = Z/DVC
        DEVA1 = -DVAL
        DEVA2 = DV/DVC + DVAL/DR(1) - DD
        DEVA3 = 0.D0
        IF (DK .LE. 0) THEN
            DBE = (DK - DFL)/DVAL
        ELSE
            DBE = DVAL/(DK + DFL)
        END IF
        DQ(10) = DQ1
        DP(10) = DBE*DQ1
    ELSE
        DVAL = DV + Z*(3.D0 - DR(1)*DR(1)/(DR(NUC)*DR(NUC))) &
          &          /(DR(NUC) + DR(NUC))
        DEVA1 = 0.D0
        DEVA2 = (DVAL - 3.D0*Z/(DR(NUC) + DR(NUC)))/DVC - DD
        DEVA3 = Z/(DR(NUC)*DR(NUC)*DR(NUC)*DSAL)
        IF (DK .LE. 0) THEN
            DP(10) = DQ1
        ELSE
            DQ(10) = DQ1
        END IF
    END IF
    DO I = 1, 5
        DP(I) = DP(10)
        DQ(I) = DQ(10)
        DEP(I) = DP(I)*DFL
        DEQ(I) = DQ(I)*DFL
    END DO
    M = 1
100 CONTINUE
    DM = M + DFL
    DSUM = DM*DM - DK*DK + DEVA1*DEVA1
    DQR = (DSAL - DEVA2)*DQ(M + 9) - DEVA3*DQ(M + 7)
    DPR = DEVA2*DP(M + 9) + DEVA3*DP(M + 7)
    DVAL = ((DM - DK)*DQR - DEVA1*DPR)/DSUM
    DSUM = ((DM + DK)*DPR + DEVA1*DQR)/DSUM
    J = -1
    DO I = 1, 5
        DPR = DR(I)**M
        DQR = DSUM*DPR
        DPR = DVAL*DPR
        IF (M .NE. 1) THEN
            IF (ABS(DPR/DP(I)) .LE. TEST .AND. ABS(DQR/DQ(I)) .LE. TEST) &
              &           J = 1
        END IF
        DP(I) = DP(I) + DPR
        DQ(I) = DQ(I) + DQR
        DEP(I) = DEP(I) + DPR*DM
        DEQ(I) = DEQ(I) + DQR*DM
    END DO
    IF (J .NE. 1) THEN
        DP(M + 10) = DVAL
        DQ(M + 10) = DSUM
        M = M + 1
        IF (M .LE. 20) GOTO 100
        NSTOP = 45
    END IF
    DO I = 1, 4
        DPNO(I, JC) = DP(I + 9)
        DQNO(I, JC) = DQ(I + 9)
    END DO
END
!*==inth.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE INTH(DP, DQ, DV, DR)
    !
    ! iHTEgPiPOBAHiE pO METOdu AdAMCA pO 5-TO~KAM bOlx{Oj KOMpOHEHTy  DP i
    ! MAlOj KOMpOHEHTy DQ   B TO~KE DR  ; dv-pOTEHciAl B |TOj TO~KE
    ! **********************************************************************
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! COMMON variables
    !
    DOUBLE PRECISION DB, DEP(5), DEQ(5), DK, DM, DSAL, DVC
    COMMON/PS1/DEP, DEQ, DB, DVC, DSAL, DK, DM
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION DP, DQ, DR, DV
    !
    ! Local variables
    !
    DOUBLE PRECISION DKOEF1, DKOEF2, DPR, DQR, DSUM
    INTEGER I
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! COMMON variables
    !
    !
    ! Dummy arguments
    !
    !
    ! Local variables
    !
    !
    !*** End of declarations rewritten by SPAG
    !
    ! DEP,DEQ pPOizBOdHyE DP i  DQ   DD=|HEPgiq/DVC    DVC CKOPOCTx CBETA
    !  B AT. Ed.        DSAL=2.*DVC   DK KBAHTOBOE ~iClO KAppA
    ! DM=|KCp.{Ag/720.
    ! DKOEF1=405./502., DKOEF2=27./502.
    ! **********************************************************************
    DATA DKOEF1/.9462151394422310D0/, DKOEF2/.5378486055776890D-1/
    DPR = DP + DM*((251.D0*DEP(1) + 2616.D0*DEP(3) + 1901.D0*DEP(5)) &
      &      - (1274.D0*DEP(2) + 2774.D0*DEP(4)))
    DQR = DQ + DM*((251.D0*DEQ(1) + 2616.D0*DEQ(3) + 1901.D0*DEQ(5)) &
      &      - (1274.D0*DEQ(2) + 2774.D0*DEQ(4)))
    DO I = 2, 5
        DEP(I - 1) = DEP(I)
        DEQ(I - 1) = DEQ(I)
    END DO
    DSUM = (DB - DV/DVC)*DR
    DEP(5) = -DK*DPR + (DSAL*DR + DSUM)*DQR
    DEQ(5) = DK*DQR - DSUM*DPR
    DP = DP + DM*((106.D0*DEP(2) + 646.D0*DEP(4) + 251.D0*DEP(5)) &
      &     - (19.D0*DEP(1) + 264.D0*DEP(3)))
    DQ = DQ + DM*((106.D0*DEQ(2) + 646.D0*DEQ(4) + 251.D0*DEQ(5)) &
      &     - (19.D0*DEQ(1) + 264.D0*DEQ(3)))
    DP = DKOEF1*DP + DKOEF2*DPR
    DQ = DKOEF1*DQ + DKOEF2*DQR
    DEP(5) = -DK*DP + (DSAL*DR + DSUM)*DQ
    DEQ(5) = DK*DQ - DSUM*DP
END
!*==dlamch.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
DOUBLE PRECISION FUNCTION DLAMCH(CMACH)
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
    DOUBLE PRECISION ONE, ZERO
    PARAMETER(ONE=1.0D+0, ZERO=0.0D+0)
!
! Dummy arguments
!
    CHARACTER CMACH
!
! Local variables
!
    DOUBLE PRECISION BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN, RND, &
      &                 SFMIN, SMALL, T
    INTEGER BETA, IMAX, IMIN, IT
    LOGICAL FIRST, LRND
    LOGICAL LSAME
    SAVE BASE, EMAX, EMIN, EPS, PREC, RMAX, RMIN, RND, SFMIN, T
    EXTERNAL DLAMC2, LSAME
!
!*** End of declarations rewritten by SPAG
!
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLAMCH determines double precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by DLAMCH:
!          = 'E' or 'e',   DLAMCH := eps
!          = 'S' or 's ,   DLAMCH := sfmin
!          = 'B' or 'b',   DLAMCH := base
!          = 'P' or 'p',   DLAMCH := eps*base
!          = 'N' or 'n',   DLAMCH := t
!          = 'R' or 'r',   DLAMCH := rnd
!          = 'M' or 'm',   DLAMCH := emin
!          = 'U' or 'u',   DLAMCH := rmin
!          = 'L' or 'l',   DLAMCH := emax
!          = 'O' or 'o',   DLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Save statement ..
!     ..
!     .. Data statements ..
    DATA FIRST/.TRUE./
!     ..
!     .. Executable Statements ..
!
    IF (FIRST) THEN
        FIRST = .FALSE.
        CALL DLAMC2(BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX)
        BASE = BETA
        T = IT
        IF (LRND) THEN
            RND = ONE
            EPS = (BASE**(1 - IT))/2
        ELSE
            RND = ZERO
            EPS = BASE**(1 - IT)
        END IF
        PREC = EPS*BASE
        EMIN = IMIN
        EMAX = IMAX
        SFMIN = RMIN
        SMALL = ONE/RMAX
        !
        !           Use SMALL plus a bit, to avoid the possibility of rounding
        !           causing overflow when computing  1/sfmin.
        !
        IF (SMALL .GE. SFMIN) SFMIN = SMALL*(ONE + EPS)
    END IF
!
    IF (LSAME(CMACH, 'E')) THEN
        RMACH = EPS
    ELSE IF (LSAME(CMACH, 'S')) THEN
        RMACH = SFMIN
    ELSE IF (LSAME(CMACH, 'B')) THEN
        RMACH = BASE
    ELSE IF (LSAME(CMACH, 'P')) THEN
        RMACH = PREC
    ELSE IF (LSAME(CMACH, 'N')) THEN
        RMACH = T
    ELSE IF (LSAME(CMACH, 'R')) THEN
        RMACH = RND
    ELSE IF (LSAME(CMACH, 'M')) THEN
        RMACH = EMIN
    ELSE IF (LSAME(CMACH, 'U')) THEN
        RMACH = RMIN
    ELSE IF (LSAME(CMACH, 'L')) THEN
        RMACH = EMAX
    ELSE IF (LSAME(CMACH, 'O')) THEN
        RMACH = RMAX
    END IF
!
    DLAMCH = RMACH
!
!     End of DLAMCH
!
END
!*==dlamc1.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
!***********************************************************************
!
SUBROUTINE DLAMC1(BETA, T, RND, IEEE1)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER BETA, T
    LOGICAL IEEE1, RND
    !
    ! Local variables
    !
    DOUBLE PRECISION A, B, C, F, ONE, QTR, SAVEC, T1, T2
    DOUBLE PRECISION DLAMC3
    LOGICAL FIRST, LIEEE1, LRND
    INTEGER LBETA, LT
    SAVE LBETA, LIEEE1, LRND, LT
    EXTERNAL DLAMC3
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    !  -- LAPACK auxiliary routine (version 2.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     October 31, 1992
    !
    !     .. Scalar Arguments ..
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DLAMC1 determines the machine parameters given by BETA, T, RND, and
    !  IEEE1.
    !
    !  Arguments
    !  =========
    !
    !  BETA    (output) INTEGER
    !          The base of the machine.
    !
    !  T       (output) INTEGER
    !          The number of ( BETA ) digits in the mantissa.
    !
    !  RND     (output) LOGICAL
    !          Specifies whether proper rounding  ( RND = .TRUE. )  or
    !          chopping  ( RND = .FALSE. )  occurs in addition. This may not
    !          be a reliable guide to the way in which the machine performs
    !          its arithmetic.
    !
    !  IEEE1   (output) LOGICAL
    !          Specifies whether rounding appears to be done in the IEEE
    !          'round to nearest' style.
    !
    !  Further Details
    !  ===============
    !
    !  The routine is based on the routine  ENVRON  by Malcolm and
    !  incorporates suggestions by Gentleman and Marovich. See
    !
    !     Malcolm M. A. (1972) Algorithms to reveal properties of
    !        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
    !
    !     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
    !        that reveal properties of floating point arithmetic units.
    !        Comms. of the ACM, 17, 276-277.
    !
    ! =====================================================================
    !
    DATA FIRST/.TRUE./
    !
    IF (FIRST) THEN
        FIRST = .FALSE.
        ONE = 1
        !
        !        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
        !        IEEE1, T and RND.
        !
        !        Throughout this routine  we use the function  DLAMC3  to ensure
        !        that relevant values are  stored and not held in registers,  or
        !        are not affected by optimizers.
        !
        !        Compute  a = 2.0**m  with the  smallest positive integer m such
        !        that
        !
        !           fl( a + 1.0 ) = a.
        !
        A = 1
        C = 1
        !
        !+       WHILE( C.EQ.ONE )LOOP
50      CONTINUE
        IF (C .EQ. ONE) THEN
            A = 2*A
            C = DLAMC3(A, ONE)
            C = DLAMC3(C, -A)
            GOTO 50
        END IF
        !+       END WHILE
        !
        !        Now compute  b = 2.0**m  with the smallest positive integer m
        !        such that
        !
        !           fl( a + b ) .gt. a.
        !
        B = 1
        C = DLAMC3(A, B)
        !
        !+       WHILE( C.EQ.A )LOOP
100     CONTINUE
        IF (C .EQ. A) THEN
            B = 2*B
            C = DLAMC3(A, B)
            GOTO 100
        END IF
        !+       END WHILE
        !
        !        Now compute the base.  a and c  are neighbouring floating point
        !        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
        !        their difference is beta. Adding 0.25 to c is to ensure that it
        !        is truncated to beta and not ( beta - 1 ).
        !
        QTR = ONE/4
        SAVEC = C
        C = DLAMC3(C, -A)
        LBETA = C + QTR
        !
        !        Now determine whether rounding or chopping occurs,  by adding a
        !        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
        !
        B = LBETA
        F = DLAMC3(B/2, -B/100)
        C = DLAMC3(F, A)
        IF (C .EQ. A) THEN
            LRND = .TRUE.
        ELSE
            LRND = .FALSE.
        END IF
        F = DLAMC3(B/2, B/100)
        C = DLAMC3(F, A)
        IF ((LRND) .AND. (C .EQ. A)) LRND = .FALSE.
        !
        !        Try and decide whether rounding is done in the  IEEE  'round to
        !        nearest' style. B/2 is half a unit in the last place of the two
        !        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
        !        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
        !        A, but adding B/2 to SAVEC should change SAVEC.
        !
        T1 = DLAMC3(B/2, A)
        T2 = DLAMC3(B/2, SAVEC)
        LIEEE1 = (T1 .EQ. A) .AND. (T2 .GT. SAVEC) .AND. LRND
        !
        !        Now find  the  mantissa, t.  It should  be the  integer part of
        !        log to the base beta of a,  however it is safer to determine  t
        !        by powering.  So we find t as the smallest positive integer for
        !        which
        !
        !           fl( beta**t + 1.0 ) = 1.0.
        !
        LT = 0
        A = 1
        C = 1
        !
        !+       WHILE( C.EQ.ONE )LOOP
150     CONTINUE
        IF (C .EQ. ONE) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3(A, ONE)
            C = DLAMC3(C, -A)
            GOTO 150
        END IF
        !+       END WHILE
        !
    END IF
    !
    BETA = LBETA
    T = LT
    RND = LRND
    IEEE1 = LIEEE1
    !
    !     End of DLAMC1
    !
END
!*==dlamc2.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
!***********************************************************************
!
SUBROUTINE DLAMC2(BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER BETA, EMAX, EMIN, T
    DOUBLE PRECISION EPS, RMAX, RMIN
    LOGICAL RND
    !
    ! Local variables
    !
    DOUBLE PRECISION A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, SIXTH, &
      &                 SMALL, THIRD, TWO, ZERO
    DOUBLE PRECISION DLAMC3
    LOGICAL FIRST, IEEE, IWARN, LIEEE1, LRND
    INTEGER GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT, NGNMIN, NGPMIN
    SAVE LBETA, LEMAX, LEMIN, LEPS, LRMAX, LRMIN, LT
    EXTERNAL DLAMC1, DLAMC3, DLAMC4, DLAMC5

    COMMON/IPRINT/IPRINT
    INTEGER IPRINT
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    !  -- LAPACK auxiliary routine (version 2.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     October 31, 1992
    !
    !     .. Scalar Arguments ..
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DLAMC2 determines the machine parameters specified in its argument
    !  list.
    !
    !  Arguments
    !  =========
    !
    !  BETA    (output) INTEGER
    !          The base of the machine.
    !
    !  T       (output) INTEGER
    !          The number of ( BETA ) digits in the mantissa.
    !
    !  RND     (output) LOGICAL
    !          Specifies whether proper rounding  ( RND = .TRUE. )  or
    !          chopping  ( RND = .FALSE. )  occurs in addition. This may not
    !          be a reliable guide to the way in which the machine performs
    !          its arithmetic.
    !
    !  EPS     (output) DOUBLE PRECISION
    !          The smallest positive number such that
    !
    !             fl( 1.0 - EPS ) .LT. 1.0,
    !
    !          where fl denotes the computed value.
    !
    !  EMIN    (output) INTEGER
    !          The minimum exponent before (gradual) underflow occurs.
    !
    !  RMIN    (output) DOUBLE PRECISION
    !          The smallest normalized number for the machine, given by
    !          BASE**( EMIN - 1 ), where  BASE  is the floating point value
    !          of BETA.
    !
    !  EMAX    (output) INTEGER
    !          The maximum exponent before overflow occurs.
    !
    !  RMAX    (output) DOUBLE PRECISION
    !          The largest positive number for the machine, given by
    !          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
    !          value of BETA.
    !
    !  Further Details
    !  ===============
    !
    !  The computation of  EPS  is based on a routine PARANOIA by
    !  W. Kahan of the University of California at Berkeley.
    !
    ! =====================================================================
    !
    !     .. Local Scalars ..
    !     ..
    !     .. External Functions ..
    !     ..
    !     .. External Subroutines ..
    !     ..
    !     .. Intrinsic Functions ..
    !     ..
    !     .. Save statement ..
    !     ..
    !     .. Data statements ..
    DATA FIRST/.TRUE./, IWARN/.FALSE./
    !     ..
    !     .. Executable Statements ..
    !
    IF (FIRST) THEN
        FIRST = .FALSE.
        ZERO = 0
        ONE = 1
        TWO = 2
        !
        !        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
        !        BETA, T, RND, EPS, EMIN and RMIN.
        !
        !        Throughout this routine  we use the function  DLAMC3  to ensure
        !        that relevant values are stored  and not held in registers,  or
        !        are not affected by optimizers.
        !
        !        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
        !
        CALL DLAMC1(LBETA, LT, LRND, LIEEE1)
        !
        !        Start to find EPS.
        !
        B = LBETA
        A = B**(-LT)
        LEPS = A
        !
        !        Try some tricks to see whether or not this is the correct  EPS.
        !
        B = TWO/3
        HALF = ONE/2
        SIXTH = DLAMC3(B, -HALF)
        THIRD = DLAMC3(SIXTH, SIXTH)
        B = DLAMC3(THIRD, -HALF)
        B = DLAMC3(B, SIXTH)
        B = ABS(B)
        IF (B .LT. LEPS) B = LEPS
        !
        LEPS = 1
        !
        !+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
50      CONTINUE
        IF ((LEPS .GT. B) .AND. (B .GT. ZERO)) THEN
            LEPS = B
            C = DLAMC3(HALF*LEPS, (TWO**5)*(LEPS**2))
            C = DLAMC3(HALF, -C)
            B = DLAMC3(HALF, C)
            C = DLAMC3(HALF, -B)
            B = DLAMC3(HALF, C)
            GOTO 50
        END IF
        !+       END WHILE
        !
        IF (A .LT. LEPS) LEPS = A
        !
        !        Computation of EPS complete.
        !
        !        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
        !        Keep dividing  A by BETA until (gradual) underflow occurs. This
        !        is detected when we cannot recover the previous A.
        !
        RBASE = ONE/LBETA
        SMALL = ONE
        DO I = 1, 3
            SMALL = DLAMC3(SMALL*RBASE, ZERO)
        END DO
        A = DLAMC3(ONE, SMALL)
        CALL DLAMC4(NGPMIN, ONE, LBETA)
        CALL DLAMC4(NGNMIN, -ONE, LBETA)
        CALL DLAMC4(GPMIN, A, LBETA)
        CALL DLAMC4(GNMIN, -A, LBETA)
        IEEE = .FALSE.
        !
        IF ((NGPMIN .EQ. NGNMIN) .AND. (GPMIN .EQ. GNMIN)) THEN
            IF (NGPMIN .EQ. GPMIN) THEN
                LEMIN = NGPMIN
                !            ( Non twos-complement machines, no gradual underflow;
                !              e.g.,  VAX )
            ELSE IF ((GPMIN - NGPMIN) .EQ. 3) THEN
                LEMIN = NGPMIN - 1 + LT
                IEEE = .TRUE.
                !            ( Non twos-complement machines, with gradual underflow;
                !              e.g., IEEE standard followers )
            ELSE
                LEMIN = MIN(NGPMIN, GPMIN)
                !            ( A guess; no known machine )
                IWARN = .TRUE.
            END IF
            !
        ELSE IF ((NGPMIN .EQ. GPMIN) .AND. (NGNMIN .EQ. GNMIN)) THEN
            IF (ABS(NGPMIN - NGNMIN) .EQ. 1) THEN
                LEMIN = MAX(NGPMIN, NGNMIN)
                !            ( Twos-complement machines, no gradual underflow;
                !              e.g., CYBER 205 )
            ELSE
                LEMIN = MIN(NGPMIN, NGNMIN)
                !            ( A guess; no known machine )
                IWARN = .TRUE.
            END IF
            !
        ELSE IF ((ABS(NGPMIN - NGNMIN) .EQ. 1) .AND. (GPMIN .EQ. GNMIN)) &
            &             THEN
            IF ((GPMIN - MIN(NGPMIN, NGNMIN)) .EQ. 3) THEN
                LEMIN = MAX(NGPMIN, NGNMIN) - 1 + LT
                !            ( Twos-complement machines with gradual underflow;
                !              no known machine )
            ELSE
                LEMIN = MIN(NGPMIN, NGNMIN)
                !            ( A guess; no known machine )
                IWARN = .TRUE.
            END IF
            !
        ELSE
            LEMIN = MIN(NGPMIN, NGNMIN, GPMIN, GNMIN)
            !         ( A guess; no known machine )
            IWARN = .TRUE.
        END IF
        !**
        ! Comment out this if block if EMIN is ok
        IF (IWARN) THEN
            FIRST = .TRUE.
            IF (IPRINT > 0) WRITE (6, FMT=99001) LEMIN
        END IF
        !**
        !
        !        Assume IEEE arithmetic if we found denormalised  numbers above,
        !        or if arithmetic seems to round in the  IEEE style,  determined
        !        in routine DLAMC1. A true IEEE machine should have both  things
        !        true; however, faulty machines may have one or the other.
        !
        IEEE = IEEE .OR. LIEEE1
        !
        !        Compute  RMIN by successive division by  BETA. We could compute
        !        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
        !        this computation.
        !
        LRMIN = 1
        DO I = 1, 1 - LEMIN
            LRMIN = DLAMC3(LRMIN*RBASE, ZERO)
        END DO
        !
        !        Finally, call DLAMC5 to compute EMAX and RMAX.
        !
        CALL DLAMC5(LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX)
    END IF
    !
    BETA = LBETA
    T = LT
    RND = LRND
    EPS = LEPS
    EMIN = LEMIN
    RMIN = LRMIN
    EMAX = LEMAX
    RMAX = LRMAX
    !
    RETURN
    !
99001 FORMAT(//' WARNING. The value EMIN may be incorrect:-', &
&        '  EMIN = ', I8, &
&        /' If, after inspection, the value EMIN looks', &
&        ' acceptable please comment out ', &
&        /' the IF block as marked within the code of routine', &
&        ' DLAMC2,', /' otherwise supply EMIN explicitly.',/)
    !
    !     End of DLAMC2
    !
END
!*==dlamc3.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
!***********************************************************************
!
DOUBLE PRECISION FUNCTION DLAMC3(A, B)
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
    DOUBLE PRECISION A, B
!
!*** End of declarations rewritten by SPAG
!
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A, B    (input) DOUBLE PRECISION
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
    DLAMC3 = A + B
!
!
!     End of DLAMC3
!
END
!*==dlamc4.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
!***********************************************************************
!
SUBROUTINE DLAMC4(EMIN, START, BASE)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER BASE, EMIN
    DOUBLE PRECISION START
    !
    ! Local variables
    !
    DOUBLE PRECISION A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
    DOUBLE PRECISION DLAMC3
    INTEGER I
    EXTERNAL DLAMC3
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    !  -- LAPACK auxiliary routine (version 2.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     October 31, 1992
    !
    !     .. Scalar Arguments ..
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DLAMC4 is a service routine for DLAMC2.
    !
    !  Arguments
    !  =========
    !
    !  EMIN    (output) EMIN
    !          The minimum exponent before (gradual) underflow, computed by
    !          setting A = START and dividing by BASE until the previous A
    !          can not be recovered.
    !
    !  START   (input) DOUBLE PRECISION
    !          The starting point for determining EMIN.
    !
    !  BASE    (input) INTEGER
    !          The base of the machine.
    !
    ! =====================================================================
    !
    !     .. Local Scalars ..
    !     ..
    !     .. External Functions ..
    !     ..
    !     .. Executable Statements ..
    !
    A = START
    ONE = 1
    RBASE = ONE/BASE
    ZERO = 0
    EMIN = 1
    B1 = DLAMC3(A*RBASE, ZERO)
    C1 = A
    C2 = A
    D1 = A
    D2 = A
    !+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
    !    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
100 CONTINUE
    IF ((C1 .EQ. A) .AND. (C2 .EQ. A) .AND. (D1 .EQ. A) .AND. (D2 .EQ. A)) &
      &     THEN
        EMIN = EMIN - 1
        A = B1
        B1 = DLAMC3(A/BASE, ZERO)
        C1 = DLAMC3(B1*BASE, ZERO)
        D1 = ZERO
        DO I = 1, BASE
            D1 = D1 + B1
        END DO
        B2 = DLAMC3(A*RBASE, ZERO)
        C2 = DLAMC3(B2/RBASE, ZERO)
        D2 = ZERO
        DO I = 1, BASE
            D2 = D2 + B2
        END DO
        GOTO 100
    END IF
!+    END WHILE
!
!
!     End of DLAMC4
!
END
!*==dlamc5.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
!***********************************************************************
!
SUBROUTINE DLAMC5(BETA, P, EMIN, IEEE, EMAX, RMAX)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! PARAMETER definitions
    !
    DOUBLE PRECISION ZERO, ONE
    PARAMETER(ZERO=0.0D0, ONE=1.0D0)
    !
    ! Dummy arguments
    !
    INTEGER BETA, EMAX, EMIN, P
    LOGICAL IEEE
    DOUBLE PRECISION RMAX
    !
    ! Local variables
    !
    DOUBLE PRECISION DLAMC3
    INTEGER EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
    DOUBLE PRECISION OLDY, RECBAS, Y, Z
    EXTERNAL DLAMC3
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    !  -- LAPACK auxiliary routine (version 2.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     October 31, 1992
    !
    !     .. Scalar Arguments ..
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DLAMC5 attempts to compute RMAX, the largest machine floating-point
    !  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
    !  approximately to a power of 2.  It will fail on machines where this
    !  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
    !  EMAX = 28718).  It will also fail if the value supplied for EMIN is
    !  too large (i.e. too close to zero), probably with overflow.
    !
    !  Arguments
    !  =========
    !
    !  BETA    (input) INTEGER
    !          The base of floating-point arithmetic.
    !
    !  P       (input) INTEGER
    !          The number of base BETA digits in the mantissa of a
    !          floating-point value.
    !
    !  EMIN    (input) INTEGER
    !          The minimum exponent before (gradual) underflow.
    !
    !  IEEE    (input) LOGICAL
    !          A logical flag specifying whether or not the arithmetic
    !          system is thought to comply with the IEEE standard.
    !
    !  EMAX    (output) INTEGER
    !          The largest exponent before overflow
    !
    !  RMAX    (output) DOUBLE PRECISION
    !          The largest machine floating-point number.
    !
    ! =====================================================================
    !
    !     .. Parameters ..
    !     ..
    !     .. Local Scalars ..
    !     ..
    !     .. External Functions ..
    !     ..
    !     .. Intrinsic Functions ..
    !     ..
    !     .. Executable Statements ..
    !
    !     First compute LEXP and UEXP, two powers of 2 that bound
    !     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
    !     approximately to the bound that is closest to abs(EMIN).
    !     (EMAX is the exponent of the required number RMAX).
    !
    LEXP = 1
    EXBITS = 1
100 CONTINUE
    TRY = LEXP*2
    IF (TRY .LE. (-EMIN)) THEN
        LEXP = TRY
        EXBITS = EXBITS + 1
        GOTO 100
    END IF
    IF (LEXP .EQ. -EMIN) THEN
        UEXP = LEXP
    ELSE
        UEXP = TRY
        EXBITS = EXBITS + 1
    END IF
    !
    !     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
    !     than or equal to EMIN. EXBITS is the number of bits needed to
    !     store the exponent.
    !
    IF ((UEXP + EMIN) .GT. (-LEXP - EMIN)) THEN
        EXPSUM = 2*LEXP
    ELSE
        EXPSUM = 2*UEXP
    END IF
    !
    !     EXPSUM is the exponent range, approximately equal to
    !     EMAX - EMIN + 1 .
    !
    EMAX = EXPSUM + EMIN - 1
    NBITS = 1 + EXBITS + P
    !
    !     NBITS is the total number of bits needed to store a
    !     floating-point number.
    !
    !
    !        Either there are an odd number of bits used to store a
    !        floating-point number, which is unlikely, or some bits are
    !        not used in the representation of numbers, which is possible,
    !        (e.g. Cray machines) or the mantissa has an implicit bit,
    !        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
    !        most likely. We have to assume the last alternative.
    !        If this is true, then we need to reduce EMAX by one because
    !        there must be some way of representing zero in an implicit-bit
    !        system. On machines like Cray, we are reducing EMAX by one
    !        unnecessarily.
    !
    IF ((MOD(NBITS, 2) .EQ. 1) .AND. (BETA .EQ. 2)) EMAX = EMAX - 1
    !
    !
    !        Assume we are on an IEEE machine which reserves one exponent
    !        for infinity and NaN.
    !
    IF (IEEE) EMAX = EMAX - 1
    !
    !     Now create RMAX, the largest machine number, which should
    !     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
    !
    !     First compute 1.0 - BETA**(-P), being careful that the
    !     result is less than 1.0 .
    !
    RECBAS = ONE/BETA
    Z = BETA - ONE
    Y = ZERO
    DO I = 1, P
        Z = Z*RECBAS
        IF (Y .LT. ONE) OLDY = Y
        Y = DLAMC3(Y, Z)
    END DO
    IF (Y .GE. ONE) Y = OLDY
    !
    !     Now multiply by BETA**EMAX to get RMAX.
    !
    DO I = 1, EMAX
        Y = DLAMC3(Y*BETA, ZERO)
    END DO
    !
    RMAX = Y
    !
    !     End of DLAMC5
    !
END
!*==lsame.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
LOGICAL FUNCTION LSAME(CA, CB)
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
    CHARACTER CA, CB
!
! Local variables
!
    INTEGER INTA, INTB, ZCODE
!
!*** End of declarations rewritten by SPAG
!
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     January 31, 1994
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
    LSAME = CA .EQ. CB
    IF (LSAME) RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
    ZCODE = ICHAR('Z')
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
    INTA = ICHAR(CA)
    INTB = ICHAR(CB)
!
    IF (ZCODE .EQ. 90 .OR. ZCODE .EQ. 122) THEN
        !
        !        ASCII is assumed - ZCODE is the ASCII code of either lower or
        !        upper case 'Z'.
        !
        IF (INTA .GE. 97 .AND. INTA .LE. 122) INTA = INTA - 32
        IF (INTB .GE. 97 .AND. INTB .LE. 122) INTB = INTB - 32
        !
    ELSE IF (ZCODE .EQ. 233 .OR. ZCODE .EQ. 169) THEN
        !
        !        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
        !        upper case 'Z'.
        !
        IF (INTA .GE. 129 .AND. INTA .LE. 137 .OR. INTA .GE. 145 .AND. &
          &        INTA .LE. 153 .OR. INTA .GE. 162 .AND. INTA .LE. 169) &
          &        INTA = INTA + 64
        IF (INTB .GE. 129 .AND. INTB .LE. 137 .OR. INTB .GE. 145 .AND. &
          &        INTB .LE. 153 .OR. INTB .GE. 162 .AND. INTB .LE. 169) &
          &        INTB = INTB + 64
        !
    ELSE IF (ZCODE .EQ. 218 .OR. ZCODE .EQ. 250) THEN
        !
        !        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
        !        plus 128 of either lower or upper case 'Z'.
        !
        IF (INTA .GE. 225 .AND. INTA .LE. 250) INTA = INTA - 32
        IF (INTB .GE. 225 .AND. INTB .LE. 250) INTB = INTB - 32
    END IF
    LSAME = INTA .EQ. INTB
!
!     RETURN
!
!     End of LSAME
!
END
!*==idamax.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
INTEGER FUNCTION IDAMAX(N, DX, INCX)
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
    INTEGER INCX, N
    DOUBLE PRECISION DX(*)
!
! Local variables
!
    DOUBLE PRECISION DMAX
    INTEGER I, IX
!
!*** End of declarations rewritten by SPAG
!
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!
    IDAMAX = 0
    IF (N .LT. 1 .OR. INCX .LE. 0) RETURN
    IDAMAX = 1
    IF (N .EQ. 1) RETURN
    IF (INCX .EQ. 1) THEN
        !
        !        code for increment equal to 1
        !
        DMAX = DABS(DX(1))
        DO I = 2, N
            IF (DABS(DX(I)) .GT. DMAX) THEN
                IDAMAX = I
                DMAX = DABS(DX(I))
            END IF
        END DO
        GOTO 99999
    END IF
!
!        code for increment not equal to 1
!
    IX = 1
    DMAX = DABS(DX(1))
    IX = IX + INCX
    DO I = 2, N
        IF (DABS(DX(IX)) .GT. DMAX) THEN
            IDAMAX = I
            DMAX = DABS(DX(IX))
        END IF
        IX = IX + INCX
    END DO
    RETURN
99999 CONTINUE
END
!*==dcopy.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE DCOPY(N, DX, INCX, DY, INCY)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER INCX, INCY, N
    DOUBLE PRECISION DX(*), DY(*)
    !
    ! Local variables
    !
    INTEGER I, IX, IY, M, MP1
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    !     copies a vector, x, to a vector, y.
    !     uses unrolled loops for increments equal to one.
    !     jack dongarra, linpack, 3/11/78.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    !
    IF (N .LE. 0) RETURN
    IF (INCX .EQ. 1 .AND. INCY .EQ. 1) THEN
        !
        !        code for both increments equal to 1
        !
        !
        !        clean-up loop
        !
        M = MOD(N, 7)
        IF (M .NE. 0) THEN
            DO I = 1, M
                DY(I) = DX(I)
            END DO
            IF (N .LT. 7) RETURN
        END IF
        MP1 = M + 1
        DO I = MP1, N, 7
            DY(I) = DX(I)
            DY(I + 1) = DX(I + 1)
            DY(I + 2) = DX(I + 2)
            DY(I + 3) = DX(I + 3)
            DY(I + 4) = DX(I + 4)
            DY(I + 5) = DX(I + 5)
            DY(I + 6) = DX(I + 6)
        END DO
    ELSE
        !
        !        code for unequal increments or equal increments
        !          not equal to 1
        !
        IX = 1
        IY = 1
        IF (INCX .LT. 0) IX = (-N + 1)*INCX + 1
        IF (INCY .LT. 0) IY = (-N + 1)*INCY + 1
        DO I = 1, N
            DY(IY) = DX(IX)
            IX = IX + INCX
            IY = IY + INCY
        END DO
        RETURN
    END IF
END
