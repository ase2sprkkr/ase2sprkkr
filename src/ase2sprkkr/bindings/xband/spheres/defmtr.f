C*==defmtr.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE DEFMTR(W,NRAD1,AV,BV,CV,IS,BAS,RO,R0,DPAS,ALAT,NATOM,
     &                  NSORT,Z,SMT,IMT,NRADMAX,WSREST)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      COMMON /IPRINT/ IPRINT
      INTEGER IPRINT
      DOUBLE PRECISION ALAT,DPAS
      INTEGER IMT,NATOM,NRAD1,NRADMAX,NSORT
      DOUBLE PRECISION AV(3),BAS(3,*),BV(3),CV(3),R0(*),RO(NRADMAX,*),
     &                 SMT(*),WSREST(*),Z(*)
      INTEGER IS(*),W(*)
C
C Local variables
C
      DOUBLE PRECISION A(:),AVW,BOUND(3,3),DV(:),ORIGIN(3),PI,PLAT(3,3),
     &                 POT(:,:),R(:),VOLCEL
      DOUBLE PRECISION DPI
      INTEGER I,J,LTMAX(3),NDEL(3),NR(:),NRMAX,OA,OBAS,ODV,ONR,OPOT,OR
      SAVE A,AVW,BOUND,DV,I,J,LTMAX,NDEL,NRMAX,OA,OBAS,ODV,ONR,OPOT,OR,
     &     ORIGIN,PI,PLAT,POT,R,VOLCEL
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE  POT,DV,R,A,NR
C
      PI = DPI()
      NRMAX = NRAD1 + 1
      CALL DEFDR(OA,NSORT)
      CALL DEFDR(OPOT,NSORT*NRMAX)
      CALL DEFI(ONR,NSORT)
      CALL DEFDR(OR,NRMAX)
      CALL DEFDR(ODV,NSORT*NRMAX)
C      call hpot(W(oa),r0,nsort,W(onr),nrad1,W(opot),Z,RO,dpas,
C     .     W(or),W(odv),nradmax)
C
      ALLOCATE (POT(NRMAX,NSORT))
      ALLOCATE (DV(NRAD1),R(NRAD1),A(NSORT),NR(NSORT))
C
      CALL HPOT(A,R0,NSORT,NR,NRAD1,POT,Z,RO,DPAS,R,DV,NRADMAX)
C
      CALL RLSE(OR)
C
      VOLCEL = AV(1)*(BV(2)*CV(3)-BV(3)*CV(2)) + AV(2)
     &         *(BV(3)*CV(1)-BV(1)*CV(3)) + AV(3)
     &         *(BV(1)*CV(2)-BV(2)*CV(1))
      VOLCEL = ABS(VOLCEL)
C    average radii:
      AVW = (VOLCEL/NATOM/(4.D0*PI)*3.D0)**(1.D0/3.D0)*ALAT
      PLAT(1,1) = AV(1)
      PLAT(1,2) = BV(1)
      PLAT(1,3) = CV(1)
      PLAT(2,1) = AV(2)
      PLAT(2,2) = BV(2)
      PLAT(2,3) = CV(2)
      PLAT(3,1) = AV(3)
      PLAT(3,2) = BV(3)
      PLAT(3,3) = CV(3)
C
      DO I = 1,3
         DO J = 1,3
            BOUND(I,J) = PLAT(I,J)
         END DO
         LTMAX(I) = 2
         NDEL(I) = 0
         ORIGIN(I) = 0
      END DO
      CALL DEFDR(OBAS,NATOM*3)
C      call defdr(owsrest,natom)
C     call hrtree(w,alat,avw,bas,bound,
C    .     is,ltmax,natom,nsort,ndel,
C    .     origin,plat,smt,z,wsrest,
C    .     W(oA),R0,W(onR),W(oPOT),NRMAX,imt)
C
      CALL HRTREE(W,ALAT,AVW,BAS,BOUND,IS,LTMAX,NATOM,NSORT,NDEL,ORIGIN,
     &            PLAT,SMT,Z,WSREST,A,R0,NR,POT,NRMAX,IMT)
C

      DEALLOCATE (POT,DV,R,A,NR)

      END
C*==hpot.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      SUBROUTINE HPOT(A,B,NCLASS,NR,NRAD1,POT,Z,RO,DPAS,R,DV,NRADMAX)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      DOUBLE PRECISION HUGE
      PARAMETER (HUGE=1.D6)
C
C Dummy arguments
C
      DOUBLE PRECISION DPAS
      INTEGER NCLASS,NRAD1,NRADMAX
      DOUBLE PRECISION A(*),B(*),DV(NRAD1),POT(NRAD1+1,*),R(NRAD1),
     &                 RO(NRADMAX,*),Z(*)
      INTEGER NR(*)
C
C Local variables
C
      INTEGER I,ICL,NP,NRMAX
      DOUBLE PRECISION ZZ
      SAVE I,ICL,NP,NRMAX,ZZ
C
C*** End of declarations rewritten by SPAG
C
C- Calculate H. potential
C ----------------------------------------------------------------------
Ci Inputs:
Ci   clabl :name of the different inequivalent atom
Ci   nclass:number of classes, atoms in same class are symmetry-related
Ci   nrmax :maximum number of mesh points
Ci   nsp   :=1 spin degenerate, =2 non-degenerate
Ci   z     :nuclear charge
Co Outputs:
Co   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Co   b     :                 -//-
Co   ierr  :0 if no error occurs, ic if read error for class ic
Co   nr    :number of mesh points
Co   pot   :spherical Hartree potential
Cw   v     :spherical potential (electronic contribution)
Cr Remarks:
C ----------------------------------------------------------------------
C
      NRMAX = NRAD1 + 1                ! - radial mesh
      DO ICL = 1,NCLASS
         DO I = 1,NRAD1
            R(I) = B(ICL)*EXP(DBLE(I-1)*DPAS)
         END DO
         NP = NRAD1
         CALL POTS(POT(2,ICL),RO(1,ICL),DV,R,DPAS,Z(ICL),NP)
         ZZ = 2*Z(ICL)
         A(ICL) = DPAS
         NR(ICL) = NRMAX
         DO I = 2,NRMAX
            POT(I,ICL) = POT(I,ICL) - ZZ/R(I-1)
         END DO
         POT(1,ICL) = -HUGE
      END DO
      END
C*==pots.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE POTS(DV,D,DP,DR,DPAS,Z,NP)
C
C iHTEgPiPOBAHiE pOTEHciAlA pO 4 TO~KAM
C DV - POTENTIAL   D CH.DENCITY  DP WORK ARRAY  DR RADIAL MESH
C DPAS EXP. PASS
C Z ATOM NUMBER     NP - NUMBER OF POINTS
C **********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION DPAS,Z
      INTEGER NP
      DOUBLE PRECISION D(*),DP(*),DR(*),DV(*)
C
C Local variables
C
      DOUBLE PRECISION DAS,DLO,DLO2
      INTEGER I,J,K
      SAVE DAS,DLO,DLO2,I,J,K
C
C*** End of declarations rewritten by SPAG
C
C$$$      double precision z
C
      DAS = DPAS/24.D0
      DO I = 1,NP
         DV(I) = D(I)*DR(I)
      END DO
      DLO = EXP(DPAS)
      DLO2 = DLO*DLO
      DP(2) = DR(1)*(D(2)-D(1)*DLO2)/(12.D0*(DLO-1.D0))
      DP(1) = DV(1)/3.D0 - DP(2)/DLO2
      DP(2) = DV(2)/3.D0 - DP(2)*DLO2
      J = NP - 1
      DO I = 3,J
         DP(I) = DP(I-1) + DAS*(13.D0*(DV(I)+DV(I-1))-(DV(I-2)+DV(I+1)))
      END DO
      DP(NP) = DP(J)
      DV(J) = DP(J)
      DV(NP) = DP(J)
      DO I = 3,J
         K = NP + 1 - I
         DV(K) = DV(K+1)
     &           /DLO + DAS*(13.D0*(DP(K+1)/DLO+DP(K))-(DP(K+2)/DLO2+
     &           DP(K-1)*DLO))
      END DO
      DV(1) = DV(3)/DLO2 + DPAS*(DP(1)+4.D0*DP(2)/DLO+DP(3)/DLO2)/3.D0
      DO I = 1,NP
         DV(I) = DV(I)/DR(I)*2
      END DO
      END
C*==fmesh.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      DOUBLE PRECISION FUNCTION FMESH(A,B,F,NR,R)
C- Computes the value of f at r for given values on a mesh
C ----------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   f     :function defined on a mesh
Ci   nr    :number of mesh points
Ci   r     :radial distance
Co Outputs:
Co   fmesh :interpolated value of f at r
Cr Remarks:
Cr   if r is insite the mesh a quadratic fit is used
Cr   if r is outside the mesh an exponential fit is used
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NP
      PARAMETER (NP=4)
C
C Dummy arguments
C
      DOUBLE PRECISION A,B,R
      INTEGER NR
      DOUBLE PRECISION F(*)
C
C Local variables
C
      DOUBLE PRECISION DELSQF,DI3INT
      INTEGER IS,NSTART
      DOUBLE PRECISION XX
      SAVE IS,NSTART,XX
      EXTERNAL DELSQF,DI3INT
C
C*** End of declarations rewritten by SPAG
C
      IF ( R.GE.B ) THEN
         XX = LOG(R/B)/A + 2
         IS = IDNINT(XX) - 1
      ELSE
         FMESH = F(2)
         RETURN
      END IF
      IF ( IS.LE.NR-2 ) THEN
         IS = MAX0(1,IS)
         FMESH = DI3INT(IS,F(IS),XX)
      ELSE
         NSTART = NR - NP + 1
         FMESH = DELSQF(NP,NSTART,F(NSTART),XX)
      END IF
      END
C*==iprint.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C- Gets last integer of print priority stack
C ----------------------------------------------------------------------
Co Outputs:
Co   iprint:defines the verbosity level
Cr Remarks:
Cr   verbosity:   0  nearly nothing is printed
Cr               10  very terse
Cr               20  terse
Cr               30  normal
Cr               40  verbose
Cr               50  very verbose
Cr               60  highest verbosity
Cr              100  low-level debugging
Cr              110  intermediate-level debugging
Cr              120  high-level debugging
C ----------------------------------------------------------------------
C
C*** Start of declarations rewritten by SPAG
C
C*** End of declarations rewritten by SPAG
C
C*==potxn.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      SUBROUTINE POTXN(A,ALAT,B,BAS,IAX,ICLASS,NPR,NR,NRMAX,PLAT,POT,
     &                 POTL,XN)
C-Calculates Hartree potential at point xn
C ----------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   alat  :length scale
Ci   b     :                 -//-
Ci   bas   :basis vectors (scaled by alat)
Ci   iax   :information about positions around a specified pair of atom
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   npr   :number of neighbors
Ci   nr    :number of mesh points
Ci   nrmax :maximum number of mesh points
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   plat  :primitive translation vectors in real space
Ci   pot   :spherical Hartree potential
Ci   xn    :cartesian coordinates where potential is calculated
Co Outputs:
Ci   potl  :ovelapping Hartree potential
Cr Remarks:
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION ALAT,POTL
      INTEGER NPR,NRMAX
      DOUBLE PRECISION A(*),B(*),BAS(3,*),PLAT(3,3),POT(NRMAX,*),XN(3)
      INTEGER IAX(5,*),ICLASS(*),NR(*)
C
C Local variables
C
      DOUBLE PRECISION D(3),R
      DOUBLE PRECISION DNRM23,FMESH
      INTEGER JBAS,JC,JPR,K
      SAVE D,JBAS,JC,JPR,K,R
      EXTERNAL DNRM23,FMESH
C
C*** End of declarations rewritten by SPAG
C
      POTL = 0.D0
      DO JPR = 1,NPR
         JBAS = IAX(2,JPR)
         JC = ICLASS(JBAS)
         DO K = 1,3
            D(K) = BAS(K,JBAS) - XN(K) + PLAT(K,1)*IAX(3,JPR)
     &             + PLAT(K,2)*IAX(4,JPR) + PLAT(K,3)*IAX(5,JPR)
         END DO
         R = SQRT(DNRM23(D))*ALAT
         POTL = POTL + FMESH(A(JC),B(JC),POT(1,JC),NR(JC),R)
      END DO
C
      END
C*==di3int.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      DOUBLE PRECISION FUNCTION DI3INT(IX,YA,X)
C- Interpolates y = f(x) for given xa and ya=f(xa)
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ix    :xa(1) is integer and equals ix
Ci          and xa(1),xa(2) and xa(3) differ exactly by 1.d0
Ci   ya    :value of f at xa
Ci   x     :x-value at which f is interpolated
Co Outputs:
Co   di3int:interpolated value
Cr Remarks:
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IX
      DOUBLE PRECISION X
      DOUBLE PRECISION YA(3)
C
C Local variables
C
      DOUBLE PRECISION XA(3)
C
C*** End of declarations rewritten by SPAG
C
      XA(1) = DBLE(IX)
      XA(2) = DBLE(IX+1)
      XA(3) = DBLE(IX+2)
C
      DI3INT = 0.5D0*(X-XA(2))*((X-XA(3))*YA(1)+(X-XA(1))*YA(3))
     &         - (X-XA(1))*(X-XA(3))*YA(2)
      END
C*==i_shell.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      SUBROUTINE I_SHELL(M,N,IARRAY)
C- shell sort of a array of integer vectors
C ----------------------------------------------------------------------
Ci Inputs:
Ci   m     :number of components in iarray
Ci   n     :number of elements in iarray
Ci   iarray:array to be sorted
Co Outputs:
Co   iarray:array to be sorted
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER M,N
      INTEGER IARRAY(M,0:N-1)
C
C Local variables
C
      INTEGER I,IT,J,K,L,LOGNB2,MM,MMM,N2,NN
C
C*** End of declarations rewritten by SPAG
C
      LOGNB2 = INT(LOG(FLOAT(N+1))*1.4426950)
      N2 = N
      DO NN = 1,LOGNB2
         N2 = N2/2
         K = N - N2
         DO J = 1,K
            I = J - 1
 20         CONTINUE
            L = I + N2
            DO MM = 1,M
               IF ( IARRAY(MM,L).LT.IARRAY(MM,I) ) THEN
                  DO MMM = 1,M
                     IT = IARRAY(MMM,I)
                     IARRAY(MMM,I) = IARRAY(MMM,L)
                     IARRAY(MMM,L) = IT
                  END DO
                  I = I - N2
                  IF ( I.LT.0 ) EXIT
                  GOTO 20
               ELSE IF ( IARRAY(MM,L).NE.IARRAY(MM,I) ) THEN
                  EXIT
               END IF
            END DO
         END DO
      END DO
C
      END
C*==hrtree.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      SUBROUTINE HRTREE(W,ALAT,AVW,BAS,BOUND,ICLASS,LTMAX,NBAS,NCLASS,
     &                  NDEL,ORIGIN,PLAT,WSR,Z,WSREST,A,B,NR,POT,NRMAX,
     &                  IMT)
C- Calculate Hartree potential on a mesh and determine muffin-tin radia
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci   avw   :average Wigner Seitz radius
Ci   bas   :basis vectors (scaled by alat)
Ci   bound :two vectors spanning the plane (scaled by alat)
Ci   clabl :name of the different inequivalent atom
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   ltmax :ltmax(i)= limit in i-dirction for unit cells considered
Ci   nbas  :number of atoms in the basis
Ci   nclass:number of classes, atoms in same class are symmetry-related
Ci   ndel  :ndel(i)=number of mesh points along the bound(i) vector
Ci   nrclas:number of atoms in the i-th class
Ci   origin:origin of the plane (scaled by alat)
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   z     :nuclear charge
Co Outputs:
Cio  wsr   :Wigner-Seitz sphere radius (in atomic units)
Cr Remarks:
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      DOUBLE PRECISION FACN1
      PARAMETER (FACN1=1.5D0)
C
C Dummy arguments
C
      DOUBLE PRECISION ALAT,AVW
      INTEGER IMT,NBAS,NCLASS,NRMAX
      DOUBLE PRECISION A(*),B(*),BAS(3,*),BOUND(3,3),ORIGIN(*),PLAT(3,3)
     &                 ,POT(NRMAX,*),WSR(*),WSREST(*),Z(*)
      INTEGER ICLASS(*),LTMAX(3),NDEL(*),NR(*),W(*)
C
C Local variables
C
      DOUBLE PRECISION DSCL(:,:),WMAX
      INTEGER I,IC,J,NEIGHM,ODSCL,OIAX1,OLOCK,ONPR1
      INTEGER IDAMAX
      COMMON /IPRINT/ IPRINT
      INTEGER IPRINT
      LOGICAL LOCK(:)
      SAVE DSCL,I,IC,J,LOCK,NEIGHM,ODSCL,OIAX1,OLOCK,ONPR1,WMAX
      EXTERNAL DEFDR,DEFI,DEFWSR,NGHBR1,POTMAX,POTSUM,RLSE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE DSCL,LOCK
C
      IF ( IPRINT.GT.90 ) THEN
         PRINT *,'******************************************'
         PRINT *,'alat',ALAT,'  avw',AVW
         PRINT *,' === bas. v. ==='
         DO J = 1,NBAS
            PRINT 99001,(BAS(I,J),I=1,3)
         END DO
         PRINT *,' === bound ==='
         DO J = 1,3
            PRINT 99001,(BOUND(I,J),I=1,3)
         END DO
         PRINT *,' === ltmax ==='
         PRINT 99002,LTMAX
         PRINT *,' === ndel ==='
         PRINT 99002,(NDEL(I),I=1,3)
         PRINT *,' === origin ==='
         PRINT 99001,(ORIGIN(I),I=1,3)
         PRINT *,' === plat ==='
         DO J = 1,3
            PRINT 99001,(PLAT(I,J),I=1,3)
         END DO
         PRINT *,'******************************************'
      END IF
C --- Calculate Hartree potential on a mesh
      IF ( NDEL(1)*NDEL(2).NE.0 ) CALL POTSUM(A,ALAT,B,BAS,BOUND,ICLASS,
     &     LTMAX,NBAS,NDEL,NR,NRMAX,ORIGIN,PLAT,POT)
C --- Determine maximum of Hartree potential
      DO IC = 1,NCLASS
         CALL DEFWSR(WSREST(IC),Z(IC))
         WSR(IC) = WSREST(IC)
      END DO
      WMAX = WSREST(IDAMAX(NCLASS,WSREST,1))/AVW
      I = 2*INT((2.D0*FACN1*WMAX+1.D0)**3)
      NEIGHM = MAX(2*INT((2.D0*FACN1*WMAX+1.D0)**3),50)
      CALL DEFI(OIAX1,5*NEIGHM*NCLASS)
      CALL DEFI(ONPR1,NCLASS)
      CALL NGHBR1(W,ALAT,BAS,FACN1,W(OIAX1),ICLASS,NBAS,NCLASS,NEIGHM,
     &            W(ONPR1),PLAT,WSREST)
C
      CALL POTMAX(W,A,ALAT,B,BAS,W(OIAX1),ICLASS,NBAS,NCLASS,W(ONPR1),
     &            NR,NRMAX,PLAT,POT,WMAX,WSR,WSREST,Z)
      CALL DEFDR(ODSCL,NCLASS*(NCLASS+2))
      CALL DEFDR(OLOCK,NCLASS)
C     call blowup(alat,bas,1.d0,0.d0,w(oiax1),iclass,nbas,
C    $           nclass,w(onpr1),0.d0,0.d0,plat,wsr,W(odscl),W(olock))
C     if(imt.ne.0)call blowup(alat,bas,1.d0,0.d0,w(oiax1),iclass,nbas,
C    $         nclass,w(onpr1),0.4d0,0.8d0,plat,wsr,W(odscl),W(olock))
C
      ALLOCATE (DSCL(NCLASS,0:NCLASS+1))
      ALLOCATE (LOCK(NCLASS))
C
      CALL BLOWUP(ALAT,BAS,1.D0,0.D0,W(OIAX1),ICLASS,NBAS,NCLASS,
     &            W(ONPR1),0.D0,0.D0,PLAT,WSR,DSCL,LOCK)
      IF ( IMT.NE.0 ) CALL BLOWUP(ALAT,BAS,1.D0,0.D0,W(OIAX1),ICLASS,
     &                            NBAS,NCLASS,W(ONPR1),0.4D0,0.8D0,PLAT,
     &                            WSR,DSCL,LOCK)
C
C
      DEALLOCATE (DSCL)
      DEALLOCATE (LOCK)

      CALL RLSE(OIAX1)
99001 FORMAT (20F21.15)
99002 FORMAT (20I10)
C
      END
C*==potsum.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      SUBROUTINE POTSUM(A,ALAT,B,BAS,BOUND,ICLASS,LTMAX,NBAS,NDEL,NR,
     &                  NRMAX,ORIGIN,PLAT,POT)
C- Calculates the overlapping Hartree potential
C ----------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   alat  :length scale
Ci   bas   :basis vectors (scaled by alat)
Ci   bound :two vectors spanning the plane (scaled by alat)
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   ltmax :ltmax(i)= limit in i-dirction for unit cells considered
Ci   nbas  :number of atoms in the basis
Ci   ndel  :ndel(i)=number of mesh points along the bound(i) vector
Ci   nr    :number of mesh points
Ci   nrmax :maximum number of mesh points
Ci   origin:origin of the plane (scaled by alat)
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   pot   :spherical Hartree potential
Co Output written to file POT
Cr Remarks:
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION ALAT
      INTEGER NBAS,NRMAX
      COMMON /IPRINT/ IPRINT
      INTEGER IPRINT
      DOUBLE PRECISION A(*),B(*),BAS(3,*),BOUND(3,3),ORIGIN(3),PLAT(3,3)
     &                 ,POT(NRMAX,*)
      INTEGER ICLASS(*),LTMAX(3),NDEL(*),NR(*)
C
C Local variables
C
      DOUBLE PRECISION D01,D02,D03,D11,D12,D13,D21,D22,D23,D31,D32,D33,
     &                 DN1,DN2,FAC1,FAC2,POTL,RAD,XN1,XN2,XN3
      DOUBLE PRECISION FMESH
      INTEGER I,I1,I2,IBAS,IC,J,J1,J2,J3
      SAVE D01,D02,D03,D11,D12,D13,D21,D22,D23,D31,D32,D33,DN1,DN2,FAC1,
     &     FAC2,I,I1,I2,IBAS,IC,J,J1,J2,J3,POTL,RAD,XN1,XN2,XN3
C
C*** End of declarations rewritten by SPAG
C
      IF ( NDEL(1).NE.1 ) FAC1 = 1.D0/DBLE(NDEL(1)-1)
      IF ( NDEL(2).NE.1 ) FAC2 = 1.D0/DBLE(NDEL(2)-1)
      IF (IPRINT > 0) THEN
        WRITE (6,99001) ORIGIN,((BOUND(I,J),I=1,3),J=1,2)
        WRITE (6,99002)
      END IF

C
      DO I1 = 0,NDEL(1) - 1
         DN1 = I1*FAC1
         DO I2 = 0,NDEL(2) - 1
            DN2 = I2*FAC2
            XN1 = ORIGIN(1) + BOUND(1,1)*DN1 + BOUND(1,2)*DN2
            XN2 = ORIGIN(2) + BOUND(2,1)*DN1 + BOUND(2,2)*DN2
            XN3 = ORIGIN(3) + BOUND(3,1)*DN1 + BOUND(3,2)*DN2
            POTL = 0.D0
            DO IBAS = 1,NBAS
               IC = ICLASS(IBAS)
C --------  ltmax(i) lattice translations in i-direction considered
               D01 = BAS(1,IBAS) - XN1
               D02 = BAS(2,IBAS) - XN2
               D03 = BAS(3,IBAS) - XN3
               DO J1 = -LTMAX(1),LTMAX(1)
                  D11 = D01 + PLAT(1,1)*J1
                  D12 = D02 + PLAT(2,1)*J1
                  D13 = D03 + PLAT(3,1)*J1
                  DO J2 = -LTMAX(2),LTMAX(2)
                     D21 = D11 + PLAT(1,2)*J2
                     D22 = D12 + PLAT(2,2)*J2
                     D23 = D13 + PLAT(3,2)*J2
                     DO J3 = -LTMAX(3),LTMAX(3)
                        D31 = D21 + PLAT(1,3)*J3
                        D32 = D22 + PLAT(2,3)*J3
                        D33 = D23 + PLAT(3,3)*J3
                        RAD = SQRT(D31*D31+D32*D32+D33*D33)*ALAT
                        POTL = POTL + FMESH(A(IC),B(IC),POT(1,IC),NR(IC)
     &                         ,RAD)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
99001 FORMAT (//25('-')//,'make plot for plane',//,'ORIGIN:',
     &        3F21.15/'R1    :',3F21.15,/'R2    :',3F21.15//,25('-'))
99002 FORMAT ('Begin to make POT ...')
      END
C*==nghbr1.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE NGHBR1(W,ALAT,BAS,FACN1,IAX,ICLASS,NBAS,NCLASS,NEIGHM,
     &                  NPR,PLAT,WSR)
C- Create a table of all neighbors within facn1*(wsr(1)+wsr(2))
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci   bas   :basis vectors (scaled by alat)
Ci   clabl :name of the different inequivalent atom
Ci   facn1 :see remarks
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   nbas  :number of atoms in the basis
Ci   nclass:number of classes, atoms in same class are symmetry-related
Ci   neighm:maximum number of neighbors
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   wsr   :Wigner-Seitz sphere radius (in atomic units)
Co Outputs:
Ci   iax   :see remarks
Ci   npr   :number of neighbors around each atom
Cr Remarks:
Cr   Creates a neighour list for a specified atom, generating iax
Cr   which contains all neigbors which fulfill:
Cr
Cr        distance(i,j) <= (wsr(i)+wsr(j))*facn1
Cr
Cr    iax(1): not used
Cr    iax(2): ibas = atom in cluster
Cr    iax(3): i
Cr    iax(4): j
Cr    iax(5): k
Cr   To be sure that at least one pair is found, fac is increased
Cr   by 1.2 until a neighbor is found.
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION ALAT,FACN1
      INTEGER NBAS,NCLASS,NEIGHM
      DOUBLE PRECISION BAS(3,*),PLAT(3,3),WSR(*)
      INTEGER IAX(5,NCLASS,*),ICLASS(*),NPR(*),W(*)
C
C Local variables
C
      DOUBLE PRECISION D2,DR(3),FAC,WI,WJ,WJPWI,WJPWI2
      DOUBLE PRECISION DRR2
      INTEGER I,I1,I2,I3,IBAS,IPR,IPR5,J,JBAS,JC,K,OIWK
      INTEGER ICLBAS
      COMMON /IPRINT/ IPRINT
      INTEGER IPRINT
      CHARACTER*72 MESSG
      SAVE D2,DR,FAC,I,I1,I2,I3,IBAS,IPR,IPR5,J,JBAS,JC,K,MESSG,OIWK,WI,
     &     WJ,WJPWI,WJPWI2
C
C*** End of declarations rewritten by SPAG
C
      CALL DEFI(OIWK,NEIGHM)
      DO JC = 1,NCLASS
         JBAS = ICLBAS(JC,ICLASS,NBAS,1)
         WJ = WSR(JC)
         IPR = 0
         FAC = FACN1
         DO WHILE ( IPR.LT.2 )
            IPR = 0
            DO IBAS = 1,NBAS
               WI = WSR(ICLASS(IBAS))
               WJPWI = (WJ+WI)/ALAT*FAC
               CALL LATLIM(PLAT,WJPWI,I1,I2,I3)
               WJPWI2 = WJPWI*WJPWI
C --------- Sweep lattice translations to find all neighbors
C --------- within wjpwi
               DO I = -I1,I1
                  DO J = -I2,I2
                     DO K = -I3,I3
                        D2 = DRR2(PLAT,BAS(1,JBAS),BAS(1,IBAS),I,J,K,DR)
                        IF ( D2.LE.WJPWI2 ) THEN
                           IF ( IPR.GE.NEIGHM ) THEN
                              WRITE (MESSG,99002) NEIGHM,I1,I2,I3
                              CALL ERRMSG(MESSG,4)
                           END IF
                           IPR5 = IPR*5
                           W(OIWK+IPR5+0) = 10000*D2
                           W(OIWK+IPR5+1) = IBAS
                           W(OIWK+IPR5+2) = I
                           W(OIWK+IPR5+3) = J
                           W(OIWK+IPR5+4) = K
                           IPR = IPR + 1
                        END IF
                     END DO
                  END DO
               END DO
            END DO
            FAC = FAC*1.2D0
         END DO
C
         CALL I_SHELL(5,IPR,W(OIWK))
         DO I = 1,IPR
            CALL INCOPY(5,W(OIWK+5*(I-1)),1,IAX(1,JC,I),1)
            IAX(1,JC,I) = JBAS
         END DO
         NPR(JC) = IPR
C
C
C ----- Printout
         IF ( IPRINT.GE.80 ) THEN
C
            DO IPR = 1,NPR(JC)
               IBAS = IAX(2,JC,IPR)
               I = IAX(3,JC,IPR)
               J = IAX(4,JC,IPR)
               K = IAX(5,JC,IPR)
               D2 = DRR2(PLAT,BAS(1,JBAS),BAS(1,IBAS),I,J,K,DR)
            END DO
            IF (IPRINT > 0) WRITE (6,*)
         END IF
C
      END DO
      IF ( IPRINT.GE.70 ) WRITE (6,99001) NEIGHM,(NPR(JC),JC=1,NCLASS)
C
99001 FORMAT (/' NGHBR1: neighm=',i4,'> npr=',50(11I4,/26x))
99002 FORMAT (' NGHBR1: too many pairs, neighm=',i3,6x,'i1,i2,i3 =',3I3,
     &        '$')
      END
C*==potmax.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      SUBROUTINE POTMAX(W,A,ALAT,B,BAS,IAX1,ICLASS,NBAS,NCLASS,NPR1,NR,
     &                  NRMAX,PLAT,POT,WMAX,WSR,WSREST,Z)
C- Calculates the overlapping Hartree potential and finds it's maximum
C ----------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   alat  :length scale
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   bas   :basis vectors (scaled by alat)
Ci   clabl :name of the different inequivalent atom
Ci   iax1  :information about positions around a specified atom
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   nbas  :number of atoms in the basis
Ci   nclass:number of classes, atoms in same class are symmetry-related
Ci   npr1  :number of neighbors
Ci   nr    :number of mesh points
Ci   nrmax :maximum number of mesh points
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   pot   :spherical Hartree potential
Ci   wmax  :largest wsr/avw
Ci   wsrest:estimation of Wigner-Seitz sphere radius (in atomic units)
Ci   z     :nuclear charge
Co Outputs:
Co   wsr   :Wigner-Seitz sphere radius (in atomic units)
Cr Remarks:
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NSTEP
      PARAMETER (NSTEP=20)
      DOUBLE PRECISION FACN2,TINY,TOL
      PARAMETER (FACN2=3.D0,TINY=1.D-5,TOL=1.D-2)
C
C Dummy arguments
C
      DOUBLE PRECISION ALAT,WMAX
      INTEGER NBAS,NCLASS,NRMAX
      DOUBLE PRECISION A(*),B(*),BAS(3,*),PLAT(3,3),POT(NRMAX,*),WSR(*),
     &                 WSREST(*),Z(*)
      INTEGER IAX1(5,NCLASS,*),ICLASS(*),NPR1(*),NR(*),W(*)
C
C Local variables
C
      DOUBLE PRECISION DLAMCH,DRR2
      DOUBLE PRECISION DR(3),POTL,RIK,RIKO,RPMX,STEP,X(:),XN(3),Y(:)
      INTEGER I,IBAS,IC,IMIN,ISTEP,K,K1,K2,K3,KBAS,KC,KCO,KPR,NEIGHM,
     &        NPR2,OIAX2
      COMMON /IPRINT/ IPRINT
      INTEGER IPRINT
      INTEGER ICLBAS
      SAVE DR,I,IBAS,IC,IMIN,ISTEP,K,K1,K2,K3,KBAS,KC,KCO,KPR,NEIGHM,
     &     NPR2,OIAX2,POTL,RIK,RIKO,RPMX,STEP,X,XN,Y
      EXTERNAL DEFI,DLAMCH,DRR2,DSCAL,ICLBAS,NGHBR2,POTXN
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATABLE x,y
      ALLOCATE (X(0:NSTEP),Y(0:NSTEP))
      NEIGHM = MAX0(INT(4.D0*(FACN2*WMAX+1.D0)**3),50)
      CALL DEFI(OIAX2,5*NEIGHM)
C
      DO IC = 1,NCLASS
         IF ( IDNINT(Z(IC)).NE.0 ) THEN
            IF ( IPRINT.GE.40 ) WRITE (6,99001)
            IBAS = ICLBAS(IC,ICLASS,NBAS,1)
C
            WSR(IC) = DLAMCH('o')
            RIKO = -1.D0
            KCO = -1
            DO KPR = 2,NPR1(IC)
               KBAS = IAX1(2,IC,KPR)
               K1 = IAX1(3,IC,KPR)
               K2 = IAX1(4,IC,KPR)
               K3 = IAX1(5,IC,KPR)
               KC = ICLASS(KBAS)
               RIK = SQRT(DRR2(PLAT,BAS(1,IBAS),BAS(1,KBAS),K1,K2,K3,DR)
     &               )*ALAT
               IF ( IDNINT(Z(KC)).NE.0 .AND. 
     &              (ABS(RIK-RIKO).GT.TINY .OR. KCO.NE.KC) ) THEN
C
                  RIKO = RIK
                  KCO = KC
                  STEP = RIK/NSTEP
                  CALL NGHBR2(ALAT,BAS,FACN2,W(OIAX2),IBAS,ICLASS,K1,K2,
     &                        K3,KBAS,NBAS,NEIGHM,NPR2,PLAT,WSREST)
                  CALL DSCAL(3,1.D0/DFLOAT(NSTEP),DR,1)
C --------  Find aproximate position of first maximun
                  DO ISTEP = 0,NSTEP
                     DO K = 1,3
                        XN(K) = BAS(K,IBAS) + ISTEP*DR(K)
                     END DO
                     CALL POTXN(A,ALAT,B,BAS,W(OIAX2),ICLASS,NPR2,NR,
     &                          NRMAX,PLAT,POT,POTL,XN)
                     X(ISTEP) = DBLE(ISTEP)*STEP
                     Y(ISTEP) = MAX(POTL,-1.D5)
                     IF ( ISTEP.NE.0 ) THEN
                        IF ( ISTEP.GE.2 .AND. Y(ISTEP).LT.Y(ISTEP-1) )
     &                       GOTO 5
                     END IF
                  END DO
                  CALL ERRMSG('bug in potmax.$',4)
C
C --------- Now a single maximum is supposed between istep-2 and istep
 5                CONTINUE
                  DO K = 1,3
                     X(K) = X(ISTEP-3+K)
                     Y(K) = Y(ISTEP-3+K)
                  END DO
                  RPMX = 0.D0
                  IMIN = (ISTEP-2)*2
                  DO WHILE ( X(3)-X(1).GT.TOL )
                     X(5) = X(3)
                     X(3) = X(2)
                     Y(5) = Y(3)
                     Y(3) = Y(2)
                     STEP = STEP*0.5D0
                     CALL DSCAL(3,0.5D0,DR,1)
                     DO I = 1,3,2
                        DO K = 1,3
                           XN(K) = BAS(K,IBAS) + (IMIN+I)*DR(K)
                        END DO
                        CALL POTXN(A,ALAT,B,BAS,W(OIAX2),ICLASS,NPR2,NR,
     &                             NRMAX,PLAT,POT,Y(I+1),XN)
                        X(I+1) = (IMIN+I)*STEP
                     END DO
                     DO I = 2,4
                        IF ( Y(I+1).LT.Y(I) ) THEN
                           DO K = 1,3
                              Y(K) = Y(I-2+K)
                              X(K) = X(I-2+K)
                           END DO
                           EXIT
                        END IF
                     END DO
                     IMIN = 2*(IMIN+I-2)
                  END DO
                  RPMX = 0.5D0*(Y(1)*(X(2)+X(3))-2.D0*Y(2)*(X(1)+X(3))
     &                   +Y(3)*(X(1)+X(2)))/(Y(1)-Y(2)-Y(2)+Y(3))
                  WSR(IC) = MIN(WSR(IC),RPMX)
               END IF
            END DO
            IF ( IPRINT.GE.40 ) WRITE (6,99002)
C        if (iprint().ge.20) write(6,303)clabl(ic),wsr(ic)
         END IF
      END DO
C
99001 FORMAT (/' POTMAX: ATOM1   ATOM2   DIST      VMAX   R(VMAX)',/,8x,
     &        42('-'))
99002 FORMAT (8x,42('-'))
C
      DEALLOCATE(x,y)
      END
C*==nghbr2.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE NGHBR2(ALAT,BAS,FACN2,IAX,IBAS,ICLASS,K1,K2,K3,KBAS,
     &                  NBAS,NEIGHM,NPR,PLAT,WSR)
C- Create a table of all neighbors around two atoms
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci   bas   :basis vectors (scaled by alat)
Ci   clabl :name of the different inequivalent atom
Ci   facn2 :see remarks
Ci   ibas  :index to first atom
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   k1-3  :second atom is translated by k1*plat1+k2*plat2+p3*plat3
Ci   kbas  :index to second atom
Ci   nbas  :number of atoms in the basis
Ci   neighm:maximum number of neighbors
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   wsr   :Wigner-Seitz sphere radius (in atomic units)
Co Outputs:
Ci   iax   :see remarks
Ci   npr   :number of neighbors
Cr Remarks:
Cr   Creates a neighour list for a pair of atoms ibas and kbas,
Cr   generating iax which contains all neigbors jbas within
Cr   facn2*wsr(jbas) around ibas or kbas , i.e. those which fulfill:
Cr
Cr        distance(i,j) <= wsr(j)*facn2
Cr    or  distance(k,j) <= wsr(j)*facn2
Cr
Cr    iax(1): not used
Cr    iax(2): ibas = atom in cluster
Cr    iax(3): i
Cr    iax(4): j
Cr    iax(5): k
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION ALAT,FACN2
      INTEGER IBAS,K1,K2,K3,KBAS,NBAS,NEIGHM,NPR
      DOUBLE PRECISION BAS(3,*),PLAT(3,3),WSR(*)
      INTEGER IAX(5,*),ICLASS(*)
C
C Local variables
C
      DOUBLE PRECISION B(3,2),D1,D2,DB11,DB12,DB21,DB22,DB31,DB32,DR11,
     &                 DR12,DR13,DR21,DR22,DR23,DR31,DR32,DR33,DRB11,
     &                 DRB12,DRB21,DRB22,DRB31,DRB32,WJF,WJF2
      INTEGER I,I1,I2,I3,II,IPR,J,JBAS,JC,K
      COMMON /IPRINT/ IPRINT
      INTEGER IPRINT

      CHARACTER*72 MESSG
      SAVE B,D1,D2,DB11,DB12,DB21,DB22,DB31,DB32,DR11,DR12,DR13,DR21,
     &     DR22,DR23,DR31,DR32,DR33,DRB11,DRB12,DRB21,DRB22,DRB31,DRB32,
     &     I,I1,I2,I3,II,IPR,J,JBAS,JC,K,MESSG,WJF,WJF2
      EXTERNAL DCOPY, LATLIM
C
C*** End of declarations rewritten by SPAG
C
      CALL DCOPY(3,BAS(1,IBAS),1,B(1,1),1)
      B(1,2) = BAS(1,KBAS) + K1*PLAT(1,1) + K2*PLAT(1,2) + K3*PLAT(1,3)
      B(2,2) = BAS(2,KBAS) + K1*PLAT(2,1) + K2*PLAT(2,2) + K3*PLAT(2,3)
      B(3,2) = BAS(3,KBAS) + K1*PLAT(3,1) + K2*PLAT(3,2) + K3*PLAT(3,3)
      IPR = 0
      DO JBAS = 1,NBAS
         JC = ICLASS(JBAS)
         WJF = WSR(JC)/ALAT*FACN2
         CALL LATLIM(PLAT,WJF,I1,I2,I3)
         WJF2 = WJF*WJF
         DB11 = BAS(1,JBAS) - B(1,1)
         DB21 = BAS(2,JBAS) - B(2,1)
         DB31 = BAS(3,JBAS) - B(3,1)
         DB12 = BAS(1,JBAS) - B(1,2)
         DB22 = BAS(2,JBAS) - B(2,2)
         DB32 = BAS(3,JBAS) - B(3,2)
         DO I = -I1,I1
            DR11 = PLAT(1,1)*I
            DR21 = PLAT(2,1)*I
            DR31 = PLAT(3,1)*I
            DO J = -I2,I2
               DR12 = DR11 + PLAT(1,2)*J
               DR22 = DR21 + PLAT(2,2)*J
               DR32 = DR31 + PLAT(3,2)*J
               DO K = -I3,I3
                  DR13 = DR12 + PLAT(1,3)*K
                  DR23 = DR22 + PLAT(2,3)*K
                  DR33 = DR32 + PLAT(3,3)*K
                  DRB11 = DR13 + DB11
                  DRB21 = DR23 + DB21
                  DRB31 = DR33 + DB31
                  DRB12 = DR13 + DB12
                  DRB22 = DR23 + DB22
                  DRB32 = DR33 + DB32
                  D1 = DRB11*DRB11 + DRB21*DRB21 + DRB31*DRB31
                  D2 = DRB12*DRB12 + DRB22*DRB22 + DRB32*DRB32
                  IF ( D1.LE.WJF2 .OR. D2.LE.WJF2 ) THEN
                     IPR = IPR + 1
                     IF ( IPR.GE.2*NEIGHM ) THEN
                        WRITE (MESSG,99002) 2*NEIGHM,I1,I2,I3
                        CALL ERRMSG(MESSG,4)
                     END IF
                     IAX(2,IPR) = JBAS
                     IAX(3,IPR) = I
                     IAX(4,IPR) = J
                     IAX(5,IPR) = K
                     II = 0
                     IF ( D1.LE.WJF2 ) II = II + 1
                     IF ( D2.LE.WJF2 ) II = II + 2
C ------------- Printout
                     IF ( IPRINT.GE.80 ) WRITE (6,99001) IPR,D1,D2,
     &                    JBAS,I,J,K,II
                  END IF
               END DO
            END DO
         END DO
      END DO
C
      NPR = IPR
C
99001 FORMAT (i3,2F21.15,4I5,2x,i3)
99002 FORMAT (' NGHBR2: too many pairs, neighm=',i3,6x,'i1,i2,i3 =',3I3,
     &        '$')
      END
C*==blowup.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      SUBROUTINE BLOWUP(ALAT,BAS,FACVOL,GAMMA,IAX1,ICLASS,NBAS,NCLASS,
     &                  NPR1,OMMAX1,OMMAX2,PLAT,WSR,DSCL,LOCK)
C- Blows up the spheres
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci   bas   :basis vectors (scaled by alat)
Ci   clabl :name of the different inequivalent atom
Ci   iax1  :information about positions around a specified atom
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   nbas  :number of atoms in the basis
Ci   nclass:number of classes, atoms in same class are symmetry-related
Ci   npr1  :number of neighbors
Ci   facvol:sum of sphere volumes = facvol * cell volume
Ci   gamma :scaling is r -> a(r+b) with a*b=gamma*(a-1)*avw
Ci   ommax1:maximum overlap divided by distance (s1+s2-d)/d  < ommax1
Ci   ommax2:maximum overlap divided by  radius  (s1+s2-d)/s1 < ommax2
Ci   plat  :primitive lattice vectors (scaled by alat)
Cio  wsr   :Wigner-Seitz sphere radius (in atomic units)
Co Outputs:
Cio  wsr   :Wigner-Seitz sphere radius (in atomic units)
Cr Remarks:
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      DOUBLE PRECISION FPI3,TINY
      PARAMETER (FPI3=4.18879020478639053D0,TINY=1.D-4)
C
C Dummy arguments
C
      DOUBLE PRECISION ALAT,FACVOL,GAMMA,OMMAX1,OMMAX2
      INTEGER NBAS,NCLASS
      DOUBLE PRECISION BAS(3,*),DSCL(NCLASS,0:NCLASS+1),PLAT(3,3),WSR(*)
      INTEGER IAX1(5,NCLASS,*),ICLASS(*),NPR1(*)
      LOGICAL LOCK(*)
C
C Local variables
C
      DOUBLE PRECISION A,AMAX,AMAX1,AMAX2,AMAX3,AVW,B,BMAX,D,DM(0:3),
     &                 DR(3),DSCLMX,GW,OMO2,OPO1,P,Q,R,RIK,RIKO,S,T,
     &                 TMO2,U,V,VOL,VOLA,VOLB,VOLSPH,WSRI,WSRK,X
      DOUBLE PRECISION DDET33,DLAMCH,DRR2
      LOGICAL FIN
      INTEGER IBAS,IC,ICMAX,ICMIN,ILOOP,K1,K2,K3,KBAS,KC,KCO,KPR,LINE,
     &        NLOOP
      INTEGER ICLBAS
      COMMON /IPRINT/ IPRINT
      INTEGER IPRINT
      CHARACTER*72 MESSG
      SAVE A,AMAX,AMAX1,AMAX2,AMAX3,AVW,B,BMAX,D,DM,DR,DSCLMX,FIN,GW,
     &     IBAS,IC,ICMAX,ICMIN,ILOOP,K1,K2,K3,KBAS,KC,KCO,KPR,LINE,
     &     MESSG,NLOOP,OMO2,OPO1,P,Q,R,RIK,RIKO,S,T,TMO2,U,V,VOL,VOLA,
     &     VOLB,VOLSPH,WSRI,WSRK,X
      EXTERNAL DLAMCH,DRR2,ERRMSG,ICLBAS
C
C*** End of declarations rewritten by SPAG
C
      OPO1 = 1.D0 + OMMAX1
      OMO2 = 1.D0 - OMMAX2
      TMO2 = 2.D0 - OMMAX2
      VOL = ABS(DDET33(PLAT))*ALAT**3
      AVW = (VOL/FPI3/NBAS)**(1.D0/3.D0)
      GW = AVW*GAMMA
C
      VOLSPH = 0.D0
      DO IBAS = 1,NBAS
         VOLSPH = VOLSPH + FPI3*WSR(ICLASS(IBAS))**3
      END DO
      IF ( IPRINT.GE.30 ) WRITE (6,99009) VOL,VOLSPH
C
C     call dcopy(nclass,1.d0,0,dscl(1,0),1)
      DSCL(1:NCLASS,0) = 1.D0
C
      DO NLOOP = 1,NCLASS + 1
         AMAX = DLAMCH('o')
C ----- lock radii of those spheres whith maximum allowed overlap
         LOCK(1:NCLASS) = .FALSE.
         DO IC = 1,NCLASS
            IBAS = ICLBAS(IC,ICLASS,NBAS,1)
            WSRI = WSR(IC)
            DO KPR = 2,NPR1(IC)
               KBAS = IAX1(2,IC,KPR)
               K1 = IAX1(3,IC,KPR)
               K2 = IAX1(4,IC,KPR)
               K3 = IAX1(5,IC,KPR)
               KC = ICLASS(KBAS)
               WSRK = WSR(KC)
               RIK = SQRT(DRR2(PLAT,BAS(1,IBAS),BAS(1,KBAS),K1,K2,K3,DR)
     &               )
               RIK = RIK*ALAT
               IF ( ABS(OPO1*RIK-WSRI-WSRK).LT.TINY ) LOCK(IC) = .TRUE.
               IF ( ABS(RIK-OMO2*WSRI-WSRK).LT.TINY ) LOCK(IC) = .TRUE.
               IF ( ABS(RIK-WSRI-OMO2*WSRK).LT.TINY ) LOCK(IC) = .TRUE.
            END DO
         END DO
C
         DO IC = 1,NCLASS
            IF ( .NOT.LOCK(IC) ) THEN
               IBAS = ICLBAS(IC,ICLASS,NBAS,1)
               RIKO = -1.D0
               KCO = -1
               WSRI = WSR(IC)
               DO KPR = 2,NPR1(IC)
                  KBAS = IAX1(2,IC,KPR)
                  K1 = IAX1(3,IC,KPR)
                  K2 = IAX1(4,IC,KPR)
                  K3 = IAX1(5,IC,KPR)
                  KC = ICLASS(KBAS)
                  RIK = SQRT(DRR2(PLAT,BAS(1,IBAS),BAS(1,KBAS),K1,K2,K3,
     &                  DR))
                  RIK = RIK*ALAT
                  IF ( ABS(RIK-RIKO).GT.TINY .OR. KCO.NE.KC ) THEN
                     WSRK = WSR(KC)
                     RIKO = RIK
                     KCO = KC
                     IF ( LOCK(KC) ) THEN
                        AMAX1 = (OPO1*RIK-WSRK+GW)/(WSRI+GW)
                        IF ( OMO2.GT.0.D0 ) AMAX2 = (RIK-WSRK+OMO2*GW)
     &                       /(WSRI+GW)/OMO2
                        AMAX3 = (RIK-OMO2*WSRK+GW)/(WSRI+GW)
                     ELSE
                        AMAX1 = (OPO1*RIK+GW+GW)/(WSRI+WSRK+GW+GW)
                        IF ( WSRK+OMO2*WSRI+TMO2*GW.GT.0.D0 )
     &                       AMAX2 = (RIK+TMO2*GW)
     &                       /(WSRK+OMO2*WSRI+TMO2*GW)
                        IF ( WSRI+OMO2*WSRK+TMO2*GW.GT.0.D0 )
     &                       AMAX3 = (RIK+TMO2*GW)
     &                       /(WSRI+OMO2*WSRK+TMO2*GW)
                     END IF
                     AMAX = MIN(AMAX,AMAX1,AMAX2,AMAX3)
                  END IF
               END DO
            END IF
         END DO
         BMAX = (1.D0-1.D0/AMAX)*GW
C
         VOLA = 0.D0
         VOLB = 0.D0
         DO IBAS = 1,NBAS
            IC = ICLASS(IBAS)
            IF ( .NOT.LOCK(IC) ) THEN
               VOLB = VOLB + (AMAX*WSR(IC)+BMAX)**3
            ELSE
               VOLA = VOLA + WSR(IC)**3
            END IF
         END DO
         VOLA = VOLA*FPI3
         VOLB = VOLB*FPI3
C
         IF ( VOL*FACVOL.LT.VOLA+VOLB ) THEN
C
            DM(:) = 0.0D0
            DO IBAS = 1,NBAS
               IC = ICLASS(IBAS)
               IF ( .NOT.LOCK(IC) ) THEN
                  A = WSR(IC) + GW
C ----------- for numerical reasons distinguish cases
                  IF ( ABS(GAMMA).GT.1.D0 ) THEN
                     B = WSR(IC)
                  ELSE
                     B = -GW
                  END IF
                  DM(0) = DM(0) + B*B*B
                  DM(1) = DM(1) + 3.D0*A*B*B
                  DM(2) = DM(2) + 3.D0*A*A*B
                  DM(3) = DM(3) + A*A*A
               END IF
            END DO
            IF ( ABS(DM(3)).GT.TINY ) THEN
               R = DM(2)/DM(3)
               S = DM(1)/DM(3)
               T = (DM(0)-(VOL*FACVOL-VOLA)/FPI3)/DM(3)
               P = S - R*R/3.D0
               Q = 2.D0*R*R*R/27.D0 - R*S/3.D0 + T
               D = P*P*P/27.D0 + Q*Q/4.D0
               U = (SQRT(D)-Q/2.D0)**(1.D0/3.D0)
               V = -P/U/3.D0
               X = U + V - R/3.D0
               IF ( ABS(GAMMA).GT.1.D0 ) THEN
                  AMAX = X + 1.D0
                  BMAX = X*GW/AMAX
               ELSE
                  AMAX = X
                  BMAX = (1.D0-1.D0/AMAX)*GW
               END IF
               IF ( IPRINT.GE.100 ) THEN
                  WRITE (6,99008) 'R S T',R,S,T
                  WRITE (6,99008) 'P Q  ',P,Q
                  WRITE (6,99008) '  D  ',D
                  WRITE (6,99008) ' U V ',U,V
                  WRITE (6,99008) ' AMAX',AMAX
                  WRITE (6,99008) ' BMAX',BMAX
                  WRITE (6,99008) ' -------------------------'
               END IF
            END IF
         END IF
C
         FIN = .TRUE.
         DSCLMX = 1.D0
         DO IC = 1,NCLASS
            DSCL(IC,NLOOP) = 1.D0
            IF ( .NOT.LOCK(IC) ) THEN
               DSCLMX = AMAX + BMAX/WSR(IC)
               WSR(IC) = DSCLMX*WSR(IC)
               DSCL(IC,0) = DSCLMX*DSCL(IC,0)
               DSCL(IC,NLOOP) = DSCLMX
               FIN = .FALSE.
            END IF
         END DO
         FIN = FIN .OR. ABS(DSCLMX-1.D0).LT.TINY/256
C
         VOLSPH = 0.D0
         DO IBAS = 1,NBAS
            VOLSPH = VOLSPH + FPI3*WSR(ICLASS(IBAS))**3
         END DO
C
C
         IF ( FIN .OR. NLOOP.EQ.NCLASS+1 ) THEN
            IF ( IPRINT.GE.30 ) THEN
               WRITE (6,99001) OMMAX1,OMMAX2
               DO LINE = 0,(NCLASS-1)/6
                  ICMIN = 1 + 6*LINE
                  ICMAX = MIN0(NCLASS,6+6*LINE)
C              write(6,301)(clabl(ic),ic=icmin,icmax)
                  WRITE (6,99007) ('============',IC=ICMIN,ICMAX+1)
                  WRITE (6,99002) (WSR(IC)/DSCL(IC,0),IC=ICMIN,ICMAX)
                  WRITE (6,99007) ('------------',IC=ICMIN,ICMAX+1)
                  DO ILOOP = 1,NLOOP
                     WRITE (6,99003) ILOOP,
     &                               (DSCL(IC,ILOOP),IC=ICMIN,ICMAX)
                  END DO
                  WRITE (6,99007) ('------------',IC=ICMIN,ICMAX+1)
                  WRITE (6,99004) (WSR(IC),IC=ICMIN,ICMAX)
                  WRITE (6,99005) (DSCL(IC,0),IC=ICMIN,ICMAX)
                  WRITE (6,99006) (WSR(IC)-WSR(IC)/DSCL(IC,0),IC=ICMIN,
     &                            ICMAX)
                  WRITE (6,99007) ('------------',IC=ICMIN,ICMAX+1)
               END DO
               WRITE (6,99009) VOL,VOLSPH
            END IF
            IF ( ABS(OMMAX1*OMMAX2).GT.TINY .AND. ABS(VOLSPH-VOL*FACVOL)
     &           .GT.TINY ) THEN
               WRITE (MESSG,99010)
               CALL ERRMSG(MESSG,0)
            END IF
            RETURN
         END IF
      END DO
C
99001 FORMAT (/' BLOWUP: scale radii;   OMMAX1=',f21.15,', OMMAX2=',
     &        f21.15)
99002 FORMAT (6x,'OLD WSR:    ',6F21.15)
99003 FORMAT (6x,'LOOP:',i3,'  * ',6F21.15)
99004 FORMAT (6x,'NEW WSR:    ',6F21.15)
99005 FORMAT (6x,'W_NEW/W_OLD:',6F21.15)
99006 FORMAT (6x,'W_NEW-W_OLD:',6F21.15)
99007 FORMAT (6x,a12,6A9)
99008 FORMAT (6x,a,3F21.15)
99009 FORMAT (/' CELL VOLUME=',f21.15,'   SUM OF SPHERE VOLUMES=',
     &        f21.15)
99010 FORMAT (' BLOWUP: impossible to reach VOL, increase OMMAX.$')
C
      END
C*==delsqf.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      DOUBLE PRECISION FUNCTION DELSQF(NDATA,IX,YA,X)
C- exponential least square fit (  a * exp(bx) )
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ndata :number of data points
Ci   ix    :xa(1) is integer and equals ix
Ci          and xa(1),xa(2),xa(3),... differ exactly by 1.d0
Ci   ya    :value of f at xa
Ci   x     :x-value at which f is interpolated
Co Outputs:
Co   delsqf:interpolated value
Cr Remarks:
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IX,NDATA
      DOUBLE PRECISION X
      DOUBLE PRECISION YA(*)
C
C Local variables
C
      DOUBLE PRECISION A,B,SS,SX,SXX,SXY,SY,XX,YY
      INTEGER IDATA
      SAVE A,B,IDATA,SS,SX,SXX,SXY,SY,XX,YY
C
C*** End of declarations rewritten by SPAG
C
      SX = 0.D0
      SY = 0.D0
      SXY = 0.D0
      SXX = 0.D0
C
      DO IDATA = 1,NDATA
         XX = DBLE(IX+IDATA-1)
         IF ( YA(IDATA).LE.0.D0 ) THEN
            DELSQF = 0.D0
            RETURN
         END IF
         YY = LOG(YA(IDATA))
         SX = SX + XX
         SY = SY + YY
         SXX = SXX + XX*XX
         SXY = SXY + XX*YY
      END DO
      SS = DBLE(NDATA)
C
      B = (SXY-SX*SY/SS)/(SXX-SX*SX/SS)
      A = EXP((SY-B*SX)/SS)
C
      DELSQF = A*EXP(B*X)
C
      END
C*==drr2.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      DOUBLE PRECISION FUNCTION DRR2(PLAT,BAS1,BAS2,I,J,K,DR)
C- Calculates the vector connecting two sites in a solid
C ----------------------------------------------------------------
Ci Inputs:
Ci   plat  :primitive lattice vectors
Ci   bas1  :basis vector of first site
Ci   bas2  :basis vector of second site
Ci   i,j,k :the number of primitive lattice vectors separating sites
Co Outputs:
Co   dr    :connecting vector bas2 - bas1
Co   drr2  :square of the length of this vector
Cr Remarks:
Cr   Using the TB package and a table of indices iax, the connecting
Cr   vector and the square of the distance is obtained by:
Cr   rsqr=drr2(plat,bas(1,iax(1)),bas(1,iax(2)),iax(3),iax(4),iax(5),dr)
C ----------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I,J,K
      DOUBLE PRECISION BAS1(3),BAS2(3),DR(3),PLAT(3,3)
C
C Local variables
C
      INTEGER IX
C
C*** End of declarations rewritten by SPAG
C
      DRR2 = 0.D0
      DO IX = 1,3
         DR(IX) = BAS2(IX) - BAS1(IX) + PLAT(IX,1)*I + PLAT(IX,2)
     &            *J + PLAT(IX,3)*K
         DRR2 = DRR2 + DR(IX)*DR(IX)
      END DO
      END
C*==iclbas.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      INTEGER FUNCTION ICLBAS(IC,ICLASS,NBAS,NRBAS)
C- Returns an index to iclbas atom in basis given class
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ic    :class label
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   nbas  :number of atoms in the basis
Ci   nrbas :the nrbas-th basis atom of class ic is seeked
Co Outputs:
Co   iclbas:the iclbas-th atom belongs to class ic
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IC,NBAS,NRBAS
      INTEGER ICLASS(*)
C
C Local variables
C
      INTEGER IBAS,N
      CHARACTER*72 MESSG
      SAVE IBAS,MESSG,N
      EXTERNAL ERRMSG
C
C*** End of declarations rewritten by SPAG
C
      ICLBAS = 1
      N = 0
      DO IBAS = 1,NBAS
         IF ( ICLASS(IBAS).EQ.IC ) N = N + 1
         IF ( N.EQ.NRBAS ) THEN
            ICLBAS = IBAS
            RETURN
         END IF
      END DO
      WRITE (MESSG,99001) IC,NRBAS,N
      CALL ERRMSG(MESSG,4)
C
99001 FORMAT (' ICLBAS: class',i3,', nrbas=',i3,', but only',i3,
     &        ' basis atoms exist.$')
      END
C*==latlim.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
      SUBROUTINE LATLIM(PQLAT,RMAXS,I1,I2,I3)
C- Set limits in X Y Z direction
C ----------------------------------------------------------------------
Ci Inputs:
Ci   pqlat :primitive lattice vectors (real or reciprocal space)
Ci   rmaxs :maximum length of connecting vector
Co Outputs:
Co   i1,i2,i3:all connecting vectors lie within these multiples of the
Co            lattice vectors.
Cr Remarks:
Cr   This routine only returns the integer part
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I1,I2,I3
      DOUBLE PRECISION RMAXS
      DOUBLE PRECISION PQLAT(3,3)
C
C Local variables
C
      DOUBLE PRECISION DET,QPLAT(3,3),X1,X2,X3
      DOUBLE PRECISION DLAMCH,DNRM23
      INTEGER IPRINT
      COMMON /IPRINT/ IPRINT
      SAVE DET,QPLAT,X1,X2,X3
      EXTERNAL DINV33,DLAMCH,DNRM23
C
C*** End of declarations rewritten by SPAG
C
      CALL DINV33(PQLAT,1,QPLAT,DET)
      X1 = RMAXS*SQRT(DNRM23(QPLAT(1,1)))
      X2 = RMAXS*SQRT(DNRM23(QPLAT(1,2)))
      X3 = RMAXS*SQRT(DNRM23(QPLAT(1,3)))
C
      I1 = 1 + INT(X1-DLAMCH('p'))
      I2 = 1 + INT(X2-DLAMCH('p'))
      I3 = 1 + INT(X3-DLAMCH('p'))
      IF ( IPRINT.GE.80 ) WRITE (6,99001) RMAXS,X1,X2,X3,I1,I2,I3
C
99001 FORMAT (/' LATLIM: rmaxs=',f21.15,', x1,x2,x3=',3F21.15,
     &        ', i1,i2,i3=',3I2)
      END
C*==defwsr.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE DEFWSR(WSR,Z)
C- Returns default value of radii for given nuclear charge
C ----------------------------------------------------------------------
Ci Inputs:
Ci   z     :nuclear charge
Co Outputs:
Co   wsr   :Wigner-Seitz sphere radius (in atomic units)
Cr Remarks:
Cr These are equilibrium average Wigner-Seitz radii for closed-packed
Cr structure of pure element.
C ----------------------------------------------------------------------
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION WSR,Z
C
C Local variables
C
      DOUBLE PRECISION CPWSR(0:100)
      INTEGER IZ
      CHARACTER*72 MESSG
      SAVE IZ,MESSG
      EXTERNAL ERRMSG
C
C*** End of declarations rewritten by SPAG
C
      DATA CPWSR/2.00D0,1.39D0,2.55D0,3.04D0,2.27D0,1.96D0,1.66D0,
     &     1.90D0,1.90D0,2.17D0,2.89D0,3.76D0,3.25D0,2.95D0,2.63D0,
     &     2.56D0,2.70D0,2.85D0,3.71D0,4.66D0,3.88D0,3.31D0,2.99D0,
     &     2.76D0,2.64D0,2.57D0,2.52D0,2.52D0,2.55D0,2.62D0,2.78D0,
     &     2.75D0,2.79D0,2.83D0,2.94D0,3.13D0,4.32D0,4.95D0,4.22D0,
     &     3.61D0,3.28D0,3.03D0,2.91D0,2.82D0,2.77D0,2.78D0,2.84D0,
     &     2.95D0,3.14D0,3.30D0,3.45D0,3.30D0,3.31D0,3.50D0,4.31D0,
     &     5.30D0,4.20D0,3.91D0,3.80D0,3.75D0,3.70D0,3.65D0,3.60D0,
     &     3.55D0,3.52D0,3.61D0,3.67D0,3.70D0,3.73D0,3.75D0,3.56D0,
     &     3.44D0,3.23D0,3.04D0,2.93D0,2.86D0,2.82D0,2.83D0,2.88D0,
     &     2.98D0,3.27D0,3.57D0,3.62D0,3.37D0,3.46D0,3.63D0,4.44D0,
     &     5.81D0,4.30D0,3.84D0,3.52D0,3.32D0,3.13D0,3.02D0,2.96D0,
     &     2.93D0,2.93D0,2.95D0,2.99D0,3.05D0,3.17D0/
C
      IZ = IDNINT(Z)
      IF ( IZ.LT.0 .OR. IZ.GT.100 ) THEN
         WRITE (MESSG,99001) IZ
         CALL ERRMSG(MESSG,0)
         WSR = CPWSR(100)
      ELSE
         WSR = CPWSR(IZ)
      END IF
99001 FORMAT ('DEFWSR: bad nuclear charge, Z=',i3)
      END
