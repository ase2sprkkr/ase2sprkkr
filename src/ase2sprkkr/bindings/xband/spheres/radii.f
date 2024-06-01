      SUBROUTINE FIND_EMPTY_SPHERES(
     >   N_OUT,   ! maximum number of outputs
     >   CENTRES, ! out: centres of the empty spheres
     >   RADII,   ! out: radii of the empty spheres, n_out:n_out+nQ
     >            ! ws-radii of atoms
     >   RMINES_, ! minimal radius (1.2)
     >   RMAXES_, ! maximal radius (2.5)
     >   ALAT,
     >   CELL,    ! 3x3 array, primitive vectors
     >   NQ,      ! number of atoms
     >   BAS,     ! atom positions, cartesian
     >   IMQ,     ! atom equivalence class
     >   NM,      ! number of classes (equivalent sites)
     >   NSORT,   ! number of core types
     >   TXTT,    ! type of the core type (chemical symbol)
     >   Z,       ! atomic numbers of the core types
     >   CONC,    ! occupations of the core types
     >   IMT,     ! the core type is used by the class nb.
     >   NG,      ! number of symmetry operations
     >   ITOP,    ! first line of point symmetry operations
     >   IQA,     ! second line of point symmetry operations
     >   VERBOSE  ! print output to the stdout
     > )
      IMPLICIT NONE

C PARAMETER definitions
C
      INTEGER NTMAX,NRMAX,NRAD1,MAXDIM,MAX2SORT
      PARAMETER (NTMAX=200,NRMAX=300,NRAD1=251,MAXDIM=1000000,
     &           MAX2SORT=NTMAX)

      INTEGER N_OUT
      DOUBLE PRECISION CENTRES(3, MAXDIM)
      DOUBLE PRECISION RADII(MAXDIM)

      DOUBLE PRECISION RMINES_, RMAXES_
      DOUBLE PRECISION CELL(3,3)
      INTEGER NQ
      DOUBLE PRECISION BAS(3,NTMAX)
      INTEGER VERBOSE
      INTEGER IPRINT
      COMMON /IPRINT/ IPRINT
C*==aa0001.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
C
C*** Start of declarations rewritten by SPAG
C
C
C COMMON variables
C
      DOUBLE PRECISION AV(3),BV(3),CV(3),QLAT(3,3),RMAXES,RMINES
      INTEGER N1ES,N2ES,N3ES,NATES
      COMMON /EMPTYS/ NATES,N1ES,N2ES,N3ES,RMINES,RMAXES
      COMMON /L2LAT / AV,BV,CV,QLAT
C
C Local variables
C
      DOUBLE PRECISION ALAT,AMTC(NRMAX,NTMAX),
     &                 BASNEW(3,MAX2SORT),CONC(NTMAX),DPAS,R0(NTMAX),
     &                 R0ACT,R0SITE(MAX2SORT),RAD,RAT(NRAD1,MAX2SORT),
     &                 RHOSITE(NRAD1,MAX2SORT),RINT,RWS(NTMAX),
     &                 SES(MAX2SORT),TAUES(3,MAX2SORT),WS,
     &                 WSREST(MAX2SORT),Z(NTMAX),ZN,ZZ(MAX2SORT)
      LOGICAL DB
      DOUBLE PRECISION DNEVMOD
      CHARACTER*100 FFF
      INTEGER I,IATOM,IDUM,IMQ(NTMAX),IMT(NTMAX),IPNT,IQA(48),
     &        ISNEW(MAX2SORT),ISORT,ISP1,ISP2,ISR,IST,IT,ITOP(48),ITYPE,
     &        J,NATOM,NATOMNEW,NDNEV,NEED_ES,NG,NM,NRADMAX,NSORT,
     &        NSORTES,NT,W(MAXDIM)
      CHARACTER*4 TXTT(NTMAX)
C
C*** End of declarations rewritten by SPAG
C
C
C COMMON variables
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
C ==================================================
C
C  parameters which can be read from input
C
      R0ACT = 1.D-2                  ! - basic translation vectors
      NRADMAX = NRMAX - 1
      DPAS = 0.05D0                     !exp. pass
      NDNEV = 14                        !interpolation
C
      N1ES = 24                         !cell division
      N2ES = 24
      N3ES = 24
C
C      rmines=1.2D0              !minimal ES
C      rmaxes=2.5D0              !maximal ES
C      need_es=0                 !0 - activate search of es
C
C ==================================================
C
C       READING INPUT DATA:
C
      NEED_ES = 0
      DB = VERBOSE .ne. 0
      RMINES = RMINES_
      RMAXES = RMAXES_
      AV = CELL(:, 1)
      BV = CELL(:, 2)
      CV = CELL(:, 3)
      NT = NSORT
      IPRINT = VERBOSE
c      DO I = 1,NSORT
c         READ (INI,*) IDUM,Z(I),TXTT(I),IDUM,CONC(I),IMT(I)
c      END DO
      IF ( DB ) THEN
         PRINT *,'NATOM:',NQ

         PRINT *,'PARAMETERS:',NT,NQ,NM,'NT,NQ,NM,-- ALAT:',ALAT
         PRINT *,'Z(i),TXTT(i),CONC(i),IMT(i)'
         DO I = 1,NT
            PRINT *,Z(I),TXTT(I),CONC(I),IMT(I)
         END DO
         PRINT *,'qx(i),qy(i),qz(i),IMQ(i)'
         DO I = 1,NQ
            PRINT *,(BAS(J,I),J=1,3),IMQ(I)
         END DO
C
         PRINT *,'Lattice vectors:'
         PRINT *,AV
         PRINT *,BV
         PRINT *,CV
      END IF
C
C
      IF ( 1.EQ.2 ) WRITE (*,*) IDUM
C     ==================================================
C
C
      CALL WKINIT(MAXDIM)
C
      CALL GETGBASIS
C
C
C amtc,dpas,r0,rom,rop,r0,rop,rom,      dummy initialisation
C
      R0 = 9999D0
      AMTC = 9999D0
      RHOSITE(:,:) = 9999D0
C
      CALL READRO(NT,AMTC,Z,R0,TXTT,RHOSITE,NRADMAX,NRAD1)
C
C Producing new radial mesh:
C
      DO IT = 1,NT
         DO I = 1,NRAD1
            RAT(I,IT) = R0(IT)*EXP(DBLE(I-1)*DPAS)
         END DO
      END DO
C
      DO I = 1,NM
         WSREST(I) = 0
         ZZ(I) = 0
         R0SITE(I) = R0ACT
      END DO
C
      RHOSITE(:,:) = 0.0D0
C
      DO ITYPE = 1,NT
C
C this is cycle over atomic types
C IMT(itype) - number of inequivalent lattice position, which is
C occupied by this atom
C
         ISORT = IMQ(IMT(ITYPE))
         ZN = Z(ITYPE)
         CALL DEFWSR(WS,ZN)
         IF ( DB ) PRINT *,'type:',ITYPE,' Z=',ZN,' Rad=',WS,' conc=',
     &                   CONC(ITYPE)
         WSREST(ISORT) = WSREST(ISORT) + WS*CONC(ITYPE)
         ZZ(ISORT) = ZZ(ISORT) + ZN*CONC(ITYPE)
         DO IPNT = 1,NRAD1
            RAD = R0ACT*EXP(DBLE(IPNT-1)*DPAS)
            IF ( RAD.LT.RAT(NRAD1-2,ITYPE) ) THEN
               RINT = MAX(0.D0,
     &                DNEVMOD(RAT(1,ITYPE),AMTC(1,ITYPE),RAD,NRAD1,
     &                NDNEV))
            ELSE
               RINT = 0
            END IF
            RHOSITE(IPNT,ISORT) = RHOSITE(IPNT,ISORT) + CONC(ITYPE)*RINT
         END DO
      END DO
C
      DO I = 1,NM
         IF ( DB ) PRINT *,' Site:',I,'  RaD:',WSREST(I),' Z:',ZZ(I)
      END DO
C
C      stop
C
      NSORT = NM
      NATOM = NQ
C
C
      IF ( DB ) THEN
         PRINT *,'av=',AV
         PRINT *,'bv=',BV
         PRINT *,'cv=',CV
         PRINT *,'imq',(IMQ(I),I=1,NATOM)
         PRINT *,'BAS'
         DO I = 1,NATOM
            PRINT *,(BAS(J,I),J=1,3)
         END DO
         DO IST = 1,NSORT
            PRINT *,'R0S',R0SITE(IST),' ZZ',ZZ(IST),WSREST(IST)
         END DO
         PRINT *,'dpas,alat,natom,nsort:',DPAS,ALAT,NATOM,NSORT
      END IF
C
      CALL DEFMTR(W,NRAD1,AV,BV,CV,IMQ,BAS,RHOSITE,R0SITE,DPAS,ALAT,
     &            NATOM,NSORT,ZZ,RWS,NEED_ES,NRAD1,WSREST)
C
      NATOMNEW = NATOM
      IF ( NEED_ES.EQ.0 ) THEN
         CALL EMPTYSPHERES(AV,BV,CV,BAS,IMQ,ALAT,NATOM,NSORT,RWS,
     &                     NSORTES,TAUES,SES,MAX2SORT,ITOP,IQA,NG,ISNEW,
     &                     BASNEW,NATOMNEW)
Cccccccccccccccccccccccccccccccccccccccccccccccccccc
      ELSE
         DO IATOM = 1,NATOM
            ISNEW(IATOM) = IMQ(IATOM)
            DO I = 1,3
               BASNEW(I,IATOM) = BAS(I,IATOM)
            END DO
         END DO
Ccccccccccccccccccccccccccccccccccccccccccccccccccc
      END IF
C
C
C ==================================================
C
C       WRITING OUTPUT DATA:
C
C
C     OPEN (22,FILE='radii.out')
C     WRITE (22,*) NATOMNEW
      IF (NATOMNEW > N_OUT) THEN
          N_OUT = -NATOMNEW
      ELSE
          N_OUT = NATOMNEW - NATOM

          DO I = 1,NATOM
             RADII(I+N_OUT) = RWS(IMQ(I))
          END DO
          DO I = 1,NATOMNEW-NATOM
             CENTRES(1:3,I) = BASNEW(1:3,I+NATOM)
             ISR = ISNEW(I+NATOM)
             RADII(I) = SES(ISR-NM)
          END DO
      END IF
      IF (DB ) THEN
        DO I = 1, NATOMNEW
          IF ( I.LE.NATOM ) THEN
            WRITE (*,99001) ISNEW(I),(BASNEW(J,I),J=1,3),RWS(IMQ(I)),
     &                       ' / A/'
          ELSE
            WRITE (*,99001) ISR,(BASNEW(J,I),J=1,3),SES(ISR-NM),' / E/'
          END IF
        END DO
      END IF
C     CLOSE (22)
      FLUSH(6)
99001 FORMAT (i4,2x,3F21.15,3x,f21.15,3x,a)

      END SUBROUTINE

C*==readro.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE READRO(NSORT,RO,Z,R0,TXT,ROM,NRADMAX,NRAD1)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NRAD1,NRADMAX,NSORT
      DOUBLE PRECISION R0(*),RO(NRADMAX+1,*),ROM(NRAD1,*),Z(*)
      CHARACTER*4 TXT(*)
      COMMON /IPRINT/ IPRINT
      INTEGER IPRINT
C
C Local variables
C
      INTEGER I,ISORT,ISR,IZC,IZN,IZNUC
C
C*** End of declarations rewritten by SPAG
C
      OPEN (8,STATUS='scratch')
      DO ISORT = 1,NSORT
         IF ( Z(ISORT).LT.0 ) THEN
            IZNUC = -10
            CALL QZZC(TXT(ISORT),IZNUC,IZC)
            Z(ISORT) = IZNUC
         END IF
         IF ( Z(ISORT).LT.0.3D0 ) THEN
C emty sphere
            DO I = 1,NRAD1
               RO(I,ISORT) = 0.D0
               ROM(I,ISORT) = 0.D0
               R0(ISORT) = 1D-5
            END DO
         ELSE
            DO ISR = 1,ISORT - 1
               IF ( NINT(Z(ISR)).EQ.NINT(Z(ISORT)) ) THEN
                  DO I = 1,NRAD1
                     RO(I,ISORT) = RO(I,ISR)
                     ROM(I,ISORT) = ROM(I,ISR)
                  END DO
                  R0(ISORT) = R0(ISR)
                  GOTO 100
               END IF
            END DO
C          write(buf,2)nint(z(isort)),nint(zc)
            if (IPRINT > 0)
     &      WRITE (6,'(/'' For atom '',a4,''    Z='',f4.0/)') TXT(ISORT)
     &             ,Z(ISORT)
            IZN = NINT(Z(ISORT))
            IZC = 0
            CALL RHFDS(IZN,IZC,RO(1,ISORT),ROM(1,ISORT),R0(ISORT))
C
         END IF
 100  END DO
      CLOSE (8)
      END
C*==qzzc.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE QZZC(TXT,IZ,IZC)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IZ,IZC
      CHARACTER*2 TXT
C
C Local variables
C
      CHARACTER*210 ELEMENTS
      INTEGER I,IZCOR(105)
      CHARACTER*2 TX
C
C*** End of declarations rewritten by SPAG
C
      DATA IZCOR/0,0,0,6*2,4,4,6*10,12,12,12*18,4*28,30,30,12*36,4*46,
     &     48,48,16*54,10*68,4*78,80,80,18*86/
      ELEMENTS = 
     &      'E H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNi'
     &      //
     &      'CuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeCsBaLa'
     &      //
     &      'CePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn'
     &      //'FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLwKu'
      IF ( IZ.GE.0 ) THEN
         IZC = IZCOR(IZ+1)
         RETURN
      END IF
      TX = TXT
      IF ( TX(2:2).EQ.'_' ) TX(2:2) = ' '
      DO I = 1,105
C        print 1,elements(i*2-1:i*2),i-1,izcor(i)
         IF ( ELEMENTS(I*2-1:I*2).EQ.TXT(1:2) ) THEN
            IZ = I - 1
            IZC = IZCOR(I)
            RETURN
         END IF
      END DO
      WRITE (16,*) TXT
      STOP 'Error : Such Chemical element is ABSENT'
      END
C*==shortn.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE SHORTN(P,P1)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1.D-7)
C
C COMMON variables
C
      DOUBLE PRECISION PLAT(3,3),QLAT(3,3)
      COMMON /L2LAT / PLAT,QLAT
C
C Dummy arguments
C
      DOUBLE PRECISION P(3),P1(3)
C
C Local variables
C
      INTEGER I
      DOUBLE PRECISION X(3)
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
C
      DO I = 1,3
         X(I) = P(1)*QLAT(1,I) + P(2)*QLAT(2,I) + P(3)*QLAT(3,I)
         X(I) = X(I) - NINT(X(I)-EPS)
C$$$        if(x(i).lt.-0.5d0+eps)x(i)=x(i)+1.d0
      END DO
      DO I = 1,3
         P1(I) = X(1)*PLAT(I,1) + X(2)*PLAT(I,2) + X(3)*PLAT(I,3)
      END DO
      END
C*==vec.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE VEC(V0,V1,V2,G1)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION G1(3,3),V0(3),V1(3),V2(3)
C
C Local variables
C
      INTEGER I,J
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,3
         V0(I) = V1(I)
         DO J = 1,3
            V0(I) = V0(I) + V2(J)*G1(J,I)
         END DO
      END DO
      CALL SHORTN(V0,V0)
      END
C*==getgbasis.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE GETGBASIS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C COMMON variables
C
      DOUBLE PRECISION AVC(3),BVC(3),CVC(3),FVC(3),HVC(3),PVC(3)
      COMMON /L2LAT / AVC,BVC,CVC,FVC,PVC,HVC
C
C Local variables
C
      INTEGER I
      DOUBLE PRECISION VOLCELL,W1
C
C*** End of declarations rewritten by SPAG
C
C
C COMMON variables
C
C
C Local variables
C
C
C*** End of declarations rewritten by SPAG
C
      CALL CROSS(FVC,BVC,CVC)
      CALL CROSS(PVC,CVC,AVC)
      CALL CROSS(HVC,AVC,BVC)
      W1 = AVC(1)*FVC(1) + AVC(2)*FVC(2) + AVC(3)*FVC(3)
      VOLCELL = ABS(W1)
C
      IF ( VOLCELL.LT.1.D-5 ) THEN
         PRINT *,AVC
         PRINT *,BVC
         PRINT *,CVC
         CALL ENDJOB(10,'GETGBASIS: Lattice vectors are complanar!')
      END IF
      DO I = 1,3
         FVC(I) = FVC(I)/W1
         PVC(I) = PVC(I)/W1
         HVC(I) = HVC(I)/W1
      END DO
      END
