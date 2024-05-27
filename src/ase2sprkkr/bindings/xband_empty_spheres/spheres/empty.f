C*==emptyspheres.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
      SUBROUTINE EMPTYSPHERES(AVC,BVC,CVC,BAS,IS,ALAT,NATOM,NSORT,S,
     &                        NSORTES,TAUES,SES,MAX2SORT,ITOP,IQA,NG,
     &                        ISNEW,BASNEW,NATOMNEW)
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
      REAL*8 AVC(3),BAS(3,*),BASNEW(*),BVC(3),CVC(3),S(*),SES(*),
     &       TAUES(3,*)
      INTEGER IQA(*),IS(*),ISNEW(*),ITOP(*)
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
      DO IGC = 1,NG
         IG = ITOP(IGC)
         IGA = IG
         IF ( IGA.GT.32 ) IGA = IGA - 32
         CALL EILANG(ANG,IGA)
         CALL TURNM(ANG,GEN(1,1,IGC))
C                              !DEF TURN MATRIX
         IF ( IG.NE.IGA ) THEN
C                              !ADD INVERSION IF NEED
            DO I = 1,3
               DO J = 1,3
                  GEN(I,J,IGC) = -GEN(I,J,IGC)
               END DO
            END DO
         END IF
         CALL VEC(VC(1,IGC),BAS(1,IQA(IGC)),VV,GEN(1,1,IGC))
C        print*,Igc,ig,(vc(i,igc),i=1,3)
      END DO
C
      DO I = 1,NG
         IG = ITOP(I)
C
      END DO
C
      CALL DEFRR(I_SNEW,MAXNPAT)
      CALL DEFRR(I_S2NEW,MAXNPAT)
      CALL DEFI(I_ISB,MAXNPATG)
      CALL DEFRR(I_TAU,MAXNPATG*3)
      CALL DEFI(I_NHSORTES,MAXNPAT)
      CALL DEFI(I_ISKIP,MAXNPAT)
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
C              print*,X,Y,Z
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
C            print*,DMIN
               IF ( DMIN.GT.DMAX ) THEN
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
