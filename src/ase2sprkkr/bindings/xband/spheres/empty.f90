!*==emptyspheres.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
INTEGER FUNCTION EMPTYSPHERES(PLAT, BAS, IS, ALAT, NATOM, NSORT, &
  &                        S, NSORTES, TAUES, SES, MAX2SORT, &
  &                        ISNEW, BASNEW, NATOMNEW, &
  &                        N_SYMMETRY_OPS, ROTATIONS, &
  &                        TRANSLATIONS)
    USE DEBUG_DUMP
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
    INTEGER NPATB1
    PARAMETER(NPATB1=6000)
!
! Dummy arguments
!
    DOUBLE PRECISION ALAT
    INTEGER MAX2SORT, NATOM, NATOMNEW, NSORT, NSORTES
    INTEGER N_SYMMETRY_OPS
    DOUBLE PRECISION BAS(3, MAX2SORT), BASNEW(3, MAX2SORT), PLAT(3, 3), &
      &       S(MAX2SORT), &
      &       SES(MAX2SORT), &
      &       TAUES(3, MAX2SORT)
    INTEGER IS(MAX2SORT), ISNEW(MAX2SORT)
    DOUBLE PRECISION, DIMENSION(3, 3, N_SYMMETRY_OPS) &
      &       :: ROTATIONS
    DOUBLE PRECISION, DIMENSION(3, N_SYMMETRY_OPS) &
      &       :: TRANSLATIONS
!
! Local variables
!
    DOUBLE PRECISION DIFF(3), S2NEW(:), SNEW(:), TAU(:, :)
    INTEGER I, IER, IG, IGA, IGC, ISB(:), ISKIP(:), I_ISB, I_ISKIP, I_NHSORTES, &
      &        I_S2NEW, I_SNEW, I_TAU, J, MAXNPAT, MAXNPATG, NHSORTES(:)

    COMMON/IPRINT/IPRINT
    INTEGER IPRINT
    INTEGER EMPTY
    EXTERNAL EMPTY
!
!*** End of declarations rewritten by SPAG
!
    ALLOCATABLE SNEW, S2NEW, TAU, NHSORTES, ISKIP, ISB
!
    MAXNPAT = MAX2SORT
    MAXNPATG = NPATB1
!
    ALLOCATE (SNEW(MAXNPAT), S2NEW(MAXNPAT), TAU(3, MAXNPATG))
    ALLOCATE (NHSORTES(MAXNPAT), ISKIP(MAXNPAT), ISB(MAXNPATG))
!
    EMPTYSPHERES = EMPTY(IER, PLAT, BAS, IS, ALAT, NATOM, NSORT, S, &
      &           BASNEW, ISNEW, &
      &           SNEW, S2NEW, ISB, TAU, MAXNPAT, MAXNPATG, NSORTES, TAUES, SES, &
      &           NHSORTES, ISKIP, ROTATIONS, TRANSLATIONS, N_SYMMETRY_OPS, &
      &           NATOMNEW)
!
!     nsortes=iter-1
    IF (IPRINT > 0) THEN
        WRITE (6, *) ' IRR. POSITIONS:'
        DO I = 1, NSORTES
            WRITE (6, '(3f21.15)') (TAUES(J, I), J=1, 3)
        END DO
    END IF
!      stop
    DEALLOCATE (SNEW, S2NEW, TAU)
    DEALLOCATE (NHSORTES, ISKIP, ISB)

END
!*==empty.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
INTEGER FUNCTION EMPTY(IER, PPLAT, BAS, IS, ALAT, NATOM, &
  &                 NSORT, S, BASNEW, &
  &                 ISNEW, SNEW, S2NEW, ISB, TAU, MAXNPAT, MAXNPATG, &
  &                 NSORTES, TAUES, SES, NHSORTES, ISKIP, GM, VC, NG, &
  &                 NATOMNEW)

    USE DEBUG_DUMP
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! COMMON variables
!
    DOUBLE PRECISION BOA, COA, RMAXES, RMINES, U(3, 3), UM1(3, 3)
    INTEGER IDATASTR, N1ES, N2ES, N3ES, NATES
    COMMON/EMPTYS/NATES, N1ES, N2ES, N3ES, RMINES, RMAXES
    COMMON/TAU2/BOA, COA
    COMMON/U/U, UM1, IDATASTR
    DOUBLE PRECISION PLAT(3, 3), QLAT(3, 3)
    COMMON/L2LAT/PLAT, QLAT
!
! Dummy arguments
!
    DOUBLE PRECISION ALAT
    INTEGER IER, MAXNPAT, MAXNPATG, NATOM, NATOMNEW, NG, NSORT, NSORTES
    DOUBLE PRECISION PPLAT(3, 3), BAS(3, *), BASNEW(3, *), GM(3, 3, NG), S(*), &
      &         S2NEW(*), SES(*), SNEW(*), TAU(3, *), TAUES(3, *), VC(3, NG)
    INTEGER IS(*), ISB(*), ISKIP(*), ISNEW(*), NHSORTES(*)
!
! Local variables
!
    DOUBLE PRECISION DIFF(3)
    DOUBLE PRECISION CS, CVOL, DEMIN, DEMAX, EPS_DIFF, EPS_DIFFMAX, PI43, RAD
    DOUBLE PRECISION RADNEW, SM, SMN, SSM, SVOLA, SVOLE
    DOUBLE PRECISION VE(3, 48), VSHIFT(3), VV(3), TEMPV(3), CORRS(3)
    DOUBLE PRECISION P(3), PSH(3)
    DOUBLE PRECISION DPI, TRNT
    INTEGER I, IATOM, IG, IIN, IN, INVMAX, IODN, ISORT, ITAV, ITER, IUN, J, K, &
      &        NATB, NATS, NSORTNEW

    INTEGER IPRINT
    COMMON/IPRINT/IPRINT
    EMPTY = 0       ! default to success

    PLAT(:, :) = PPLAT(:, :)
! just wrap the atoms to -0.5, 0.5
   
    basnew(:,1:NATOM) = MATMUL(transpose(qlat), bas(:,1:NATOM))
    basnew(:,1:NATOM) = basnew(:,1:NATOM) - DNINT(basnew(:,1:NATOM) - 1.0D-7)
    bas(:,1:NATOM) = MATMUL(plat, basnew(:,1:NATOM))  

    VSHIFT(:) = 0

    BOA = 1
    COA = 1
!---------------------------------------
    IUN = 6
    NATES = 0
    SSM = 0
    DO I = 1, NSORT
        SNEW(I) = S(I)/ALAT
        SSM = MAX(SSM, SNEW(I))
        S2NEW(I) = SNEW(I)**2
    END DO
    INVMAX = 0
    IF (RMAXES .LT. 0) INVMAX = 1
    IF (RMINES .LT. 0) INVMAX = -INVMAX
    RMINES = ABS(RMINES)
    RMAXES = ABS(RMAXES)
    SM = RMAXES/ALAT + SSM
    SMN = RMINES/ALAT
    DO IATOM = 1, NATOM
        ISNEW(IATOM) = IS(IATOM)
        DO I = 1, 3
            BASNEW(I, IATOM) = BAS(I, IATOM)
        END DO
    END DO
    NATOMNEW = NATOM
    NSORTNEW = NSORT
    ITER = 0
    IF (IPRINT > 0) THEN
        PRINT *, 'BOA,COA:', BOA, COA, ALAT
        WRITE (IUN, *) '  Searching of Empty spheres  '
        WRITE (*, *) '  Searching of Empty spheres  '
        WRITE (IUN, '(a,3i3)') ' Mesh: ', N1ES, N2ES, N3ES
        WRITE (IUN, 99001) 'Smin, Smax (a.u.): ', RMINES, RMAXES
    END IF
100 CONTINUE
    ITER = ITER + 1
    IF (IPRINT > 0) THEN
        WRITE (IUN, *) '  Iteration ', ITER
        WRITE (*, *) '  Iteration ', ITER
    END IF
!
    CALL GETEPOS(MAXNPATG, PLAT, BASNEW, NATOMNEW, ISNEW, TAU, NATB, &
      &             ISB, SNEW, S2NEW, N1ES, N2ES, N3ES, P, RAD, IER, SM)
!
    IF (IPRINT > 0) THEN
        WRITE (IUN, 99001) ' Maximum possible sphere is ', RAD*ALAT
        WRITE (IUN, 99001) ' Position in a,b,c:', P(1), P(2), P(3)
    END IF
    IF (RAD .LT. SMN) THEN
        IF (IPRINT > 0) THEN
            WRITE (IUN, *) 'skipped'
            PRINT *, 'RAD=', RAD*ALAT, ' RMINES=', RMINES
        END IF
        EMPTY = IER
        GOTO 300
    END IF
!
    IF (IER == 0) THEN
        ISKIP(ITER) = 0
        IF (RAD .GT. RMAXES/ALAT) THEN
            IF (INVMAX .EQ. 0) THEN
                IF (IPRINT > 0) WRITE (IUN, 99001) 'Contracted to be ', RMAXES
                RAD = RMAXES/ALAT        ! rad is used later
            ELSE
                ISKIP(ITER) = INVMAX
                IF (IPRINT > 0) THEN
                    WRITE (IUN, 99001) 'Sphere is too large', RAD*ALAT
                    WRITE (IUN, 99001) &
                      &                      'Small auxiliary sphere is inserted instead'
                END IF
            END IF
        END IF
        !
        ! multiplication :
        !        Convert VSHIFT (Cartesian) to fractional shift via QLAT, then apply to P
        PSH(:) = MATMUL(VSHIFT, QLAT)
        P(:) = P(:) - PSH(:)
        P(:) = P(:) - DNINT(P(:) -1D-7)
        !         Build VE from shifted fractional P (no direct VSHIFT subtraction)
        VE(:, 1) = MATMUL(PLAT, P)
        ITAV = 0

150     CONTINUE

        ITAV = ITAV + 1
        IN = 1

        CORRS(:) = 0
        IODN = 0
        DEMIN = 1.D10
        DEMAX = 0D0
        
        groups: DO IG = 1, NG
            CALL VEC(TEMPV, VC(1, IG), P, GM(1, 1, IG))
            VV(:) = MATMUL(PLAT, TEMPV)

            DO IIN = IN, 1, -1
                DIFF(:) = VV(:) - VE(:, IIN)
                CALL SHORTN(DIFF, DIFF)
                EPS_DIFF = sum(DIFF**2)
                IF (EPS_DIFF .LT. 1.D-5) THEN  !The same point
                    CYCLE groups
                END IF
            END DO

            IN = IN + 1
            VE(:, IN) = VV(:)

            EPS_DIFF = SQRT(EPS_DIFF)
            DEMIN = MIN(EPS_DIFF, DEMIN)
            IF (EPS_DIFF .LT. RAD*1.5D0) THEN
                CORRS(:) = CORRS(:) + DIFF(:)
                IODN = IODN + 1
            ELSE
                DEMAX = MAX(DEMAX, EPS_DIFF)
            END IF
        END DO groups

        IF (IODN .NE. 0) THEN
            IF (IPRINT > 0) THEN
                WRITE (IUN, '(a,i3,f21.15, 3f21.15)') &
                  &             'a bit shifted position' &
                  &             , IODN, REAL(RAD), CORRS
                WRITE (IUN, 99001) 'OLD:', VE(1, 1), VE(2, 1)/BOA, VE(3, 1)/COA
            END IF

            VE(:, 1) = VE(:, 1) + CORRS/(IODN + 1)
            !           Update fractional coordinates P from finalized VE(:,1)
            P(:) = MATMUL(VE(:, 1), QLAT)
            P(:) = P(:) - DNINT(P(:) -1D-7)
            VE(:, 1) = MATMUL(PLAT, P)
            
            IF (IPRINT > 0) WRITE (IUN, 99001) &
              &              'new:', VE(1, 1), VE(2, 1)/BOA, VE(3, 1)/COA
            CALL REGET_RE(TAU, ISB, SNEW, NATB, VSHIFT, VE(1, 1), RADNEW)

            IF (IPRINT > 0) WRITE (IUN, 99001) &
              &              'rad,radtst', RAD*ALAT, RADNEW*ALAT
            RAD = RADNEW
            IF (ITAV .LE. 5) GOTO 150   ! to avoide infinite cycle
            IF (IPRINT > 0) WRITE (IUN, *) 'can not find averaged position'
            EMPTY = -1
            RETURN
        END IF

        IF (DEMIN .LT. 2.D0*RAD) THEN
            ! overlapping spheres
            IF (IPRINT > 0) &
              &          WRITE (IUN, *) 'radius is decreased from', RAD, ' to', &
              &                    DEMIN/2.D0
            RAD = DEMIN/2.D0            ! decrease radius to touching
            IF (RAD .LT. SMN) THEN
                IF (IPRINT > 0) THEN
                    WRITE (IUN, *) 'sphere is less then rmines'
                    WRITE (IUN, *) 'rad=', RAD*ALAT, ' rmines=', RMINES
                END IF
                ISKIP(ITER) = -1         ! will be skipped later
            END IF
        END IF
        IF (IPRINT > 0) THEN
            WRITE (IUN, '(a,i4,a)') 'it will be ', IN, ' ES per unit cell:'
        END IF
        IF (NATOMNEW + IN .GT. MAXNPAT) THEN
            PRINT *, 'EMPTY: too many atoms', NATOMNEW + IN, ' maxnpat=', &
              &            MAXNPAT
            PRINT *, 'try to increase array sizes in empty.f'
            EMPTY = -2
            RETURN
        END IF
        DO IATOM = NATOMNEW + 1, NATOMNEW + IN
            ISNEW(IATOM) = NSORTNEW + 1
            DO I = 1, 3
                BASNEW(I, IATOM) = VE(I, IATOM - NATOMNEW) + VSHIFT(I)
            END DO
            IF (IPRINT > 0) &
              &          WRITE (IUN, '(2x,3f15.8)') (BASNEW(I, IATOM), I=1, 3)
        END DO
        !
        TAUES(:, ITER) = VE(:, 1)
        !      write(iun,'(i3,3f15.8)')iter,vv
        NSORTNEW = NSORTNEW + 1
        NATOMNEW = NATOMNEW + IN
        NHSORTES(ITER) = IN
        IF (ISKIP(ITER) .EQ. 0) THEN
            SNEW(NSORTNEW) = MIN(RAD, RMAXES/ALAT)
        ELSE
            SNEW(NSORTNEW) = RMINES/ALAT
        END IF
        S2NEW(NSORTNEW) = SNEW(NSORTNEW)**2
        SES(ITER) = SNEW(NSORTNEW)*ALAT
        IF (IPRINT > 0) &
          &      PRINT *, IN, ' Empty spheres added, S=', SNEW(NSORTNEW)*ALAT
        NATES = NATES + IN
        GOTO 100
    END IF
!
300 CONTINUE
    NSORTES = ITER - 1
    IF (INVMAX .GE. 0) THEN
        ! check if some spheres should be skipped due to small radius
        DO ISORT = 1, NSORTES
            IF (ISKIP(ISORT) .LT. 0) INVMAX = ISKIP(ISORT)
        END DO
    END IF
    IF (INVMAX .LT. 0) THEN
        ! skip auxiliary spheres
        !$$$        write(iun,*)'nsortes,natomnew',nsortes,natomnew
        !$$$        do isort=1,nsortes
        !$$$          write(iun,'(i3,3f21.15,2i4,f21.15)')isort,(taues(k,isort),k=1,3)
        !$$$     $         ,nhsortes(isort),iskip(isort),ses(isort)
        !$$$        enddo
        !$$$        do iatom=1,natomnew
        !$$$          write(iun,'(i3,3f21.15)')iatom,(basnew(k,iatom),k=1,3)
        !$$$        enddo
        !$$$        write(iun,*)'rmines',rmines
        IATOM = 1
        ISORT = 1
        DO WHILE (ISORT .LE. NSORTES)
            IF (IPRINT > 0) &
              &         WRITE (IUN, *) 'isort,nsortes,s', ISORT, NSORTES, SES(ISORT)
            NATS = NHSORTES(ISORT)
            IF (ISKIP(ISORT) .LT. 0) THEN
                IF (IPRINT > 0) &
                  &            WRITE (IUN, 99001) 'sphere is skipped', &
                  &                           (TAUES(K, ISORT), K=1, 3)
                DO I = ISORT + 1, NSORTES
                    NHSORTES(I - 1) = NHSORTES(I)
                    SES(I - 1) = SES(I)
                    ISKIP(I - 1) = ISKIP(I)
                    DO K = 1, 3
                        TAUES(K, I - 1) = TAUES(K, I)
                    END DO
                END DO
                DO I = IATOM + NATS, NATOMNEW
                    DO K = 1, 3
                        BASNEW(K, I - NATS) = BASNEW(K, I)
                    END DO
                END DO
                NSORTES = NSORTES - 1
                NATOMNEW = NATOMNEW - NATS
            ELSE
                IATOM = IATOM + NATS
                ISORT = ISORT + 1
            END IF
        END DO
        !$$$        write(iun,*)'nsortes,natomnew',nsortes,natomnew
        !$$$        do isort=1,nsortes
        !$$$          write(iun,'(i3,3f21.15,i4,f21.15)')isort,(taues(k,isort),k=1,3)
        !$$$     $         ,nhsortes(isort),ses(isort)
        !$$$        enddo
        !$$$        do iatom=1,natomnew
        !$$$          write(iun,'(i3,3f21.15)')iatom,(basnew(k,iatom),k=1,3)
        !$$$        enddo
    END IF
!$$$      if(nsortes.gt.0)then
!$$$        open(24,status='scratch',form='unformatted')
!$$$        do iatom=natom+1,natomnew
!$$$          write(24)(basnew(i,iatom),i=1,3)
!$$$c$$$        write(iun,'(2x,3f15.8,i4)')(basnew(i,iatom),i=1,3),isnew(iatom)
!$$$        enddo
!$$$      endif
! MT -> ASA spheres
    PI43 = 4.D0*DPI()/3.D0
    SVOLA = 0.D0
    DO IATOM = 1, NATOM
        SVOLA = SVOLA + PI43*S(IS(IATOM))**3
    END DO
    SVOLE = 0.D0
    DO ISORT = 1, NSORTES
        SVOLE = SVOLE + NHSORTES(ISORT)*PI43*SES(ISORT)**3
        !        print*,'SES::',ses(isort),isort
    END DO
    CVOL = TRNT(PLAT(:, 1), PLAT(:, 2), PLAT(:, 3))*ALAT**3
    CS = CVOL/(SVOLA + SVOLE)
    SVOLA = SVOLA*CS
    SVOLE = SVOLE*CS
    CS = CS**(1.D0/3.D0)
    DO ISORT = 1, NSORT
        S(ISORT) = CS*S(ISORT)
    END DO
    DO ISORT = 1, NSORTES
        SES(ISORT) = CS*SES(ISORT)
    END DO
    IF (IPRINT > 0) THEN
        WRITE (IUN, 99001) 'Vat, Ves, V', SVOLA, SVOLE, SVOLA + SVOLE
        WRITE (IUN, 99001) 'Vat/V, Ves/V', SVOLA/CVOL, SVOLE/CVOL
    END IF
99001 FORMAT(1x, a, 3F21.15)
END
!*==getepos.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
SUBROUTINE GETEPOS(NPAT, PLAT, BAS, NATOM, IS, TAU, NATB, ISB, S, S2, N1, &
    &                   N2, N3, RES, RAD, IER, SM)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER IER, N1, N2, N3, NATB, NATOM, NPAT
    DOUBLE PRECISION RAD, SM
    DOUBLE PRECISION RES(3)
    DOUBLE PRECISION PLAT(3, 3), BAS(3, *), S(*), S2(*), TAU(3, NPAT),SHIFT(3)
    INTEGER IS(NPAT), ISB(NPAT)
    
    COMMON/IPRINT/IPRINT
    INTEGER IPRINT   
    !
    ! Local variables
    !
    DOUBLE PRECISION ANA, ANB, ANC, DD, DDMAX, DIAG, DIST2, DMAX, DMIN, D(3), P(3), &
      &       XYZ(3), V(3)
    INTEGER I, IA, IATOM, IB, IC, J, K
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    DIAG = maxval([ &
                  sum((PLAT(:, 1) + PLAT(:, 2) + PLAT(:, 3))**2), &
                  sum((PLAT(:, 1) - PLAT(:, 2) + PLAT(:, 3))**2), &
                  sum((PLAT(:, 1) + PLAT(:, 2) - PLAT(:, 3))**2), &
                  sum((PLAT(:, 1) - PLAT(:, 2) - PLAT(:, 3))**2) &
                  ])

    !
    IER = 0
    DIAG = SQRT(DIAG)/2
    DDMAX = (SM + DIAG)**2
    TAU = 0.0D0
    !
    NATB = 0

    DO I = -3, 3
        V(1) = I
        DO J = -3, 3
            V(2) = J
            DO K = -3, 3
                V(3) = K
                D = matmul(PLAT, V)
                DO IATOM = 1, NATOM
                    XYZ = D + BAS(:, IATOM)
                    DD = dot_product(XYZ, XYZ)
                    IF (DD .LT. DDMAX) THEN
                        NATB = NATB + 1
                        TAU(:, NATB) = XYZ
                        IF (NATB .GT. NPAT) THEN
                            IER = -3
                            IF ( IPRINT > 0 ) THEN
                                PRINT *, 'IER=', IER
                                PRINT *, &
                                &                        'GETEPOS: too many atoms are generated ', &
                                &                        NATB, ' npat=', NPAT
                                PRINT *, 'try to increase npatb1 in empty.f'
                            END IF
                            RETURN
                        END IF
                        ISB(NATB) = IS(IATOM)
                    END IF
                END DO                   ! iatom
            END DO                      ! k
        END DO                         ! j
    END DO                            ! i
    !
    P(:) = 0
    !
    ANA = 0.5D0/N1
    ANB = 0.5D0/N2
    ANC = 0.5D0/N3
    DMAX = 0
    !      npnt=0

    SHIFT(:) = MATMUL(PLAT, (/0.D0,0.D0,ANC/))
    
    DO IA = -N1, N1 - 1, 2
        P(1) = ANA*IA
        DO IB = -N2, N2 - 1, 2
            P(2) = ANB*IB
            P(3) = -N3 * ANC
            XYZ = matmul(PLAT, P)

            DO IC = -N3, N3 - 1, 2
                DMIN = 1.D10
                DO I = 1, NATB
                    DIST2 = sum((XYZ - TAU(:, I))**2)
                    IF (DIST2 .LT. S2(ISB(I))) GOTO 20
                    DD = SQRT(DIST2) - S(ISB(I))
                    IF (DD .LT. DMAX) GOTO 20
                    DMIN = MIN(DD, DMIN)                    
                END DO                   ! i
                IF (DMIN .GT. DMAX) THEN
                    DMAX = DMIN
                    RES(:) = P(:)
                END IF
                XYZ = XYZ + SHIFT
20          END DO                        ! ic
        END DO                          ! ib
    END DO                            ! ia
    !      print*,iter,'MAX is:',DMAX,iae,ibe,ice,pae,pbe,pce
    !      print*,iter,'npnt=',npnt
    RAD = DMAX
END
!*==reget_re.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
SUBROUTINE REGET_RE(TAU, ISB, S, NATB, VSHIFT, VS, RAD)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER NATB
    DOUBLE PRECISION RAD
    INTEGER ISB(*)
    DOUBLE PRECISION S(*), TAU(3, *), VS(*), VSHIFT(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION D, DMIN, V(3)
    INTEGER I, K
    !
    !*** End of declarations rewritten by SPAG
    !
    DO K = 1, 3
        V(K) = VS(K) + VSHIFT(K)
    END DO
    DMIN = 1.D10
    DO I = 1, NATB
        D = (V(1) - TAU(1, I))**2 + (V(2) - TAU(2, I))**2 + (V(3) - TAU(3, I)) &
           &       **2
        D = SQRT(D) - S(ISB(I))
        IF (D .LT. 0.D0) THEN
            PRINT *, 'REGET_RE: error'
            PRINT *, I, D + S(ISB(I)), ISB(I), S(ISB(I))
            PRINT *, (V(K), K=1, 3)
            PRINT *, TAU(1, I), TAU(2, I), TAU(3, I)
            RETURN
        END IF
        DMIN = MIN(D, DMIN)
    END DO                            ! natb
    RAD = DMIN
END
