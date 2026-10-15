INTEGER FUNCTION FIND_EMPTY_SPHERES( &
  &   N_OUT, &
  &   CENTRES, &
  &   RADII, &
  &   RMINES_, &
  &   RMAXES_, &
  &   ALAT, &
  &   CELL, &
  &   NQ, &
  &   BAS, &
  &   IMQ, &
  &   NM, &
  &   NSORT, &
  &   TXTT, &
  &   Z, &
  &   CONC, &
  &   IMT, &
  &   MESH, &
  &   N_SYMMETRY_OPS, &
  &   ROTATIONS, &
  &   TRANSLATIONS, &
  &   VERBOSE &
  & )
    IMPLICIT NONE
    INTEGER EMPTYSPHERES
    EXTERNAL EMPTYSPHERES

! PARAMETER definitions
!
    INTEGER NTMAX, NRMAX, NRAD1, MAXDIM, MAX2SORT
    PARAMETER(NTMAX=200, NRMAX=300, NRAD1=251, MAXDIM=1000000, &
      &           MAX2SORT=NTMAX)

    INTEGER N_OUT
    DOUBLE PRECISION CENTRES(3, MAXDIM)
    DOUBLE PRECISION RADII(MAXDIM)

    DOUBLE PRECISION RMINES_, RMAXES_
    DOUBLE PRECISION CELL(3, 3)
    INTEGER NQ
    DOUBLE PRECISION BAS(3, NTMAX)
    INTEGER MESH(3)
    INTEGER N_SYMMETRY_OPS
    DOUBLE PRECISION, DIMENSION(3, 3, N_SYMMETRY_OPS) &
      &       :: ROTATIONS
    DOUBLE PRECISION, DIMENSION(3, N_SYMMETRY_OPS) &
      &       :: TRANSLATIONS
    INTEGER VERBOSE
    INTEGER IPRINT
    COMMON/IPRINT/IPRINT

!*==aa0001.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
!*** Start of declarations rewritten by SPAG
!
!
! COMMON variables
!
    DOUBLE PRECISION PLAT(3, 3), QLAT(3, 3), RMAXES, RMINES
    INTEGER N1ES, N2ES, N3ES, NATES
    COMMON/EMPTYS/NATES, N1ES, N2ES, N3ES, RMINES, RMAXES
    COMMON/L2LAT/PLAT, QLAT
!
! Local variables
!
    DOUBLE PRECISION ALAT, AMTC(NRMAX, NTMAX), &
      &                 BASNEW(3, MAX2SORT), CONC(NTMAX), DPAS, R0(NTMAX), &
      &                 R0ACT, R0SITE(MAX2SORT), RAD, RAT(NRAD1, MAX2SORT), &
      &                 RHOSITE(NRAD1, MAX2SORT), RINT, RWS(NTMAX), &
      &                 SES(MAX2SORT), TAUES(3, MAX2SORT), WS, &
      &                 WSREST(MAX2SORT), Z(NTMAX), ZN, ZZ(MAX2SORT)
    LOGICAL DB
    DOUBLE PRECISION DNEVMOD
    CHARACTER*100 FFF
    INTEGER I, IATOM, IDUM, IMQ(NQ), IMT(NQ), IPNT, &
      &        ISNEW(MAX2SORT), ISORT, ISP1, ISP2, ISR, IST, IT, ITYPE, &
      &        J, NATOM, NATOMNEW, NDNEV, NEED_ES, NG, NM, NRADMAX, NSORT, &
      &        NSORTES, NT, WORK(MAXDIM)
    CHARACTER*4 TXTT(NSORT)

    R0ACT = 1.D-2                  ! - basic translation vectors
    NRADMAX = NRMAX - 1
    DPAS = 0.05D0                     !exp. pass
    NDNEV = 14                        !interpolation
!
    N1ES = MESH(1)                    !cell division
    N2ES = MESH(2)
    N3ES = MESH(3)
!
!      rmines=1.2D0              !minimal ES
!      rmaxes=2.5D0              !maximal ES
!      need_es=0                 !0 - activate search of es
!
! ==================================================
!
!       READING INPUT DATA:
!
    NEED_ES = 0
    DB = VERBOSE .ne. 0
    RMINES = RMINES_
    RMAXES = RMAXES_
    NT = NSORT
    IPRINT = VERBOSE
!      DO I = 1,NSORT
!         READ (INI,*) IDUM,Z(I),TXTT(I),IDUM,CONC(I),IMT(I)
!      END DO
    IF (DB) THEN
        PRINT *, 'NATOM:', NQ

        PRINT *, 'PARAMETERS:', NT, NQ, NM, 'NT,NQ,NM,-- ALAT:', ALAT
        PRINT *, 'Z(i),TXTT(i),CONC(i),IMT(i)'
        DO I = 1, NT
            PRINT *, Z(I), TXTT(I), CONC(I), IMT(I)
        END DO
        PRINT *, 'qx(i),qy(i),qz(i),IMQ(i)'
        DO I = 1, NQ
            PRINT *, (BAS(J, I), J=1, 3), IMQ(I)
        END DO
        !
        PRINT *, 'Lattice vectors:'
        PRINT *, CELL(:, 1)
        PRINT *, CELL(:, 2)
        PRINT *, CELL(:, 3)
    END IF
!
!
    IF (1 .EQ. 2) WRITE (*, *) IDUM
!     ==================================================
!
!
    CALL WKINIT(MAXDIM)
    PLAT(:, :) = CELL(:, :)
    CALL GETGBASIS
!
!
! amtc,dpas,r0,rom,rop,r0,rop,rom,      dummy initialisation
!
    R0 = 9999D0
    AMTC = 9999D0
    RHOSITE(:, :) = 9999D0
!
    CALL READRO(NT, AMTC, Z, R0, TXTT, RHOSITE, NRADMAX, NRAD1)
!
! Producing new radial mesh:
!
    DO IT = 1, NT
        DO I = 1, NRAD1
            RAT(I, IT) = R0(IT)*EXP(DBLE(I - 1)*DPAS)
        END DO
    END DO
!
    DO I = 1, NM
        WSREST(I) = 0
        ZZ(I) = 0
        R0SITE(I) = R0ACT
    END DO
!
    RHOSITE(:, :) = 0.0D0
!
    DO ITYPE = 1, NT
        !
        ! this is cycle over atomic types
        ! IMT(itype) - number of inequivalent lattice position, which is
        ! occupied by this atom
        !
        ISORT = IMQ(IMT(ITYPE))
        ZN = Z(ITYPE)
        CALL DEFWSR(WS, ZN)
        IF (DB) PRINT *, 'type:', ITYPE, ' Z=', ZN, ' Rad=', WS, ' conc=', &
          &                   CONC(ITYPE)
        WSREST(ISORT) = WSREST(ISORT) + WS*CONC(ITYPE)
        ZZ(ISORT) = ZZ(ISORT) + ZN*CONC(ITYPE)
        DO IPNT = 1, NRAD1
            RAD = R0ACT*EXP(DBLE(IPNT - 1)*DPAS)
            IF (RAD .LT. RAT(NRAD1 - 2, ITYPE)) THEN
                RINT = MAX(0.D0, &
                  &                DNEVMOD(RAT(1, ITYPE), AMTC(1, ITYPE), RAD, NRAD1, &
                  &                NDNEV))
            ELSE
                RINT = 0
            END IF
            RHOSITE(IPNT, ISORT) = RHOSITE(IPNT, ISORT) + CONC(ITYPE)*RINT
        END DO
    END DO
!
    DO I = 1, NM
        IF (DB) PRINT *, ' Site:', I, '  RaD:', WSREST(I), ' Z:', ZZ(I)
    END DO
!
!      stop
!
    NSORT = NM
    NATOM = NQ
!
!
    IF (DB) THEN
        PRINT *, 'av=', CELL(:, 1)
        PRINT *, 'bv=', CELL(:, 2)
        PRINT *, 'cv=', CELL(:, 3)
        PRINT *, 'imq', (IMQ(I), I=1, NATOM)
        PRINT *, 'BAS'
        DO I = 1, NATOM
            PRINT *, (BAS(J, I), J=1, 3)
        END DO
        DO IST = 1, NSORT
            PRINT *, 'R0S', R0SITE(IST), ' ZZ', ZZ(IST), WSREST(IST)
        END DO
        PRINT *, 'dpas,alat,natom,nsort:', DPAS, ALAT, NATOM, NSORT
    END IF
!
    CALL DEFMTR(WORK, NRAD1, CELL(:, 1), CELL(:, 2), CELL(:, 3), IMQ, BAS, RHOSITE, R0SITE, DPAS, ALAT, &
      &            NATOM, NSORT, ZZ, RWS, NEED_ES, NRAD1, WSREST)
!
    NATOMNEW = NATOM
    IF (NEED_ES .EQ. 0) THEN
        FIND_EMPTY_SPHERES = EMPTYSPHERES(CELL, BAS, IMQ, &
          &                     ALAT, NATOM, NSORT, RWS, &
          &                     NSORTES, TAUES, SES, MAX2SORT, &
          &                     ISNEW, BASNEW, NATOMNEW, &
          &                     N_SYMMETRY_OPS, ROTATIONS, TRANSLATIONS)
        !ccccccccccccccccccccccccccccccccccccccccccccccccccc
        IF (FIND_EMPTY_SPHERES .ne. 0) THEN
            RETURN
        END IF
    ELSE
        FIND_EMPTY_SPHERES = 1
        DO IATOM = 1, NATOM
            ISNEW(IATOM) = IMQ(IATOM)
            DO I = 1, 3
                BASNEW(I, IATOM) = BAS(I, IATOM)
            END DO
        END DO
        !cccccccccccccccccccccccccccccccccccccccccccccccccc
    END IF
!
!
! ==================================================
!
!       WRITING OUTPUT DATA:
!
!
!     OPEN (22,FILE='radii.out')
!     WRITE (22,*) NATOMNEW
    IF (NATOMNEW > N_OUT) THEN
        N_OUT = -NATOMNEW
    ELSE
        N_OUT = NATOMNEW - NATOM

        DO I = 1, NATOM
            RADII(I + N_OUT) = RWS(IMQ(I))
        END DO
        DO I = 1, NATOMNEW - NATOM
            CENTRES(1:3, I) = BASNEW(1:3, I + NATOM)
            ISR = ISNEW(I + NATOM)
            RADII(I) = SES(ISR - NM)
        END DO
    END IF
    IF (DB) THEN
        DO I = 1, NATOMNEW
            IF (I .LE. NATOM) THEN
                WRITE (*, 99001) ISNEW(I), (BASNEW(J, I), J=1, 3), RWS(IMQ(I)), &
                  &                       ' / A/'
            ELSE
                WRITE (*, 99001) ISR, (BASNEW(J, I), J=1, 3), SES(ISR - NM), ' / E/'
            END IF
        END DO
    END IF
!     CLOSE (22)
    FLUSH (6)
99001 FORMAT(i4, 2x, 3F21.15, 3x, f21.15, 3x, a)

END FUNCTION

!*==readro.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE READRO(NSORT, RO, Z, R0, TXT, ROM, NRADMAX, NRAD1)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER NRAD1, NRADMAX, NSORT
    DOUBLE PRECISION R0(*), RO(NRADMAX + 1, *), ROM(NRAD1, *), Z(*)
    CHARACTER*4 TXT(*)
    COMMON/IPRINT/IPRINT
    INTEGER IPRINT
    !
    ! Local variables
    !
    INTEGER I, ISORT, ISR, IZC, IZN, IZNUC
    !
    !*** End of declarations rewritten by SPAG
    !
    OPEN (8, STATUS='scratch')
    DO ISORT = 1, NSORT
        IF (Z(ISORT) .LT. 0) THEN
            IZNUC = -10
            CALL QZZC(TXT(ISORT), IZNUC, IZC)
            Z(ISORT) = IZNUC
        END IF
        IF (Z(ISORT) .LT. 0.3D0) THEN
            ! emty sphere
            DO I = 1, NRAD1
                RO(I, ISORT) = 0.D0
                ROM(I, ISORT) = 0.D0
                R0(ISORT) = 1D-5
            END DO
        ELSE
            DO ISR = 1, ISORT - 1
                IF (NINT(Z(ISR)) .EQ. NINT(Z(ISORT))) THEN
                    DO I = 1, NRAD1
                        RO(I, ISORT) = RO(I, ISR)
                        ROM(I, ISORT) = ROM(I, ISR)
                    END DO
                    R0(ISORT) = R0(ISR)
                    GOTO 100
                END IF
            END DO
            !          write(buf,2)nint(z(isort)),nint(zc)
            if (IPRINT > 0) &
              &      WRITE (6, '(/'' For atom '',a4,''    Z='',f4.0/)') TXT(ISORT) &
              &             , Z(ISORT)
            IZN = NINT(Z(ISORT))
            IZC = 0
            CALL RHFDS(IZN, IZC, RO(1, ISORT), ROM(1, ISORT), R0(ISORT))
            !
        END IF
100 END DO
    CLOSE (8)
END
!*==qzzc.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE QZZC(TXT, IZ, IZC)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER IZ, IZC
    CHARACTER*2 TXT
    !
    ! Local variables
    !
    CHARACTER*210 ELEMENTS
    INTEGER I, IZCOR(105)
    CHARACTER*2 TX
    !
    !*** End of declarations rewritten by SPAG
    !
    DATA IZCOR/0, 0, 0, 6*2, 4, 4, 6*10, 12, 12, 12*18, 4*28, 30, 30, 12*36, 4*46, &
      &     48, 48, 16*54, 10*68, 4*78, 80, 80, 18*86/
    ELEMENTS = &
      &      'E H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNi' &
      &      // &
      &      'CuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeCsBaLa' &
      &      // &
      &      'CePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn' &
      &      //'FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLwKu'
    IF (IZ .GE. 0) THEN
        IZC = IZCOR(IZ + 1)
        RETURN
    END IF
    TX = TXT
    IF (TX(2:2) .EQ. '_') TX(2:2) = ' '
    DO I = 1, 105
        !        print 1,elements(i*2-1:i*2),i-1,izcor(i)
        IF (ELEMENTS(I*2 - 1:I*2) .EQ. TXT(1:2)) THEN
            IZ = I - 1
            IZC = IZCOR(I)
            RETURN
        END IF
    END DO
    WRITE (16, *) TXT
    STOP 'Error : Such Chemical element is ABSENT'
END
!*==shortn.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE SHORTN(P, P1)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! PARAMETER definitions
    !
    DOUBLE PRECISION EPS
    PARAMETER(EPS=1.D-7)
    !
    ! COMMON variables
    !
    DOUBLE PRECISION PLAT(3, 3), QLAT(3, 3)
    COMMON/L2LAT/PLAT, QLAT
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION P(3), P1(3)
    !
    ! Local variables
    !
    INTEGER I
    DOUBLE PRECISION X(3)
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
    !
    DO I = 1, 3
        X(I) = P(1)*QLAT(1, I) + P(2)*QLAT(2, I) + P(3)*QLAT(3, I)
        X(I) = X(I) - NINT(X(I) - EPS)
        !$$$        if(x(i).lt.-0.5d0+eps)x(i)=x(i)+1.d0
    END DO
    DO I = 1, 3
        P1(I) = X(1)*PLAT(I, 1) + X(2)*PLAT(I, 2) + X(3)*PLAT(I, 3)
    END DO
END
!*==vec.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE VEC(V0, V1, V2, G1)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION G1(3, 3), V0(3), V1(3), V2(3)
    !
    ! Local variables
    !
    INTEGER I, J
    !
    !*** End of declarations rewritten by SPAG
    !
    V0 = V1 + MATMUL(V2, G1)
    V0 = V0 - DNINT(V0 - 1.0D-7)  
END
!*==getgbasis.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE GETGBASIS
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! COMMON variables
    !
    DOUBLE PRECISION AVC(3), BVC(3), CVC(3), FVC(3), HVC(3), PVC(3)
    COMMON/L2LAT/AVC, BVC, CVC, FVC, PVC, HVC
    !
    ! Local variables
    !
    INTEGER I
    DOUBLE PRECISION VOLCELL, W1
    !
    !*** End of declarations rewritten by SPAG
    !
    !
    ! COMMON variables
    !
    !
    ! Local variables
    !
    !
    !*** End of declarations rewritten by SPAG
    !
    CALL CROSS(FVC, BVC, CVC)
    CALL CROSS(PVC, CVC, AVC)
    CALL CROSS(HVC, AVC, BVC)
    W1 = AVC(1)*FVC(1) + AVC(2)*FVC(2) + AVC(3)*FVC(3)
    VOLCELL = ABS(W1)
    !
    IF (VOLCELL .LT. 1.D-5) THEN
        PRINT *, AVC
        PRINT *, BVC
        PRINT *, CVC
        CALL ENDJOB(10, 'GETGBASIS: Lattice vectors are complanar!')
    END IF
    DO I = 1, 3
        FVC(I) = FVC(I)/W1
        PVC(I) = PVC(I)/W1
        HVC(I) = HVC(I)/W1
    END DO
END
