!*==defmtr.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE DEFMTR(W, NRAD1, AV, BV, CV, IS, BAS, RO, R0, DPAS, ALAT, NATOM, &
    &                  NSORT, Z, SMT, IMT, NRADMAX, WSREST)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    COMMON/IPRINT/IPRINT
    INTEGER IPRINT
    DOUBLE PRECISION ALAT, DPAS
    INTEGER IMT, NATOM, NRAD1, NRADMAX, NSORT
    DOUBLE PRECISION AV(3), BAS(3, *), BV(3), CV(3), R0(*), RO(NRADMAX, *), &
      &                 SMT(*), WSREST(*), Z(*)
    INTEGER IS(*), W(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION A(:), AVW, BOUND(3, 3), DV(:), ORIGIN(3), PI, PLAT(3, 3), &
      &                 POT(:, :), R(:), VOLCEL
    DOUBLE PRECISION DPI
    INTEGER I, J, LTMAX(3), NDEL(3), NR(:), NRMAX, OA, OBAS, ODV, ONR, OPOT, OR
    SAVE A, AVW, BOUND, DV, I, J, LTMAX, NDEL, NRMAX, OA, OBAS, ODV, ONR, OPOT, OR, &
      &     ORIGIN, PI, PLAT, POT, R, VOLCEL
    !
    !*** End of declarations rewritten by SPAG
    !
    ALLOCATABLE POT, DV, R, A, NR
    !
    PI = DPI()
    NRMAX = NRAD1 + 1
    CALL DEFDR(OA, NSORT)
    CALL DEFDR(OPOT, NSORT*NRMAX)
    CALL DEFI(ONR, NSORT)
    CALL DEFDR(OR, NRMAX)
    CALL DEFDR(ODV, NSORT*NRMAX)
    !      call hpot(W(oa),r0,nsort,W(onr),nrad1,W(opot),Z,RO,dpas,
    !     .     W(or),W(odv),nradmax)
    !
    ALLOCATE (POT(NRMAX, NSORT))
    ALLOCATE (DV(NRAD1), R(NRAD1), A(NSORT), NR(NSORT))
    !
    CALL HPOT(A, R0, NSORT, NR, NRAD1, POT, Z, RO, DPAS, R, DV, NRADMAX)
    !
    CALL RLSE(OR)
    !
    VOLCEL = AV(1)*(BV(2)*CV(3) - BV(3)*CV(2)) + AV(2) &
            &         *(BV(3)*CV(1) - BV(1)*CV(3)) + AV(3) &
            &         *(BV(1)*CV(2) - BV(2)*CV(1))
    VOLCEL = ABS(VOLCEL)
    !    average radii:
    AVW = (VOLCEL/NATOM/(4.D0*PI)*3.D0)**(1.D0/3.D0)*ALAT
    PLAT(1, 1) = AV(1)
    PLAT(1, 2) = BV(1)
    PLAT(1, 3) = CV(1)
    PLAT(2, 1) = AV(2)
    PLAT(2, 2) = BV(2)
    PLAT(2, 3) = CV(2)
    PLAT(3, 1) = AV(3)
    PLAT(3, 2) = BV(3)
    PLAT(3, 3) = CV(3)
    !
    DO I = 1, 3
        DO J = 1, 3
            BOUND(I, J) = PLAT(I, J)
        END DO
        LTMAX(I) = 2
        NDEL(I) = 0
        ORIGIN(I) = 0
    END DO
    CALL DEFDR(OBAS, NATOM*3)
    !      call defdr(owsrest,natom)
    !     call hrtree(w,alat,avw,bas,bound,
    !    .     is,ltmax,natom,nsort,ndel,
    !    .     origin,plat,smt,z,wsrest,
    !    .     W(oA),R0,W(onR),W(oPOT),NRMAX,imt)
    !
    CALL HRTREE(W, ALAT, AVW, BAS, BOUND, IS, LTMAX, NATOM, NSORT, NDEL, ORIGIN, &
      &            PLAT, SMT, Z, WSREST, A, R0, NR, POT, NRMAX, IMT)
    !

    DEALLOCATE (POT, DV, R, A, NR)

END
!*==hpot.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
SUBROUTINE HPOT(A, B, NCLASS, NR, NRAD1, POT, Z, RO, DPAS, R, DV, NRADMAX)
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! PARAMETER definitions
    !
    DOUBLE PRECISION HUGE
    PARAMETER(HUGE=1.D6)
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION DPAS
    INTEGER NCLASS, NRAD1, NRADMAX
    DOUBLE PRECISION A(*), B(*), DV(NRAD1), POT(NRAD1 + 1, *), R(NRAD1), &
      &                 RO(NRADMAX, *), Z(*)
    INTEGER NR(*)
    !
    ! Local variables
    !
    INTEGER I, ICL, NP, NRMAX
    DOUBLE PRECISION ZZ
    SAVE I, ICL, NP, NRMAX, ZZ
    !
    !*** End of declarations rewritten by SPAG
    !
    !- Calculate H. potential
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   clabl :name of the different inequivalent atom
    !i   nclass:number of classes, atoms in same class are symmetry-related
    !i   nrmax :maximum number of mesh points
    !i   nsp   :=1 spin degenerate, =2 non-degenerate
    !i   z     :nuclear charge
    !o Outputs:
    !o   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !o   b     :                 -//-
    !o   ierr  :0 if no error occurs, ic if read error for class ic
    !o   nr    :number of mesh points
    !o   pot   :spherical Hartree potential
    !w   v     :spherical potential (electronic contribution)
    !r Remarks:
    ! ----------------------------------------------------------------------
    !
    NRMAX = NRAD1 + 1                ! - radial mesh
    DO ICL = 1, NCLASS
        DO I = 1, NRAD1
            R(I) = B(ICL)*EXP(DBLE(I - 1)*DPAS)
        END DO
        NP = NRAD1
        CALL POTS(POT(2, ICL), RO(1, ICL), DV, R, DPAS, Z(ICL), NP)
        ZZ = 2*Z(ICL)
        A(ICL) = DPAS
        NR(ICL) = NRMAX
        DO I = 2, NRMAX
            POT(I, ICL) = POT(I, ICL) - ZZ/R(I - 1)
        END DO
        POT(1, ICL) = -HUGE
    END DO
END
!*==pots.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE POTS(DV, D, DP, DR, DPAS, Z, NP)
    !
    ! iHTEgPiPOBAHiE pOTEHciAlA pO 4 TO~KAM
    ! DV - POTENTIAL   D CH.DENCITY  DP WORK ARRAY  DR RADIAL MESH
    ! DPAS EXP. PASS
    ! Z ATOM NUMBER     NP - NUMBER OF POINTS
    ! **********************************************************************
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION DPAS, Z
    INTEGER NP
    DOUBLE PRECISION D(*), DP(*), DR(*), DV(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION DAS, DLO, DLO2
    INTEGER I, J, K
    SAVE DAS, DLO, DLO2, I, J, K
    !
    !*** End of declarations rewritten by SPAG
    !
    !$$$      double precision z
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
    DO I = 1, NP
        DV(I) = DV(I)/DR(I)*2
    END DO
END
!*==fmesh.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
DOUBLE PRECISION FUNCTION FMESH(A, B, F, NR, R)
!- Computes the value of f at r for given values on a mesh
! ----------------------------------------------------------------------
!i Inputs:
!i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
!i   b     :                 -//-
!i   f     :function defined on a mesh
!i   nr    :number of mesh points
!i   r     :radial distance
!o Outputs:
!o   fmesh :interpolated value of f at r
!r Remarks:
!r   if r is insite the mesh a quadratic fit is used
!r   if r is outside the mesh an exponential fit is used
! ----------------------------------------------------------------------
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
    INTEGER NP
    PARAMETER(NP=4)
!
! Dummy arguments
!
    DOUBLE PRECISION A, B, R
    INTEGER NR
    DOUBLE PRECISION F(*)
!
! Local variables
!
    DOUBLE PRECISION DELSQF, DI3INT
    INTEGER IS, NSTART
    DOUBLE PRECISION XX
    SAVE IS, NSTART, XX
    EXTERNAL DELSQF, DI3INT
!
!*** End of declarations rewritten by SPAG
!
    IF (R .GE. B) THEN
        XX = LOG(R/B)/A + 2
        IS = IDNINT(XX) - 1
    ELSE
        FMESH = F(2)
        RETURN
    END IF
    IF (IS .LE. NR - 2) THEN
        IS = MAX0(1, IS)
        FMESH = DI3INT(IS, F(IS), XX)
    ELSE
        NSTART = NR - NP + 1
        FMESH = DELSQF(NP, NSTART, F(NSTART), XX)
    END IF
END
!*==iprint.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!- Gets last integer of print priority stack
! ----------------------------------------------------------------------
!o Outputs:
!o   iprint:defines the verbosity level
!r Remarks:
!r   verbosity:   0  nearly nothing is printed
!r               10  very terse
!r               20  terse
!r               30  normal
!r               40  verbose
!r               50  very verbose
!r               60  highest verbosity
!r              100  low-level debugging
!r              110  intermediate-level debugging
!r              120  high-level debugging
! ----------------------------------------------------------------------
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
!*==potxn.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
SUBROUTINE POTXN(A, ALAT, B, BAS, IAX, ICLASS, NPR, NR, NRMAX, PLAT, POT, &
    &                 POTL, XN)
    !-Calculates Hartree potential at point xn
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   alat  :length scale
    !i   b     :                 -//-
    !i   bas   :basis vectors (scaled by alat)
    !i   iax   :information about positions around a specified pair of atom
    !i   iclass:the jth atom belongs to class iclass(j)
    !i   npr   :number of neighbors
    !i   nr    :number of mesh points
    !i   nrmax :maximum number of mesh points
    !i   plat  :primitive lattice vectors (scaled by alat)
    !i   plat  :primitive translation vectors in real space
    !i   pot   :spherical Hartree potential
    !i   xn    :cartesian coordinates where potential is calculated
    !o Outputs:
    !i   potl  :ovelapping Hartree potential
    !r Remarks:
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION ALAT, POTL
    INTEGER NPR, NRMAX
    DOUBLE PRECISION A(*), B(*), BAS(3, *), PLAT(3, 3), POT(NRMAX, *), XN(3)
    INTEGER IAX(5, *), ICLASS(*), NR(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION D(3), R
    DOUBLE PRECISION DNRM23, FMESH
    INTEGER JBAS, JC, JPR, K
    SAVE D, JBAS, JC, JPR, K, R
    EXTERNAL DNRM23, FMESH
    !
    !*** End of declarations rewritten by SPAG
    !
    POTL = 0.D0
    DO JPR = 1, NPR
        JBAS = IAX(2, JPR)
        JC = ICLASS(JBAS)
        DO K = 1, 3
            D(K) = BAS(K, JBAS) - XN(K) + PLAT(K, 1)*IAX(3, JPR) &
              &             + PLAT(K, 2)*IAX(4, JPR) + PLAT(K, 3)*IAX(5, JPR)
        END DO
        R = SQRT(DNRM23(D))*ALAT
        POTL = POTL + FMESH(A(JC), B(JC), POT(1, JC), NR(JC), R)
    END DO
    !
END
!*==di3int.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
DOUBLE PRECISION FUNCTION DI3INT(IX, YA, X)
!- Interpolates y = f(x) for given xa and ya=f(xa)
! ----------------------------------------------------------------------
!i Inputs:
!i   ix    :xa(1) is integer and equals ix
!i          and xa(1),xa(2) and xa(3) differ exactly by 1.d0
!i   ya    :value of f at xa
!i   x     :x-value at which f is interpolated
!o Outputs:
!o   di3int:interpolated value
!r Remarks:
! ----------------------------------------------------------------------
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
    INTEGER IX
    DOUBLE PRECISION X
    DOUBLE PRECISION YA(3)
!
! Local variables
!
    DOUBLE PRECISION XA(3)
!
!*** End of declarations rewritten by SPAG
!
    XA(1) = DBLE(IX)
    XA(2) = DBLE(IX + 1)
    XA(3) = DBLE(IX + 2)
!
    DI3INT = 0.5D0*(X - XA(2))*((X - XA(3))*YA(1) + (X - XA(1))*YA(3)) &
      &         - (X - XA(1))*(X - XA(3))*YA(2)
END
!*==i_shell.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
SUBROUTINE I_SHELL(M, N, IARRAY)
    !- shell sort of a array of integer vectors
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   m     :number of components in iarray
    !i   n     :number of elements in iarray
    !i   iarray:array to be sorted
    !o Outputs:
    !o   iarray:array to be sorted
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER M, N
    INTEGER IARRAY(M, 0:N - 1)
    !
    ! Local variables
    !
    INTEGER I, IT, J, K, L, LOGNB2, MM, MMM, N2, NN
    !
    !*** End of declarations rewritten by SPAG
    !
    LOGNB2 = INT(LOG(FLOAT(N + 1))*1.4426950)
    N2 = N
    DO NN = 1, LOGNB2
        N2 = N2/2
        K = N - N2
        DO J = 1, K
            I = J - 1
20          CONTINUE
            L = I + N2
            DO MM = 1, M
                IF (IARRAY(MM, L) .LT. IARRAY(MM, I)) THEN
                    DO MMM = 1, M
                        IT = IARRAY(MMM, I)
                        IARRAY(MMM, I) = IARRAY(MMM, L)
                        IARRAY(MMM, L) = IT
                    END DO
                    I = I - N2
                    IF (I .LT. 0) EXIT
                    GOTO 20
                ELSE IF (IARRAY(MM, L) .NE. IARRAY(MM, I)) THEN
                    EXIT
                END IF
            END DO
        END DO
    END DO
    !
END
!*==hrtree.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
SUBROUTINE HRTREE(W, ALAT, AVW, BAS, BOUND, ICLASS, LTMAX, NBAS, NCLASS, &
    &                  NDEL, ORIGIN, PLAT, WSR, Z, WSREST, A, B, NR, POT, NRMAX, &
    &                  IMT)
    !- Calculate Hartree potential on a mesh and determine muffin-tin radia
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   alat  :length scale
    !i   avw   :average Wigner Seitz radius
    !i   bas   :basis vectors (scaled by alat)
    !i   bound :two vectors spanning the plane (scaled by alat)
    !i   clabl :name of the different inequivalent atom
    !i   iclass:the jth atom belongs to class iclass(j)
    !i   ltmax :ltmax(i)= limit in i-dirction for unit cells considered
    !i   nbas  :number of atoms in the basis
    !i   nclass:number of classes, atoms in same class are symmetry-related
    !i   ndel  :ndel(i)=number of mesh points along the bound(i) vector
    !i   nrclas:number of atoms in the i-th class
    !i   origin:origin of the plane (scaled by alat)
    !i   plat  :primitive lattice vectors (scaled by alat)
    !i   z     :nuclear charge
    !o Outputs:
    !io  wsr   :Wigner-Seitz sphere radius (in atomic units)
    !r Remarks:
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! PARAMETER definitions
    !
    DOUBLE PRECISION FACN1
    PARAMETER(FACN1=1.5D0)
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION ALAT, AVW
    INTEGER IMT, NBAS, NCLASS, NRMAX
    DOUBLE PRECISION A(*), B(*), BAS(3, *), BOUND(3, 3), ORIGIN(*), PLAT(3, 3) &
      &                 , POT(NRMAX, *), WSR(*), WSREST(*), Z(*)
    INTEGER ICLASS(*), LTMAX(3), NDEL(*), NR(*), W(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION DSCL(:, :), WMAX
    INTEGER I, IC, J, NEIGHM, ODSCL, OIAX1, OLOCK, ONPR1
    INTEGER IDAMAX
    COMMON/IPRINT/IPRINT
    INTEGER IPRINT
    LOGICAL LOCK(:)
    SAVE DSCL, I, IC, J, LOCK, NEIGHM, ODSCL, OIAX1, OLOCK, ONPR1, WMAX
    EXTERNAL DEFDR, DEFI, DEFWSR, NGHBR1, POTMAX, POTSUM, RLSE
    !
    !*** End of declarations rewritten by SPAG
    !
    ALLOCATABLE DSCL, LOCK
    !
    IF (IPRINT .GT. 90) THEN
        PRINT *, '******************************************'
        PRINT *, 'alat', ALAT, '  avw', AVW
        PRINT *, ' === bas. v. ==='
        DO J = 1, NBAS
            PRINT 99001, (BAS(I, J), I=1, 3)
        END DO
        PRINT *, ' === bound ==='
        DO J = 1, 3
            PRINT 99001, (BOUND(I, J), I=1, 3)
        END DO
        PRINT *, ' === ltmax ==='
        PRINT 99002, LTMAX
        PRINT *, ' === ndel ==='
        PRINT 99002, (NDEL(I), I=1, 3)
        PRINT *, ' === origin ==='
        PRINT 99001, (ORIGIN(I), I=1, 3)
        PRINT *, ' === plat ==='
        DO J = 1, 3
            PRINT 99001, (PLAT(I, J), I=1, 3)
        END DO
        PRINT *, '******************************************'
    END IF
    ! --- Calculate Hartree potential on a mesh
    IF (NDEL(1)*NDEL(2) .NE. 0) CALL POTSUM(A, ALAT, B, BAS, BOUND, ICLASS, &
      &     LTMAX, NBAS, NDEL, NR, NRMAX, ORIGIN, PLAT, POT)
    ! --- Determine maximum of Hartree potential
    DO IC = 1, NCLASS
        CALL DEFWSR(WSREST(IC), Z(IC))
        WSR(IC) = WSREST(IC)
    END DO
    WMAX = WSREST(IDAMAX(NCLASS, WSREST, 1))/AVW
    I = 2*INT((2.D0*FACN1*WMAX + 1.D0)**3)
    NEIGHM = MAX(2*INT((2.D0*FACN1*WMAX + 1.D0)**3), 50)
    CALL DEFI(OIAX1, 5*NEIGHM*NCLASS)
    CALL DEFI(ONPR1, NCLASS)
    CALL NGHBR1(W, ALAT, BAS, FACN1, W(OIAX1), ICLASS, NBAS, NCLASS, NEIGHM, &
      &            W(ONPR1), PLAT, WSREST)
    !
    CALL POTMAX(W, A, ALAT, B, BAS, W(OIAX1), ICLASS, NBAS, NCLASS, W(ONPR1), &
      &            NR, NRMAX, PLAT, POT, WMAX, WSR, WSREST, Z)
    CALL DEFDR(ODSCL, NCLASS*(NCLASS + 2))
    CALL DEFDR(OLOCK, NCLASS)
    !     call blowup(alat,bas,1.d0,0.d0,w(oiax1),iclass,nbas,
    !    $           nclass,w(onpr1),0.d0,0.d0,plat,wsr,W(odscl),W(olock))
    !     if(imt.ne.0)call blowup(alat,bas,1.d0,0.d0,w(oiax1),iclass,nbas,
    !    $         nclass,w(onpr1),0.4d0,0.8d0,plat,wsr,W(odscl),W(olock))
    !
    ALLOCATE (DSCL(NCLASS, 0:NCLASS + 1))
    ALLOCATE (LOCK(NCLASS))
    !
    CALL BLOWUP(ALAT, BAS, 1.D0, 0.D0, W(OIAX1), ICLASS, NBAS, NCLASS, &
      &            W(ONPR1), 0.D0, 0.D0, PLAT, WSR, DSCL, LOCK)
    IF (IMT .NE. 0) CALL BLOWUP(ALAT, BAS, 1.D0, 0.D0, W(OIAX1), ICLASS, &
      &                            NBAS, NCLASS, W(ONPR1), 0.4D0, 0.8D0, PLAT, &
      &                            WSR, DSCL, LOCK)
    !
    !
    DEALLOCATE (DSCL)
    DEALLOCATE (LOCK)

    CALL RLSE(OIAX1)
99001 FORMAT(20F21.15)
99002 FORMAT(20I10)
    !
END
!*==potsum.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
SUBROUTINE POTSUM(A, ALAT, B, BAS, BOUND, ICLASS, LTMAX, NBAS, NDEL, NR, &
    &                  NRMAX, ORIGIN, PLAT, POT)
    !- Calculates the overlapping Hartree potential
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :                 -//-
    !i   alat  :length scale
    !i   bas   :basis vectors (scaled by alat)
    !i   bound :two vectors spanning the plane (scaled by alat)
    !i   iclass:the jth atom belongs to class iclass(j)
    !i   ltmax :ltmax(i)= limit in i-dirction for unit cells considered
    !i   nbas  :number of atoms in the basis
    !i   ndel  :ndel(i)=number of mesh points along the bound(i) vector
    !i   nr    :number of mesh points
    !i   nrmax :maximum number of mesh points
    !i   origin:origin of the plane (scaled by alat)
    !i   plat  :primitive lattice vectors (scaled by alat)
    !i   pot   :spherical Hartree potential
    !o Output written to file POT
    !r Remarks:
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION ALAT
    INTEGER NBAS, NRMAX
    COMMON/IPRINT/IPRINT
    INTEGER IPRINT
    DOUBLE PRECISION A(*), B(*), BAS(3, *), BOUND(3, 3), ORIGIN(3), PLAT(3, 3) &
      &                 , POT(NRMAX, *)
    INTEGER ICLASS(*), LTMAX(3), NDEL(*), NR(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION D01, D02, D03, D11, D12, D13, D21, D22, D23, D31, D32, D33, &
      &                 DN1, DN2, FAC1, FAC2, POTL, RAD, XN1, XN2, XN3
    DOUBLE PRECISION FMESH
    INTEGER I, I1, I2, IBAS, IC, J, J1, J2, J3
    SAVE D01, D02, D03, D11, D12, D13, D21, D22, D23, D31, D32, D33, DN1, DN2, FAC1, &
      &     FAC2, I, I1, I2, IBAS, IC, J, J1, J2, J3, POTL, RAD, XN1, XN2, XN3
    !
    !*** End of declarations rewritten by SPAG
    !
    IF (NDEL(1) .NE. 1) FAC1 = 1.D0/DBLE(NDEL(1) - 1)
    IF (NDEL(2) .NE. 1) FAC2 = 1.D0/DBLE(NDEL(2) - 1)
    IF (IPRINT > 0) THEN
        WRITE (6, 99001) ORIGIN, ((BOUND(I, J), I=1, 3), J=1, 2)
        WRITE (6, 99002)
    END IF

    !
    DO I1 = 0, NDEL(1) - 1
        DN1 = I1*FAC1
        DO I2 = 0, NDEL(2) - 1
            DN2 = I2*FAC2
            XN1 = ORIGIN(1) + BOUND(1, 1)*DN1 + BOUND(1, 2)*DN2
            XN2 = ORIGIN(2) + BOUND(2, 1)*DN1 + BOUND(2, 2)*DN2
            XN3 = ORIGIN(3) + BOUND(3, 1)*DN1 + BOUND(3, 2)*DN2
            POTL = 0.D0
            DO IBAS = 1, NBAS
                IC = ICLASS(IBAS)
                ! --------  ltmax(i) lattice translations in i-direction considered
                D01 = BAS(1, IBAS) - XN1
                D02 = BAS(2, IBAS) - XN2
                D03 = BAS(3, IBAS) - XN3
                DO J1 = -LTMAX(1), LTMAX(1)
                    D11 = D01 + PLAT(1, 1)*J1
                    D12 = D02 + PLAT(2, 1)*J1
                    D13 = D03 + PLAT(3, 1)*J1
                    DO J2 = -LTMAX(2), LTMAX(2)
                        D21 = D11 + PLAT(1, 2)*J2
                        D22 = D12 + PLAT(2, 2)*J2
                        D23 = D13 + PLAT(3, 2)*J2
                        DO J3 = -LTMAX(3), LTMAX(3)
                            D31 = D21 + PLAT(1, 3)*J3
                            D32 = D22 + PLAT(2, 3)*J3
                            D33 = D23 + PLAT(3, 3)*J3
                            RAD = SQRT(D31*D31 + D32*D32 + D33*D33)*ALAT
                            POTL = POTL + FMESH(A(IC), B(IC), POT(1, IC), NR(IC) &
                              &                         , RAD)
                        END DO
                    END DO
                END DO
            END DO
        END DO
    END DO
    !
99001 FORMAT(//25('-')//, 'make plot for plane', //, 'ORIGIN:', &
&        3F21.15/'R1    :', 3F21.15, /'R2    :', 3F21.15//, 25('-'))
99002 FORMAT('Begin to make POT ...')
END
!*==nghbr1.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE NGHBR1(W, ALAT, BAS, FACN1, IAX, ICLASS, NBAS, NCLASS, NEIGHM, &
    &                  NPR, PLAT, WSR)
    !- Create a table of all neighbors within facn1*(wsr(1)+wsr(2))
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   alat  :length scale
    !i   bas   :basis vectors (scaled by alat)
    !i   clabl :name of the different inequivalent atom
    !i   facn1 :see remarks
    !i   iclass:the jth atom belongs to class iclass(j)
    !i   nbas  :number of atoms in the basis
    !i   nclass:number of classes, atoms in same class are symmetry-related
    !i   neighm:maximum number of neighbors
    !i   plat  :primitive lattice vectors (scaled by alat)
    !i   wsr   :Wigner-Seitz sphere radius (in atomic units)
    !o Outputs:
    !i   iax   :see remarks
    !i   npr   :number of neighbors around each atom
    !r Remarks:
    !r   Creates a neighour list for a specified atom, generating iax
    !r   which contains all neigbors which fulfill:
    !r
    !r        distance(i,j) <= (wsr(i)+wsr(j))*facn1
    !r
    !r    iax(1): not used
    !r    iax(2): ibas = atom in cluster
    !r    iax(3): i
    !r    iax(4): j
    !r    iax(5): k
    !r   To be sure that at least one pair is found, fac is increased
    !r   by 1.2 until a neighbor is found.
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION ALAT, FACN1
    INTEGER NBAS, NCLASS, NEIGHM
    DOUBLE PRECISION BAS(3, *), PLAT(3, 3), WSR(*)
    INTEGER IAX(5, NCLASS, *), ICLASS(*), NPR(*), W(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION D2, DR(3), FAC, WI, WJ, WJPWI, WJPWI2
    DOUBLE PRECISION DRR2
    INTEGER I, I1, I2, I3, IBAS, IPR, IPR5, J, JBAS, JC, K, OIWK
    INTEGER ICLBAS
    COMMON/IPRINT/IPRINT
    INTEGER IPRINT
    CHARACTER*72 MESSG
    SAVE D2, DR, FAC, I, I1, I2, I3, IBAS, IPR, IPR5, J, JBAS, JC, K, MESSG, OIWK, WI, &
      &     WJ, WJPWI, WJPWI2
    !
    !*** End of declarations rewritten by SPAG
    !
    CALL DEFI(OIWK, NEIGHM)
    DO JC = 1, NCLASS
        JBAS = ICLBAS(JC, ICLASS, NBAS, 1)
        WJ = WSR(JC)
        IPR = 0
        FAC = FACN1
        DO WHILE (IPR .LT. 2)
            IPR = 0
            DO IBAS = 1, NBAS
                WI = WSR(ICLASS(IBAS))
                WJPWI = (WJ + WI)/ALAT*FAC
                CALL LATLIM(PLAT, WJPWI, I1, I2, I3)
                WJPWI2 = WJPWI*WJPWI
                ! --------- Sweep lattice translations to find all neighbors
                ! --------- within wjpwi
                DO I = -I1, I1
                    DO J = -I2, I2
                        DO K = -I3, I3
                            D2 = DRR2(PLAT, BAS(1, JBAS), BAS(1, IBAS), I, J, K, DR)
                            IF (D2 .LE. WJPWI2) THEN
                                IF (IPR .GE. NEIGHM) THEN
                                    WRITE (MESSG, 99002) NEIGHM, I1, I2, I3
                                    CALL ERRMSG(MESSG, 4)
                                END IF
                                IPR5 = IPR*5
                                W(OIWK + IPR5 + 0) = 10000*D2
                                W(OIWK + IPR5 + 1) = IBAS
                                W(OIWK + IPR5 + 2) = I
                                W(OIWK + IPR5 + 3) = J
                                W(OIWK + IPR5 + 4) = K
                                IPR = IPR + 1
                            END IF
                        END DO
                    END DO
                END DO
            END DO
            FAC = FAC*1.2D0
        END DO
        !
        CALL I_SHELL(5, IPR, W(OIWK))
        DO I = 1, IPR
            CALL INCOPY(5, W(OIWK + 5*(I - 1)), 1, IAX(1, JC, I), 1)
            IAX(1, JC, I) = JBAS
        END DO
        NPR(JC) = IPR
        !
        !
        ! ----- Printout
        IF (IPRINT .GE. 80) THEN
            !
            DO IPR = 1, NPR(JC)
                IBAS = IAX(2, JC, IPR)
                I = IAX(3, JC, IPR)
                J = IAX(4, JC, IPR)
                K = IAX(5, JC, IPR)
                D2 = DRR2(PLAT, BAS(1, JBAS), BAS(1, IBAS), I, J, K, DR)
            END DO
            IF (IPRINT > 0) WRITE (6, *)
        END IF
        !
    END DO
    IF (IPRINT .GE. 70) WRITE (6, 99001) NEIGHM, (NPR(JC), JC=1, NCLASS)
    !
99001 FORMAT(/' NGHBR1: neighm=', i4, '> npr=', 50(11I4, /26x))
99002 FORMAT(' NGHBR1: too many pairs, neighm=', i3, 6x, 'i1,i2,i3 =', 3I3, &
&        '$')
END
!*==potmax.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
SUBROUTINE POTMAX(W, A, ALAT, B, BAS, IAX1, ICLASS, NBAS, NCLASS, NPR1, NR, &
    &                  NRMAX, PLAT, POT, WMAX, WSR, WSREST, Z)
    !- Calculates the overlapping Hartree potential and finds it's maximum
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   alat  :length scale
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   bas   :basis vectors (scaled by alat)
    !i   clabl :name of the different inequivalent atom
    !i   iax1  :information about positions around a specified atom
    !i   iclass:the jth atom belongs to class iclass(j)
    !i   nbas  :number of atoms in the basis
    !i   nclass:number of classes, atoms in same class are symmetry-related
    !i   npr1  :number of neighbors
    !i   nr    :number of mesh points
    !i   nrmax :maximum number of mesh points
    !i   plat  :primitive lattice vectors (scaled by alat)
    !i   pot   :spherical Hartree potential
    !i   wmax  :largest wsr/avw
    !i   wsrest:estimation of Wigner-Seitz sphere radius (in atomic units)
    !i   z     :nuclear charge
    !o Outputs:
    !o   wsr   :Wigner-Seitz sphere radius (in atomic units)
    !r Remarks:
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! PARAMETER definitions
    !
    INTEGER NSTEP
    PARAMETER(NSTEP=20)
    DOUBLE PRECISION FACN2, TINY, TOL
    PARAMETER(FACN2=3.D0, TINY=1.D-5, TOL=1.D-2)
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION ALAT, WMAX
    INTEGER NBAS, NCLASS, NRMAX
    DOUBLE PRECISION A(*), B(*), BAS(3, *), PLAT(3, 3), POT(NRMAX, *), WSR(*), &
      &                 WSREST(*), Z(*)
    INTEGER IAX1(5, NCLASS, *), ICLASS(*), NPR1(*), NR(*), W(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION DLAMCH, DRR2
    DOUBLE PRECISION DR(3), POTL, RIK, RIKO, RPMX, STEP, X(:), XN(3), Y(:)
    INTEGER I, IBAS, IC, IMIN, ISTEP, K, K1, K2, K3, KBAS, KC, KCO, KPR, NEIGHM, &
      &        NPR2, OIAX2
    COMMON/IPRINT/IPRINT
    INTEGER IPRINT
    INTEGER ICLBAS
    SAVE DR, I, IBAS, IC, IMIN, ISTEP, K, K1, K2, K3, KBAS, KC, KCO, KPR, NEIGHM, &
      &     NPR2, OIAX2, POTL, RIK, RIKO, RPMX, STEP, X, XN, Y
    EXTERNAL DEFI, DLAMCH, DRR2, DSCAL, ICLBAS, NGHBR2, POTXN
    !
    !*** End of declarations rewritten by SPAG
    !
    ALLOCATABLE x, y
    ALLOCATE (X(0:NSTEP), Y(0:NSTEP))
    NEIGHM = MAX0(INT(4.D0*(FACN2*WMAX + 1.D0)**3), 50)
    CALL DEFI(OIAX2, 5*NEIGHM)
    !
    DO IC = 1, NCLASS
        IF (IDNINT(Z(IC)) .NE. 0) THEN
            IF (IPRINT .GE. 40) WRITE (6, 99001)
            IBAS = ICLBAS(IC, ICLASS, NBAS, 1)
            !
            WSR(IC) = DLAMCH('o')
            RIKO = -1.D0
            KCO = -1
            DO KPR = 2, NPR1(IC)
                KBAS = IAX1(2, IC, KPR)
                K1 = IAX1(3, IC, KPR)
                K2 = IAX1(4, IC, KPR)
                K3 = IAX1(5, IC, KPR)
                KC = ICLASS(KBAS)
                RIK = SQRT(DRR2(PLAT, BAS(1, IBAS), BAS(1, KBAS), K1, K2, K3, DR) &
                  &               )*ALAT
                IF (IDNINT(Z(KC)) .NE. 0 .AND. &
                  &              (ABS(RIK - RIKO) .GT. TINY .OR. KCO .NE. KC)) THEN
                    !
                    RIKO = RIK
                    KCO = KC
                    STEP = RIK/NSTEP
                    CALL NGHBR2(ALAT, BAS, FACN2, W(OIAX2), IBAS, ICLASS, K1, K2, &
                      &                        K3, KBAS, NBAS, NEIGHM, NPR2, PLAT, WSREST)
                    CALL DSCAL(3, 1.D0/DFLOAT(NSTEP), DR, 1)
                    ! --------  Find aproximate position of first maximun
                    DO ISTEP = 0, NSTEP
                        DO K = 1, 3
                            XN(K) = BAS(K, IBAS) + ISTEP*DR(K)
                        END DO
                        CALL POTXN(A, ALAT, B, BAS, W(OIAX2), ICLASS, NPR2, NR, &
                          &                          NRMAX, PLAT, POT, POTL, XN)
                        X(ISTEP) = DBLE(ISTEP)*STEP
                        Y(ISTEP) = MAX(POTL, -1.D5)
                        IF (ISTEP .NE. 0) THEN
                            IF (ISTEP .GE. 2 .AND. Y(ISTEP) .LT. Y(ISTEP - 1)) &
                              &                       GOTO 5
                        END IF
                    END DO
                    CALL ERRMSG('bug in potmax.$', 4)
                    !
                    ! --------- Now a single maximum is supposed between istep-2 and istep
5                   CONTINUE
                    DO K = 1, 3
                        X(K) = X(ISTEP - 3 + K)
                        Y(K) = Y(ISTEP - 3 + K)
                    END DO
                    RPMX = 0.D0
                    IMIN = (ISTEP - 2)*2
                    DO WHILE (X(3) - X(1) .GT. TOL)
                        X(5) = X(3)
                        X(3) = X(2)
                        Y(5) = Y(3)
                        Y(3) = Y(2)
                        STEP = STEP*0.5D0
                        CALL DSCAL(3, 0.5D0, DR, 1)
                        DO I = 1, 3, 2
                            DO K = 1, 3
                                XN(K) = BAS(K, IBAS) + (IMIN + I)*DR(K)
                            END DO
                            CALL POTXN(A, ALAT, B, BAS, W(OIAX2), ICLASS, NPR2, NR, &
                              &                             NRMAX, PLAT, POT, Y(I + 1), XN)
                            X(I + 1) = (IMIN + I)*STEP
                        END DO
                        DO I = 2, 4
                            IF (Y(I + 1) .LT. Y(I)) THEN
                                DO K = 1, 3
                                    Y(K) = Y(I - 2 + K)
                                    X(K) = X(I - 2 + K)
                                END DO
                                EXIT
                            END IF
                        END DO
                        IMIN = 2*(IMIN + I - 2)
                    END DO
                    RPMX = 0.5D0*(Y(1)*(X(2) + X(3)) - 2.D0*Y(2)*(X(1) + X(3)) &
                      &                   + Y(3)*(X(1) + X(2)))/(Y(1) - Y(2) - Y(2) + Y(3))
                    WSR(IC) = MIN(WSR(IC), RPMX)
                END IF
            END DO
            IF (IPRINT .GE. 40) WRITE (6, 99002)
            !        if (iprint().ge.20) write(6,303)clabl(ic),wsr(ic)
        END IF
    END DO
!
99001 FORMAT(/' POTMAX: ATOM1   ATOM2   DIST      VMAX   R(VMAX)', /, 8x, &
&        42('-'))
99002 FORMAT(8x, 42('-'))
!
    DEALLOCATE (x, y)
END
!*==nghbr2.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE NGHBR2(ALAT, BAS, FACN2, IAX, IBAS, ICLASS, K1, K2, K3, KBAS, &
    &                  NBAS, NEIGHM, NPR, PLAT, WSR)
    !- Create a table of all neighbors around two atoms
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   alat  :length scale
    !i   bas   :basis vectors (scaled by alat)
    !i   clabl :name of the different inequivalent atom
    !i   facn2 :see remarks
    !i   ibas  :index to first atom
    !i   iclass:the jth atom belongs to class iclass(j)
    !i   k1-3  :second atom is translated by k1*plat1+k2*plat2+p3*plat3
    !i   kbas  :index to second atom
    !i   nbas  :number of atoms in the basis
    !i   neighm:maximum number of neighbors
    !i   plat  :primitive lattice vectors (scaled by alat)
    !i   wsr   :Wigner-Seitz sphere radius (in atomic units)
    !o Outputs:
    !i   iax   :see remarks
    !i   npr   :number of neighbors
    !r Remarks:
    !r   Creates a neighour list for a pair of atoms ibas and kbas,
    !r   generating iax which contains all neigbors jbas within
    !r   facn2*wsr(jbas) around ibas or kbas , i.e. those which fulfill:
    !r
    !r        distance(i,j) <= wsr(j)*facn2
    !r    or  distance(k,j) <= wsr(j)*facn2
    !r
    !r    iax(1): not used
    !r    iax(2): ibas = atom in cluster
    !r    iax(3): i
    !r    iax(4): j
    !r    iax(5): k
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION ALAT, FACN2
    INTEGER IBAS, K1, K2, K3, KBAS, NBAS, NEIGHM, NPR
    DOUBLE PRECISION BAS(3, *), PLAT(3, 3), WSR(*)
    INTEGER IAX(5, *), ICLASS(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION B(3, 2), D1, D2, DB11, DB12, DB21, DB22, DB31, DB32, DR11, &
      &                 DR12, DR13, DR21, DR22, DR23, DR31, DR32, DR33, DRB11, &
      &                 DRB12, DRB21, DRB22, DRB31, DRB32, WJF, WJF2
    INTEGER I, I1, I2, I3, II, IPR, J, JBAS, JC, K
    COMMON/IPRINT/IPRINT
    INTEGER IPRINT

    CHARACTER*72 MESSG
    SAVE B, D1, D2, DB11, DB12, DB21, DB22, DB31, DB32, DR11, DR12, DR13, DR21, &
      &     DR22, DR23, DR31, DR32, DR33, DRB11, DRB12, DRB21, DRB22, DRB31, DRB32, &
      &     I, I1, I2, I3, II, IPR, J, JBAS, JC, K, MESSG, WJF, WJF2
    EXTERNAL DCOPY, LATLIM
    !
    !*** End of declarations rewritten by SPAG
    !
    CALL DCOPY(3, BAS(1, IBAS), 1, B(1, 1), 1)
    B(1, 2) = BAS(1, KBAS) + K1*PLAT(1, 1) + K2*PLAT(1, 2) + K3*PLAT(1, 3)
    B(2, 2) = BAS(2, KBAS) + K1*PLAT(2, 1) + K2*PLAT(2, 2) + K3*PLAT(2, 3)
    B(3, 2) = BAS(3, KBAS) + K1*PLAT(3, 1) + K2*PLAT(3, 2) + K3*PLAT(3, 3)
    IPR = 0
    DO JBAS = 1, NBAS
        JC = ICLASS(JBAS)
        WJF = WSR(JC)/ALAT*FACN2
        CALL LATLIM(PLAT, WJF, I1, I2, I3)
        WJF2 = WJF*WJF
        DB11 = BAS(1, JBAS) - B(1, 1)
        DB21 = BAS(2, JBAS) - B(2, 1)
        DB31 = BAS(3, JBAS) - B(3, 1)
        DB12 = BAS(1, JBAS) - B(1, 2)
        DB22 = BAS(2, JBAS) - B(2, 2)
        DB32 = BAS(3, JBAS) - B(3, 2)
        DO I = -I1, I1
            DR11 = PLAT(1, 1)*I
            DR21 = PLAT(2, 1)*I
            DR31 = PLAT(3, 1)*I
            DO J = -I2, I2
                DR12 = DR11 + PLAT(1, 2)*J
                DR22 = DR21 + PLAT(2, 2)*J
                DR32 = DR31 + PLAT(3, 2)*J
                DO K = -I3, I3
                    DR13 = DR12 + PLAT(1, 3)*K
                    DR23 = DR22 + PLAT(2, 3)*K
                    DR33 = DR32 + PLAT(3, 3)*K
                    DRB11 = DR13 + DB11
                    DRB21 = DR23 + DB21
                    DRB31 = DR33 + DB31
                    DRB12 = DR13 + DB12
                    DRB22 = DR23 + DB22
                    DRB32 = DR33 + DB32
                    D1 = DRB11*DRB11 + DRB21*DRB21 + DRB31*DRB31
                    D2 = DRB12*DRB12 + DRB22*DRB22 + DRB32*DRB32
                    IF (D1 .LE. WJF2 .OR. D2 .LE. WJF2) THEN
                        IPR = IPR + 1
                        IF (IPR .GE. 2*NEIGHM) THEN
                            WRITE (MESSG, 99002) 2*NEIGHM, I1, I2, I3
                            CALL ERRMSG(MESSG, 4)
                        END IF
                        IAX(2, IPR) = JBAS
                        IAX(3, IPR) = I
                        IAX(4, IPR) = J
                        IAX(5, IPR) = K
                        II = 0
                        IF (D1 .LE. WJF2) II = II + 1
                        IF (D2 .LE. WJF2) II = II + 2
                        ! ------------- Printout
                        IF (IPRINT .GE. 80) WRITE (6, 99001) IPR, D1, D2, &
                          &                    JBAS, I, J, K, II
                    END IF
                END DO
            END DO
        END DO
    END DO
    !
    NPR = IPR
    !
99001 FORMAT(i3, 2F21.15, 4I5, 2x, i3)
99002 FORMAT(' NGHBR2: too many pairs, neighm=', i3, 6x, 'i1,i2,i3 =', 3I3, &
&        '$')
END
!*==blowup.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
SUBROUTINE BLOWUP(ALAT, BAS, FACVOL, GAMMA, IAX1, ICLASS, NBAS, NCLASS, &
    &                  NPR1, OMMAX1, OMMAX2, PLAT, WSR, DSCL, LOCK)
    !- Blows up the spheres
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   alat  :length scale
    !i   bas   :basis vectors (scaled by alat)
    !i   clabl :name of the different inequivalent atom
    !i   iax1  :information about positions around a specified atom
    !i   iclass:the jth atom belongs to class iclass(j)
    !i   nbas  :number of atoms in the basis
    !i   nclass:number of classes, atoms in same class are symmetry-related
    !i   npr1  :number of neighbors
    !i   facvol:sum of sphere volumes = facvol * cell volume
    !i   gamma :scaling is r -> a(r+b) with a*b=gamma*(a-1)*avw
    !i   ommax1:maximum overlap divided by distance (s1+s2-d)/d  < ommax1
    !i   ommax2:maximum overlap divided by  radius  (s1+s2-d)/s1 < ommax2
    !i   plat  :primitive lattice vectors (scaled by alat)
    !io  wsr   :Wigner-Seitz sphere radius (in atomic units)
    !o Outputs:
    !io  wsr   :Wigner-Seitz sphere radius (in atomic units)
    !r Remarks:
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! PARAMETER definitions
    !
    DOUBLE PRECISION FPI3, TINY
    PARAMETER(FPI3=4.18879020478639053D0, TINY=1.D-4)
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION ALAT, FACVOL, GAMMA, OMMAX1, OMMAX2
    INTEGER NBAS, NCLASS
    DOUBLE PRECISION BAS(3, *), DSCL(NCLASS, 0:NCLASS + 1), PLAT(3, 3), WSR(*)
    INTEGER IAX1(5, NCLASS, *), ICLASS(*), NPR1(*)
    LOGICAL LOCK(*)
    !
    ! Local variables
    !
    DOUBLE PRECISION A, AMAX, AMAX1, AMAX2, AMAX3, AVW, B, BMAX, D, DM(0:3), &
      &                 DR(3), DSCLMX, GW, OMO2, OPO1, P, Q, R, RIK, RIKO, S, T, &
      &                 TMO2, U, V, VOL, VOLA, VOLB, VOLSPH, WSRI, WSRK, X
    DOUBLE PRECISION DDET33, DLAMCH, DRR2
    LOGICAL FIN
    INTEGER IBAS, IC, ICMAX, ICMIN, ILOOP, K1, K2, K3, KBAS, KC, KCO, KPR, LINE, &
      &        NLOOP
    INTEGER ICLBAS
    COMMON/IPRINT/IPRINT
    INTEGER IPRINT
    CHARACTER*72 MESSG
    SAVE A, AMAX, AMAX1, AMAX2, AMAX3, AVW, B, BMAX, D, DM, DR, DSCLMX, FIN, GW, &
      &     IBAS, IC, ICMAX, ICMIN, ILOOP, K1, K2, K3, KBAS, KC, KCO, KPR, LINE, &
      &     MESSG, NLOOP, OMO2, OPO1, P, Q, R, RIK, RIKO, S, T, TMO2, U, V, VOL, VOLA, &
      &     VOLB, VOLSPH, WSRI, WSRK, X
    EXTERNAL DLAMCH, DRR2, ERRMSG, ICLBAS
    !
    !*** End of declarations rewritten by SPAG
    !
    OPO1 = 1.D0 + OMMAX1
    OMO2 = 1.D0 - OMMAX2
    TMO2 = 2.D0 - OMMAX2
    VOL = ABS(DDET33(PLAT))*ALAT**3
    AVW = (VOL/FPI3/NBAS)**(1.D0/3.D0)
    GW = AVW*GAMMA
    !
    VOLSPH = 0.D0
    DO IBAS = 1, NBAS
        VOLSPH = VOLSPH + FPI3*WSR(ICLASS(IBAS))**3
    END DO
    IF (IPRINT .GE. 30) WRITE (6, 99009) VOL, VOLSPH
    !
    !     call dcopy(nclass,1.d0,0,dscl(1,0),1)
    DSCL(1:NCLASS, 0) = 1.D0
    !
    DO NLOOP = 1, NCLASS + 1
        AMAX = DLAMCH('o')
        ! ----- lock radii of those spheres whith maximum allowed overlap
        LOCK(1:NCLASS) = .FALSE.
        DO IC = 1, NCLASS
            IBAS = ICLBAS(IC, ICLASS, NBAS, 1)
            WSRI = WSR(IC)
            DO KPR = 2, NPR1(IC)
                KBAS = IAX1(2, IC, KPR)
                K1 = IAX1(3, IC, KPR)
                K2 = IAX1(4, IC, KPR)
                K3 = IAX1(5, IC, KPR)
                KC = ICLASS(KBAS)
                WSRK = WSR(KC)
                RIK = SQRT(DRR2(PLAT, BAS(1, IBAS), BAS(1, KBAS), K1, K2, K3, DR) &
                  &               )
                RIK = RIK*ALAT
                IF (ABS(OPO1*RIK - WSRI - WSRK) .LT. TINY) LOCK(IC) = .TRUE.
                IF (ABS(RIK - OMO2*WSRI - WSRK) .LT. TINY) LOCK(IC) = .TRUE.
                IF (ABS(RIK - WSRI - OMO2*WSRK) .LT. TINY) LOCK(IC) = .TRUE.
            END DO
        END DO
        !
        DO IC = 1, NCLASS
            IF (.NOT. LOCK(IC)) THEN
                IBAS = ICLBAS(IC, ICLASS, NBAS, 1)
                RIKO = -1.D0
                KCO = -1
                WSRI = WSR(IC)
                DO KPR = 2, NPR1(IC)
                    KBAS = IAX1(2, IC, KPR)
                    K1 = IAX1(3, IC, KPR)
                    K2 = IAX1(4, IC, KPR)
                    K3 = IAX1(5, IC, KPR)
                    KC = ICLASS(KBAS)
                    RIK = SQRT(DRR2(PLAT, BAS(1, IBAS), BAS(1, KBAS), K1, K2, K3, &
                      &                  DR))
                    RIK = RIK*ALAT
                    IF (ABS(RIK - RIKO) .GT. TINY .OR. KCO .NE. KC) THEN
                        WSRK = WSR(KC)
                        RIKO = RIK
                        KCO = KC
                        IF (LOCK(KC)) THEN
                            AMAX1 = (OPO1*RIK - WSRK + GW)/(WSRI + GW)
                            IF (OMO2 .GT. 0.D0) AMAX2 = (RIK - WSRK + OMO2*GW) &
                              &                       /(WSRI + GW)/OMO2
                            AMAX3 = (RIK - OMO2*WSRK + GW)/(WSRI + GW)
                        ELSE
                            AMAX1 = (OPO1*RIK + GW + GW)/(WSRI + WSRK + GW + GW)
                            IF (WSRK + OMO2*WSRI + TMO2*GW .GT. 0.D0) &
                              &                       AMAX2 = (RIK + TMO2*GW) &
                              &                       /(WSRK + OMO2*WSRI + TMO2*GW)
                            IF (WSRI + OMO2*WSRK + TMO2*GW .GT. 0.D0) &
                              &                       AMAX3 = (RIK + TMO2*GW) &
                              &                       /(WSRI + OMO2*WSRK + TMO2*GW)
                        END IF
                        AMAX = MIN(AMAX, AMAX1, AMAX2, AMAX3)
                    END IF
                END DO
            END IF
        END DO
        BMAX = (1.D0 - 1.D0/AMAX)*GW
        !
        VOLA = 0.D0
        VOLB = 0.D0
        DO IBAS = 1, NBAS
            IC = ICLASS(IBAS)
            IF (.NOT. LOCK(IC)) THEN
                VOLB = VOLB + (AMAX*WSR(IC) + BMAX)**3
            ELSE
                VOLA = VOLA + WSR(IC)**3
            END IF
        END DO
        VOLA = VOLA*FPI3
        VOLB = VOLB*FPI3
        !
        IF (VOL*FACVOL .LT. VOLA + VOLB) THEN
            !
            DM(:) = 0.0D0
            DO IBAS = 1, NBAS
                IC = ICLASS(IBAS)
                IF (.NOT. LOCK(IC)) THEN
                    A = WSR(IC) + GW
                    ! ----------- for numerical reasons distinguish cases
                    IF (ABS(GAMMA) .GT. 1.D0) THEN
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
            IF (ABS(DM(3)) .GT. TINY) THEN
                R = DM(2)/DM(3)
                S = DM(1)/DM(3)
                T = (DM(0) - (VOL*FACVOL - VOLA)/FPI3)/DM(3)
                P = S - R*R/3.D0
                Q = 2.D0*R*R*R/27.D0 - R*S/3.D0 + T
                D = P*P*P/27.D0 + Q*Q/4.D0
                U = (SQRT(D) - Q/2.D0)**(1.D0/3.D0)
                V = -P/U/3.D0
                X = U + V - R/3.D0
                IF (ABS(GAMMA) .GT. 1.D0) THEN
                    AMAX = X + 1.D0
                    BMAX = X*GW/AMAX
                ELSE
                    AMAX = X
                    BMAX = (1.D0 - 1.D0/AMAX)*GW
                END IF
                IF (IPRINT .GE. 100) THEN
                    WRITE (6, 99008) 'R S T', R, S, T
                    WRITE (6, 99008) 'P Q  ', P, Q
                    WRITE (6, 99008) '  D  ', D
                    WRITE (6, 99008) ' U V ', U, V
                    WRITE (6, 99008) ' AMAX', AMAX
                    WRITE (6, 99008) ' BMAX', BMAX
                    WRITE (6, 99008) ' -------------------------'
                END IF
            END IF
        END IF
        !
        FIN = .TRUE.
        DSCLMX = 1.D0
        DO IC = 1, NCLASS
            DSCL(IC, NLOOP) = 1.D0
            IF (.NOT. LOCK(IC)) THEN
                DSCLMX = AMAX + BMAX/WSR(IC)
                WSR(IC) = DSCLMX*WSR(IC)
                DSCL(IC, 0) = DSCLMX*DSCL(IC, 0)
                DSCL(IC, NLOOP) = DSCLMX
                FIN = .FALSE.
            END IF
        END DO
        FIN = FIN .OR. ABS(DSCLMX - 1.D0) .LT. TINY/256
        !
        VOLSPH = 0.D0
        DO IBAS = 1, NBAS
            VOLSPH = VOLSPH + FPI3*WSR(ICLASS(IBAS))**3
        END DO
        !
        !
        IF (FIN .OR. NLOOP .EQ. NCLASS + 1) THEN
            IF (IPRINT .GE. 30) THEN
                WRITE (6, 99001) OMMAX1, OMMAX2
                DO LINE = 0, (NCLASS - 1)/6
                    ICMIN = 1 + 6*LINE
                    ICMAX = MIN0(NCLASS, 6 + 6*LINE)
                    !              write(6,301)(clabl(ic),ic=icmin,icmax)
                    WRITE (6, 99007) ('============', IC=ICMIN, ICMAX + 1)
                    WRITE (6, 99002) (WSR(IC)/DSCL(IC, 0), IC=ICMIN, ICMAX)
                    WRITE (6, 99007) ('------------', IC=ICMIN, ICMAX + 1)
                    DO ILOOP = 1, NLOOP
                        WRITE (6, 99003) ILOOP, &
                          &                               (DSCL(IC, ILOOP), IC=ICMIN, ICMAX)
                    END DO
                    WRITE (6, 99007) ('------------', IC=ICMIN, ICMAX + 1)
                    WRITE (6, 99004) (WSR(IC), IC=ICMIN, ICMAX)
                    WRITE (6, 99005) (DSCL(IC, 0), IC=ICMIN, ICMAX)
                    WRITE (6, 99006) (WSR(IC) - WSR(IC)/DSCL(IC, 0), IC=ICMIN, &
                      &                            ICMAX)
                    WRITE (6, 99007) ('------------', IC=ICMIN, ICMAX + 1)
                END DO
                WRITE (6, 99009) VOL, VOLSPH
            END IF
            IF (ABS(OMMAX1*OMMAX2) .GT. TINY .AND. ABS(VOLSPH - VOL*FACVOL) &
              &           .GT. TINY) THEN
                WRITE (MESSG, 99010)
                CALL ERRMSG(MESSG, 0)
            END IF
            RETURN
        END IF
    END DO
!
99001 FORMAT(/' BLOWUP: scale radii;   OMMAX1=', f21.15, ', OMMAX2=', &
&        f21.15)
99002 FORMAT(6x, 'OLD WSR:    ', 6F21.15)
99003 FORMAT(6x, 'LOOP:', i3, '  * ', 6F21.15)
99004 FORMAT(6x, 'NEW WSR:    ', 6F21.15)
99005 FORMAT(6x, 'W_NEW/W_OLD:', 6F21.15)
99006 FORMAT(6x, 'W_NEW-W_OLD:', 6F21.15)
99007 FORMAT(6x, a12, 6A9)
99008 FORMAT(6x, a, 3F21.15)
99009 FORMAT(/' CELL VOLUME=', f21.15, '   SUM OF SPHERE VOLUMES=', &
&        f21.15)
99010 FORMAT(' BLOWUP: impossible to reach VOL, increase OMMAX.$')
!
END
!*==delsqf.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
DOUBLE PRECISION FUNCTION DELSQF(NDATA, IX, YA, X)
!- exponential least square fit (  a * exp(bx) )
! ----------------------------------------------------------------------
!i Inputs:
!i   ndata :number of data points
!i   ix    :xa(1) is integer and equals ix
!i          and xa(1),xa(2),xa(3),... differ exactly by 1.d0
!i   ya    :value of f at xa
!i   x     :x-value at which f is interpolated
!o Outputs:
!o   delsqf:interpolated value
!r Remarks:
! ----------------------------------------------------------------------
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
    INTEGER IX, NDATA
    DOUBLE PRECISION X
    DOUBLE PRECISION YA(*)
!
! Local variables
!
    DOUBLE PRECISION A, B, SS, SX, SXX, SXY, SY, XX, YY
    INTEGER IDATA
    SAVE A, B, IDATA, SS, SX, SXX, SXY, SY, XX, YY
!
!*** End of declarations rewritten by SPAG
!
    SX = 0.D0
    SY = 0.D0
    SXY = 0.D0
    SXX = 0.D0
!
    DO IDATA = 1, NDATA
        XX = DBLE(IX + IDATA - 1)
        IF (YA(IDATA) .LE. 0.D0) THEN
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
!
    B = (SXY - SX*SY/SS)/(SXX - SX*SX/SS)
    A = EXP((SY - B*SX)/SS)
!
    DELSQF = A*EXP(B*X)
!
END
!*==drr2.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
DOUBLE PRECISION FUNCTION DRR2(PLAT, BAS1, BAS2, I, J, K, DR)
!- Calculates the vector connecting two sites in a solid
! ----------------------------------------------------------------
!i Inputs:
!i   plat  :primitive lattice vectors
!i   bas1  :basis vector of first site
!i   bas2  :basis vector of second site
!i   i,j,k :the number of primitive lattice vectors separating sites
!o Outputs:
!o   dr    :connecting vector bas2 - bas1
!o   drr2  :square of the length of this vector
!r Remarks:
!r   Using the TB package and a table of indices iax, the connecting
!r   vector and the square of the distance is obtained by:
!r   rsqr=drr2(plat,bas(1,iax(1)),bas(1,iax(2)),iax(3),iax(4),iax(5),dr)
! ----------------------------------------------------------------
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
    INTEGER I, J, K
    DOUBLE PRECISION BAS1(3), BAS2(3), DR(3), PLAT(3, 3)
!
! Local variables
!
    INTEGER IX
!
!*** End of declarations rewritten by SPAG
!
    DRR2 = 0.D0
    DO IX = 1, 3
        DR(IX) = BAS2(IX) - BAS1(IX) + PLAT(IX, 1)*I + PLAT(IX, 2) &
                &            *J + PLAT(IX, 3)*K
        DRR2 = DRR2 + DR(IX)*DR(IX)
    END DO
END
!*==iclbas.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
INTEGER FUNCTION ICLBAS(IC, ICLASS, NBAS, NRBAS)
!- Returns an index to iclbas atom in basis given class
! ----------------------------------------------------------------------
!i Inputs:
!i   ic    :class label
!i   iclass:the jth atom belongs to class iclass(j)
!i   nbas  :number of atoms in the basis
!i   nrbas :the nrbas-th basis atom of class ic is seeked
!o Outputs:
!o   iclbas:the iclbas-th atom belongs to class ic
! ----------------------------------------------------------------------
    IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
    INTEGER IC, NBAS, NRBAS
    INTEGER ICLASS(*)
!
! Local variables
!
    INTEGER IBAS, N
    CHARACTER*72 MESSG
    SAVE IBAS, MESSG, N
    EXTERNAL ERRMSG
!
!*** End of declarations rewritten by SPAG
!
    ICLBAS = 1
    N = 0
    DO IBAS = 1, NBAS
        IF (ICLASS(IBAS) .EQ. IC) N = N + 1
        IF (N .EQ. NRBAS) THEN
            ICLBAS = IBAS
            RETURN
        END IF
    END DO
    WRITE (MESSG, 99001) IC, NRBAS, N
    CALL ERRMSG(MESSG, 4)
!
99001 FORMAT(' ICLBAS: class', i3, ', nrbas=', i3, ', but only', i3, &
&        ' basis atoms exist.$')
END
!*==latlim.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
!
SUBROUTINE LATLIM(PQLAT, RMAXS, I1, I2, I3)
    !- Set limits in X Y Z direction
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   pqlat :primitive lattice vectors (real or reciprocal space)
    !i   rmaxs :maximum length of connecting vector
    !o Outputs:
    !o   i1,i2,i3:all connecting vectors lie within these multiples of the
    !o            lattice vectors.
    !r Remarks:
    !r   This routine only returns the integer part
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    INTEGER I1, I2, I3
    DOUBLE PRECISION RMAXS
    DOUBLE PRECISION PQLAT(3, 3)
    !
    ! Local variables
    !
    DOUBLE PRECISION DET, QPLAT(3, 3), X1, X2, X3
    DOUBLE PRECISION DLAMCH, DNRM23
    INTEGER IPRINT
    COMMON/IPRINT/IPRINT
    SAVE DET, QPLAT, X1, X2, X3
    EXTERNAL DINV33, DLAMCH, DNRM23
    !
    !*** End of declarations rewritten by SPAG
    !
    CALL DINV33(PQLAT, 1, QPLAT, DET)
    X1 = RMAXS*SQRT(DNRM23(QPLAT(1, 1)))
    X2 = RMAXS*SQRT(DNRM23(QPLAT(1, 2)))
    X3 = RMAXS*SQRT(DNRM23(QPLAT(1, 3)))
    !
    I1 = 1 + INT(X1 - DLAMCH('p'))
    I2 = 1 + INT(X2 - DLAMCH('p'))
    I3 = 1 + INT(X3 - DLAMCH('p'))
    IF (IPRINT .GE. 80) WRITE (6, 99001) RMAXS, X1, X2, X3, I1, I2, I3
    !
99001 FORMAT(/' LATLIM: rmaxs=', f21.15, ', x1,x2,x3=', 3F21.15, &
&        ', i1,i2,i3=', 3I2)
END
!*==defwsr.f    processed by SPAG 7.10RU at 07:46 on 26 Mar 2020
SUBROUTINE DEFWSR(WSR, Z)
    !- Returns default value of radii for given nuclear charge
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   z     :nuclear charge
    !o Outputs:
    !o   wsr   :Wigner-Seitz sphere radius (in atomic units)
    !r Remarks:
    !r These are equilibrium average Wigner-Seitz radii for closed-packed
    !r structure of pure element.
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    !
    !*** Start of declarations rewritten by SPAG
    !
    ! Dummy arguments
    !
    DOUBLE PRECISION WSR, Z
    !
    ! Local variables
    !
    DOUBLE PRECISION CPWSR(0:100)
    INTEGER IZ
    CHARACTER*72 MESSG
    SAVE IZ, MESSG
    EXTERNAL ERRMSG
    !
    !*** End of declarations rewritten by SPAG
    !
    DATA CPWSR/2.00D0, 1.39D0, 2.55D0, 3.04D0, 2.27D0, 1.96D0, 1.66D0, &
      &     1.90D0, 1.90D0, 2.17D0, 2.89D0, 3.76D0, 3.25D0, 2.95D0, 2.63D0, &
      &     2.56D0, 2.70D0, 2.85D0, 3.71D0, 4.66D0, 3.88D0, 3.31D0, 2.99D0, &
      &     2.76D0, 2.64D0, 2.57D0, 2.52D0, 2.52D0, 2.55D0, 2.62D0, 2.78D0, &
      &     2.75D0, 2.79D0, 2.83D0, 2.94D0, 3.13D0, 4.32D0, 4.95D0, 4.22D0, &
      &     3.61D0, 3.28D0, 3.03D0, 2.91D0, 2.82D0, 2.77D0, 2.78D0, 2.84D0, &
      &     2.95D0, 3.14D0, 3.30D0, 3.45D0, 3.30D0, 3.31D0, 3.50D0, 4.31D0, &
      &     5.30D0, 4.20D0, 3.91D0, 3.80D0, 3.75D0, 3.70D0, 3.65D0, 3.60D0, &
      &     3.55D0, 3.52D0, 3.61D0, 3.67D0, 3.70D0, 3.73D0, 3.75D0, 3.56D0, &
      &     3.44D0, 3.23D0, 3.04D0, 2.93D0, 2.86D0, 2.82D0, 2.83D0, 2.88D0, &
      &     2.98D0, 3.27D0, 3.57D0, 3.62D0, 3.37D0, 3.46D0, 3.63D0, 4.44D0, &
      &     5.81D0, 4.30D0, 3.84D0, 3.52D0, 3.32D0, 3.13D0, 3.02D0, 2.96D0, &
      &     2.93D0, 2.93D0, 2.95D0, 2.99D0, 3.05D0, 3.17D0/
    !
    IZ = IDNINT(Z)
    IF (IZ .LT. 0 .OR. IZ .GT. 100) THEN
        WRITE (MESSG, 99001) IZ
        CALL ERRMSG(MESSG, 0)
        WSR = CPWSR(100)
    ELSE
        WSR = CPWSR(IZ)
    END IF
99001 FORMAT('DEFWSR: bad nuclear charge, Z=', i3)
END
