      subroutine find_symmetry(n_operations, operations, slen,
     >   spacegrp, cell, angles, latvec,
     >   n, cpositions, natoms, types, positions, align, Magnetic,
     >   verbose)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      INTEGER :: n_operations
      INTEGER :: operations(64, 2)
c     spacegroup, either number or identifier, see the symdata.f
      integer slen
      byte, dimension(slen) :: spacegrp
c     primitive cell lengths
      double precision :: cell(3)
c     angles
      double precision :: angles(3)
c     cell vectors
      double precision :: latvec(3,3)
c     number of the atoms (only the equivalent)
      integer :: n
c     positions of the atoms ...lattice vectors
      double precision :: cpositions(3, n)
c     ...cartesian
      integer natoms
      double precision :: positions(3, natoms)
c     types
      integer types(natoms)
c     correct the atom positions by aligning them to a grid
      integer :: align
c     magnetic field direction
      double precision :: Magnetic(3)
c     do an output to stdout
      integer :: verbose

      integer num

      parameter (MAXDIM=1 000 000)
      integer W(MAXDIM)
      include 'param.fi'
C      include 'W.H'
      DIMENSION G(3,3),IN(64,64),it(65)
     >  ,GEN(3,3,64),ANG(3),v(3,64)
      LOGICAL EQW,exist
      character buf*78
      COMMON /TAU/ boa,coa,alat,alfa,beta,gamma,nsort
      COMMON /output/ output
      integer output
      call wkinit(MAXDIM)

      if (verbose .ne. 0 ) then
          output = 6
      else
#ifdef _WIN32
          open(unit=20, file='NUL', action='write')
#else
          open(unit=20, file='/dev/null', action='write')
#endif
          output = 20
      end if
c
c Definitions of all possible turn matrixes.
c
      DO 20 IG=1,64
        IGA=IG
        IF(IGA.GT.32)IGA=IGA-32
        CALL EILANG(ANG,IGA)
        CALL TURNM(ANG,GEN(1,1,IG))
C                              !DEF TURN MATRIX
        IF(IG.NE.IGA)THEN
C                              !ADD INVERSION IF NEED
          DO I=1,3
           DO J=1,3
              GEN(I,J,IG)=-GEN(I,J,IG)
           END DO
          END DO
        ENDIF
20    CONTINUE

      call defrr(itau,3*NSMAX)
      call multable(in)
      !READ(11,*)BOA,COA,alfa,beta,gamma
!     call readin(num,w(itau),buf)

      buf=""
      buf = transfer(spacegrp(:MIN(78,slen)), buf(:MIN(78,slen)))

      boa = cell(2) / cell(1)
      boa = cell(3) / cell(1)
      alpha = angles(1)
      beta = angles(2)
      gamma = angles(3)
      nsort = n

      call init_in(num,w(itau),buf, cpositions, align .ne. 0)
      if(num.eq.0)then
        call read_dst(latvec, natoms, positions, types,
     $                magnetic, W,GEN, n_operations, operations)
      else
        call genvec
        call gener(W,in,it,n,gen,v,inv,iad,num,buf)
        CALL QAT(MAGNETIC, W,GEN,V,IT,N,inv,iad,w(itau),
     >           n_operations, operations)
        call rlse(itau)
      endif

      end subroutine


      subroutine init_in(num,tau,buf, positions, align)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer num
      double precision :: positions(3, nsort)
      logical :: align
      include 'param.fi'
      character*2 F*1,ac*9,bc*9,cc*9
     1,sc*6 ,lat(7),lattic*8,buf*78,TMP*10
      integer numa
      DIMENSION  A(3),B(3),C(3)
      DIMENSION  A0(3),B0(3),C0(3)
      dimension TAU(3,NSMAX)
      common /u/u(3,3),a,b,c,um1(3,3),IDATASTR
      common /lat/ lattic
      COMMON /TAU/ boa,coa,alat,alfa,beta,gamma,nsort
      data lat/'C ','H ','T ','TG','O ','M ','TC'/
      DATA EPS/1.d-5/
      PI=DPI()
      alat=1d0
c 1 - cubic
c 2 - hexagonal
c 3 - tetragonal
c 4 - trigonal
c 5 - orthorombic
c 6 - monoclinic
c 7 - triclinic

      IDATASTR=1
      if(ichar(buf(1:1)).lt.ichar('9'))then
        read(buf,*) NUMA
        NUM=numa
        if(num.eq.0)return
        iorig=0
        if(numa.gt.1000.and.numa.lt.1231)then
          NUM=numa-1000
          iorig=1
        endif
        call spacegroup(buf,num,iorig)
      else
        iorig=1
        !spacegroup(11:)=' '
        !do i=1,10
        !  ic=ichar(spacegroup(i:i))
        !  if(ic.lt.32.or.ic.gt.122)goto 201 ! if not digit or letter
        !  if(ic.ge.97) buf(i:i)=char(ic-32) ! to upper case
        !  num=-i
        !enddo
        ! 201    !continue
        call spacegroup(buf,num,iorig)
      endif
c      iorig=2
c      if(buf(4:4).eq.'*')iorig=1
      write(ilun(11),*)NUM,' / International group :',
     $     buf(5:13),' ORIGIN SET:',iorig
      write(ilun(6), *) NUM,' ',buf(5:13)
      numa=num
      parlat=1.D0
      alfa=0
      beta=90
      gamma=90
      F=buf(5:5)
      itrig=0
      READ(buf,16)NUMG,INVERSION,IADDPARAM,Icub
      !READ(11,*)BOA,COA,alfa,beta,gamma
      if(boa.lt.0.d0)boa=abs(boa/alat)
      if(coa.lt.0.d0)coa=abs(coa/alat)
      IF( BOA.eq.0.or.icub.lt.4) BOA = 1 !***************************
      IF( COA.eq.0.or.icub.eq.1) COA = 1
      if(icub.eq.-4)then
c      if(F.eq.'R'.and.
c     $     abs(COA-1d0).lt.1.e-8.and.abs(alfa).gt.1.e-10)then
        write(ilun(6), *) 'Trigonal axis'
        itrig=1
        ca=cos(PI/180*alfa)
        coa=sqrt((3*ca+1.5d0)/(1d0-ca))
        alfa=90
        parlat=1d0/sqrt(1d0/3+coa**2/9)
        F='r'
        buf(5:5)='r'
c        if(num.eq.161)buf(47:)='666'
c        if(num.eq.167)buf(41:43)='666'
        icub=4
      endif
      if(alfa.eq.0.or.icub.ne.7) alfa=90
      if(beta.eq.0.or.icub.lt.6) beta=90
      if(gamma.eq.0.)gamma=90
      if(icub.eq.4.or.icub.eq.2)then
        gamma=120                       ! hexagonal axis (hexagonal,trigonal)
      else if(icub.lt.6)then
        gamma=90                        ! cubic,tetragonal,orthorombic!
      endif
      if(icub.eq.6)then
c
c    Idiotic treatment of the uniq axis
c
        if(buf(39:40).eq.'02'.or.buf(39:40).eq.'34'
     $       .or.buf(45:46).eq.'34')then
c  UNIQ AXIS B
          write(ilun(6),*),'UNIQ AXIS B'
          if(abs(beta-90.).lt.0.001)then
            print*,' ***** WARNING ******'
            print*,' You did not put correct beta'
            print*,' Data for gamma will be accepted as beta:',gamma
            beta=gamma
            gamma=90
          endif
        else
c  UNIQ AXIS C
          write(ilun(6),*),'UNIQ AXIS C'
          beta=90
        endif
      endif
cccccc1     format(3f21.15)  
1     format(3f3.1)   
      bkoef=1
      if(F.eq.'P')then
        aC='1. 0. 0. '
        bC='0. 1. 0. '
        cC='0. 0. 1. '
        SC='S*'
      else if(F.eq.'I')then
        aC='-.5 .5 .5'
        bC=' .5-.5 .5'
        cC=' .5 .5-.5'
        SC='BC*'
      else if(F.eq.'F')then
        aC=' .0 .5 .5'
        bC=' .5 .0 .5'
        cC=' .5 .5 .0'
        SC='FC*'
      else if(F.eq.'C')then             !.or.F.eq.'A')then
        aC=' .5 .5 .0'
        bC='-.5 .5 .0'
        cC='0. 0. 1. '
        SC='CFC*'
      else if(F.eq.'A')then
        aC='1. 0. 0. '
        bC='0.  .5-.5'
        cC='0.  .5 .5'
        SC='CFC*'                       !'AFC*'
c$$$      else if(F.eq.'B')then
c$$$        aC='1. 0. 0. '
c$$$        bC='0.  .5-.5'
c$$$        cC='0.  .5 .5'
c$$$        SC='CFC*'
      else if(F.eq.'B')then
        aC='.5 0. -.5'
        bC='0. 1. 0.0'
        cC='.5 0.  .5'
        SC='CFC*'                       !'BFC*'

c$$$      else if(F.eq.'C'.or.F.eq.'A')then
c$$$        aC=' .5 .5 .0'
c$$$        bC='-.5 .5 .0'
c$$$        cC='0. 0. 1. '
c$$$        SC='CFC*'
c$$$      else if(F.eq.'B')then
c$$$        aC='1. 0. 0. '
c$$$        bC='0.  .5-.5'
c$$$        cC='0.  .5 .5'
c$$$        SC='CFC*'
      ELSE if(F.eq.'R')then
        aC=' 2. 1. 1.'
        bC='-1. 1. 1.'
        cC='-1.-2. 1.'
c        aC=' 1. 0. 1.'
c        bC=' 0. 1. 1.'
c        cC='-1.-1. 1.'
         bkoef=1.d0/3
        SC='R*'
      ELSE if(F.eq.'r')then
        aC='1. 0. 0. '
        bC='0. 1. 0. '
        cC='0. 0. 1. '
        bkoef=1.d0!/3
        SC='R*'
      else
        call endjob(10,'error in GROUP NAME')
      endif
      if(numa.ne.0)then
        read(aC,1)a
        read(bC,1)b
        read(cC,1)c
        do i=1,3
          a0(i)=a(i)*bkoef
          b0(i)=b(i)*bkoef
          c0(i)=c(i)*bkoef
        enddo  
      endif
 16   format(i3,24x,i1,i3,i2,i4,3(i3,i1,i1,i1))
      lattic=sc
      lattic(index(lattic,'*'):)=lat(icub)
      write(ilun(16),*)' LAT:',lattic
c      REWIND(14)
c      write(14)LATTIC,buf
c!
c!      matrix:
c!
        write(ilun(16),*)'alpha,beta,gamma:',alfa,beta,gamma!,boa,coa

        alfa=alfa*PI/180
        beta=beta*PI/180
        gamma=gamma*PI/180
        if(icub.eq.1.or.icub.eq.3.or.icub.eq.5.or.icub.eq.6.
     &    or.icub.eq.7.or.numa.eq.0)then
          U(1,1)=1
          U(1,2)=0
          U(1,3)=0
          U(2,1)=boa*cos(gamma)
          U(2,2)=boa*sin(gamma)
          U(2,3)=0
          U(3,1)=coa*(cos(beta)-cos(alfa)*cos(gamma))/sin(gamma)
          U(3,2)=-coa*cos(alfa)
          U(3,3)=sqrt(coa**2-u(3,1)**2-u(3,2)**2)
c$$$
c$$$          U(1,1)=1
c$$$          U(1,2)=0
c$$$          U(1,3)=0
c$$$          U(2,1)=boa*cos(gamma)
c$$$          U(2,2)=boa*sin(gamma)
c$$$          U(2,3)=0
c$$$          U(3,1)=0
c$$$          U(3,2)=0
c$$$          U(3,3)=coa
        elseif (F.eq.'r')then
c          coa=1-ca/(ca+0.5)

          U(1,1)=sqrt(3d0)/6*parlat
          U(2,1)=sqrt(3d0)/6*parlat
          U(3,1)=-sqrt(3d0)/3*parlat
          U(1,2)=-0.5*parlat
          U(2,2)=0.5*parlat
          U(3,2)=0.
          U(1,3)=1d0/3*coa*parlat
          U(2,3)=1d0/3*coa*parlat
          U(3,3)=1d0/3*coa*parlat
        else

          U(2,1)=0
          U(2,2)=1
          U(2,3)=0
          U(1,1)=sin(gamma)
          U(1,2)=cos(gamma)
          U(1,3)=0
          U(3,1)=0
          U(3,2)=0
          U(3,3)=coa

c$$$          U(1,1)=0
c$$$          U(1,2)=-boa
c$$$          U(1,3)=0
c$$$          U(2,1)=sin(gamma)
c$$$          U(2,2)=-cos(gamma)
c$$$          U(2,3)=0
c$$$          U(3,1)=coa*(cos(beta)-cos(alfa)*cos(gamma))/sin(gamma)
c$$$          U(3,2)=-coa*cos(alfa)
c$$$          U(3,3)=sqrt(coa**2-u(3,1)**2-u(3,2)**2)
        endif
        write(ilun(16),*)'    U matrix                  U^(-1)'
        call DINV33(U,0,UM1,DETRM)
        do i=1,3
          write(ilun(16),48)
     $         (U(i,j),j=1,3),(UM1(i,j),j=1,3)
 48       format(3f21.15, ' : ',3f21.15)
        enddo
c! 
c!
c!
      do k=1,3
        a(k)=0
        b(k)=0
        c(k)=0
        do j=1,3
          a(k)=a(k)+a0(j)*u(j,k)!*parlat
          b(k)=b(k)+b0(j)*u(j,k)!*parlat
          c(k)=c(k)+c0(j)*u(j,k)!*parlat
        enddo
      enddo
      write(ilun(16),*)' Primitive lattice vectors:'
      write(ilun(16),11)'A',a,a0
      write(ilun(16),11)'B',b,b0
      write(ilun(16),11)'C',c,c0
11     format(' ',a,' ',3f21.15,' <--> ',3f21.15)
!       write(ilun(11),17)a,b,c
 17    format(3f21.15,' / A in lat.par'/3f21.15,' / B in lat.par'/
     $      3f21.15,' / C in lat.par')
c      write(14)a(1),a(2)/boa,a(3)/coa,b(1),b(2)/boa,b(3)/coa,
c     1 c(1),c(2)/boa,c(3)/coa,u,um1

c
c Atomic positions:
c

c      READ(11,*)nsort
c      READ(11,*)
      DO  I=1,nsort
c         READ(11,*)vv
        if( align )then
          do k=1,3
            DumA=positions(k,i)*24
            IDumA=nint(DumA)
            if(abs(DumA-IDumA).lt.1.e-3) then
              positions(k,i)=IDumA/24.d0
            end if
          enddo
        end if
        do k=1,3
          tau(k,i)=0
            do j=1,3
              tau(k,i)=tau(k,i)+positions(j, i)*u(j,k)
            enddo
        enddo
      enddo
      end
      subroutine gener(W,m,it,n,g,v,inversion,IADDPARAM,num,buf)
      IMPLICIT DOUBLE PRECISION
     $ (A-H,O-Z)
C      include 'W.H'
        INTEGER W(*)
      include 'param.fi'
      character*2 op*64,ti*3,filename*11,file*60,F*1,ac*9,bc*9,cc*9
     1, FILEN*100,sc*6 ,lat(7),lattic*8,buf*78,names*24
      dimension m(64,64),G(3,3,64),v(3,64),ivv(3,3),vv(3)
      logical hex
      DIMENSION  A(3),B(3),C(3)
      INTEGER IT(65),P
      common /lat/ lattic
      COMMON /TAU/ boa,coa,alat,alfa,beta,gamma,nsort
      common /u/u(3,3),a,b,c,UM1(3,3),IDATASTR
      data lat/'C ','H ','T ','TG','O ','M ','TC'/
      DATA EPS/1.d-5/
c     call spacegroup(buf,num)
      do i=1,65
         it(i)=-100
      enddo
      READ(buf,11)NUMG,INVERSION,IADDPARAM,Icub,n,
     1 (it(i),(ivv(j,i),j=1,3),i=1,n)
c      print*,NUMG,INVERSION,IADDPARAM,Icub,n,
c     1 (it(i),(ivv(j,i),j=1,3),i=1,n)
      write(ilun(6),*) buf

11    format(i3,24x,i1,i3,i2,i4,3(i3,i1,i1,i1))
      do i=1,n
        do k=1,3
          v(k,it(i))=0
          do j=1,3
            v(k,it(i))=v(k,it(i))+ivv(j,i)*u(j,k)/12!*parlat
          enddo
        enddo
      enddo
      N1=N
 5    n=n1
      DO I=1,N
        do j=1,n
          DO K=1,N1
            IF(M(it(i),it(j)).EQ.IT(K))GOTO 10
          ENDDO
          N1=N1+1
          IT(N1)=m(it(i),it(j))
          call vec(v(1,it(n1)),v(1,it(j)),v(1,it(i)),g(1,1,it(j)) )
c          print*,IT(N1),(v(kk,it(n1)),kk=1,3)
c          print*,'New:',n1,':',it(i),'*',it(j),'=',it(n1)
c          print*,'(1):',v(1,it(i)),v(2,it(i)),v(3,it(i))
c          print*,'(2):',v(1,it(j)),v(2,it(j)),v(3,it(j))
c          print*,'(=):',v(1,it(n1)),v(2,it(n1)),v(3,it(n1))
          do l=1,3
            do mw=1,3
              glm=0
              do k=1,3
                glm=glm+g(l,k,it(i))*g(k,mw,it(j))
              enddo
              if(abs(glm-g(l,mw,it(n1))).gt.1.e-5)then
                write(ilun(6),*) i,j,n1,'(',it(i),it(j),it(n1),')'
                write(ilun(6),*) l,mw,glm,g(l,mw,it(n1))
              endif
            enddo
          enddo
 10       CONTINUE
        enddo
      ENDDO
      IF(N1.NE.N)GOTO5
      DO 300 II=2,N
        I=II-1
        K=I
        P=IT(I)
C     
        DO 260 J=II,N
          IF(IT(J).GE.P) GO TO 260
          K=J
          P=IT(J)
 260    CONTINUE
C     
        IF(K.EQ.I) GO TO 300
        IT(K)=IT(I)
        IT(I)=P
C     
 300  CONTINUE
      do i=1,64
        op(i:i)='0'
      enddo
      do i=1,n
        op(it(i):it(i))='1'
      enddo
c      call PNTGRP(0,IPG,N,TI,OP)
c      write(ilun(16),*) 'Point group is :', Ti,' N=',N
      icb=0
      if(numg.gt.194)icb=1 
      write(ilun(16),*)' Symmetry operations:',n!icub
      do i=1,n
          do k=1,3
             vv(k)=0
            do j=1,3
              vv(k)=vv(k)+v(j,it(i))*um1(j,k)
            enddo
          enddo
          call namesym(it(i),icub,names,vv)
          write(ilun(16),116)it(i),names!vv
      enddo
 116  format(i3,' :   ',a,2x,3f21.15)
      END
c
      subroutine vec(v0,v1,v2,g1)
      IMPLICIT DOUBLE PRECISION
     $ (A-H,O-Z)
      dimension v0(3),v1(3),v2(3),g1(3,3)
      do i=1,3
        v0(i)=v1(i)
        do j=1,3
          v0(i)=v0(i)+v2(j)*g1(j,i)
        enddo
      enddo
      CALL SHORTN(V0,V0)
      end
      SUBROUTINE genvec
      IMPLICIT DOUBLE PRECISION
     $ (A-H,O-Z)
       INCLUDE 'param.fi'
      INTEGER DTMAX
C
      DIMENSION DLAT(3,NVCMAX0),D(NVCMAX0)
      DIMENSION  A(3),B(3),C(3)
      COMMON /VG/ DLAT
      common /u/ u(9),a,b,c,UM1(3,3),IDATASTR
C
C
C*** GENERATE STR. VECTORS OF THE DIRECT AND RECIPR. SPACE ***
C
      I=0
      DO III=-3,3
        DO JJJ=-3,3
          DO KKK=-3,3
            I=I+1
            DLAT(1,I)=III*A(1)+JJJ*B(1)+KKK*C(1)
            DLAT(2,I)=III*A(2)+JJJ*B(2)+KKK*C(2)
            DLAT(3,I)=III*A(3)+JJJ*B(3)+KKK*C(3)
            d(I)=DLAT(1,I)**2+DLAT(2,I)**2+DLAT(3,I)**2
          ENDDO
        ENDDO
      ENDDO
      N=I
      DO 300 II=2,N
      I=II-1
      K=I
      P=D(I)
C
      DO 260 J=II,N
        IF(D(J).GE.P) GO TO 260
        K=J
        P=D(J)
  260 CONTINUE
C
      IF(K.EQ.I) GO TO 300
      D(K)=D(I)
      D(I)=P
C
      DO 280 J=1,3
        P=DLAT(J,I)
        DLAT(J,I)=DLAT(J,K)
        DLAT(J,K)=P
  280 CONTINUE
C
  300 CONTINUE

      END

c
      SUBROUTINE QAT(Magnetic, W,Gm,VC,IT,NG,inversion,iadop,tau,
     $ ng_all, operations)
      IMPLICIT DOUBLE PRECISION
     $ (A-H,O-Z)
      DOUBLE PRECISION :: Magnetic(3)
C      INCLUDE 'W.H'
      INTEGER W(*)
      INCLUDE 'param.fi'
c      parameter(a00=0.529177d0)
      INTEGER :: ng_all
      INTEGER, DIMENSION(64, 2) :: operations
      parameter(max2sort=1000)
      parameter(EPSMIS=1d-6)
      character *80 tnam(max2sort)
      DIMENSION it(65),Gm(3,3,64),VC(3,64),ioldnum(max2sort)
      DIMENSION TAU(3,NSMAX),VV(3)
      COMMON /TAU/ boa,coa,alat,alfa,beta,gamma,NSORT
      COMMON /VG/ D(3,NVCMAX0)
c      dimension  TAUX(NAMAX),TAUY(NAMAX),
c     1           TAUZ(NAMAX),ISS(NAMAX)
      common /lat/ lattic
      common /u/ u(9),a0(3),b0(3),c0(3),UM1(3,3),IDATASTR
      DIMENSION VE(3,48),ANG(3),f(3),g(3),h(3)
      character*72 a,b,yes,chem*8,name*8,struc*8,lattic*8
      LOGICAL EQW
      nnn=0
c
c reciprocal vectors
c
      do i=1,max2sort
        tnam(i)=' '
        ioldnum(i)=i
      enddo
      write(ilun(16),*) 'Length of  basic vectors:'

      a0l=sqrt(a0(1)*a0(1)+a0(2)*a0(2)+a0(3)*a0(3))
      b0l=sqrt(b0(1)*b0(1)+b0(2)*b0(2)+b0(3)*b0(3))
      c0l=sqrt(c0(1)*c0(1)+c0(2)*c0(2)+c0(3)*c0(3))
      write(ilun(16),*) '  A          B         C  '
      write(ilun(16),12)a0l,b0l,c0l
      write(ilun(16),*) 'Angles between basic vectors:'
      write(ilun(16),*) '  AB         BC        CA '
      write(ilun(16),12)
     $ acos((a0(1)*b0(1)+a0(2)*b0(2)+a0(3)*b0(3))/(a0l*b0l))*180/DPI(),
     $ acos((c0(1)*b0(1)+c0(2)*b0(2)+c0(3)*b0(3))/(c0l*b0l))*180/DPI(),
     $ acos((a0(1)*c0(1)+a0(2)*c0(2)+a0(3)*c0(3))/(a0l*c0l))*180/DPI()
      CALL CROSS(F,B0,C0)
      CALL CROSS(G,C0,A0)
      CALL CROSS(H,A0,B0)
      WL=A0(1)*F(1)+A0(2)*F(2)+A0(3)*F(3)
      W1=1.d0/WL
      OMEGA=ABS(WL)
C     
      DO  I=1,3
        F(I)=F(I)*W1
        G(I)=G(I)*W1
        H(I)=H(I)*W1
      enddo

      open(3,status='scratch',form='unformatted')
      DO IS=1,NSORT
        eps_angmax=0.d0
        IN=1
        VE(1,1)=TAU(1,IS)
        VE(2,1)=TAU(2,IS)
        VE(3,1)=TAU(3,IS)
        DO IG=1,NG
          CALL VEC(VV,VC(1,IT(IG)),TAU(1,IS),Gm(1,1,IT(IG)))
          DO IIN=1,IN
            ang(1)=vE(1,iIN)-VV(1)
            ang(2)=vE(2,iIN)-VV(2)
            ang(3)=vE(3,iIN)-VV(3)
            call shortn(ang,ang)
            eps_ang=ang(1)**2+ang(2)**2+ang(3)**2
            if(eps_ang.Lt.EPSMIS)then
              eps_angmax=max(eps_ang,eps_angmax)
              GOTO 800
            endif
          ENDDO

          IN=IN+1
          VE(1,IN)=VV(1)
          VE(2,IN)=VV(2)
          VE(3,IN)=VV(3)
800       CONTINUE
        ENDDO
        nnn=nnn+in
        write(ilun(16),*) 'all:',IN, '  Max. mismatch:',sqrt(eps_angmax)
        write(3)in
        DO I=1,IN
          iinv=0
c$$$          if(inversion.eq.0)then
c$$$            DO IIN=1,IN
c$$$              ang(1)=vE(1,iIN)+vE(1,i)
c$$$              ang(2)=vE(2,iIN)+vE(2,i)
c$$$              ang(3)=vE(3,iIN)+vE(3,i)
c$$$              call shortn(ang,ang)
c$$$              if(ang(1)**2+ang(2)**2+ang(3)**2.Lt.1.e-5)then
c$$$                iinv=iin
c$$$                GOTO 801
c$$$              endif
c$$$            ENDDO
c$$$801         continue
c$$$            if(i.ne.iinv)then
c$$$c!05-19-94 Perlov - clearifying of INVERSION:
c$$$              do jj=1,3
c$$$                ve(jj,iinv)=-ve(jj,i)
c$$$              enddo
c$$$            endif
c$$$            if(iinv.eq.0)call endjob(10,'inv error')
c$$$          endif
          do k=1,3
            vv(k)=0
            do j=1,3
              vv(k)=vv(k)+ve(j,i)*um1(j,k)
            enddo
          enddo

c          write(ilun(16),12) (VE(JJ,I),JJ=1,3),IS,iinv,
          write(ilun(16),12) vv,IS,iinv,
     $         VE(1,I)*f(1)+VE(2,I)*f(2)+VE(3,I)*f(3),
     $         VE(1,I)*g(1)+VE(2,I)*g(2)+VE(3,I)*g(3),
     $         VE(1,I)*h(1)+VE(2,I)*h(2)+VE(3,I)*h(3)
          write(3) (VE(JJ,I),JJ=1,3),IS,iinv
        ENDDO
      ENDDO
12    FORMAT(3F21.15,I3,i2,' :',3f21.15)
      rewind(3)
c      write(14)nnn,nsort
c      if(nnn.gt.NAMAX)then
c        print*,' MAX number of atoms:',NAMAX,' Here:',nnn
c        call endjob(10,'A lot of atoms ')
c      endif
      NAMAX=nnn
c      write(14) alat,boa,coa,alfa,beta,gamma

      call defrr(itaux,NAMAX)
      call defrr(itauy,NAMAX)
      call defrr(itauz,NAMAX)
      call defi(iss,NAMAX)
      write(ilun(16),*)' Atomic positions in lat.par: '
      write(ilun(16),*)NAMAX,'    -   NATOMS'
      call defrr(i_alcs,NAMAX*3)
      CALL DEFTAU(Magnetic, W(iTAUX),W(iTAUY),W(iTAUZ),W(ISS),NSORT,
     $W(i_alcs))
      close(unit=3,status='delete')
      call defi(i_KK,NAMAX*NG)
      call KTO_KYDA(W(iTAUX),W(iTAUY),W(iTAUZ),
     $     Gm,It,VC,namax,ng,W(i_KK))
      NATOM=nnn
      call defi(i_iatpos,natom)
      call defi(i_newis,natom)
      call defi(i_nwt,natom)
      call defi(i_ntw,natom)
      call defi(i_isnew,natom)
      call defi(i_itg,64)
      call check_sym(NG,NG,NATOM,w(i_kk),w(i_alcs),w(i_isnew)
     $       ,it,nsort,w(i_iatpos),w(i_newis),w(iss),w(i_nwt),w(i_ntw),
     $     w(i_itg),ng_all,nsortnew,ier,W(iTAUX),W(iTAUY),W(iTAUZ),Gm
     $     ,tnam,ioldnum, operations)
      END

      SUBROUTINE DEFTAU(Magnetic, TAUX,TAUY,TAUZ,ISS,NSORT,euler)
      IMPLICIT DOUBLE PRECISION
     $ (A-H,O-Z)
      double precision, dimension(3) :: magnetic
      dimension taux(*)     ! - x - sites of atoms
      dimension tauy(*)     ! - y - sites of atoms
      dimension tauz(*)     ! - z - sites of atoms
      dimension isS(*)      ! - number of sort of the atom
      dimension euler(3,*)  ! - Euler angles of field for all atoms
      COMMON /TAU/ boa,coa,alat,alfa,beta,gamma,nsrt
      ijatom=0
      Hx=0
      Hy=0
      Hz=0
      !read(11,*)HX,HY,HZ
      Hx=Magnetic(1)
      Hy=Magnetic(2)
      Hz=Magnetic(3)
      if(Hx**2+Hy**2+Hz**2.gt.1.d-7)then
        write(ilun(11),17)HX,HY,HZ
 17     format(/3f21.15,' / Magnetic field direction'/)
        anorm=sqrt(hx*hx+hy*hy+hz*hz)
        beth=acos(hz/anorm)
        if(hx*hy.lt.1.d-10)then
          alfh=0.d0
        else
          alfh=atan2(hy,hx)
        endif
        gamh=0.d0
      else
        alfh=0d0
        beth=0d0
        gamh=150.d0
      endif
      do in=1,nsort
        ijsort=ijatom
        read(3)na
        DO IATOM = 1 , NA
         ijatom=ijatom+1
          read(3) TAUXx, TAUYx, TAUZx, IS,iinv
          write(ilun(16),12) TAUXx, TAUYx, TAUZx, IS,ijatom
 12       FORMAT(3F21.15,I3,' : ',i3)
c          write(14) TAUXx, TAUYx/boa, TAUZx/coa, IS,iinv+ijsort
          taux(ijatom)=tauxx
          tauy(ijatom)=tauyx
          tauz(ijatom)=tauzx
          iss(ijatom)=is
          euler(1,ijatom)=alfh
          euler(2,ijatom)=beth
          euler(3,ijatom)=gamh
        ENDdo
      ENDdo
      END

      SUBROUTINE SHORTN(P,P1)
      IMPLICIT DOUBLE PRECISION
     $ (A-H,P-Z), INTEGER(O)
      INCLUDE 'param.fi'
      COMMON /VG/ DLAT(3,NVCMAX0)
      DIMENSION P(3),P1(3)
      ANRM2(X,Y,Z)=X*X+Y*Y+Z*Z
      NKD=NVCMAX0
      P1(1)=P(1)
      P1(2)=P(2)
      P1(3)=P(3)
      CRIT0=10000.d0
      DO 10 K=1,NKD
        CRIT=ANRM2(P1(1)+DLAT(1,K),P1(2)+DLAT(2,K),P1(3)+DLAT(3,K))
        IF(CRIT.LT.CRIT0)THEN
        CRIT0=CRIT
        JS=K
      ENDIF
10    CONTINUE
      P1(1)=P1(1)+DLAT(1,JS)
      P1(2)=P1(2)+DLAT(2,JS)
      P1(3)=P1(3)+DLAT(3,JS)
      END


c

      SUBROUTINE RESORTI3(N,IT,ia1,ia2)
      DOUBLE PRECISION
     $ IT(N),P
      dimension ia1(*),ia2(*)
      DO II=2,N
        I=II-1
        K=I
        P=IT(I)
C     
        DO J=II,N
          IF(IT(J).LT.P) THEN
            K=J
            P=IT(J)
          ENDIF
        ENDDO
C     
        IF(K.NE.I)THEN
          ia1s=ia1(k)
          ia2s=ia2(k)
          ia1(k)=ia1(i)
          ia2(k)=ia2(i)
          IT(K)=IT(I)
          ia1(i)=ia1s
          ia2(i)=ia2s
          IT(I)=P
        ENDIF
C     
      ENDDO
      END

      subroutine KTO_KYDA(taux,tauy,tauz,Gm,It,VC,natom,ng,KK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter (EPSMIS=1D-6)
      dimension taux(*)                 ! - x - sites of atoms
      dimension tauy(*)                 ! - y - sites of atoms
      dimension tauz(*)                 ! - z - sites of atoms
      DIMENSION it(65),Gm(3,3,64),VC(3,64),KK(ng,natom)
     $     ,ang(3),vv(3),vat(3)
      write(ilun(16),*)' Transformations:'
      eps_angmax=0.d0
      DO IG=1,NG
        do iatom=1,natom
          vat(1)=taux(iatom)
          vat(2)=tauy(iatom)
          vat(3)=tauz(iatom)
          CALL VEC(VV,VC(1,IT(IG)),vat,Gm(1,1,IT(IG)))
          DO JATOM=1,NATOM
            ang(1)=taux(Jatom)-VV(1)
            ang(2)=tauy(Jatom)-VV(2)
            ang(3)=tauz(Jatom)-VV(3)
            call shortn(ang,ang)
            eps_ang=ang(1)**2+ang(2)**2+ang(3)**2
c           write(*,*) 'EPSMIS',EPSMIS
            if(eps_ang.Lt.EPSMIS)then
              eps_angmax=max(eps_ang,eps_angmax)
              GOTO 800
            endif
          ENDDO
          call endjob(10,' Strange... very strange')
 800      KK(IG,iatom)=jatom
        enddo
        write(ilun(16),11)IT(iG),(KK(IG,iatom),iatom=1,natom)
      ENDDO
c      write(14)(IT(iG),(KK(IG,iatom),iatom=1,natom),IG=1,ng)
 11   format(I2,':',48I3)
      end




      subroutine check_sym(NGD,NG,NATOM,kto_kyda,a,is,
     $     it,nsort0,iat_pos,newis,is0,ntw,nwt,
     $     itg,NG_ALL,NSORT,ier,TAUX,TAUY,TAUZ,Gm,tnam,ioldnum,
     $     operations)
c
C      INCLUDE 'PREC.H'
        implicit double precision (a-h,o-z)
        integer :: operations(64, 2)

      dimension taux(*)     ! - x - sites of atoms
      dimension tauy(*)     ! - y - sites of atoms
      dimension tauz(*)     ! - z - sites of atoms
      character *80 tnam(*)

c
      dimension it(64),ang(3),h0(3),h(3),a(3,*),KTO_KYDA(NGD,NATOM),
     $     is(*),gm(3,3,64),itg(64),iat_pos(*),newis(*),
     $     is0(*),nwt(*),ntw(*),g(3,3),ioldnum(*)
      character*1 znak
      character*20 angs
c
ci nsort0 - old number of sorts
ci is0    - old sort for each atom
c     
      hx(ia)=sin(a(2,ia))*cos(a(1,ia)-a(3,ia))
      hy(ia)=sin(a(2,ia))*sin(a(1,ia)-a(3,ia))
      hz(ia)=cos(a(2,ia))
c

c$$$      print*,'INPUT:',nsort0
c$$$      print*,'IT:',it
c$$$      print*,'NGD,NG,NATOM',NGD,NG,NATOM
c$$$      do i=1,natom
c$$$        print*,i,is0(i),(a(ii,i),ii=1,3)
c$$$        print*,(kto_kyda(iop,iatom),iop=1,ng)
c$$$      enddo
      iun=ilun(16)
      ier=0
c      read(14)(IT(iG),(KTO_KYDA(IG,iatom),iatom=1,natom),IG=1,ng)
      NG_ALL=NG
      do i=1,64
        itg(i)=i
      enddo
      if(a(3,1).gt.100d0)goto 145
      ign=0
      do i0=1,nsort0
        newis(i0)=0
      enddo
      DO IG=1,NG
        iznak=0
        igg=it(iG)
        if(igg.gt.32)then
          igg=igg-32
        endif
        call eilang(ang,igg)
        call turnm(ang,g)
        DO IA=1,NATOM
          if(newis(is0(ia)).eq.0)then
            sum=0.d0
            h0(1)=hx(ia)        !sin(a(2,ia))*cos(a(3,ia)-a(1,ia))
            h0(2)=hy(ia)        !sin(a(2,ia))*sin(a(3,ia)-a(1,ia))
            h0(3)=hz(ia)        !cos(a(2,ia))
            do  i=1,3
              h(i)=0.d0
              do  j=1,3
                h(i)=h(i)+g(i,j)*h0(j)
              enddo
            enddo
c            if(ia.eq.1)print*,it(IG),'h0:',h0,'h:',h
            j=KTO_KYDA(IG,ia)
            h0(1)=hx(j)         !sin(a(2,j))*cos(a(3,j)-a(1,j))
            h0(2)=hy(j)         !sin(a(2,j))*sin(a(3,j)-a(1,j))
            h0(3)=hz(j)         !cos(a(2,j))
            if (abs(h0(1)-h(1)).lt.1.e-2.and.
     $           abs(h0(2)-h(2)).lt.1.e-2.and.
     $           abs(h0(3)-h(3)).lt.1.e-2.and.iznak.ge.0)then
              iznak=1
            else if (abs(h0(1)+h(1)).lt.1.e-2.and.
     $             abs(h0(2)+h(2)).lt.1.e-2.and.
     $             abs(h0(3)+h(3)).lt.1.e-2.and.iznak.le.0)then
              iznak=-1
            else
c$$$  print*,'IG:',IG,'   ia,j',ia,j
c$$$  print*,hx(ia),hy(ia),hz(ia)
c$$$  print*,h
c$$$  print*,hx(j),hy(j),hz(j)
              it(ig)=0
              NG_ALL=NG_ALL-1
              goto 100
            endif
          endif
        enddo
 100    continue
        it(ig)=it(ig)*iznak
        if(it(ig).ne.0)then
          ign=ign+1
          itg(ign)=ig
        endif
      ENDDO
      
      do iatom=1,natom
        is(iatom)=0
      enddo
      isn=0
      do iatom=1,natom
        if(is(iatom).eq.0)then
          isn=isn+1
          iat_pos(isn)=iatom
          is(iatom)=isn
        endif
        do ig=1,ng
          if(iT(ig).ne.0)then
            do jatom=iatom+1,natom
              if(KTO_KYDA(IG,jatom).eq.iatom)then
                is(jatom)=is(iatom)
              endif
            enddo
          endif
        enddo
      enddo

c$$$        do ist=1,isn
c$$$          print*,'SORTS:',ist,iat_pos(ist)
c$$$        enddo

      write(iun,*)' With these magnetic field directions we are with'
      write(iun,*)NG_ALL,' symmetry operatios instead of ', NG,' : '
      write(iun,13)(it(itg(i)),i=1,NG_ALL)
 13   format(24i3)
c
 145  continue
      if(NG_ALL.ne.NG.and.nsort0.ne.isn)then
        write(iun,*)' Sort redefinitions:'
        write(iun,*)' New number of sorts:',isn
        nsort=isn
        do iatom=1,natom
          write(iun,12)iatom, is(iatom), (a(i,iatom)*180/DPI(),i=1,3)
     $         ,hx(iatom),hy(iatom),hz(iatom)
 12       format(2i3,':',3f21.15,' h:',3f21.15)
        enddo
          
        do i0=1,nsort0
          newis(i0)=0
        enddo
        do isort=1,nsort
          i0=is0(iat_pos(isort))
          newis(i0)=newis(i0)+1
        enddo
c     atoms resortings:
        ia0=0
        do isort=1,nsort
          iatom=iat_pos(isort)
          na=0
          do ig=1,NG
            if(it(ig).ne.0)then
              newatom=KTO_KYDA(IG,iatom)
              do ia=1,na
                if(newatom.eq.nwt(ia+ia0))goto 45
              enddo
              na=na+1
              nwt(na+ia0)=newatom
              ntw(newatom)=na+ia0
 45           continue
            endif           
          enddo          
          ia0=ia0+na
        enddo
        ier=1
c$$$        do i0=1,nsort0
c$$$          print*,'NEW_SORTS:',i0,newis(i0)
c$$$        enddo
c$$$        print*,'ATOMS,all:',na, 'were:',natom
c       do iatom=1,natom
c          print*,'ATOMS:',iatom,nwt(iatom),ntw(iatom)
c       enddo
      else                              ! copy old is0 to is
        nsort=nsort0
        do iatom=1,natom
          is(iatom)=is0(iatom)
        enddo
      endif
      write(ilun(11),'(2I5,3A)')nsort,natom,' / nsort,natom      ',
     $     '    taus         ','          ICL  IQ(inp) CL(inp) '
      do i=1,nsort
        do iatom=1,natom
          if(is(iatom).eq.i)then
            !write 17
            write(ilun(17),17)taux(iatom),tauy(iatom),tauz(iatom),i
     $           ,ioldnum(iatom),tnam(is0(iatom))!,is0(iatom)
          endif
        enddo
 17     format(3f21.15,i6,2x,i6,'  ',a12)
c#'  olst #:',i4,'-',i4)
      enddo
      !write 17
      write(ilun(17),*)NG_ALL,
     $     ' / Point symmetry operations for reciprocal sums'
      write(ilun(17),22)(it(itg(i)),i=1,ng_all)
      write(ilun(17),22) (KTO_KYDA(itg(I),1),i=1,ng_all)
      do i=1, ng_all
          operations(i,1) = it(itg(i))
          operations(i,2) = KTO_KYDA(itg(I),1)
      end do
 22   format(24i3)
      do i=1,NG_ALL
       iop=it(itg(i))
c$$$       write(17,*)iop,' / #OP'
c$$$       if(iop.gt.0)then
c$$$         write(17,177)((Gm(i1,i2,iop),i1=1,3),i2=1,3)
c$$$       else
c$$$         write(17,177)((-Gm(i1,i2,-iop),i1=1,3),i2=1,3)
c$$$       endif
       znak='+'
       if(iop.lt.0)then
         iop=-iop
         znak='-'
       endif
       if(iop.gt.32)then
         iop=iop-32
         if(znak.eq.'-')then
           znak='+'
         else
           znak='-'
         endif
       endif
       call EILNAM(angs,Iop)
       !write 17
       write(ilun(17),*)znak,' ',angs,' / ',
     >       it(itg(i)),KTO_KYDA(itg(I),1)
 177   format(3f21.15)
      enddo
      end

      subroutine namesym(ig,icub,nam,vv)
      IMPLICIT DOUBLE PRECISION
     $ (A-H,O-Z)
      character*6 frac(12)
      character*12 hex(32),cub(24),name,nam*24,rhomb(9)
      dimension vv(3)
      data cub/
     $      '+x   +y   +z','-x   +y   -z','+x   -y   -z','-x   -y   +z'
     $     ,'-y   +z   -x','-y   -z   +x','+y   -z   -x','+y   +z   +x'
     $     ,'-z   -x   +y','+z   +x   +y','+z   -x   -y','-z   +x   -y'
     $     ,'+y   +x   -z','+y   -x   +z','-y   +x   +z','-y   -x   -z'
     $     ,'+z   -y   +x','-z   -y   -x','-z   +y   +x','+z   +y   -x'
     $     ,'-x   -z   -y','+x   +z   -y','-x   +z   +y','+x   -z   +y'/
      data rhomb/
     $      '+x   +y   +z','-y   -x   -z','+y   +z   +x', '+z   +x  +y',
     $      ' ','-x   -z   -y',' ',' ','-z   -y   -x'/
      data hex/
     $      '+x   +y   +z','-x   +y-x -z','+x   +x-y -z','-x   -y   +z',
     $     20*' '
     $     ,'+y   +y-x +z','+y-x -x   +z','-y   +x-y +z','+x-y +x   +z'
     $     ,'+y   +x   -z','+y-x +y   -z','-y   -x   -z','+x-y -y   -z'/
      data frac
     $     /' ','+1/12','+1/6','+1/4','+1/3','+5/12','+1/2','+7/12',
     $     '+2/3','+3/4','+5/6','+11/12'/
      iq=mod(ig-1,32)+1
      if(icub.eq.-4)then
        if(iq.lt.7)then
          name=rhomb(iq)
        else
          name=rhomb(iq-23)
        endif
      else if(icub.ne.4.and.icub.ne.2)then
        name=cub(iq)
      else
        name=hex(iq)
      endif

      if(ig.gt.32)then
        do i=1,12
          if(name(i:i).eq.'+')then
            name(i:i)='-'
          else if(name(i:i).eq.'-')then
            name(i:i)='+'
          endif
        enddo
      endif

c     print*,'IG,ICUB,name:',ig,icub,iq,name
c      print1,name,vv
 1    format(a,3f21.15,' : ',$)
      ilet=-4
      do i=1,3
      l=(i-1)*8+1
        ilet=ilet+5
        if(vv(i).lt.0)vv(i)=vv(i)+1
        in=nint(vv(i)*12)
        in=mod(in,12)
        do k=ilet,12
          if(name(k:k).eq.' ')exit
          nam(l:l)=name(k:k)
          l=l+1
        enddo
        nam(l:)=frac(in+1)
c        print2,frac(in+1)
 2      format(a,$)
      enddo
c      print*,' '
      end


      subroutine vecn(v0,v1,v2,g1)
      IMPLICIT DOUBLE PRECISION
     $ (A-H,O-Z)
      dimension v0(3),v1(3),v2(3),g1(3,3)
      do i=1,3
        v0(i)=v1(i)
        do j=1,3
          v0(i)=v0(i)+v2(j)*g1(j,i)
        enddo
      enddo
c      CALL SHORTN(V0,V0)
      end



