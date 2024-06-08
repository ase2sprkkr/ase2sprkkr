      subroutine read_dst(vectors, natom,
     >                    positions, types, magnetic,
     >                     W,GEN, n_operations, operations)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION vectors(3,3)
      INTEGER natom
      DOUBLE PRECISION positions(3, natom)
      INTEGER types(natom)
      DOUBLE PRECISION magnetic(3)
C      include 'W.H'
        INTEGER W(*)
      parameter(max2sort=1000)
C      include 'LEWD.H'
      character*2 lattic*8,buf*78,text*60
      dimension  a(3),b(3),c(3)
      dimension gen(1),ioldnum(max2sort)
      common /u/ u(3,3),a,b,c,um1(3,3),idatastr
      common/l2lat/ fvc(3),pvc(3),hvc(3)
      dimension  it(64)
      common /lat/ lattic
      character *40 tnam(max2sort)

c
      do i=1,3
        do j=1,3
          u(i,j)=0
        enddo
        u(i,i)=1d0
      enddo
      pi=dpi()
      lattic='-AUTO- +'
      buf='**************************************'
cccccccccccc      REWIND( ilun(17))
c     REWIND( 11 )
c      read(ilun(17),'(A)')text
c      read(ilun(17),*)alat,boa,coa,alfa,beta,gamma
c      if(boa.lt.0.d0)boa=abs(boa/alat)
c      if(coa.lt.0.d0)coa=abs(coa/alat)
c      read(11,*)a
c      read(11,*)b
c      read(11,*)c
c      read(11,*)natom
      a = vectors(:,1)
      b = vectors(:,2)
      c = vectors(:,3)
      NAMAX=natom
      CALL GENVEC

      call defrr(itau1,3*NATOM)
      call defrr(itau0,3*NATOM)
      call defrr(itau00,3*NATOM)
      call defrr(itaux,NATOM)
      call defrr(itauy,NATOM)
      call defrr(itauz,NATOM)
      call defi(is,NATOM)
      call defi(is1,NATOM)
      call defi(iss,64*NATOM)
      call defrr(i_alcs,NATOM*3)
      call def_sorts(natom,types,positions,magnetic,nsort,a,b,c
     $     ,W(itau0),W(itau1)
     $     ,w(is),W(is1),W(iss),GEN,W(itaux),W(itauy),W(itauz),
     $     w(i_alcs),it,ng,tnam,ioldnum,w(itau00))
      call defi(i_iatpos,natom)
      call defi(i_newis,natom)
      call defi(i_nwt,natom)
      call defi(i_ntw,natom)
      call defi(i_isnew,natom)
      call defi(i_itg,64)
      call check_sym(64,NG,NATOM,w(iss),w(i_alcs),w(i_isnew)
     $       ,it,nsort,w(i_iatpos),w(i_newis),w(is),w(i_nwt),w(i_ntw),
     $       w(i_itg),n_operations,nsortnew,ier,W(iTAUX),W(iTAUY),
     $       W(iTAUZ)
     $     ,Gen,tnam,ioldnum, operations)


      end

      subroutine def_sorts(natom,types,positions,magnetic, nsort,a,b,c
     $     ,tau0,tau1,is,is1,iss,
     $     GEN,taux,tauy,tauz,euler,it,ng,name,ioldnum,tau00)
C
C     This subroutine checks all possible crystal symmetry operations
c     and looking for the sutable ones. Besides it defines new classes.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(max2sort=1000)
C
C     Adjustable arrays:
C
      integer types(natom)
      double precision positions(3,natom)
      double precision magnetic(3)

      dimension tau0(3,natom),tau1(3,natom),tau00(3,natom),
     $     is1(natom),is(natom),iss(64,natom),
     $     taux(*),tauy(*),tauz(*),euler(3,*)
C
C     other oter arrays
C
      dimension a(3),b(3),c(3),gen(3,3,64),orig(3),a0(3),b0(3),c0(3)
C
C     Inner arrays:
C
      dimension  it(64),D1(3),f(3),g(3),h(3),ang(3),ioldnum(*)


      character*3 ti
      character*40  name(max2sort)
      character *40 tnam(max2sort)                    !,op*64
      character lattic*8

C------------------------------------------------

      CALL CROSS(F,B,C)
      CALL CROSS(G,C,A)
      CALL CROSS(H,A,B)
C
C     strarting assumptions:
C
      inversion=0               ! no inversion
      iaddparam=0               ! we will not use it
      icb=0                     ! not cubic
      write(ilun(17),*)0,' / No group information'
      write(ilun(6),*)' Primitive lattice vectors:'
       write(ilun(17),17)a,b,c
 17    format(3f21.15,' / A in lat.par'/3f21.15,' / B in lat.par'/
     $      3f21.15,' / C in lat.par')

      write(ilun(6),11)'A',A
      write(ilun(6),11)'B',B
      write(ilun(6),11)'C',C
 11   format(' ',a,' ',3f21.15)

      W=A(1)*F(1)+A(2)*F(2)+A(3)*F(3)
      W1=1/W
C
      DO  I=1,3
        orig(i)=0
        F(I)=F(I)*W1
        G(I)=G(I)*W1
        H(I)=H(I)*W1
      ENDDO

      write(ilun(6),*)' Primitive reciprocal vectors:'
      write(ilun(6),11)'F',F
      write(ilun(6),11)'G',G
      write(ilun(6),11)'H',H
      natom_input=natom
      do i=1,natom
        !read(11,*)(tau0(k,i),k=1,3),tnam(i)
        tau0(:,i) = positions(:,i)
        write(tnam(i)(:20),*) types(i)
        write(tnam(i)(21:),'(i4.4)')i
        tau00(1,i)=tau0(1,i)
        tau00(2,i)=tau0(2,i)
        tau00(3,i)=tau0(3,i)
        call tofirst(tau0(1,i),a,b,c,f,g,h)
        ioldnum(i)=i
        is(i)=-1
      enddo
      Hx=magnetic(1)
      Hy=magnetic(2)
      Hz=magnetic(3)
      !read(11,*)HX,HY,HZ
      if(Hx**2+Hy**2+Hz**2.gt.1.d-7)then
c        read(11,*)HX,HY,HZ
        write(ilun(17),172)HX,HY,HZ
 172    format(/3f21.15,' / Magnetic field direction'/)
        anorm=sqrt(hx*hx+hy*hy+hz*hz)
        beth=acos(hz/anorm)
        if(hx*hy.lt.1.d-10)then
          alfh=0.d0
        else
          alfh=atan2(hy,hx)
        endif
        gamh=0d0
      else
        alfh=0d0
        beth=0d0
        gamh=150.d0
      endif
      nsort=0
      do iatom=1,natom
        euler(1,iatom)=alfh
        euler(2,iatom)=beth
        euler(3,iatom)=gamh

        if(is(iatom).lt.0)then
          nsort=nsort+1
          is(iatom)=nsort
          do jatom=iatom+1,natom
            if(tnam(jatom)(:20) .eq. tnam(iatom)(:20))is(jatom)=nsort
          enddo
        endif
      enddo
      call check_uniq(is,tau0,a,b,c,f,g,h,natom,ioldnum)
      if(natom.ne.natom_input)then
        write(ilun(6),*) '*** NEW BASIS: ***'
        write(ilun(6),11)'A',A
        write(ilun(6),11)'B',B
        write(ilun(6),11)'C',C
        write(ilun(6),*) '*** NEW ATOMS: ***',natom
        do iatom=1,natom
          print18,(tau0(i,iatom),i=1,3),is(iatom)
        enddo
      endif

C
C-----------------------------------------------------
C
      ismtr=0
      igid=0                            ! group index
      DO  IG=1,64
C
C     Checking all possible crystal rotations
C
c$$$        op(ig:ig)='0'
        if(ICHECK(A,B,C,F,G,H,gen(1,1,IG)).eq.1)then

C
C     lattice vectrors rotation OK, checking basis:
C
          do iat=1,natom

C
C     generating new basis set in the rotated coordinat system:
C
            do i=1,3
              tau1(i,iat)=0
              do l=1,3
                tau1(i,iat)=tau1(i,iat)+tau0(l,iat)*gen(l,i,IG)
              enddo
            enddo
            call tofirst(tau1(1,iat),a,b,c,f,g,h)
      enddo
c
c     In the original coordin system the origin is in the first atom.
c
          do i1=1,natom
C
C     Now we are checking all possible origins, that have to coinside with
C     some of the atom of the same sort as the first one
C
            if(is(i1).eq.is(1))then

              do iat=1,natom
C
C     For each atom in the original basis
C
                IOP=0
                do jat=1,natom
                  if(is(jat).eq.is(iat))then
C
C     We are looking for the correponding atom in the rotated basis:
C
                    DD=0
                    do i=1,3
                      D1(i)=tau1(i,iat)-tau1(i,i1)
     $                     -tau0(i,jat)+tau0(i,1)
                    enddo
                    call tofirst(D1,a,b,c,f,g,h)
                    DD=D1(1)**2+D1(2)**2+D1(3)**2
                    if(DD.lt.1.d-5)then
C
C     We found it!:
C
                      IOP=1
                      is1(iat)=jat
c     exit
                      goto 101
                    endif
                  endif

                enddo           ! end of rotated atom cycle
 101            continue
                if(IOP.eq.0)then
C     exit
                  goto 102
                endif
              enddo             ! end of original set cycle
 102          continue
              if(IOP.eq.1)then
c     exit
                goto 103
              endif
            endif
          enddo                 ! end of new origin cycle
 103      continue
          if(IOP.ne.0)then
            ismtr=ismtr+1
            do iatom=1,natom
              iss(ismtr,iatom)=is1(iatom)
            enddo
c$$$            op(ig:ig)='1'
            it(ismtr)=ig
            igid=igid+iopind(ig)
            if(ig.eq.33)then
              do iat=1,natom
                iat1=is1(iat)
                orig(1)=(tau0(1,iat)+tau0(1,iat1))/2
                orig(2)=(tau0(2,iat)+tau0(2,iat1))/2
                orig(3)=(tau0(3,iat)+tau0(3,iat1))/2
                if(is1(iat).eq.iat)then
c     exit
                  goto 104
                endif
              enddo
 104          continue
              inversion=ismtr
              write(ilun(6),*)'Inversion center:'
              write(ilun(6),17)orig
            endif
          endif
        ENDIF
      ENDDO
C
C     End of checking!
C
C     Trying to define point group (Just for fun)
C
c$$$      TI='=-='
c$$$      call pntgrp(0,IPG,N,TI,OP)
c$$$      ng=n

      call group_by_ind(igid,ti,6)
      write(ilun(6),*)'============================================'
      write(ilun(6),*)'Group:',TI
      write(ilun(6),*)'============================================'
      write(ilun(6),*)
      ng=ismtr
C     print*,OP
      do iat=1,natom
        is1(iat)=iat
        if(natom.eq.natom_input)then
          do i=1,3
            tau0(i,iat)=tau0(i,iat)-orig(i)
          enddo
        else
          inversion=0
        endif
      enddo

      if(inversion.ne.0)then
        do iat=1,natom
          invatom=iss(inversion,iat)
c             print*,invatom,iat
c             print*,(tau0(k,invatom),k=1,3)
c             print*,(tau0(k,iat),k=1,3)
          if(invatom.ne.iat)then
            tau0(1,invatom)=-tau0(1,iat)
            tau0(2,invatom)=-tau0(2,iat)
            tau0(3,invatom)=-tau0(3,iat)
          endif
        enddo
      endif

      do iat=1,natom

        do iop=1,ismtr
          NN=iss(iop,iat)
          if(is1(NN).gt.iat)is1(NN)=iat
        enddo
      enddo
c      print*,'*********'
c      do iop=1,ismtr
c        print*,(iss(iop,iat),iat=1,natom)
c      enddo

c      print*,'IS1:',(is1(II),II=1,natom)
      do iatom=1,natom
        ioldnum(iatom)=is1(iatom)*1000+ioldnum(iatom)
      enddo
c      call RESORTI4(Natom,IS1,tau0,tnam)
      call RESORTI4(Natom,ioldnum,tau0,tnam)
      isort=1
      do iatom=1,natom
        is1(iatom)=ioldnum(iatom)/1000
        ioldnum(iatom)=ioldnum(iatom)-is1(iatom)*1000
      enddo
      do iat=1,natom
c        print*,tnam(iat)(21:)
        if(natom.eq.natom_input)read(tnam(iat)(21:),*) ioldnum(iat)
c        print*,iat,ioldnum(iat)
        if(iat.gt.1)then
          if(is1(iat).ne.is1(iat-1))then
            tau1(1,isort)=tau0(1,iat-1)
            tau1(2,isort)=tau0(2,iat-1)
            tau1(3,isort)=tau0(3,iat-1)
            name(isort)=tnam(iat-1)
            isort=isort+1
          endif
        endif
        write(ilun(6),174)tau0(1,iat),tau0(2,iat),tau0(3,iat)
     $       ,isort,iat,is1(iat)
        is(iat)=isort
 174    format(3g14.6,4i6)
      enddo
      tau1(1,isort)=tau0(1,iat-1)
      tau1(2,isort)=tau0(2,iat-1)
      tau1(3,isort)=tau0(3,iat-1)
      name(isort)=tnam(iat-1)

      inv=0
      if(inversion.eq.0)inv=1
      write(ilun(6),*)' Symmetry operations:',ismtr
      write(ilun(6),116)(it(i),i=1,ismtr)

 116  format(24i3)
      nsort=isort
      if(natom.eq.natom_input)then
        do iat=1,natom
          niat=iat
          call tofirst01(tau0(1,iat),a,b,c,f,g,h)
          taux(iat)=tau00(1,ioldnum(iat))
          tauy(iat)=tau00(2,ioldnum(iat))
          tauz(iat)=tau00(3,ioldnum(iat))
c          print*,'-->',taux(iat),tauy(iat),tauz(iat),is(iat),is1(iat)
        enddo
      else
        write(ilun(17),'(2I5,3A)')nsort,natom_input,
     $        ' / nsort,natom      '
     $       ,'    taus         ','          ICL  IQ(inp) CL(inp) '

        do iat=1,natom_input
          taux(iat)=tau00(1,iat)
          tauy(iat)=tau00(2,iat)
          tauz(iat)=tau00(3,iat)
          call tofirst(tau00(1,iat),a,b,c,f,g,h)
          do j=1,natom
            if(dist_s(tau00(1,iat),tau0(1,j)).lt.1d-5)is1(iat)=is(j)
          enddo
          write(ilun(17),317)taux(iat),tauy(iat),tauz(iat),is1(iat)
     $         ,iat,name(is1(iat)) !,is0(iatom)
 317      format(3f21.15,i6,2x,i6,'  ',a12)
c          print*,'-->',taux(iat),tauy(iat),tauz(iat),is1(iat)
        enddo
        write(ilun(17),*)NG,
     $       ' / Point symmetry operations for reciprocal sums'
        write(ilun(17),22)(it(i),i=1,ng)
 22   format(24i3)
        write(ilun(17),*)' !!! NON-PRIMITIVE SET !!! '
        write(ilun(17),17)a,b,c
        write(ilun(17),172)0.,0.,0.
        write(ilun(17),*)nsort,natom,' / nsort,natom '
        do j=1,natom
          iat=ioldnum(j)
          write(ilun(17),317)taux(iat),tauy(iat),tauz(iat),is(j)
     $         ,iat,name(is(j)) !,is0(iatom)
        enddo       
      endif

      write(ilun(6),*)'==========================='
      write(ilun(6),*)nsort,'          /nsort'
      do isort=1,nsort
 18     format(3f21.15,'   ''',i4,'''  /')
        write(ilun(6),18)tau1(1,isort),tau1(2,isort),
     $       tau1(3,isort),isort
      enddo
c
c     Now all tau1 contain the sort coordinats and tau0 and taux..
c     contain the atomic coordinats
      if(natom.ne.natom_input)stop
      end
c
      subroutine tofirst(t,a,b,c,f,g,h)
      implicit double precision (a-h,o-z)
C
C     This routine reduses atom. corr. in the first unit cell: i.e.
c     t=alpha*A+beta*B+gamma*C, where -0.5<{alpha,beta,gamma}<0.5,
c     (A,B,C - primitive translations, t -at. coordinats)
C
      dimension t(3), a(3),b(3),c(3),f(3),g(3),h(3)
      al=t(1)*f(1)+t(2)*f(2)+t(3)*f(3)+1000
      bl=t(1)*g(1)+t(2)*g(2)+t(3)*g(3)+1000
      cl=t(1)*h(1)+t(2)*h(2)+t(3)*h(3)+1000
      al=al-nint(al+1d-5)
      bl=bl-nint(bl+1d-5)
      cl=cl-nint(cl+1d-5)
      t(1)=al*a(1)+bl*b(1)+cl*c(1)
      t(2)=al*a(2)+bl*b(2)+cl*c(2)
      t(3)=al*a(3)+bl*b(3)+cl*c(3)
      return
c$$$        print*,'abcfgh'
c$$$        print*,a
c$$$        print*,b
c$$$        print*,c
c$$$        print*,f
c$$$        print*,g
c$$$        print*,h
c$$$        print*,al,bl,cl, '--'
      end

      subroutine tofirst01(t,a,b,c,f,g,h)
      implicit double precision (a-h,o-z)
C
C     This routine reduses atom. corr. in the first unit cell: i.e.
c     t=alpha*A+beta*B+gamma*C, where 0<{alpha,beta,gamma}<1.,
c     (A,B,C - primitive translations, t -at. coordinats)
C
      dimension t(3), a(3),b(3),c(3),f(3),g(3),h(3)
      al=t(1)*f(1)+t(2)*f(2)+t(3)*f(3)+1000
      bl=t(1)*g(1)+t(2)*g(2)+t(3)*g(3)+1000
      cl=t(1)*h(1)+t(2)*h(2)+t(3)*h(3)+1000
      al=al-int(al)
      bl=bl-int(bl)
      cl=cl-int(cl)
      t(1)=al*a(1)+bl*b(1)+cl*c(1)
      t(2)=al*a(2)+bl*b(2)+cl*c(2)
      t(3)=al*a(3)+bl*b(3)+cl*c(3)
      return
      end

      SUBROUTINE RESORTI4(N,IT,a1,tnam)
      implicit double precision (a-h,o-z)
      integer IT(N),P
      dimension a1(3,*)
      character*40 tnam(*),tn
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
          tn=tnam(k)
          tnam(k)=tnam(i)
          tnam(i)=tn
          do l=1,3
            a1s=a1(l,k)
            a1(l,k)=a1(l,i)
            a1(l,i)=a1s
          enddo
          IT(K)=IT(I)
          IT(I)=P
        ENDIF
C
      ENDDO
      END

      function ICHECK(A,B,C,F,G,H,gen)
C
C     Thsi function cheks if the matrix gen(3,3) leaves lattice
c     vectors unchangable.
C
      implicit double precision (a-h,o-z)
      dimension a(3),b(3),c(3),gen(3,3)
      dimension f(3),g(3),h(3),ang(3)
      ICHECK=1
      do i=1,3
        ang(i)=0
        do l=1,3
          ang(i)=ang(i)+A(l)*gen(l,i)
        enddo
      enddo
      call tofirst(ang,a,b,c,f,g,h)
      if(ang(1)**2+ang(2)**2+ang(3)**2.gt.1.d-5)then
        ICHECK=0
        return
      endif

      do i=1,3
        ang(i)=0
        do l=1,3
          ang(i)=ang(i)+B(l)*gen(l,i)
        enddo
      enddo
      call tofirst(ang,a,b,c,f,g,h)
      if(ang(1)**2+ang(2)**2+ang(3)**2.gt.1.d-5)then
        ICHECK=0
        return
      endif
      
      do i=1,3
        ang(i)=0
        do l=1,3
          ang(i)=ang(i)+C(l)*gen(l,i)
        enddo
      enddo
      call tofirst(ang,a,b,c,f,g,h)
      if(ang(1)**2+ang(2)**2+ang(3)**2.gt.1.d-5)then
        ICHECK=0
        return
      endif
      end
c
      subroutine check_uniq(is,tau0,a,b,c,f,g,h,natom,ioldnum)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension tau0(3,natom),is(natom),ioldnum(natom)
      dimension tran(3)
C     other oter arrays
C
      dimension a(3),b(3),c(3),f(3),g(3),h(3),ang(3)

 100  continue
      DLENGTHMIN=1d20
      ntran=0
      inew=0
      do i=1,natom
c       print*,'i=',i,is(i),'   ',tnam(i)
        if(is(i).eq.is(1))then
          do j=i+1,natom
            if(is(j).eq.is(1))then
              do k=1,natom
                iok=0
                do ii=1,3
                  ang(ii)=tau0(ii,k)+tau0(ii,i)-tau0(ii,j)
                enddo
                call tofirst(ang,a,b,c,f,g,h)
                do kk=1,natom
                  if(is(kk).eq.is(k).and.
     $                 abs(tau0(1,kk)-ang(1))+
     $                 abs(tau0(2,kk)-ang(2))+
     $                 abs(tau0(3,kk)-ang(3)).lt.1.d-5)then
                    iok=1
                    goto 91
                  endif
                enddo
 91             continue
                if(iok.eq.0)goto 92
              enddo
 92           continue
              if(iok.eq.1)then
                inew=1
c$$$                write(ilun(6),*)'*********Warning!*********'
c$$$                write(ilun(6),*)'There is an extra trans. vector:'
c$$$                write(ilun(6),*) (tau0(ii,i)-tau0(ii,j),ii=1,3)
Cabc                DLENGTH=dist_s(tau0(ii,i),tau0(ii,j))
                DLENGTH=dist_s(tau0(1,i),tau0(1,j))         ! OS fix
                if(DLENGTH.lt.DLENGTHMIN)then
                  do ii=1,3
                    tran(ii)=tau0(ii,i)-tau0(ii,j)
                    DLENGTHMIN=DLENGTH
                  enddo                  
                endif
c                call newabc(a,b,c,tran(1,ntran)


c                goto 100
              endif
            endif
          enddo
        endif
 93     continue
      enddo
      if(inew.eq.0)return
c      print*,'***** MINIMAL VECTOR****:',tran
      call newabc(a,b,c,tran)
c      print*,'NEW VECTRORS:'
c      print*,a
c      print*,b
c      print*,c
      CALL CROSS(F,B,C)
      CALL CROSS(G,C,A)
      CALL CROSS(H,A,B)
      W=A(1)*F(1)+A(2)*F(2)+A(3)*F(3)
      W1=1/W
C
      DO  I=1,3
        F(I)=F(I)*W1
        G(I)=G(I)*W1
        H(I)=H(I)*W1
      ENDDO
      do i=1,natom
c        print111,'   do',(tau0(k,i),k=1,3)
        call tofirst(tau0(1,i),a,b,c,f,g,h)
c        print111,'posle',(tau0(k,i),k=1,3)
 111   format(' ',a,' ',3f21.15)
      enddo
      do i=1,natom
        do j=i+1,natom
          if(dist_s(tau0(1,i),tau0(1,j)).lt.1d-5)is(j)=-abs(is(j))
        enddo
      enddo
      iat=0
      do i=1,natom
        if(is(i).gt.0)then
          iat=iat+1
          is(iat)=is(i)
          ioldnum(iat)=ioldnum(i)
          call cp_v(tau0(1,iat),tau0(1,i))
        endif
      enddo
      natom=iat
c      print*,'NEW NATOM=',natom
      goto 100
      end
      subroutine newabc(a,b,c,d)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension a(3),b(3),c(3),d(3)
      ic=0
      vol0=spat_s(a,b,c)             
      if(spat_s(a,b,d).lt.vol0)then
        ic=3
Cabc        vol0=vol                    ! OS fix
      endif
      if(spat_s(b,c,d).lt.vol0)then
        ic=1
Cabc        vol0=vol                    ! OS fix
      endif
      if(spat_s(a,c,d).lt.vol0)then
        ic=2
      endif
      if(ic.eq.1)call cp_v(a,d)
      if(ic.eq.2)call cp_v(b,d)
      if(ic.eq.3)call cp_v(c,d)
      end
      subroutine cp_v(d,a)
      real*8 a(3),d(3)
      d(1) = a(1)
      d(2) = a(2)
      d(3) = a(3)
      end

      function dist_s(d,a)
      real*8 a(3),d(3),dist_s
      dist_s=(d(1)- a(1))**2+(d(2) - a(2))**2+(d(3) - a(3))**2
      dist_s=sqrt(dist_s)
      end

      function spat_s(a,b,c)

      real*8 a(3),b(3),c(3),d(3),s,spat_s
C
      d(1) = a(2)*b(3) - b(2)*a(3)
      d(2) = a(3)*b(1) - b(3)*a(1)
      d(3) = a(1)*b(2) - b(1)*a(2)
      s=c(1)*d(1)+c(2)*d(2)+c(3)*d(3)

      spat_s= abs(s)
      if(spat_s.lt.1d-7)spat_s=1d+7
      return
      end
