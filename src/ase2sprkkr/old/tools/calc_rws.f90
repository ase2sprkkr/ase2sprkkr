program calc_rws 

implicit none

integer :: i,j,N,Z,JWS,var,var_rbas,var_zinput,irmtalg,l1,n1,m1,new_rws_scale,jrmt,calc_rmt,variante
integer :: kmax,kx,ky,kz,n_step_vol,n_step_vol_fix,var_anz,var_scale_RWS,test_sprkkr,scale_by_hand
integer,dimension(:),allocatable :: Z_input
real*8 :: VRWS,pi,ALAT,scale_RWS,R1,RWS_tab,x0,y0,z0,d,VRMT,VRMT_sum,vol_cell_rws,dummy,SWS
real*8 :: RWS_old,VRBAS,VRWS_sum,DX,DIJ,AUX,sx,sy,sz,x1,y1,z1,scale_RMT_RWS,var_constant
real*8 :: x_max,y_max,z_max,step_x,step_y,step_z,x2,y2,z2,alpha_read,delta_x,delta_y,delta_z
real*8 :: a,b,c,d2,phi,RMT_min,q,p,u,v,u_help,v_help,D3,vol_ks1,vol_ks2,d1,VRWS1,VRWS2
real*8,dimension(3) :: RBAS1,RBAS2,RBAS3,pos_nn
real*8,dimension(:),allocatable :: RWS,RMT_out,RWS_new,RMT
real*8,dimension(:,:),allocatable :: atom_pos,atoms_diff,tess_eq,tess_bound
real*8,dimension(:,:),allocatable :: x_bound,y_bound,z_bound,alpha
logical :: irmtalg_exp,cubic,check_rmt_scaled
complex*8 :: z_imag,x11,x21,x31,y11,y21,y31

pi = 3.141592653589793238462643D0
VRWS=0.0D0
JWS=721
R1=0.000001D0
VRWS_sum=0.0D0
irmtalg_exp=.true.
z_imag=(0.0,1.0)

open(335,file='atompositionen')
open(337,file='rmt_rws_data')
open(338,file='overlapp')
open(339,file='overlapp_spheres') !for drawing circles

write(*,*) ''
write(*,*) 'Anz. Atome pro EZ (via file (1))'
read(5,*) var_anz

if(var_anz.ne.1) then

 write(*,*) 'Anzahl - Eingabe von Hand:'
 read(5,*) N

else if(var_anz.eq.1) then

N=0

do
 read(335,*,end=111) dummy
 N=N+1
end do

111 rewind(335)

N=N-1

write(*,*) 'Anz. Atome pro EZ'
write(*,*) N

end if

allocate(Z_input(N))
allocate(RWS(N))
allocate(atom_pos(N,1:3))
allocate(atoms_diff(1:N,1:N))
allocate(RWS_new(N))
allocate(RMT(1:N))
allocate(RMT_out(1:N))
allocate(tess_eq(1:N,1:N))
allocate(tess_bound(1:N,1:N))
allocate(x_bound(1:N,1:N))
allocate(y_bound(1:N,1:N))
allocate(z_bound(1:N,1:N))
allocate(alpha(1:N,1:N))

!------------------------------ Einlesen Kernladungszahl
write(*,*) 'Kernladungszahl via file (1)'
read(5,*) var_zinput

if(var_zinput.ne.1) then

do i=1,N
  write(*,*) 'Atom',i
  read(5,*) Z_input(i)
end do

else if(var_zinput.eq.1) then

!------------------------------ Einlesen Kernladungszahl automatisch

do i=1,N+1
  if(i.eq.1) then
    read(335,*) dummy
  else
    read(335,*) Z_input(i-1)
    write(*,*) i-1,'Z',Z_input(i-1)
  end if
end do

rewind(335)

end if

!------------------------------ Einlesen Gitterkonstante
write(*,*) 'Gitterkonstante [a.u.]'
write(*,*) 'via file (1)'
read(5,*) var_constant
if(var_constant.ne.1) then
  write(*,*) 'Gitterkonstante [a.u.] - Eingabe von Hand:'
  read(5,*) ALAT
else if(var_constant.eq.1) then
  open(335,file='atompositionen')
  read(335,*) ALAT
  write(*,*) 'ALAT:',ALAT
end if

rewind(335)

!------------------------------ Einlesen Basisvektoren
write(*,*) 'Basisvektoren:'
write(*,*) 'via file (1)'
read(5,*) var_rbas

if (var_rbas.ne.1) then

do i=1,3

   if(i.eq.1) then

      write(*,*) 'RBAS1'
      read(5,*) RBAS1(:)

   else if(i.eq.2) then

      write(*,*) 'RBAS2'
      read(5,*) RBAS2(:)

   else if(i.eq.3) then

      write(*,*) 'RBAS3'
      read(5,*) RBAS3(:)

   end if
   
end do

else if(var_rbas.eq.1) then

!------------------------------ automatisches Einlesen Basisvektoren
open(336,file='basisvektoren')

do i=1,3

   if(i.eq.1) then

      read(336,*) RBAS1(1),RBAS1(2),RBAS1(3)
      write(*,*) RBAS1(1),RBAS1(2),RBAS1(3)

   else if(i.eq.2) then

      read(336,*) RBAS2(1),RBAS2(2),RBAS2(3)
      write(*,*) RBAS2(1),RBAS2(2),RBAS2(3)

   else if(i.eq.3) then

      read(336,*) RBAS3(1),RBAS3(2),RBAS3(3)
      write(*,*) RBAS3(1),RBAS3(2),RBAS3(3)

   end if
   
end do

close(336)

end if

!------------------------------ Berechnung Volumen Einheitszelle
  VRBAS = RBAS1(1)*(RBAS2(2)*RBAS3(3)-RBAS2(3)*RBAS3(2))
  VRBAS = VRBAS + RBAS1(2)*(RBAS2(3)*RBAS3(1)-RBAS2(1)*RBAS3(3))
  VRBAS = VRBAS + RBAS1(3)*(RBAS2(1)*RBAS3(2)-RBAS2(2)*RBAS3(1))

!    write(*,*) RBAS1(:)
!    write(*,*) RBAS2(:)
!    write(*,*) RBAS3(:)
!    write(*,*) VRBAS
!    stop

  VRBAS = abs(VRBAS)*ALAT**3.D0

!    write(*,*) VRBAS
!    stop

write(*,*) ''
write(*,*) 'Variante Berechnung RMT und RWS:'
write(*,*) 'SPRKKR (1), HUTSEPOT (2), Tesellierung (experimentell) (3)'
write(*,*) ''
read(5,*) variante

if(variante.eq.1) then
!------------------------------ Einlesen RWS aus Tabelle unde Berechnung des Gesamtvolumens der ASA-Kugeln
do i = 1,N
 
  Z = Z_input(i)

!  write(*,*) i,Z_input(i)
!  stop

  call readrws(Z,N,RWS_tab)

  RWS(i) = RWS_tab

  VRWS = (4.D0/3.D0)*pi*RWS(i)**3.D0
  VRWS_sum = VRWS_sum+VRWS

end do

!------------------------------ Skalierung der RWS
!write(*,*) VRWS_sum,VRBAS
!write(*,*) VRWS_sum-VRBAS
!stop
if(abs(VRWS_sum-VRBAS).gt.1.D-15) then

  scale_RWS = (VRBAS/VRWS_sum)**(1.D0/3.D0)

  VRWS = 0.0D0
  VRWS_sum = 0.0D0

  do i =1,N

     RWS_old = RWS(i)
     RWS_new(i) = scale_RWS*RWS(i)
    
     DX = log(RWS(i)/R1)/(JWS-1.0D0) !siehe xband

     if(i.eq.1) write(*,*) ''
!     if(mod(i,2).eq.1) then !for potentials using symmetry option (more than one pot on a site)
!     if(i.eq.1) write(*,*) ''
!     if(i.eq.1) write(*,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
!     if(i.eq.1) write(337,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
!     if(calc_rmt.eq.1) then
!       write(*,99001) i,R1,DX,jrmt,RMT_out(i),721,RWS_new(i)
!       write(337,99001) i,R1,DX,jrmt,RMT_out(i),721,RWS_new(i)
!     else
!       write(*,99001) i,R1,DX,0,RMT_out(i),721,RWS_new(i)
!       write(337,99001) i,R1,DX,0,RMT_out(i),721,RWS_new(i)
!     end if
!      end if

     VRWS = (4.D0/3.D0)*pi*RWS_new(i)**3.D0
     VRWS_sum = VRWS_sum+VRWS
  
  end do

!    write(*,*) 'VRWS_sum',VRWS_sum
!    write(*,*) 'VRBAS   ',VRBAS

else

  do i = 1,N

     RWS_old = RWS(i)
     RWS_new(i) = RWS_old

     DX = log(RWS(i)/R1)/(JWS-1.0D0) !siehe xband

     if(i.eq.1) write(*,*) 'RWS NOT scaled'
!     if(mod(i,2).eq.1) then !for potentials using symmetry option (more than one pot on a site)
!     if(i.eq.1) write(*,*) ''
!     if(i.eq.1) write(*,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
!     if(i.eq.1) write(337,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
!     if(calc_rmt.eq.1) then
!       write(*,99001) i,R1,DX,jrmt,RMT_out(i),721,RWS_new(i)
!       write(337,99001) i,R1,DX,jrmt,RMT_out(i),721,RWS_new(i)
!     else
!       write(*,99001) i,R1,DX,0,RMT_out(i),721,RWS_new(i)
!       write(337,99001) i,R1,DX,0,RMT_out(i),721,RWS_new(i)
!     end if
!      end if

  end do

end if

!------------------------------ Einlesen Atompositionen fuer Berechnung RMT
 write(*,*) 'Atompositionen:'
 write(*,*) 'via file (1)'
 read(5,*) var

if (var.ne.1) then
 
  do i=1,N
    write(*,*) 'Atom',i
    read(5,*) atom_pos(i,1:3)
  end do

else if(var.eq.1) then

!------------------------------ automatisches Einlesen Atompositionen fuer Berechnung RMT

  do i=1,N+1
   if(i.eq.1) then
     read(335,*) dummy
   else
    read(335,*) dummy,atom_pos(i-1,1),atom_pos(i-1,2),atom_pos(i-1,3)
    write(*,*) i-1,atom_pos(i-1,1),atom_pos(i-1,2),atom_pos(i-1,3)
   end if
  end do

end if

close(335)

!------------------------------ Einlesen Algortihmus zur Berechnung von RMT
 write(*,*) ''
 write(*,99002) !'IRMTALG (1) - eigene Routine'
 write(*,99003) !'IRMTALG (2) - Variante aus potio.f, nicht vollstaendig'
 write(*,99004) !'IRMTALG (3) - Variante aus potio.f'
 write(*,99005) !'IRMTALG (4) - experimentelle Variante'
 read(5,*) irmtalg

!------------------------------ Berechnung RMT
if(N.eq.1) then
 write(*,*) 'relative Postition naechster Nachbar:'
 read(5,*) pos_nn(:)
 RMT_out(1) = sqrt((pos_nn(1))**2.D0+(pos_nn(2))**2.D0+(pos_nn(3))**2.D0)*ALAT
 RMT_out(1) =  RMT_out(1)*0.5D0
 write(*,*) ''
 write(*,*) 'RMT half atomic distance next neighbor!'
else
 do i=1,N

    RMT_out(i) = 9999D0

  do j=1,N

    if(i.ne.j) then

      atoms_diff(i,j)=(atom_pos(i,1)-atom_pos(j,1))**2.D0
      atoms_diff(i,j)=atoms_diff(i,j)+(atom_pos(i,2)-atom_pos(j,2))**2.D0
      atoms_diff(i,j)=atoms_diff(i,j)+(atom_pos(i,3)-atom_pos(j,3))**2.D0
      atoms_diff(i,j)=sqrt(atoms_diff(i,j))*ALAT

      !write(*,*) atoms_diff(i,j)

      RMT(i) = (1.D0/(1.D0 + RWS_new(j)/RWS_new(i)))*atoms_diff(i,j)

 !------------------------------ RMT Algorithmus 1 (wird fuer Algortihmus 2 und 3 auch verwendet)     
      if(RMT(i).lt.RMT_out(i)) then
         RMT_out(i) = RMT(i)
      end if
      !write(*,*) i,j,RMT(i),RMT_out(i)
      !top

      jrmt = int(log(RMT_out(i)/R1)/DX) + 1 !Berechnung JRMT
      !write(*,*) jrmt
!------------------------------ RMT Algorithmus 2
     if(irmtalg.eq.2) then !siehe potio.f (?DNRM -> Hubert fragen)

       !CALL DNRM2(3,X,1)
       !DIJ = DNRM2(3,DQBAS,1)*ALAT ?

       AUX = DIJ/(RWS_new(i)+RWS_new(j))
       RMT_out(i) = MIN(RMT_out(i),AUX*RWS_new(i))
       RMT_out(j) = MIN(RMT_out(j),AUX*RWS_new(j))

      jrmt = int(log(RMT_out(i)/R1)/DX) + 1 !Berechnung JRMT

!------------------------------ RMT Algorithmus 3
     else if(irmtalg.eq.3) then !siehe potio.f

       DIJ = ALAT
       AUX = DIJ/2.D0
       RMT_out(i) = MIN(RMT_out(i),AUX)
       RMT_out(j) = MIN(RMT_out(j),AUX)

     if(ABS(RMT_out(i)-RWS_new(i)).LT.1D-6 ) then !siehe potio.f
        RMT(i) = 0.97D0*RWS(i)
        write(*,*) 'RMT scaled'
     end if

      jrmt = int(log(RMT_out(i)/R1)/DX) + 1 !Berechnung JRMT

     end if
    end if
  end do
 !write(*,*) RMT_out(i)
 end do
end if
 !do i=1,N
 ! write(*,*) RMT_out(i)
 !end do
 !stop
    !if(i.eq.1.and.irmtalg.eq.1) write(*,*) 'IRMTALG = 1'
    !if(i.eq.1.and.irmtalg.eq.2) write(*,*) 'IRMTALG = 2'
    !if(i.eq.1.and.irmtalg.eq.3) write(*,*) 'IRMTALG = 3'
    !if(i.eq.1.and.irmtalg.eq.4) write(*,*) 'IRMTALG = 4'
!------------------------------ RMT Algorithmus (experimentelle Variante)
  if(irmtalg.eq.4) then

  do i=1,N
    do j=1,N

       if(i.ne.j.and.RMT_out(i).gt.RWS_new(i).and.RMT_out(i).gt.RMT_out(j)) then
         !write(*,*) '->',i,RMT_out(i).gt.RWS_new(i),RWS_new(j)
         RMT_out(j) = (RMT_out(i)/abs(RMT_out(i)-RWS(i)))*RWS_new(j)
         !write(*,*) '-->',i,RMT_out(j)
         RMT_out(i) = RWS_new(i)
         !write(*,*) '--->',i,RMT_out(i)
       end if
     
    end do

    if(i.eq.1.and.irmtalg_exp) write(*,*) '!!! IRMTALG EXP !!!'  

    write(*,*) 'Atom',i,'RMT:',RMT_out(i)

  end do

  end if
!------------------------
!if calculation or JRMT is necessary
    if(i.eq.1) then
      write(*,*) ''
      write(*,*) 'calculation of JRMT (1)'
      write(*,*) ''
      read(5,*) calc_rmt
    end if
!------------------------
  do i=1,N
    if(irmtalg.le.4) then 
!      if(mod(i,2).eq.1) then !for potentials using symmetry option (more than one pot on a site)
      if(i.eq.1) then
        write(*,*) ''
        write(*,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
        write(337,*) 'Variante SPRKKR'
        if(irmtalg.eq.1) write(337,*) 'IRMTALG (1) - eigene Routine'
        if(irmtalg.eq.2) write(337,*) 'IRMTALG (2) - Variante aus potio.f, nicht vollstaendig'
        if(irmtalg.eq.3) write(337,*) 'IRMTALG (3) - Variante aus potio.f'
        if(irmtalg.eq.4) write(337,*) 'IRMTALG (4) - experimentelle Variante'
        write(337,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
      end if
      if(calc_rmt.eq.1) then
        write(*,99001) i,R1,DX,jrmt,RMT_out(i),721,RWS_new(i)
        write(337,99001) i,R1,DX,jrmt,RMT_out(i),721,RWS_new(i)
      else
        write(*,99001) i,R1,DX,0,RMT_out(i),721,RWS_new(i)
        write(337,99001) i,R1,DX,0,RMT_out(i),721,RWS_new(i)
      end if
!       end if
    end if
  end do
!---------------- testing new RMT after scaling by hand
  check_rmt_scaled=.false.
  do i=1,N
     do j=1,N
        if(i.ne.j) then
          atoms_diff(i,j)=(atom_pos(i,1)-atom_pos(j,1))**2.D0
          atoms_diff(i,j)=atoms_diff(i,j)+(atom_pos(i,2)-atom_pos(j,2))**2.D0
          atoms_diff(i,j)=atoms_diff(i,j)+(atom_pos(i,3)-atom_pos(j,3))**2.D0
          atoms_diff(i,j)=sqrt(atoms_diff(i,j))*ALAT
        if(RMT_out(i)+RMT_out(j).gt.atoms_diff(i,j)) then
          write(*,*) ''
          write(*,*) 'check RMT radii for atom',i,'and neighbor',j
          write(*,*) 'sum RMT',i,j,':',RMT_out(i)+RMT_out(j)
          write(*,*) 'distance atoms:',atoms_diff(i,j)
          write(*,*) 'overlap:',abs(RMT_out(i)+RMT_out(j)-atoms_diff(i,j))/atoms_diff(i,j),'%'
          write(*,*) 'routine stopped!!!'
          stop
        else
          check_rmt_scaled=.true.
        end if
        end if
      enddo
  end do
  if(check_rmt_scaled) then
    write(*,*) ''
    write(*,*) 'check scaled RMT radii for all atoms ok!'
  end if
!-----------------
else if(variante.eq.2) then
!------------------------------ RMT Algorithmus (hutsepot -> symmetry/find_rmt3d.F90)
!------------------------------ Einlesen Atompositionen fuer Berechnung RMT
 write(*,*) 'Atompositionen:'
 write(*,*) 'via file (1)'
 read(5,*) var

if (var.ne.1) then
 
  do i=1,N
    write(*,*) 'Atom',i
    read(5,*) atom_pos(i,1:3)
  end do

else if(var.eq.1) then

!------------------------------ automatisches Einlesen Atompositionen fuer Berechnung RMT

  do i=1,N+1
   if(i.eq.1) then
     read(335,*) dummy
   else
    read(335,*) dummy,atom_pos(i-1,1),atom_pos(i-1,2),atom_pos(i-1,3)
    write(*,*) i-1,atom_pos(i-1,1),atom_pos(i-1,2),atom_pos(i-1,3)
   end if
  end do

end if

close(335)
!----------------------------- Start Algorithmus HUTSEPOT
 if(N.eq.1) then
   write(*,*) ''
   write(*,*) 'for one atom not implemented!'
   write(*,*) 'use SPRKKR variant'
   write(*,*) 'program stopped!'
   stop
 end if
  VRMT=0.0D0
  VRMT_sum=0.0D0
  VRWS=0.0D0
  VRWS_sum=0.0D0
  write(*,*) ''
  write(*,*) 'scaling RMT (hutsepot) by hand? (1)'
  read(5,*) scale_by_hand
  do i=1,N
     x0=atom_pos(i,1)
     y0=atom_pos(i,2)
     z0=atom_pos(i,3)
     RMT_out(i)=1.D10
     do j=1,N
        sx=x0-atom_pos(j,1)
        sy=y0-atom_pos(j,2)
        sz=z0-atom_pos(j,3)
        do l1=-3,3
           do n1=-3,3
              do m1=-3,3
                 x1=l1*RBAS1(1)+n1*RBAS2(1)+m1*RBAS3(1)+sx
                 y1=l1*RBAS1(2)+n1*RBAS2(2)+m1*RBAS3(2)+sy
                 z1=l1*RBAS1(3)+n1*RBAS2(3)+m1*RBAS3(3)+sz
                 d=sqrt(x1*x1+y1*y1+z1*z1)
                 if(d.gt.1.d-3) RMT_out(i)=min(RMT_out(i),d)
              enddo
           enddo
        enddo
     enddo
  if(scale_by_hand.ne.1) then
     RMT_out(i) = RMT_out(i)*ALAT/2.D0
     if(i.eq.1.and.irmtalg_exp) write(*,*) 'IRMTALG HUTSEPOT'  
     VRMT = (4.D0/3.D0)*pi*RMT_out(i)**3.D0
     VRMT_sum = VRMT_sum+VRMT
  end if
  end do
!--------------- test scaling by hand
  if(scale_by_hand.eq.1) then
  do i=1,N
     if(Z_input(i).eq.38) RMT_out(i)=RMT_out(i)*1.19
     if(Z_input(i).eq.8) RMT_out(i)=RMT_out(i)/1.19
     if(Z_input(i).eq.57) RMT_out(i)=RMT_out(i)*1.19
     if(Z_input(i).eq.13) RMT_out(i)=RMT_out(i)/1.19
     RMT_out(i) = RMT_out(i)*ALAT/2.D0
     if(i.eq.1.and.irmtalg_exp) write(*,*) 'SCALING by hand!!!'
     VRMT = (4.D0/3.D0)*pi*RMT_out(i)**3.D0
     VRMT_sum = VRMT_sum+VRMT
  end do
!---------------- testing new RMT after scaling by hand
  check_rmt_scaled=.false.
  do i=1,N
     do j=1,N
        if(i.ne.j) then
          atoms_diff(i,j)=(atom_pos(i,1)-atom_pos(j,1))**2.D0
          atoms_diff(i,j)=atoms_diff(i,j)+(atom_pos(i,2)-atom_pos(j,2))**2.D0
          atoms_diff(i,j)=atoms_diff(i,j)+(atom_pos(i,3)-atom_pos(j,3))**2.D0
          atoms_diff(i,j)=sqrt(atoms_diff(i,j))*ALAT
        if(RMT_out(i)+RMT_out(j).gt.atoms_diff(i,j)) then
          write(*,*) ''
          write(*,*) 'check RMT radii for atom',i,'and neighbor',j
          write(*,*) 'sum RMT',i,j,':',RMT_out(i)+RMT_out(j)
          write(*,*) 'distance atoms:',atoms_diff(i,j)
          write(*,*) 'overlap:',abs(RMT_out(i)+RMT_out(j)-atoms_diff(i,j))/atoms_diff(i,j),'%'
          write(*,*) 'routine stopped!!!'
          stop
        else
          check_rmt_scaled=.true.
        end if
        end if
      enddo
  end do
  if(check_rmt_scaled) then
    write(*,*) ''
    write(*,*) 'check scaled RMT radii for all atoms ok!'
  end if
  end if
!-----------------
!change between linear or cubic scaling of RWS
    write(*,*) ''
    write(*,*) 'linear(1) or cubic scaling of RWS'
    write(*,*) ''
    read(5,*) var_scale_RWS
    if(var_scale_RWS.eq.1) cubic=.false.
    if(var_scale_RWS.ne.1) cubic=.true.
!if calculation or JRMT is necessary
    write(*,*) ''
    write(*,*) 'calculation of JRMT (1)'
    write(*,*) ''
    read(5,*) calc_rmt
!--------------------------------------------
    if(.not.cubic) then
     write(*,*) 'linear scaling RWS'
     scale_RMT_RWS = (VRBAS/VRMT_sum)**(1.D0/3.D0)
     write(*,*) ''
     do i=1,N
        RWS_new(i) = scale_RMT_RWS*RMT_out(i)
        DX = log(RWS_new(i)/R1)/(JWS-1.0D0)
        jrmt = int(log(RMT_out(i)/R1)/DX) + 1
        VRWS = (4.D0/3.D0)*pi*RWS_new(i)**3.D0
        VRWS_sum = VRWS_sum + VRWS

!--------------------------------------------
!        if(mod(i,2).eq.1) then !for potentials using symmetry option (more than one pot on a site)
        if(i.eq.1) write(*,*) ''
        if(calc_rmt.eq.1) then
          if(i.eq.1) write(*,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
          if(i.eq.1) write(337,*) 'Variante HUTSEPOT - linear scaling'
          if(i.eq.1) write(337,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
          write(*,99001) i,R1,DX,jrmt,RMT_out(i),721,RWS_new(i)
          write(337,99001) i,R1,DX,jrmt,RMT_out(i),721,RWS_new(i)
        else
          if(i.eq.1) write(*,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
          if(i.eq.1) write(337,*) 'Variante HUTSEPOT - linear scaling'
          if(i.eq.1) write(337,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
          write(*,99001) i,R1,DX,0,RMT_out(i),721,RWS_new(i)
          write(337,99001) i,R1,DX,0,RMT_out(i),721,RWS_new(i)
        end if
!         end if
     end do   
     write(*,*) 'VRWS_sum',VRWS_sum
     write(*,*) 'VRBAS   ',VRBAS
    end if
!----------------------------------------- scale RWS from RMT by cubic formula (see PHYSICAL REVIEW B 80, 125123 2009)
    if(cubic) then
      write(*,*) 'cubic scaling RWS'
      RMT_min=1.d+10
      a = 0.D0
      b = 0.D0
      c = 0.D0
      d2 = 0.D0
      do i=1,N
        RMT_min = min(RMT_out(i),RMT_min)     
      end do
!        write(*,*) RMT_min
!        stop
!see http://archive.today/x4cJS#selection-233.3-267.106
!d2*x^3 + c*x^2 + b*x + a-vbas = 0
!a*x^3 + b*x^2 + c*x + d2-vbas = 0
      do i=1,N
         d2 = d2 + (4.D0*pi/3.D0)*RMT_min**3.D0*(RMT_out(i)/RMT_min)**3.0D0
         c = c + (4.D0*pi/1.D0)*RMT_min**2.D0*(RMT_out(i)/RMT_min)
         b = b + (4.D0*pi/1.D0)*RMT_min*(RMT_min/RMT_out(i))
         a = a + (4.D0*pi/3.D0)*(RMT_min/RMT_out(i))**3.0D0
      end do
         d2 = d2 - VRBAS
!      write(*,*) a,b,c,d2
!      stop
!test: roots 3,1,2 (real part), 0,0,0 (imag part)
!         a=1.
!         b=-6.
!         c=11.
!         d2=-6.
!         a=1.  !D>0: z^3+z+10=(z-1-2i)(z-1+2i)(z+2)
!         b=0.
!         c=1.
!         d2=10.
         p = ( 3.0D0*c/a - (b/a)**2.0D0 )/3.0D0
         q = ( 2.0D0*(b/a)**3.D0 - 9.0D0*b*c/a/a + 27.D0*d2/a )/27.D0
         D3 = (p/3.0D0)**3.D0 + (q/2.D0)**2.0D0
!Depending on the sign of D, you follow different strategy.
!          If D3<0, three distinct real roots.
!          If D3=0, three real roots of which at least two are equal.
!          If D3>0, one real and two complex roots.
!      write(*,*) p,q,D3
!      stop
       if(sign(1.D0,D3)+1.0D0.lt.1.0D-16) then
          u_help = -q/2.D0 - sqrt(abs(D3))
          v_help = -q/2.D0 - sqrt(abs(D3))
       else
          u_help = -q/2.D0 + sqrt(abs(D3))
          v_help = -q/2.D0 - sqrt(abs(D3))
       end if
!       write(*,*) u_help,v_help
!       stop
       if(sign(1.D0,u_help)+1.0D0.lt.1.0D-16) then
          u = -(abs(u_help))**(1.D0/3.D0)
       else
          u = (u_help)**(1.D0/3.D0)
       end if
       if(sign(1.D0,v_help)+1.0D-16.lt.1.0D-16) then
          v = -(abs(v_help))**(1.D0/3.D0)
       else
          v = (v_help)**(1.D0/3.D0)
       end if
!          write(*,*) u,v
!          stop
        if(D3.ge.1.D-16) then
          y11 = u + v
          y21 = -(u+v)/2.D0 + z_imag*(u-v)*sqrt(3.D0)/2.D0
          y31 = -(u+v)/2.D0 - z_imag*(u-v)*sqrt(3.D0)/2.D0
        else
          write(*,*) -q/2.D0/sqrt(abs(p)**3.D0/27.D0)
          phi = acos(-q/2.D0/sqrt(abs(p)**3.D0/27.D0))
          y11 =  2.D0 * sqrt(abs(p)/3.D0) * cos(phi/3.D0)
          y21 = -2.D0 * sqrt(abs(p)/3.D0) * cos((phi+pi)/3.D0)
          y31 = -2.D0 * sqrt(abs(p)/3.D0) * cos((phi-pi)/3.D0)
        end if
!          write(*,*) phi
!          write(*,*) y11,y21,y31
!          stop
!Finally, find the three roots
          x11 = y11 - b/a/3.D0
          x21 = y21 - b/a/3.D0
          x31 = y31 - b/a/3.D0
!          write(*,*) x11,x21,x31
!          stop
      SWS=0.0D0
      do i=1,N
          RWS_new(i) = RMT_out(i) + RMT_min*real(x11)/RMT_out(i)
          DX = log(RWS_new(i)/R1)/(JWS-1.0D0)
          jrmt = int(log(RMT_out(i)/R1)/DX) + 1
          VRWS = (4.D0/3.D0)*pi*RWS_new(i)**3.D0
          VRWS_sum = VRWS_sum + VRWS
!          if(mod(i,2).eq.1) then !for potentials using symmetry option (more than one pot on a site)
          if(i.eq.1) write(*,*) ''
          if(i.eq.1) write(*,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
          if(i.eq.1) write(337,*) 'Variante HUTSEPOT - cubic scaling'
          if(i.eq.1) write(337,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
          if(calc_rmt.eq.1) then
            write(*,99001) i,R1,DX,jrmt,RMT_out(i),721,RWS_new(i)
            write(337,99001) i,R1,DX,jrmt,RMT_out(i),721,RWS_new(i)
          else
            write(*,99001) i,R1,DX,0,RMT_out(i),721,RWS_new(i)
            write(337,99001) i,R1,DX,0,RMT_out(i),721,RWS_new(i)
          end if
!           end if
      end do   
!---------------------------------- rescale RWS determined by cubic RMT
         write(*,*) 'VRWS_sum',VRWS_sum
         write(*,*) 'VRBAS   ',VRBAS
         write(*,*) 'new RWS scaling ? (1)'
         read(5,*) new_rws_scale
         if(new_rws_scale.eq.1) then
             rewind(337)
             scale_RWS = (VRBAS/VRWS_sum)**(1.D0/3.D0)
             VRWS = 0.0D0
             VRWS_sum = 0.0D0
              do i =1,N
                RWS_new(i) = scale_RWS*RWS_new(i)
                DX = log(RWS_new(i)/R1)/(JWS-1.0D0) !siehe xband
                jrmt = int(log(RMT_out(i)/R1)/DX) + 1
!                if(mod(i,2).eq.1) then !for potentials using symmetry option (more than one pot on a site)
                if(i.eq.1) write(*,*) ''
                if(i.eq.1) write(*,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
                if(i.eq.1) write(337,*) 'Variante HUTSEPOT - cubic scaling, additional rescaling'
                if(i.eq.1) write(337,*) 'i,R1,DX,JRMT,RMT,JRWS,RWS(new)'
                if(calc_rmt.eq.1) then
                  write(*,99001) i,R1,DX,jrmt,RMT_out(i),721,RWS_new(i)
                  write(337,99001) i,R1,DX,jrmt,RMT_out(i),721,RWS_new(i)
                else
                  write(*,99001) i,R1,DX,0,RMT_out(i),721,RWS_new(i)
                  write(337,99001) i,R1,DX,0,RMT_out(i),721,RWS_new(i)
                end if
!                 end if
                VRWS = (4.D0/3.D0)*pi*RWS_new(i)**3.D0
                VRWS_sum = VRWS_sum+VRWS
              end do
              write(*,*) 'VRWS_sum',VRWS_sum
              write(*,*) 'VRBAS   ',VRBAS
         end if
    end if
!--------------------------------- testing of VRWS_sum and VRBAS determined by SPRKKR
        write(*,*) 'SPRKKR test for VRWS ? (1)'
        read(5,*) test_sprkkr
        if(test_sprkkr.eq.1) then
          do i=1,N
             SWS = SWS + RWS_new(i)**3 !siehe potio.f
          end do
          SWS = (SWS/real(N))**(1.0D0/3.0D0) !siehe potio.f
          VRWS_sum = SWS**3*4D0*PI/3D0
          write(*,*) 'VRWS_sum (SPRKKR)',real(N)*VRWS_sum
          write(*,*) 'VRBAS            ',VRBAS
        end if
  !end if
!-------------------------------------------- calculation of overlapp of spheres for HUTSEPOT cubic and linear
     do i=1,N
       do j=1,N 
          if(i.ne.j) then 
            !write(*,*) i,j,RWS_new(i),RWS_new(j),'!!!'
            !call kugelschnitt(atoms_diff(i,j),RWS_new(i),RWS_new(j),N,i,j) !Problem Uebergabe RWS_new(j)
            atoms_diff(i,j)=(atom_pos(i,1)-atom_pos(j,1))**2.D0
            atoms_diff(i,j)=atoms_diff(i,j)+(atom_pos(i,2)-atom_pos(j,2))**2.D0
            atoms_diff(i,j)=atoms_diff(i,j)+(atom_pos(i,3)-atom_pos(j,3))**2.D0
            atoms_diff(i,j)=sqrt(atoms_diff(i,j))*ALAT
            if(atoms_diff(i,j).lt.1.0D-16) then
              d1=1.0D-16
            else
              d1=(RWS_new(i)**2.D0-RWS_new(j)**2.D0+atoms_diff(i,j)**2.D0)/(2.D0*atoms_diff(i,j)) !Abstand Kugeln 1 vom Mittelpkt Kreisegment
            end if
            !write(*,*) atoms_diff(i,j),d1,'!!!!!'
            if(RWS_new(i)-d1.lt.(-1.0D-16)) then
              a=0.D0
            else
              a=sqrt(RWS_new(i)**2.D0-d1**2.D0) !Radius Kreisegment Schnittflaeche Kugeln
            end if

            vol_ks1 = (RWS_new(i)-sqrt(RWS_new(i)**2.D0-a**2.D0))**2.D0 !Vol Kreissegment
            vol_ks1 = vol_ks1*(2.D0*RWS_new(i)+sqrt(RWS_new(i)**2.D0-a**2.D0))
            vol_ks1 = pi/3.D0*vol_ks1

            vol_ks2 = (RWS_new(j)-sqrt(RWS_new(j)**2.D0-a**2.D0))**2.D0 !Vol Kreissegment
            vol_ks2 = vol_ks2*(2.D0*RWS_new(j)+sqrt(RWS_new(j)**2.D0-a**2.D0))
            vol_ks2 = pi/3.D0*vol_ks2

            !write(*,*) 'vol_ks_1',vol_ks1
            !write(*,*) 'vol_ks_2',vol_ks2
            VRWS1 = (4.D0/3.D0)*pi*RWS_new(i)**3.D0
            VRWS2 = (4.D0/3.D0)*pi*RWS_new(j)**3.D0
            VRWS1 = VRWS1 + VRWS2
            if(i.eq.1) write(338,*) 'i,j,vol_ks1+vol_ks2,%'
            write(338,99008) i,j,vol_ks1+vol_ks2,(vol_ks1+vol_ks2)*100./VRWS1
            !write(338,*) i,j,vol_ks1+vol_ks2,(vol_ks1+vol_ks2)*100./VRWS1
!------------------------------ plot of circles
            if(i.ge.1.and.i.le.2) then     !number of atoms to be plottet
              write(339,99009) 'set object circle at',(i-1)*atoms_diff(i,j),',',0,'size',RWS_new(i)
              write(339,99009) 'set object circle at',(j-1)*atoms_diff(i,j),',',0,'size',RWS_new(j)
            end if
!------------------------------
          end if
        end do
      end do
else if(variante.eq.3) then
!------------------------------ RWS Tesselierung
  kmax=100

  if(irmtalg.eq.6) then
  write(*,*) 'alpha'
  read(5,*) alpha_read

  do i=1,1!N

    n_step_vol=0
    n_step_vol_fix=0

    do j=1,N

       if(i.ne.j) then

         alpha(i,j) = alpha_read
        
         x_max = atom_pos(j,1)-atom_pos(i,1)
         y_max = atom_pos(j,2)-atom_pos(i,2)
         z_max = atom_pos(j,3)-atom_pos(i,3)
        
         step_x=x_max/kmax
         if(step_x.gt.1.0D-16) delta_x=step_x
         step_y=y_max/kmax
         if(step_y.gt.1.0D-16) delta_y=step_y
         step_z=z_max/kmax
         if(step_z.gt.1.0D-16) delta_z=step_z

         x2 = 0.0
         y2 = 0.0
         z2 = 0.0
        
             do kz=1,kmax
                z2 = z2 + step_z
                x2 = x2 + step_x
                y2 = y2 + step_y
       
                tess_eq(i,j) = (x2-atom_pos(i,1))*(atom_pos(j,1)-atom_pos(i,1))
                tess_eq(i,j) = tess_eq(i,j) + (y2-atom_pos(i,2))*(atom_pos(j,2)-atom_pos(i,2))
                tess_eq(i,j) = tess_eq(i,j) + (z2-atom_pos(i,3))*(atom_pos(j,3)-atom_pos(i,3))
       
                tess_bound(i,j) = (atom_pos(j,1)-atom_pos(i,1))**2.D0
                tess_bound(i,j) = tess_bound(i,j) + (atom_pos(j,2)-atom_pos(i,2))**2.D0
                tess_bound(i,j) = tess_bound(i,j) + (atom_pos(j,3)-atom_pos(i,3))**2.D0
                tess_bound(i,j) = tess_bound(i,j)*alpha(i,j)

                n_step_vol = n_step_vol + 1

                if(abs(tess_eq(i,j)-tess_bound(i,j)).lt.1.D-3) then
                  x_bound(i,j) = x2
                  y_bound(i,j) = y2
                  z_bound(i,j) = z2

                  n_step_vol_fix = n_step_vol
                  write(*,*) n_step_vol_fix
                  n_step_vol = 0
                  n_step_vol_fix = 4.0D0*n_step_vol_fix**3.D0
                end if
             end do

         if(i.eq.1) then
           write(*,*) '!',j,x_bound(i,j),y_bound(i,j),z_bound(i,j),'!'
         end if
       end if
    end do
    vol_cell_rws = n_step_vol_fix*delta_x*delta_y*delta_z
    RWS_new(i) = (3.D0*vol_cell_rws/(4.0D0*pi))**(1.D0/3.D0)
    RWS_new(i) = RWS_new(i)*ALAT
    write(*,*) VRBAS,vol_cell_rws*ALAT**3.D0
    write(*,*) 'Atom',i,'RWS:',RWS_new(i)
  end do
  end if

else
  write(*,*) ''
  write(*,*) 'keine der Varianten gewaehlt'
  write(*,*) ''
  stop
end if

close(337)
close(338)
close(339)

!99001 format (1X,I4,2(F16.10),2(I5,E16.10))
99001 format (1X,I3,2x,ES16.10,2x,ES16.10,2x,I3,2x,ES16.10,2x,I3,2X,ES16.10)
99002 format ('IRMTALG (1) - eigene Routine')
99003 format ('IRMTALG (2) - Variante aus potio.f, nicht vollstaendig')
99004 format ('IRMTALG (3) - Variante aus potio.f')
99005 format ('IRMTALG (4) - experimentelle Variante')
99006 format ('IRMTALG (5) - Variante Hutsepot')
99007 format ('IRMTALG (6) - RWS Tesselierung')
99008 format (1x,I2,1x,I2,1x,f7.3,1x,f7.3)
99009 format (1x,A20,1x,f7.3,A1,f7.3,1x,A4,f7.3)

end program

subroutine readrws(Z,N,RWS_tab)

!implicit none

!read tabulated RWS (see xband rws.vst)

!#   ********************************************************************
!#   *                                                                  *
!#   *     the Wigner-Seitz-radius  RWS  for atomic number  Z           *
!#   *                                                                  *
!#   *     column 1: atomic number  Z                                   *
!#   *     column 2: experimental values  (incomplete)                  *
!#   *               taken from H.L.Skrivers's book                     *
!#   *     column 3: data used in the program spheres                   *
!#   *                                                                  *
!#   *     if the entry in column 2 is missing the value is taken       *
!#   *     from column 3                                                *
!#   *                                                                  *
!#   ********************************************************************

integer :: Z,Z_file
real*8 :: RWS_tab,RWS1,RWS2

open(333,file='tabulated_rws')

do
    
   read(333,*) Z_file,RWS1,RWS2
!   write(*,*) Z,Z_file,RWS1,RWS2

   if(Z_file.eq.Z) then
      
      RWS_tab = RWS1

      if(RWS_tab.lt.1.0D0) then
        
         RWS_tab = RWS2

      end if

   rewind(333)
   exit

   end if

end do

end

function DNRM2(N,X,INCX)

!implicit none

real*8 ONE,ZERO
parameter (ONE=1.0D+0,ZERO=0.0D+0)
integer INCX,N
real*8 X(*)

real*8 ABSXI,NORM,SCALEVAR,SSQ
integer IX

!  DNRM2 returns the euclidean norm of a vector via the function
!  name, so that DNRM2 := sqrt( x'*x )

      IF ( N.LT.1 .OR. INCX.LT.1 ) THEN
         NORM = ZERO
      ELSE IF ( N.EQ.1 ) THEN
         NORM = ABS(X(1))
      ELSE
         SCALEVAR = ZERO
         SSQ = ONE

         DO IX = 1,1 + (N-1)*INCX,INCX
            IF ( X(IX).NE.ZERO ) THEN
               ABSXI = ABS(X(IX))
               IF ( SCALEVAR.LT.ABSXI ) THEN
                  SSQ = ONE + SSQ*(SCALEVAR/ABSXI)**2
                  SCALEVAR = ABSXI
               ELSE
                  SSQ = SSQ + (ABSXI/SCALEVAR)**2
               END IF
            END IF
         END DO
         NORM = SCALEVAR*SQRT(SSQ)
      END IF

      DNRM2 = NORM

end

subroutine kugelschnitt(d,r1,r2,N,i,j)

implicit none

integer :: i,j,N
real*8 :: pi,a,vol_ks1,vol_ks2,d1
real*8 :: r1(N),r2(N)
real*8 :: d(N,N)

!allocate(r1(N),r2(N),d(N,N))

pi = 3.141592653589793238462643D0

write(*,*) i,j

!write(*,*) 'r1'
!read(5,*) r1

!write(*,*) 'r2'
!read(5,*) r2

!write(*,*) 'd'
!read(5,*) d

!Vol KS Kreis 1

write(*,*) r1(i),r2(j)
stop

if(d(i,j).lt.1.0D-16) then
  d1=1.0D-16
else
  d1=(r1(i)**2.D0-r2(j)**2.D0+d(i,j)**2.D0)/(2.D0*d(i,j)) !Abstand Kugeln 1 vom Mittelpkt Kreisegment
end if

if(r1(i)-d1.lt.(-1.0D-16)) then
  a=0.D0
else
  a=sqrt(r1(i)**2.D0-d1**2.D0) !Radius Kreisegment Schnittflaeche Kugeln
end if

write(*,*) 'a',a,r1(i),r2(j),d1,'!!!'
vol_ks1 = (r1(i)-sqrt(r1(i)**2.D0-a**2.D0))**2.D0 !Vol Kreissegment
vol_ks1 = vol_ks1*(2.D0*r1(i)+sqrt(r1(i)**2.D0-a**2.D0))
vol_ks1 = pi/3.D0*vol_ks1

vol_ks2 = (r2(j)-sqrt(r2(j)**2.D0-a**2.D0))**2.D0 !Vol Kreissegment
vol_ks2 = vol_ks2*(2.D0*r2(j)+sqrt(r2(j)**2.D0-a**2.D0))
vol_ks2 = pi/3.D0*vol_ks2

!write(*,*) 'vol_ks_1',vol_ks1
!write(*,*) 'vol_ks_2',vol_ks2
write(338,*) 'i,j,vol_ks_1+vol_ks_2',i,j,vol_ks1+vol_ks2

end
