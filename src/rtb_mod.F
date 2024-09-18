      module rtb
c
c Rotational Translational Block method for large scale normal mode analysis
c
c Copyright (c) 2000-2005 Florence Tama, Yves-Henri Sanejouand
c Developed at Paul Sabatier University, Toulouse, France
c
c Copyright (c) 2024 Florence Tama, Osamu Miyashita (block cutoff, f90)
c Developed at Nagoya University, Nagoya, Japan
c              RIKEN, Kobe, Japan
c
c References:
c 
c RTB method
c 
c Durand, P, Trinquier, G & Sanejouand, YH. New Approach for Determining
c Low-Frequency Normal-Modes in Macromolecules. Biopolymers 34: 759-71
c (1994).
c 
c Tama, F, Gadea, FX, Marques, O & Sanejouand, YH. Building-block
c approach for determining low-frequency normal modes of
c macromolecules. Proteins 41: 1-7 (2000).
c 
c Elastic network model
c 
c Tirion, MM. Large amplitude elastic motions in proteins from a single-
c parameter, atomic analysis. Phys Rev Lett 77: 1905-8 (1996). 
c
      contains

c -----------------------------------------------
      SUBROUTINE SCHMIDT(M,N,C,NMAX)
c ----------------------------------------------- 
C
C     ORTHOGONALISATION PAR GRAM-SCHMIDT
C
      IMPLICIT REAL*8(A-H,O-Z)
c      PARAMETER (NMAX=NDDPERB)   
      DIMENSION C(6,NMAX),REC(6,6)
C
      ANORM = 0.D0
      DO I=1,N
       ANORM = ANORM + C(1,I)*C(1,I)
      ENDDO
      ANORM = 1.D0 / ( DSQRT(ANORM) )
C
      DO I=1,N
       C(1,I) = ANORM * C(1,I)
      ENDDO
C 
      DO  I=2,M      
       DO J=1,I-1
        REC(J,I) = 0.D0
        DO K=1,N
         REC(J,I) = REC(J,I) + C(J,K)*C(I,K)
        ENDDO
       ENDDO
       DO K=1,N
        AAA = 0.D0
        DO J=1,I-1
         AAA = AAA + C(J,K)*REC(J,I)
        ENDDO
        C(I,K) = C(I,K) - AAA
       ENDDO
C
       ANORM = 0.D0
       DO K=1,N
        ANORM = ANORM + C(I,K)*C(I,K)
       ENDDO 
       ANORM = 1.D0 / ( DSQRT(ANORM) )
C
       DO K=1,N
        C(I,K) = ANORM*C(I,K)
       ENDDO
C
      ENDDO

      RETURN
      END  
C
      subroutine time_stamp(arg)
      character*20 arg
      real dta(2),dt
      
      call dtime(dta,dt)
      write(6,*) arg
      write(6,'(1X,A5,1X,F12.2,1X,F12.2,1X,F12.2)') 'time>'
     $     ,dta(1), dta(2), dta(1)+dta(2)
      return
      end


      subroutine readpdb_count_atoms(natom)
         integer, intent(out) :: natom

         integer :: i
         character*80 lign80
         integer ios

         open (unit=75, file='pdb', status='old', form='formatted')
         i = 0
         do
            read (75, '(A)', iostat=ios) lign80
            if (ios /= 0) exit
            if (lign80(1:4) == 'ATOM') then
               i = i + 1
            end if
          end do
         close(75)
         natom = i
      end subroutine


      subroutine readpdb(at,atname,resname,resi,
     .     xat,yat,zat,tot,massat,sid,natom,nb)

      implicit none
      
      integer natom, nb
      character*4 atname(natom),resname(natom)
      integer sid(natom),at(natom),resi(natom)
      real*8 massat(natom),xat(natom),yat(natom),zat(natom),tot(natom)

      integer i
      character*80 lign80

      open (unit=75, file='pdb', status='old', form='formatted')

      i=0
 25   continue
      read(75,'(A)',end=30) lign80
      if (lign80(1:4).ne.'ATOM') goto 25
      i=i+1

      read(lign80,
     .     '(6X,I5,1X,A4,1X,A4,1X,I4,4X,F8.3,F8.3,
     .F8.3,1X,F5.2,F6.2,6X,I4)')
     .     at(i),atname(i),resname(i),resi(i),xat(i),yat(i),
     .     zat(i),tot(i),massat(i),sid(i)

      goto 25
 30   continue

      close(25)

      natom=i
      nb=sid(natom)

      return
      end

      subroutine setlgbloc(sid,n,lgbloc,natom,nb)

      integer, intent(in) :: sid(natom)   ! segids in pdb
      integer, intent(out) :: n(nb)   ! dimension of each block
      integer, intent(out) :: lgbloc   ! largest block
      integer, intent(in) :: natom, nb

      integer lgblocmax, i, nn, test

      test=sid(1)
      nn=1
      lgblocmax = -999
      lgbloc = 0
      do i=1,natom

         if (sid(i).eq.test) then
            lgbloc=lgbloc+1
            if(lgbloc .gt. lgblocmax) lgblocmax = lgbloc
         else
            n(nn)=3*lgbloc
            nn=nn+1
            if ( nn .gt. nb) then
               write(6,*) 'error number of block larger than nb:' ,nb 
               stop
            endif
            lgbloc=1
            test=sid(i)
         endif
      enddo

      n(nn)=lgbloc*3

      lgbloc = lgblocmax*3

      return
      end

      subroutine calcbloc(amass, xat,yat,zat,
     .     n, RT,
     .     nbrt,nb,natom,lgbloc,cutoff,kij, si, sj, sz)

      implicit none

      integer, intent(in) :: natom   ! number of atoms
      integer, intent(in) :: nb      ! number of blocks
      integer, intent(in) :: lgbloc  ! largest block size
      integer, intent(in) :: nbrt    ! degree of freedom of one block == 6
      integer, intent(in) :: n(nb)   ! 3 * number of atoms in each block
      real*8, intent(in) :: xat(natom),yat(natom),zat(natom)
      real*8, intent(in) :: kij      ! Tirion potential force constant
      real*8, intent(in) :: cutoff   ! Tirion potential cutoff
      
      real*8, intent(out) :: RT(nbrt,lgbloc,nb)
      real*8, intent(out) :: amass(3*natom)
      integer, intent(out) :: si(:), sj(:)
      real*8, intent(out) :: sz(:)
      
      real*8 :: corlin(3*natom)
      real*8 :: corlinm(3*natom)
      real*8 :: HT(3*natom,9)
      real*8 :: HH(lgbloc,lgbloc)
      real*8 :: HD(lgbloc,nbrt)
      real*8 :: D(nbrt,lgbloc)

#ifdef BLOCK_CUTOFF
      real*8 :: rmax_block(1:nb)
      real*8 :: cx_block(1:nb)
      real*8 :: cy_block(1:nb)
      real*8 :: cz_block(1:nb)
#endif

      real*8 :: s(nbrt,nbrt)
      real*8 :: xg, yg, zg
      real*8 :: rx, ry, rz

      integer ntoti,ntotj,itot,jtot,ntot
      integer i, j, ii, jj, mmm
      integer ni, nj, indxi, indxj, ibloc, jbloc

      integer  k
      integer  comp
      integer  ct

      real*8 cutoff_block
      real*8 dist,dist2,amasstot

      integer max_sn

      write(6,*) 'starting calcbloc'

      max_sn = size(sz)

c    change xyz in corlin

      do i=1,natom
        ii=3*i-2
        corlin(ii)=xat(i)
        corlin(ii+1)=yat(i)
        corlin(ii+2)=zat(i)
        amass(ii)=1.00
        amass(ii+1)=1.00
        amass(ii+2)=1.00
      enddo

C     
C     CHANGEMENT D'ORIGINE DES COORDONNEES CARTESIENNES 
C     
      AMASSTOT = 0.D0 
      DO I=1,natom
        AMASSTOT = AMASSTOT + AMASS(I)
      ENDDO      
      AMASSTOT = AMASSTOT / 3.D0 
C     
      XG = 0.D0
      YG = 0.D0  
      ZG = 0.D0
      DO I=0,(natom-1 )
         XG = XG + AMASS(3*I+1)*CORLIN(3*I+1)
         YG = YG + AMASS(3*I+2)*CORLIN(3*I+2)
         ZG = ZG + AMASS(3*I+3)*CORLIN(3*I+3)
      ENDDO
      XG = XG / AMASSTOT
      YG = YG / AMASSTOT
      ZG = ZG / AMASSTOT
C     
      DO I=0,(natom-1)      
         CORLINM(3*I+1)=CORLIN(3*I+1) - XG
         CORLINM(3*I+2)=CORLIN(3*I+2) - YG
         CORLINM(3*I+3)=CORLIN(3*I+3) - ZG
c	 write(94,*)CORLINM(3*I+1),CORLINM(3*I+2),CORLINM(3*I+3) 
      ENDDO
      write(6,*) 'eng of origin shift'
c      call time_stamp(" ")

C     Calcul des elements de matrice pour i = j
C
      write(6,*) 'begin of diagonal of hessian'
      do i=1,natom
          ii=3*i-2
          HT(ii,1)=0.d0
          HT(ii,2)=0.d0
          HT(ii,3)=0.d0
          HT(ii+1,4)=0.d0
          HT(ii+1,5)=0.d0
          HT(ii+1,6)=0.d0
          HT(ii+2,7)=0.d0
          HT(ii+2,8)=0.d0
          HT(ii+2,9)=0.d0
       enddo

       comp=0
       do i=1,natom-1
          ii=3*i-2
          do j=i+1, natom
             jj=3*j-2
                dist2=0.d0
                dist=0.d0
                rx=0.d0
                ry=0.d0
                rz=0.d0
                rx=xat(i)-xat(j)
                if (abs(rx).le.cutoff) then
                   ry=yat(i)-yat(j)
                   if (abs(ry).le.cutoff) then
                      rz=zat(i)-zat(j)
                      if (abs(rz).le.cutoff) then
                      dist2=rx*rx+ry*ry+rz*rz
                      if(dist2.le.(cutoff**2)) then
                         comp=comp+1
                         HT(ii,1)=HT(ii,1)+rx*rx*kij/dist2
                         HT(ii,2)=HT(ii,2)+rx*ry*kij/dist2
                         HT(ii,3)=HT(ii,3)+rx*rz*kij/dist2
                         HT(ii+1,4)=HT(ii+1,4)+rx*ry*kij/dist2
                         HT(ii+1,5)=HT(ii+1,5)+ry*ry*kij/dist2
                         HT(ii+1,6)=HT(ii+1,6)+ry*rz*kij/dist2
                         HT(ii+2,7)=HT(ii+2,7)+rz*rx*kij/dist2
                         HT(ii+2,8)=HT(ii+2,8)+rz*ry*kij/dist2
                         HT(ii+2,9)=HT(ii+2,9)+rz*rz*kij/dist2
                         HT(jj,1)=HT(jj,1)+rx*rx*kij/dist2
                         HT(jj,2)=HT(jj,2)+rx*ry*kij/dist2
                         HT(jj,3)=HT(jj,3)+rx*rz*kij/dist2
                         HT(jj+1,4)=HT(jj+1,4)+rx*ry*kij/dist2
                         HT(jj+1,5)=HT(jj+1,5)+ry*ry*kij/dist2
                         HT(jj+1,6)=HT(jj+1,6)+ry*rz*kij/dist2
                         HT(jj+2,7)=HT(jj+2,7)+rz*rx*kij/dist2
                         HT(jj+2,8)=HT(jj+2,8)+rz*ry*kij/dist2
                         HT(jj+2,9)=HT(jj+2,9)+rz*rz*kij/dist2

                      endif
                     endif
                   endif
                endif
             enddo
       enddo
       write(6,*) 'number of interactions =' , comp
       write(6,*) 'end of diagonal of hessian'
c       call time_stamp(" ")
c      write(31,*) HT

C     <------------------------------------------ DEBUT DE BOUCLE SUR LES BLOCS
C
       NTOT=0
       
       write(6,*) 'begin making RT'
       DO IBLOC=1,NB  
C     
C     REMPLISSAGE DES BLOCS

c     <------------------------------------------Initialisation du bloc

          do i=1,N(IBLOC)
             do j=1,N(IBLOC)
                HH(i,j)=0.d0
             enddo
          enddo
          
C 
C     CHANGEMENT D'ORIGINE DANS CHAQUE BLOC (CDG DU BLOC)
c     ---------------------------------------------------------

          AMASSTOT = 0.D0
          DO I=1,N(IBLOC)
             AMASSTOT = AMASSTOT + AMASS(NTOT+I)
          ENDDO      
          AMASSTOT=AMASSTOT/3.d0

C     
          XG = 0.D0
          YG = 0.D0  
          ZG = 0.D0
          
          DO I=0,( (N(IBLOC)/3)-1 )
             XG = XG + AMASS(NTOT+3*I+1)*CORLINM(NTOT+3*I+1)
             YG = YG + AMASS(NTOT+3*I+2)*CORLINM(NTOT+3*I+2)
             ZG = ZG + AMASS(NTOT+3*I+3)*CORLINM(NTOT+3*I+3)
          ENDDO
          
          XG = XG / AMASSTOT
          YG = YG / AMASSTOT
          ZG = ZG / AMASSTOT
C     
#ifdef BLOCK_CUTOFF
          cx_block(ibloc)=xg
          cy_block(ibloc)=yg
          cz_block(ibloc)=zg
          
          dist=0d0
#endif
          DO I=0,( (N(IBLOC)/3)-1)      
             CORLINM(NTOT+3*I+1)= CORLINM(NTOT+3*I+1) - XG
             CORLINM(NTOT+3*I+2)= CORLINM(NTOT+3*I+2) - YG
             CORLINM(NTOT+3*I+3)= CORLINM(NTOT+3*I+3) - ZG
c             write(93,*) CORLINM(NTOT+3*I+1),
c     .      CORLINM(NTOT+3*I+2),CORLINM(NTOT+3*I+3)
#ifdef BLOCK_CUTOFF 
             rx = corlinm(ntot+3*i+1)
             ry = corlinm(ntot+3*i+2)
             rz = corlinm(ntot+3*i+3)
             if( dist**2 .lt. rx**2 + ry**2 + rz**2 ) then
                dist = dsqrt( rx**2 + ry**2 + rz**2 )
             endif
#endif

          ENDDO

#ifdef BLOCK_CUTOFF          
          rmax_block(ibloc) = dist
#endif
          
C     
C     PROJECTION DES TRANSLATIONS ET ROTATIONS --------------------------
C     

          DO I=1,6
             DO J=1,N(IBLOC)
                D(I,J) = 0.D0 
                RT(I,J, ibloc) = 0.D0
             ENDDO
          ENDDO
C     
          DO I=1, N(IBLOC)/3
             II=NTOT 
             D(1,1+3*(I-1))= DSQRT(AMASS(II+1+3*(I-1)))
             D(2,2+3*(I-1))= DSQRT(AMASS(II+2+3*(I-1)))
             D(3,3+3*(I-1))= DSQRT(AMASS(II+3+3*(I-1)))
             D(4,2+3*(I-1))=
     $            -DSQRT(AMASS(II+2+3*(I-1)))*CORLINM(II+3+3*(I-1))
c             write(95,*) AMASS(II+2+3*(I-1)),CORLINM(II+3+3*(I-1))
             D(4,3+3*(I-1))=
     $            DSQRT(AMASS(II+3+3*(I-1)))*CORLINM(II+2+3*(I-1))
             D(5,1+3*(I-1))=
     $            DSQRT(AMASS(II+1+3*(I-1)))*CORLINM(II+3+3*(I-1))
c             write(95,*) AMASS(II+1+3*(I-1)),CORLINM(II+3+3*(I-1))
             D(5,3+3*(I-1))=
     $            -DSQRT(AMASS(II+3+3*(I-1)))*CORLINM(II+1+3*(I-1))
             D(6,1+3*(I-1))=
     $            -DSQRT(AMASS(II+1+3*(I-1)))*CORLINM(II+2+3*(I-1))
             D(6,2+3*(I-1))=
     $            DSQRT(AMASS(II+2+3*(I-1)))*CORLINM(II+1+3*(I-1))
          ENDDO
C
          MMM=6

          CALL SCHMIDT(MMM,N(IBLOC),D,lgbloc) 
C
          do i=1,6
             do j=1,N(IBLOC)
                RT(i,j, ibloc)=D(i,j)
             enddo
          enddo

c          write(6, *) ibloc, ( (rt(i,j),j=1,n(ibloc)),i=1,6)

          NTOT=NTOT+N(IBLOC)

       ENDDO
 
       write(6,*) 'end of making RT'
c      call time_stamp(" ")
C       
C
C   <-------------------------------------------FIN DE BOUCLE SUR LES BLOCS

      write(6,*) 'begin final step'

      indxi=0
      indxj=0
      NTOTI=0
      NTOTJ=0
      si(:)=0
      sj(:)=0
      sz(:)=0.0
      ct = 0
      do IBLOC=1,NB
         do JBLOC=IBLOC,NB

#ifdef BLOCK_CUTOFF
            cutoff_block =
     $           rmax_block(ibloc) + rmax_block(jbloc) + cutoff

            rx = cx_block(ibloc) - cx_block(jbloc)
            if( abs(rx) .gt. cutoff_block ) then
               NTOTJ=NTOTJ+N(JBLOC)
               indxj=indxj+6
               cycle
c     don't forget to increment index
            endif
            
            ry = cy_block(ibloc) - cy_block(jbloc)
            if( abs(ry) .gt. cutoff_block ) then
               NTOTJ=NTOTJ+N(JBLOC)
               indxj=indxj+6
               cycle
            endif
            
            rz = cz_block(ibloc) - cz_block(jbloc)
            if( abs(rz) .gt. cutoff_block ) then
               NTOTJ=NTOTJ+N(JBLOC)
               indxj=indxj+6
               cycle
            endif
            
            if( rx**2 + ry**2 + rz**2 .gt. cutoff_block**2 ) then
               NTOTJ=NTOTJ+N(JBLOC)
               indxj=indxj+6
               cycle
            endif

#ifdef DEBUG
            if( (ibloc .eq. 1) .and. (jbloc.eq.10) ) then
               write(6,*) 'blocks',ibloc,
     $              cx_block(ibloc),cy_block(ibloc),cz_block(ibloc)
               write(6,*) 'blocks',jbloc,
     $              cx_block(jbloc),cy_block(jbloc),cz_block(jbloc)
               rx = cx_block(ibloc) - cx_block(jbloc)
               ry = cy_block(ibloc) - cy_block(jbloc)
               rz = cz_block(ibloc) - cz_block(jbloc)
               write(6,*) rx, ry, rz
               write(6,*) sqrt(rx**2 + ry**2 + rz**2)
               write(6,*) cutoff_block
               write(6,*) rmax_block(ibloc), rmax_block(jbloc), cutoff
            endif
#endif
               
#endif

ccc   initialize
            do i=1,N(IBLOC)
               do j=1,N(JBLOC)
                  HH(i,j)=0.d0
               enddo
            enddo
            
c         call time_stamp("b hessian")
            do i=1,N(IBLOC)/3
               ii=3*i-2
               do j=1,N(JBLOC)/3
                  jj=3*j-2
                  
                  itot=0
                  jtot=0
                  
                  itot=i+ntoti/3 
                  jtot=j+ntotj/3
#ifdef DEBUG
                  if( (ibloc .eq. 1) .and. (jbloc.eq.10) ) then
                     write(6,*) 'block_atomid',ibloc, jbloc, itot,jtot
                  endif
#endif
                  if (itot.ne.jtot) then
                     rx=0.d0
                     ry=0.d0
                     rz=0.d0
                     dist2=0.d0
                     dist=0.d0
                     rx=xat(itot)-xat(jtot)
                     if(abs(rx).le.cutoff) then
                        ry=yat(itot)-yat(jtot)
                        if(abs(ry).le.cutoff) then
                           rz=zat(itot)-zat(jtot)
                           if(abs(rz).le.cutoff) then
                              dist2=rx*rx+ry*ry+rz*rz
                              dist=dsqrt(dist2)
                              if (dist.le.cutoff) then
#ifdef DEBUG
                              if( (ibloc .eq. 1) .and. 
     .                                (jbloc .eq. 10 ) ) then
                                 write(6,*) itot, xat(itot), 
     .                                yat(itot), zat(itot)
                                 write(6,*) jtot, xat(jtot), 
     .                                yat(jtot), zat(jtot)
                                 write(6,*) dist, cutoff
                                 write(6,*) ibloc
                                 write(6,*) jbloc
                              endif
#endif
                                 HH(ii,jj)=-rx*rx*kij/dist2
                                 HH(ii+1,jj)=-rx*ry*kij/dist2
                                 HH(ii+1,jj+2)=-ry*rz*kij/dist2
                                 HH(ii+2,jj)=-rz*rx*kij/dist2
                                 HH(ii+2,jj+1)=-rz*ry*kij/dist2
                                 HH(ii,jj+1)=-rx*ry*kij/dist2
                                 HH(ii,jj+2)=-rx*rz*kij/dist2
                                 HH(ii+1,jj+1)=-ry*ry*kij/dist2
                                 HH(ii+2,jj+2)=-rz*rz*kij/dist2
                              endif
                           endif
                        endif
                     endif   
                  endif
                  
                  if (itot.eq.jtot) then
                     HH(ii,jj)=1e-8+HT(3*itot-2,1)
                     HH(ii,jj+1)=HT(3*itot-2,2)
                     HH(ii,jj+2)=HT(3*itot-2,3)
                     HH(ii+1,jj)=HT(3*itot-1,4)
                     HH(ii+1,jj+1)=1e-8+HT(3*itot-1,5)
                     HH(ii+1,jj+2)=HT(3*itot-1,6)
                     HH(ii+2,jj)=HT(3*itot,7)
                     HH(ii+2,jj+1)=HT(3*itot,8)
                     HH(ii+2,jj+2)=1e-8+HT(3*itot,9)
                  endif
                  
               enddo
            enddo
            
c         call time_stamp("through loop")

c         call time_stamp("first read")
                        
#ifdef BLAS
C     SOMETIMES FAST
            call DGEMM('N','T',N(IBLOC),6,N(JBLOC), 1d0, HH(1,1),lgbloc,
     $           RT(1,1),nbrt,0.0d0,HD,lgbloc)
#else
            do i=1,6
               do j=1,N(IBLOC)
                  HD(j,i)=0.d0
                  do k=1,N(JBLOC)
                     HD(j,i)=HD(j,i)+HH(j,k)*RT(i,k, jbloc)
#ifdef DEBUG
                     write(*,*) 'j i HD(j,i)', j,i,HD(j,i)
                     write(*,*) 'HH(j,i)',HH(j,k)
                     write(*,*) 'RT(i,k)',RT(i,k)
#endif
                  enddo
               enddo
            enddo
#endif
            
#ifdef BLAS
            call DGEMM('N','N',6,6,N(IBLOC),1d0,RT(1,1),nbrt,
     $           HD,lgbloc,0d0,s,nbrt)
#else 
            do i=1,6
               do j=1,6
                  s(i,j)=0.d0
                  do k=1,N(IBLOC)
                     s(i,j)=s(i,j)+RT(i,k, ibloc)*HD(k,j)
                  enddo
               enddo
            enddo
#endif
            do i=1, 6
               do j=1,6
                  ni=i+indxi
                  nj=j+indxj
                  if (nj.ge.ni) then
                     if (s(i,j).ne.0.d0) then
                        ct = ct + 1
                        if( ct .gt. max_sn ) then
                           print *, 'inrease max_sn', max_sn
                           stop
                        endif
                        si(ct) = ni
                        sj(ct) = nj
                        sz(ct) = s(i,j)
#ifdef DEBUG      
                        write(*,*) 'ni nj s(i,j)',ni,nj,s(i,j)
                        write(*,*) 'ni nj iblock jbloc', ni,nj,
     .                       ibloc,jbloc
#endif
                     endif
                  endif
               enddo
            enddo
            
            NTOTJ=NTOTJ+N(JBLOC)
            indxj=indxj+6
            
         ENDDO
         
         indxi=indxi+6
         indxj=indxi
         
         NTOTI=NTOTI+N(IBLOC)
         NTOTJ=NTOTI
         
         
      ENDDO
#ifdef DEBUG      
      write(6,*) s(1,1),s(1,2),s(2,1),s(2,3)
#endif
      
      write(6,*) 'end final step'
      write(6,*) 'number of hessian elements =', ct

      return
      end



      subroutine rtb2mode(n,matvec,rt,natom,nb,
     .     lblocg,nbrt,nvectot,nstep)

      integer, intent(in) ::  nb,natom,nbrt,lblocg
      integer, intent(in) :: n(nb)   ! size of each block
      real*8, intent(in) :: rt(nbrt,lblocg,nb)
      real*8, intent(in) :: matvec(6*nb,nvectot)
      integer, intent(in) :: nstep   ! label in output mode file name
      integer, intent(in) :: nvectot   ! total number of modes to covert

      real*8 covec(3*natom)

      integer nvec, ntot
      integer ibloc, i, j, ii

      real*8 norm
      character*80 filename

      write(6,*)' nvectot =',nvectot

      do nvec=1,nvectot
         
         write(filename, '(A, I3.3, A, I3.3)') 
     &         'mov', nstep, '.mod', nvec
         
         open(unit=50,file=filename,status='unknown',form='formatted')
                  
         ntot=0
         ii=0         
         do ibloc=1,nb
            do j=1,n(ibloc)
               covec(j+ntot)=0.d0
               do i=1,6
                  covec(j+ntot)=covec(j+ntot)+
     .                 rt(i,j,ibloc)*matvec(i+ii,nvec)
               enddo
            enddo
            ii=ii+6
            ntot=ntot+n(ibloc)
         
         enddo
         
         do i=1,natom
            ii=3*i-2
            write(50,*) covec(ii) ,covec(ii+1) ,covec(ii+2)
         enddo
         close(50)
         
         norm = 0
         do i=1,natom
            ii=3*i-2
            norm=norm+covec(ii)**2+covec(ii+1)**2
     .           +covec(ii+2)**2
         enddo
         
         norm=dsqrt(norm)
         write(6,*) 'vec, norm =' , nvec, norm
         norm=0
         
         
      enddo
      
      close(51)

      return
      end

      end module rtb
