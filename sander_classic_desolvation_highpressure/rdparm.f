#ifdef MEM_ALLOC
      SUBROUTINE RDPARM(xptr,ixptr,ihptr,ipptr,KFORM,NF)
#else
      SUBROUTINE RDPARM(x,ix,ih,ipairs,KFORM,NF)
#endif
C
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995                      **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
C
#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
#ifdef MEM_ALLOC
      pointer(xptr,  x)
      pointer(ixptr, ix)
      pointer(ihptr, ih)
      pointer(ipptr, ipairs)
      dimension x(1)
      integer   ix(1), ih(1), ipairs(1)
#else
      dimension x(*)
      integer   ix(*), ih(*), ipairs(*)
#endif
C
#include "md.h"
#include "memory.h"
#include "parms.h"
#include "files.h"
#include "box.h"
#include "nmr.h"
#include "pol.h"
#include "iewald.h"
#ifdef LES
#  include "les.h"
#endif
C
C     ----- READ THE MOLECULAR TOPOLOGY -----
C
      IPTRES = 0
      IPTATM = 0
      NSPSOL = 0
      NSPSTR = 0
      IF(KFORM.LE.0) GO TO 500
C
C     ----- FORMATTED INPUT -----
C
      READ(NF,9108) (ITITL(I),I=1,20)
      READ(NF,9118) NATOM,NTYPES,NBONH,MBONA,NTHETH,MTHETA,NPHIH,MPHIA,
     $        NHPARM,NPARM,NNB,NRES,NBONA,NTHETA,NPHIA,
     $        NUMBND,NUMANG,NPTRA,NATYP,NPHB,ifpert,IDUM,IDUM,IDUM,IDUM,
     $        IDUM,IDUM,IFBOX,NMXRS,IFCAP
      if (ifpert.ne.0) then
          write(6,*) ' *** PARM TOPOLOGY FILE SET UP FOR PERTURBATION'
          write(6,*) '     - ONLY USABLE WITH GIBBS'
          call mexit(6,1)
      endif
#ifdef LES
      if (nparm .ne. 1) then
         write(6,*) ' *** THIS VERSION ONLY ACCEPTS TOPOLOGY FILES'
         write(6,*) '     THAT WERE CREATED BY ADDLES, WITH NPARM=1'
         write(6,*) '     USE A VERSION COMPILED WITHOUT -DLES '
         call mexit(6,1)
      end if
#else
      if (nparm .eq. 1) then
         write(6,*) ' *** THIS VERSION WILL NOT ACCEPT TOPOLOGY FILES'
         write(6,*) '     THAT WERE CREATED BY ADDLES, WITH NPARM=1'
         write(6,*) '     USE A VERSION COMPILED WITH -DLES '
         call mexit(6,1)
      end if
#endif
#ifndef ARCHIVE
      write(6,8118) NATOM,NTYPES,NBONH,MBONA,NTHETH,MTHETA,NPHIH,MPHIA,
     $              NHPARM,NPARM,NNB,NRES,NBONA,NTHETA,NPHIA,
     $ NUMBND,NUMANG,NPTRA,NATYP,NPHB,IFBOX,NMXRS,IFCAP
#endif
 8118 format(t2,
     $  'NATOM  = ',i7,' NTYPES = ',i7,' NBONH = ',i7,' MBONA  = ',i7,
     $/' NTHETH = ',i7,' MTHETA = ',i7,' NPHIH = ',i7,' MPHIA  = ',i7,
     $/' NHPARM = ',i7,' NPARM  = ',i7,' NNB   = ',i7,' NRES   = ',i7,
     $/' NBONA  = ',i7,' NTHETA = ',i7,' NPHIA = ',i7,' NUMBND = ',i7,
     $/' NUMANG = ',i7,' NPTRA  = ',i7,' NATYP = ',i7,' NPHB   = ',i7,
     $/' IFBOX  = ',i7,' NMXRS  = ',i7,' IFCAP = ',i7/)
c
c     --- make sure we do not exceed memory limits in commons ---
c
      NTTYP = NTYPES*(NTYPES+1)/2
c
#ifdef MEM_ALLOC
C     TODO - check the implications of the arbitrary limits
#endif
c jianfa changed nttyp 2.20, 2014
      if (numbnd.gt.5000 .or. numang.gt.9000 .or. nptra.gt.2000 .or.
     $     nphb.gt.95000 .or. natyp.gt.10000 .or. nttyp.gt.500000) 
     $ go to 1000
C
C     ----- NOW ALL THE DATA IS AVAILABLE FOR MEMORY ALLOCATION.
C           DO THE PARTITION AND LOAD THE ARRAYS -----
C
#ifdef MEM_ALLOC
      call LOCMEM(xptr,ixptr,ihptr,ipptr)
#else
      call LOCMEM(x,ix,ih,ipairs)
#endif
C
      NTYPE = NTYPES*NTYPES
C
C     ----- READ THE SYMBOLS AND THE CHARGES AND THE MASSES -----
C
      READ(NF,9108) (IH(m04+I-1),I = 1,NATOM)
      READ(NF,9128) (X(L15+I-1),I = 1,NATOM)
      READ(NF,9128) (X(L20+I-1),I = 1,NATOM)
      READ(NF,9118) (IX(I04+I-1),I = 1,NATOM)
      READ(NF,9118) (IX(I+I08-1),I = 1,NATOM)
      READ(NF,9118) (IX(I+I06-1),I = 1,NTYPE)
      READ(NF,9108) (IH(I+m02-1),I=1,NRES)
      READ(NF,9118) (IX(I+I02-1),I=1,NRES)
      IX(I02+NRES) = NATOM+1
C
C     ----- READ THE PARAMETERS -----
C
      READ(NF,9128) (RK(I),    I = 1,NUMBND)
      READ(NF,9128) (REQ(I),   I = 1,NUMBND)
      READ(NF,9128) (TK(I),    I = 1,NUMANG)
      READ(NF,9128) (TEQ(I),   I = 1,NUMANG)
      READ(NF,9128) (PK(I),    I = 1,NPTRA)
      READ(NF,9128) (PN(I),    I = 1,NPTRA)
      READ(NF,9128) (PHASE(I), I = 1,NPTRA)
      READ(NF,9128) (SOLTY(I), I = 1,NATYP)
      READ(NF,9128) (CN1(I),   I = 1,NTTYP)
      READ(NF,9128) (CN2(I),   I = 1,NTTYP)
C
C     ----- READ THE BONDING INFORMATIONS -----
C
      READ(NF,9118) (IX(I+I12-1),IX(I+I14-1),IX(I+I16-1),I = 1,NBONH)
      READ(NF,9118) (IX(I+I18-1),IX(I+I20-1),IX(I+I22-1),I = 1,NBONA)
      READ(NF,9118) (IX(I+I24-1),IX(I+I26-1),IX(I+I28-1),IX(I+I30-1),
     $               I = 1,NTHETH)
      READ(NF,9118) (IX(I+I32-1),IX(I+I34-1),IX(I+I36-1),IX(I+I38-1),
     $               I = 1,NTHETA)
      READ(NF,9118) (IX(I+I40-1),IX(I+I42-1),IX(I+I44-1),IX(I+I46-1),
     $               IX(I+I48-1),I = 1,NPHIH)
      READ(NF,9118) (IX(I+I50-1),IX(I+I52-1),IX(I+I54-1),IX(I+I56-1),
     $               IX(I+I58-1),I = 1,NPHIA)
      READ(NF,9118) (IX(I+I10-1),I=1,NNB)
C
C     ----- READ THE H-BOND PARAMETERS -----
C
      READ(NF,9128) (ASOL(I),I=1,NPHB)
      READ(NF,9128) (BSOL(I),I=1,NPHB)
      READ(NF,9128) (HBCUT(I),I=1,NPHB)

c     If ewald, and if doing an accuracy test, need to zero vdw & 
c     hbond terms...
      if ( iewald .eq. 1 )
     $     call zero_if_acc_test(asol,bsol,nphb,cn1,cn2,nttyp)
c
C
C     ----- READ ISYMBL,ITREE,JOIN AND IROTAT ARRAYS -----
C
      READ(NF,9108) (IX(I+j01-1),I=1,NATOM)
      READ(NF,9108) (IH(I+m08-1),I=1,NATOM)
      READ(NF,9118) (IX(I+I64-1),I=1,NATOM)
      READ(NF,9118) (IX(I+I66-1),I=1,NATOM)
C
C     ----- READ THE BOUNDARY CONDITION STUFF -----
C
      NSPM = 1
      IX(I70) = NATOM
      IF(IFBOX.GT.0) THEN
        READ(NF,9118) IPTRES,NSPM,NSPSOL
        READ(NF,9118) (IX(I+I70-1),I=1,NSPM)
        READ(NF,9128) BETA,(BOX(I),I=1,3)
        IF(IPTRES.GT.0) IPTATM = IX(I02+IPTRES-1+1)-1
      END IF
C
C     ----- LOAD THE CAP INFORMATION IF NEEDED -----
C
      IF(IFCAP.GT.0) THEN
        READ(NF,9118) NATCAP
        READ(NF,9128) CUTCAP,XCAP,YCAP,ZCAP
        write(6,8128) natcap,cutcap,xcap,ycap,zcap
 8128   format(t2,'cap option on:',
     $   /,'NATCAP = ',I5,'  CUTCAP = ',F10.3,' XYZCAP = ',3F8.3)
      ELSE
         NATCAP = 0
         CUTCAP = 0.0D0
         XCAP = 0.0D0
         YCAP = 0.0D0
         ZCAP = 0.0D0
      ENDIF
c
      if(ipol.gt.0)then
      write(6,3333)nf
 3333   format(t2,
     .    '%MD/POL-I-POLREAD, reading polarizabilities from topology',
     $     i3/)
        read(nf,9128) (x(l25+i-1),i=1,natom)
      endif
#ifdef LES
      if (nparm.eq.1) then
        read (nf,9118) nlesty
        lestmp=nlesty*nlesty
c
c check the array sizes to make sure we don't overflow
c
c LES types
c
        if (nlesty.gt.maxlestyp) then
          write (6,*) 'Exceeded MAXLESTYP',nlesty
          stop
        endif
c
c LES atoms
c
        if (natom.gt.maxles) then
          write (6,*) 'Exceeded MAXLES',natom
          stop
        endif

        read (nf,9118) (lestyp(i),i=1,natom)
        read (nf,9128) (lesfac(i),i=1,lestmp)
        read (nf,9118) (cnum(i), i=1,natom)
        read (nf,9118) (subsp(i), i=1,natom)
        write (6,*) 'LES parameters were found'
c       write (6,9118) (lestyp(i),i=1,natom)
c       write (6,9128) (lesfac(i),i=1,lestmp)
      endif
#endif
      GO TO 520
C
C     ----- UNFORMATTED INPUT -----
C
  500 CONTINUE
      READ(NF) (ITITL(I),I=1,20)
      READ(NF) NATOM,NTYPES,NBONH,MBONA,NTHETH,MTHETA,NPHIH,MPHIA,
     $         NHPARM,NPARM,NNB,NRES,NBONA,NTHETA,NPHIA,
     $         NUMBND,NUMANG,NPTRA,NATYP,NPHB,ifpert,IDUM,IDUM,IDUM,
     $         IDUM,IDUM,IDUM,IFBOX,NMXRS,IFCAP
      if (ifpert.ne.0) then
          write(6,*) ' *** PARM TOPOLOGY FILE SET UP FOR PERTURBATION'
          write(6,*) '     - ONLY USABLE WITH GIBBS'
          call mexit(6,1)
      endif
      write(6,8118) NATOM,NTYPES,NBONH,MBONA,NTHETH,MTHETA,NPHIH,MPHIA,
     $              NHPARM,NPARM,NNB,NRES,NBONA,NTHETA,NPHIA,
     $              NUMBND,NUMANG,NPTRA,NATYP,NPHB,IFBOX,NMXRS,IFCAP
c
c     --- make sure we do not exceed memory limits in commons ---
c
      NTTYP = NTYPES*(NTYPES+1)/2
c
c    jianfa changed nttyp 2.20, 2014
 
      if (numbnd.gt.5000 .or. numang.gt.9000 .or. nptra.gt.2000 .or.
     $     nphb.gt.95000 .or. natyp.gt.10000 .or. nttyp.gt.500000)
     $ go to 1000
C
C
C     ----- NOW ALL THE DATA IS AVAILABLE FOR MEMORY ALLOCATION.
C           DO THE PARTITION AND LOAD THE ARRAYS -----
C
#ifdef MEM_ALLOC
      call LOCMEM(xptr,ixptr,ihptr,ipptr)
#else
      call LOCMEM(x,ix,ih,ipairs)
#endif
C
      NTYPE = NTYPES*NTYPES
C
C     ----- READ THE SYMBOLS AND THE CHARGES AND THE MASSES -----
C
      READ(NF) (ih(m04+I-1),I = 1,NATOM)
      READ(NF) (X(L15+I-1),I = 1,NATOM)
      READ(NF) (X(L20+I-1),I = 1,NATOM)
      READ(NF) (IX(I04+I-1),I = 1,NATOM)
      READ(NF) (IX(I+I08-1),I = 1,NATOM)
      READ(NF) (IX(I+I06-1),I = 1,NTYPE)
      READ(NF) (ih(I+m02-1),I=1,NRES)
      READ(NF) (IX(I+I02-1),I=1,NRES)
      IX(I02+NRES) = NATOM+1
C
C     ----- READ THE PARAMETERS -----
C
      READ(NF) (RK(I),    I = 1,NUMBND)
      READ(NF) (REQ(I),   I = 1,NUMBND)
      READ(NF) (TK(I),    I = 1,NUMANG)
      READ(NF) (TEQ(I),   I = 1,NUMANG)
      READ(NF) (PK(I),    I = 1,NPTRA)
      READ(NF) (PN(I),    I = 1,NPTRA)
      READ(NF) (PHASE(I), I = 1,NPTRA)
      READ(NF) (SOLTY(I), I = 1,NATYP)
      READ(NF) (CN1(I),   I = 1,NTTYP)
      READ(NF) (CN2(I),   I = 1,NTTYP)
C
C     ----- READ THE BONDING INFORMATIONS -----
C
      READ(NF) (IX(I+I12-1),IX(I+I14-1),IX(I+I16-1),I = 1,NBONH)
      READ(NF) (IX(I+I18-1),IX(I+I20-1),IX(I+I22-1),I = 1,NBONA)
      READ(NF) (IX(I+I24-1),IX(I+I26-1),IX(I+I28-1),IX(I+I30-1),
     $          I = 1,NTHETH)
      READ(NF) (IX(I+I32-1),IX(I+I34-1),IX(I+I36-1),IX(I+I38-1),
     $          I = 1,NTHETA)
      READ(NF) (IX(I+I40-1),IX(I+I42-1),IX(I+I44-1),IX(I+I46-1),
     $          IX(I+I48-1),I = 1,NPHIH)
      READ(NF) (IX(I+I50-1),IX(I+I52-1),IX(I+I54-1),IX(I+I56-1),
     $          IX(I+I58-1),I = 1,NPHIA)
      READ(NF) (IX(I+I10-1),I=1,NNB)
C
C     ----- READ THE H-BOND PARAMETERS -----
C
      READ(NF) (ASOL(I),I=1,NPHB)
      READ(NF) (BSOL(I),I=1,NPHB)
      READ(NF) (HBCUT(I),I=1,NPHB)
C
C     ----- READ ISYMBL,ITREE,JOIN AND IROTAT ARRAYS -----
C
      READ(NF) (ix(I+j01-1),I=1,NATOM)
      READ(NF) (ih(I+m08-1),I=1,NATOM)
      READ(NF) (IX(I+I64-1),I=1,NATOM)
      READ(NF) (IX(I+I66-1),I=1,NATOM)
C
C     ----- READ THE BOUNDARY STUFF -----
C
      NSPM = 1
      IX(I70) = NATOM
      IF(IFBOX.GT.0) THEN
        READ(NF) IPTRES,NSPM,NSPSOL
        READ(NF) (IX(I+I70-1),I=1,NSPM)
        READ(NF) BETA,(BOX(I),I=1,3)
        IF(IPTRES.GT.0) IPTATM = IX(I02+IPTRES-1+1)-1
      END IF
C
C     ----- LOAD THE CAP INFORMATION IF NEEDED -----
C
      IF(IFCAP.GT.0) THEN
        READ(NF) NATCAP,CUTCAP,XCAP,YCAP,ZCAP
        write(6,8128) natcap,cutcap,xcap,ycap,zcap
      ELSE
         NATCAP = 0
         CUTCAP = 0.0D0
         XCAP = 0.0D0
         YCAP = 0.0D0
         ZCAP = 0.0D0
      ENDIF
c
      if(ipol.gt.0)then
        write(6,3333)nf
        read(nf) (x(l25+i-1),i=1,natom)
      endif
  520 CONTINUE
      if(ipol.ne.0) then
        ip14 = 0
        do 511 i = 1,nphih
          i1 = IX(I+I40-1)/3+1
          i2 = IX(I+I42-1)/3+1
          i3 = IX(I+I44-1)/3+1
          i4 = IX(I+I46-1)/3+1
          if(ix(i+i44-1).ge.0.and.ix(i+i46-1).ge.0)then
               ih(i1+m12-1) = ih(i1+m12-1) + 1
               call istuff(ih(i1+m12-1),i1,ih(m14),i4)
               ip14 = ip14+1
          endif
  511   continue
        do 512 i = 1,nphia
          i1 = IX(I+I50-1)/3+1
          i2 = IX(I+I52-1)/3+1
          i3 = abs(IX(I+I54-1))/3+1
          i4 = abs(IX(I+I56-1))/3+1
          if(ix(i+i54-1).ge.0.and.ix(i+i56-1).ge.0) then
               ih(i1+m12-1) = ih(i1+m12-1) + 1
               if(ih(i1+m12-1).gt.15) then
                  print *,'too many dihedrals '
                  call mexit(6,1)
               endif
               call istuff(ih(i1+m12-1),i1,ih(m14),i4)
               ip14 = ip14+1
          endif
  512   continue

        write(6,3377)ip14
 3377   format(t2,i4,' pairs for pol 1-4.'/)
      endif
C
C     ----- CALCULATE INVERSE, TOTAL MASSES -----
c
c       -- save the masses for removal of mass weighted velocity,
c          leaving the inverse masses in the legacy, L20 area
c
C
      tmass = 0.0d0
c     -- index over molecules
      j = L75-1
      jj = I70-1
c     -- index over mass->invmass
      k = L20-1
c     -- index over saved mass
      l = L70-1
      DO N = 1,NSPM
        j = j + 1
        jj = jj + 1
        X(j) = 0.0d0
        NATSM = IX(jj)
        DO NN = 1,NATSM
          k = k+1
          l = l+1
c
c         -- sum molecule
c
          X(j) = X(j) + X(k)
c
c         -- save mass in "new" L70 area
c
          X(l) = X(k)
c
c         -- make inverse in "old" L20 area
c
          X(k) = 1.0d0 / X(k)
        enddo
        tmass = tmass + X(j)
      enddo
      tmassinv = 1.0d0 / tmass
C
C     ----- SCALE THE CHARGES IF DIELC.GT.1.0E0 -----
C
      IF (DIELC .ne. 1.0E0) THEN
           DUMD = SQRT(DIELC)
           DO 140 I = 1,NATOM
                X(I+L15-1) = X(I+L15-1)/DUMD
  140      CONTINUE
      ENDIF
C
C     ----- GENERATE THE SOLVENT BEGINNING INDEX -----
C
      IF (NSPSOL .GT. 1) THEN
           DO 170 I = 1,NSPSOL-1
                NSPSTR = NSPSTR+IX(I70+I-1)
  170      CONTINUE
      ENDIF
C
C     ----- INVERT THE HBCUT ARRAY -----
C
      DO 180 I = 1,NPHB
           IF(HBCUT(I).LE.0.001E0) HBCUT(I) = 1.0E-10
           HBCUT(I) = 1.0E0/HBCUT(I)
  180 CONTINUE
c
c     ----- duplicate dihedral pointers for vector ephi -----
c
      idum = nphih
      call dihdup(nphih,idum,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48),pn)
      call dihdup(mphia,nphia,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58),
     $            pn)
c
c     --- pre-calculate some parameters for vector ephi ---
c
      CALL DIHPAR(NPTRA,PK,PN,PHASE,GAMC,GAMS,IPN,FMN)
C
#ifdef CHARMM
c
c    ---read in   1-4 parameters at end of parm file:
c
      if (kform.le.0) then
         read(NF) (CN114(I),   I = 1,NTTYP)
         read(NF) (CN214(I),   I = 1,NTTYP)
      else
         read(NF,9128) (CN114(I),   I = 1,NTTYP)
         read(NF,9128) (CN214(I),   I = 1,NTTYP)
      end if
c
c    --- read in   urey-bradley parameters:
c
      if (kform.le.0) then
        read(nf) (rkub(i),i=1,numang)
        read(nf) (rub(i),i=1,numang)
      else
        read(nf,9128) (rkub(i),i=1,numang)
        read(nf,9128) (rub(i),i=1,numang)
      end if
#endif
      RETURN
 1000 continue
      write(6,'(/,5x,a)') 'rdparm: a parameter array overflowed'
      write(6,'(/,5x,a)') '        (e.g. the table of dihedral params)'
      call mexit(6, 1)
 9108 FORMAT(20A4)
 9118 FORMAT(12I6)
 9128 FORMAT(5E16.8)
      END
C
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995                      **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
      subroutine istuff(i,j,iarray,k)
C
c routine to correctly load a strange shaped 2 dim matrix 
c
      dimension iarray(15,*)
      iarray(i,j) = k
      return
      end
