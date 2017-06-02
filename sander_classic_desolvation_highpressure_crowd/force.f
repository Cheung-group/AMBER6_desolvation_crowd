#ifdef SGI_MP
# define SHARED_MEMORY
#endif
#ifdef CNX_MP
# define SHARED_MEMORY
#endif
#ifdef CRAY_MP
# define SHARED_MEMORY
#endif
#ifdef SHARED_MEMORY
c     ========================= SHARED MEMORY =========================
c     
c     SGI_MP: Rather than clutter up the standard version of FORCE with
c     revisions necessary for the shared memory implementation of AMBER,
c     an additional copy of force is placed into the file "forcemp.f"
c     and included here.  When the #ifdef SGI_MP define token is turned
c     on at CPP time, "forcemp.f" is included and compiled instead of
c     force defined below.  
c
c     NOTE: If you change force.f, you may also need to change the 
c     code in forcemp.f.  
c
#include "forcemp.f"
c
c     ======================= END SHARED MEMORY =======================
#else
#include "vector.h"
c-----------------------------------------------------------------------
# ifdef MEM_ALLOC
      SUBROUTINE FORCE(xx,ix,ih,ipptr,X,F,ENER,VIR)
# else
      SUBROUTINE FORCE(xx,ix,ih,ipairs,X,F,ENER,VIR)
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
c     bugfixes
c       31
c
#ifdef DPREC 
      implicit double precision (a-h,o-z)
      parameter  (nwdbyt = 8)
#else
      parameter  (nwdbyt = 4)
#endif
#ifdef MEM_ALLOC
      pointer(ipptr, ipairs)
      integer   ipairs(1)
#else
      integer   ipairs(*)
#endif
      dimension xx(*)
      integer   ix(*), ih(*)

#include "extra.h"
#ifdef MPI
# include "parallel.h"
# include "ew_parallel.h"
# include "mpif.h"
#endif
      parameter (NT_LIM=5000)
#ifdef ISTAR2
      integer*2 itrp
#endif
      LOGICAL BELLY,NOCRST,FSTCLL
      DATA FSTCLL/.TRUE./
      SAVE FSTCLL
c     
c     EWALD: bypass the standard AMBER RESNBA and NONBON and 
c     call special pairlist and nonbond force and energy
c     evaluation routines...
c
# ifdef DPREC
      double precision atvir(6), molvir(6), subvir(6)
# else
      real atvir(6), molvir(6), subvir(6)
# endif
c
#include "md.h"
#include "box.h"
#include "iewald.h"
#include "nmr.h"
#include "memory.h"
#include "parms.h"
#include "pol.h"
#include "files.h"
#include "lgcom.h"
#include "avfrc.h"
c
      DIMENSION ENMR(3),DEVDIS(4),DEVANG(4),DEVTOR(4),VT(4)
C
C Save variables for secondary cutoff:
C
      SAVE VIR2ND, ENB2ND, EEL2ND, EHB2ND
c
      dimension itrp(3,NT_LIM)
      DIMENSION X(*),F(*),ENE(30),ENER(*),VIR(*),VIR2ND(4)
      save itrp
C
      aveper = 0.0d0
      aveind = 0.0d0
      avetot = 0.0d0
      ptot = 0.0d0
      epolar = 0.0d0
      emcur = 0.0d0
      do 897 jn = 1,30
        ene(jn) = 0.0d0
  897 continue
c
      DUM = 0.0D0
      ZERO = 0.0D0
      NREP = 15
      NREPC = 5
      BELLY = IBELLY.GT.0
      NOCRST = .false.
      KBOND = NBONA-MBONA
      NTTYP = NTYPES*(NTYPES+1)/2
#ifdef MPI
      if (mpi_orig) then
c     =========================== AMBER/MPI ===========================
c
c     Check to see if we are done yet in mpi_orig case (tec3).
c     This is done by monitoring the status of an integer notdone.
c     If notdone .eq. 1 then we keep going.  notdone is set to zero
c     when we no longer want to call force().  This perhaps is not the
c     most efficient means to implement this check...
c
         call mpi_bcast(notdone,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         if (notdone .ne. 1) return
c
c     JV Send xyz coords to all nodes
c
        call second(time0)
        call sendxyz(xx,vir)
        CALL SECOND(TIME1)
        TIMSTS(21) = TIMSTS(21) + (TIME1-TIME0)
c
c     ========================= END AMBER/MPI =========================
      end if
#endif /* MPI */
C
C     ----- SET THE BOX RELATED VARIABLES -----
C
      CALL MDBOX
C
C     ----- ZERO OUT THE ENERGIES AND FORCES -----
C
      call zerout(3,enmr)
      call zerout(4,vir)
      call zerout(4,vt)
      CALL ZEROUT(NREP+NREPC,ENE)
      CALL ZEROUT(3*NATOM,F)
      CALL ZEROUT(3*NATOM,FAF)
      DUM = 0.0D+00
C
C     ----- PARTITION THE WORK ARRAY FOR NNBOND -----
C
c     L95  is start of 3*natom real scratch
c     L100 is start of 4*natom real scratch
c     m06  is start of 4*natom Integer scratch
c
c     Core ptr     Size   Type   NNBOND          RESNBA
c      L95       3*Natom  Real  XWIJ(3,*)       XWIJ(3,*),RRW(*)
c      LW10        Natom  Real  RW(*)           
c      LW11        Natom  Real  RW2(*)
c      LW12        Natom  Real  RW3(*)
c      LW15      3*Natom  Real  FW(3,*)
c      m06         Natom  Int   JPW(*)          IWH(*)
c      m08         Natom  Int                   IWA(*)
c      m10         Natom  Int                   IEXW(*)
c      iw10        Natom  Int                   IRP(*)
c      
      LRRW = L95 + 3*NMXRS*NMXRS
      LW10 = L100
      LW11 = LW10 + NATOM
      LW12 = LW11 + NATOM
      LW15 = LW12 + NATOM
      iw10 = m10  + natom
C
C ----------------------------------------------------------------
C Do weight changes, if requested
C ----------------------------------------------------------------
C
      IF (NMROPT.GT.0)
     *  CALL NMRCAL(X,F,IH(M04),IH(M02),IX(I02),XX(L20),ENMR,DEVDIS,
     *              DEVANG,DEVTOR,TEMP0,TAUTP,CUT,NTB,XX(LNMR01),
     *              IX(INMR02),XX(L95),31,6,RK,TK,PK,CN1,
     *              CN2,ASOL,BSOL,XX(L15),NUMBND,NUMANG,NPTRA-NIMPRP,
     *              NIMPRP,NTTYP,NPHB,NATOM,NATOM,NTYPES,NRES,RAD,
     *              WEL,RADHB,WELHB,RWELL,ISFTRP,-1,'WEIT')
C
C ----------------------------------------------------------------
C Re-calculate the pairlist whenever NTNB.NE.0:
C ----------------------------------------------------------------
C
      IF (NTNB.NE.0) THEN
         NPAIR = 0
         NHB = 0
C
C If CUT2ND.GT.0, then we are carrying out dual cutoffs. The energy and forces
C for interactions between CUT and CUT2ND are updated only when the
C pairlist is updated. The procedure used is: 1) Call RESNBA to determine
C the pairlist for the secondary cutoff range; 2) call NNBOND to determine
C the energy (ENB2ND,EEL2ND,EHB2ND) and forces (F2ND == XX(LVM01)) 
C due to these pairs; 3) Call RESNBA again to determine the 
C primary pairlist. -- dap (12/92)
C
C Pairlist and energy for secondary region when CUT2ND>0:
C
         IF (CUT2ND.GT.0.01D0) THEN
            CALL SECOND(TIME0)
            CALL RESNBA(MAXNB,   NATOM,   NRES,    NPAIR,    NHB,
     +                  IX(I02), IX(I66), IX(I04), IX(I06),  IX(I08),
     +                  IX(I10), IX(I78), 
#ifdef MEM_ALLOC
     +                  ipptr,  X, 
#else
     +                  ipairs,  X, 
#endif
     +                  IX(I62), IH(M06), IH(M08), ADUM,     IH(IW10), 
     +                  NTYPES,  XX(L95), XX(LRRW),IH(M10),  CUT2ND, 
     +                  NWDVAR
#ifndef MPI
     + )
#else
     +                  ,IX(IVM03),nsolpr,ix(i68))
c
c     JV Balance number of pairs across all nodes
c           call balance(IX(I78),ipairs,IX(IVM03),npair,nsolw,nsolpr)
#endif
            CALL SECOND(TIME1)
            TIMSTS(1) = TIMSTS(1) + (TIME1-TIME0)
C
            CALL ZEROUT(3*NATOM,XX(LVM01))
            CALL SECOND(TIME0)
            CALL NONBON(NATOM,   IX(I78),  ipairs , IX(I04),  IX(I06),
     +                  X,       XX(LVM01),cn1,     cn2,      asol,
     +                  bsol,    XX(L15),  ENE(2),  ENE(4),   ENE(3),
     +                  XX(L95), XX(LW10), XX(L45), XX(LW15), IH(M06),
     +                  VIR2ND,  NTYPES,   IX(I01), NSOLW,    XX(LW11),
     +                  XX(LW12),IX(I02)
#ifndef MPI
     +                  )
#else
c     =========================== AMBER/MPI ===========================
c
c     JV add force, ene, vir, npair, nhb copies from all nodes
c
     +                 ,IX(IVM03))
c
       call fdist(XX(LVM01),XX(LW10),ene,vir2nd,npair,nhb)
c
c     ========================= END AMBER/MPI =========================
#endif
            CALL SECOND(TIME1)
            TIMSTS(2) = TIMSTS(2) + (TIME1-TIME0)
            ENB2ND = ENE(2)
            EEL2ND = ENE(3)
            EHB2ND = ENE(4)
            ESECND = ENE(2) + ENE(3) + ENE(4)
c            if (master) WRITE(6,449) NPAIR,NHB,ESECND
         END IF
C
C Regular pairlist:
C
C
C If iprr.ne.0, then pairs will be read from a pairlist file. Allows
C one to "freeze" the pairlist, if desired. -- dac (4/90)
C
C Set IPRR=0 after reading the pairlist the first time -- dap (6/93)
C
         IF (IPRR.NE.0) THEN
#ifdef MEM_ALLOC
            CALL PAIRR(NPAIR,NHB,IX(I78),ipptr)
#else
            CALL PAIRR(NPAIR,NHB,IX(I78),ipairs)
#endif
            IPRR = 0
         ELSE IF (iewald .EQ. 1) THEN
            CALL SECOND(TIME0)
               call ewald_list(x,ix(i04),ix(i06),ix(i08),ix(i10),
     .           ntypes,natom,xx,ix,
#ifdef MEM_ALLOC
     +           ipptr,xx(l15))
#else
     +           ipairs,xx(l15))
#endif
         ELSE
            CALL SECOND(TIME0)
            CALL RESNBA(MAXNB,   NATOM,   NRES,    NPAIR,    NHB,
     +                  IX(I02), IX(I66), IX(I04), IX(I06),  IX(I08),
     +                  IX(I10), IX(I78), 
#ifdef MEM_ALLOC
     +                  ipptr,  X,
#else
     +                  ipairs,  X,
#endif
     +                  IX(I62), IH(M06), IH(M08), ADUM,     IH(IW10),
     +                  NTYPES,  XX(L95), XX(LRRW),IH(M10),  ZERO,
     +                  NWDVAR
#ifndef MPI
     + )
#else
     +                  ,IX(IVM03),nsolpr,ix(i68))
c
            CALL SECOND(TIME1)
            TIMSTS(1) = TIMSTS(1)+(TIME1-TIME0)
            time0 = time1
c
c           JV Balance number of pairs across all nodes
c
            if ( ilbnob .eq. 0 ) then
              call balance(IX(I78),ipairs,IX(IVM03),npair,nsolw,nsolpr)
            else
               if (nsolw .gt. 0) then
                  doqik = .true.
               else
                  doqik = .false.
               endif
            endif
c
            CALL SECOND(TIME1)
            TIMSTS(21) = TIMSTS(21) + (TIME1-TIME0)
#endif
#ifndef MPI
c            WRITE(6,478) NPAIR,NHB
            CALL SECOND(TIME1)
            TIMSTS(1) = TIMSTS(1)+(TIME1-TIME0)
#endif
         END IF
         IF (IPRW.NE.0) CALL PAIRW(NATOM,NPAIR,NHB,IX(I78),ipairs)
C
C Call to recalculate three body list (only for polarization)
C
         if (n3b.gt.0) then
            CALL SECOND(TIME0)
            call tripl(natom,ix(i78),ipairs,ix(j01),itrp,nt_lim)
            CALL SECOND(TIME1)
            TIMSTS(10) = TIMSTS(10)+(TIME1-TIME0)
         end if
C
C For ligand grid calculations, need to determine additional pairs
C for interactions to ligand atoms
C
         IF (LGRON.eq.1) THEN
            IF (FSTCLL .OR. (LGCLLS.GT.LGSKP-1 .OR. LGSKP.LE.0)) THEN
C           TINE LGNBAD(IAR1      ,IPAIRS    ,JPW       ,NPHIH     ,
C    *                  IPH       ,KPH       ,LPH       ,NPHIA     ,
C    *                  IPA       ,KPA       ,LPA)
            CALL LGNBAD(IX(I78)   ,IPAIRS    ,IH(M06)   ,NPHIH     ,
     *                  IX(I40)   ,IX(I44)   ,IX(I46)   ,NPHIA     ,
     *                  IX(I50)   ,IX(I54)   ,IX(I56)   )      
            FSTCLL = .FALSE.
            END IF
         END IF
      END IF
C
      emtot = 0.0d0
C
C ----------------------------------------------------------------
C Calls required only if polarization calculation being performed:
C ----------------------------------------------------------------
C
      if(ipol.gt.0) then
          call politr(natom,x,f,xx(l15),xx(l25),ix(i78),ipairs,
     $                xx(l05),xx(l10),xx(l65),epolar,xx(l95),xx(l100),
     $                ix(i02),aveper,aveind,avetot,emtot,instep,nres,
     $                ih(m12),ih(m14),ih(m16),ix(i01),
     $                ih(m06),xx(lw15))
          ene(21) = epolar
          call second(time0)
c
c     subroutine polder(natom,x,f,p,    q,      iara,   iarb,
c         xrc,    xij,    r2,     fw,      vt, n14,    ni14,
c         iarx,   mark,   jpw)
c
          call polder(natom,x,f,xx(l10),xx(l15),ix(i78),ipairs,
     $       xx(l45),xx(l85),xx(l95),xx(l100),vt,ih(m12),ih(m14),
     $       ih(m16),ix(i01),ih(m06))
c
          emcur = emtot
          call second(time1)
          timsts(7) = timsts(7) + time1 - time0
      else
c
c ----dipole momement stuff for non pol systems
c
         CALL DIPLE(X,XX(L15),IX(I02),NRES,PTOT,EMSQ,EMTOT,EMCUR,AVEPER,
     *              AVETOT)
      end if
C
C ----------------------------------------------------------------
C Calculate the non-bonded contributions 
C ----------------------------------------------------------------
C
      IF(NTF.LT.8) THEN
           CALL SECOND(TIME0)
C
           if (iewald .eq. 1) then
#ifndef NEW_NB_VIRIAL
              call ewald_force(X,natom,ix(i04),ix(i06),ntypes,
     $             XX(l15),cn1,cn2,asol,bsol,
     $             eelt,evdw,ehb,F,XX,IX,ipairs,
     $             atvir,molvir,subvir,XX(l45),virvsene)
#else
              call ewald_force(X,natom,ix(i04),ix(i06),ntypes,
     $             XX(l15),cn1,cn2,asol,bsol,
     $             eelt,F,XX,IX,ipairs,
     $             XX(l45),virvsene)
#endif
              ene(2) = evdw
              ene(3) = eelt
              ene(4) = ehb
#ifdef MPI
            if(mytaskid.eq.0)then
#endif
              ene(2) = evdw
              ene(3) = eelt
              ene(4) = ehb
c for now use molecular virial. For atomic need to put in
c corrections for intramolecular forces  i.e. add fi.ri for intramolecular
c forces fi on atom i. The molecular virial is calculated using the
c coords wrt center of molecule in XX(L45).  
#ifndef NEW_NB_VIRIAL
              vir(1) = 0.5d0*molvir(1)
              vir(2) = 0.5d0*molvir(4)
              vir(3) = 0.5d0*molvir(6)
#endif
#ifdef MPI
            else
               ene(2) = 0.
               ene(3) = 0.
               ene(4) = 0.

#ifndef NEW_NB_VIRIAL
               vir(1) = 0.
               vir(2) = 0.
               vir(3) = 0.
#endif
            endif
#endif
c MEA CULPA kluge to slip in runtime error estimate
#ifndef NEW_NB_VIRIAL
              ener(21) = virvsene
#endif
           else
c
           CALL AFRCMK(F,0,1,NATOM)
           CALL Nonbon(NATOM,   IX(I78), ipairs,  IX(I04),  IX(I06),
     +                 X,       F,       cn1,     cn2,      asol,
     +                 bsol,    XX(L15), ENE(2),  ENE(4),   ENE(3),
     +                 XX(L95), XX(LW10),XX(L45), XX(LW15), IH(M06),
     +                 VIR,     NTYPES,  IX(I01), NSOLW,    XX(LW11),
     +                 XX(LW12),IX(I02) 
#ifndef MPI
     +                 )
#else
     +                 ,IX(IVM03))
#endif
           endif
c
           CALL SECOND(TIME1)
           TIMSTS(2) = TIMSTS(2)+(TIME1-TIME0)
           CALL AFRDIF(F,0,1,NATOM)
           CALL AFRADD(0,1,NATOM)
      END IF
C
C ----------------------------------------------------------------
C Calculate the other contributions
C ----------------------------------------------------------------
C
#ifndef NEW_NB_VIRIAL
      vir(1) = vir(1) + vt(1)
      vir(2) = vir(2) + vt(2)
      vir(3) = vir(3) + vt(3)
      vir(4) = vir(4) + vt(4)
#endif
      CALL SECOND(TIME0)

c
      e3bod = 0.0d0
      if(n3b.gt.0) then
           CALL SECOND(TIME0)
           call threeb(x,f,itrp,e3bod,nt_lim)
           CALL SECOND(TIME1)
             TIMSTS(10) = TIMSTS(10)+(TIME1-TIME0)
      endif
c
      GOTO (41,42,43,44,45,46,150,165),NTF
C
C     ----- BOND ENERGY CONTRIBUTION -----
C
      ebdev = 0.0
   41 IF(NBONH.GT.0) THEN
      CALL BOND(NBONH,IX(I12),IX(I14),IX(I16),X,F,ENE(6),NOCRST)
      END IF
   42 IF(MBONA.GT.0) THEN
      CALL BOND(MBONA,IX(I18),IX(I20),IX(I22),X,F,ENE(7),NOCRST)
      END IF
      if (nbonh+mbona .gt. 0) ebdev = sqrt( ebdev/(nbonh+mbona) )
      CALL SECOND(TIME1)
      TIMSTS(3) = TIMSTS(3)+(TIME1-TIME0)
      TIME0 = TIME1
C
C     ----- CALCULATE THE CONSTRAINT BOND CONTRIBUTION SINCE
C           CONSTRAINED BONDS ARE NOT INCLUDED FOR SHAKE -----
C
   43 CONTINUE
      IF(KBOND.GT.0) THEN
      CALL CBOND(KBOND,IX(I18+MBONA),IX(I20+MBONA),
     +           IX(I22+MBONA),X,F,ENE(17),NOCRST)
      END IF
C
C     ----- ANGLE ENERGY CONTRIBUTION -----
C
      eadev = 0.0
      IF(NTHETH.GT.0) THEN
      CALL ANGL(NTHETH,IX(I24),IX(I26),IX(I28),IX(I30),X,F,
     +          ENE(8),NTHETH,DUM,NOCRST)
      END IF
   44 IF(NTHETA.GT.0) THEN
      CALL ANGL(NTHETA,IX(I32),IX(I34),IX(I36),IX(I38),X,F,
     +          ENE(9),MTHETA,ENE(18),NOCRST)
      END IF
      if (ntheth+ntheta .gt. 0)
     .  eadev = 57.296*sqrt( eadev/(ntheth+ntheta) )
      CALL SECOND(TIME1)
      TIMSTS(4) = TIMSTS(4)+(TIME1-TIME0)
      TIME0 = TIME1
C
C     ----- DIHEDRAL ENERGY CONTRIBUTION -----
C
   45 IF(NPHIH.GT.0) THEN
      CALL EPHI(NPHIH,IX(I40),IX(I42),IX(I44),IX(I46),IX(I48),
     +          XX(L15),IX(I04),X,F,
     +     ENE(10),ENE(11),ENE(12),NPHIH,DUM,NOCRST)
      END IF
   46 IF(NPHIA.GT.0) THEN
      CALL EPHI(NPHIA,IX(I50),IX(I52),IX(I54),IX(I56),IX(I58),
     +          XX(L15),IX(I04),X,F,
     +     ENE(13),ENE(14),ENE(15),MPHIA,ENE(19),NOCRST)
      END IF
      CALL SECOND(TIME1)
      TIMSTS(5) = TIMSTS(5)+(TIME1-TIME0)
  150 CONTINUE
C
C     ----- CALCULATE THE POSITION CONSTRAINT ENERGY -----
C
C Put the necessary calls in to add the postitional restraint forces
C into the averaged forces, if they're required. Note that this
C will probably not work properly with MPI; need to check into this. -- dap
C
      IF (NATC.GT.0 .OR. IFCAP.GT.0) CALL AFRCMK(F,0,1,NATOM)

      IF(NATC.GT.0) THEN 
      CALL XCONST(NATC,ENE(20),IX(I60),X,F,XX(L55),XX(L60))
      END IF
C
C     ----- CALCULATE THE CAP FORCE IF NECESSARY -----
C
      IF(IFCAP.GT.0) CALL CAPWAT(NATOM,X,F)

      IF (NATC.GT.0 .OR. IFCAP.GT.0) CALL AFRDIF(F,0,1,NATOM)
      IF (NATC.GT.0 .OR. IFCAP.GT.0) CALL AFRADD(0,1,NATOM)

#ifdef MPI
# ifdef NEW_NB_VIRIAL
            if(mytaskid.eq.0)then
              vir(1) = vir(1)+0.5d0*molvir(1)
              vir(2) = vir(2)+0.5d0*molvir(4)
              vir(3) = vir(3)+0.5d0*molvir(6)
            endif
              ener(21) = virvsene
      vir(1) = vir(1) + vt(1)
      vir(2) = vir(2) + vt(2)
      vir(3) = vir(3) + vt(3)
      vir(4) = vir(4) + vt(4)
# endif
c     =========================== AMBER/MPI ===========================
c
c     JV add force, ene, vir, npair, nhb copies from all nodes
c
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call second(time0)
      call fdist(F,XX(LW10),ene,vir,npair,nhb)
      CALL SECOND(TIME1)
      TIMSTS(22) = TIMSTS(22) + (TIME1-TIME0)
c
c     ========================= END AMBER/MPI =========================
#endif
c
C
C If doing averaged forces and NMRAFA > 0, then we add the NMR
C forces into the averaged forces. The forces before calling the NMR
C routines are stored in AFTMP and the differences after calling
C the nmr routines are added into FAF afterwards
C
      CALL AFRCMK(F,0,1,NATOM)
      eshf = 0.0
      epcshf = 0.0
      enoe = 0.0
      if (iredir(4) .ne. 0) call noecalc(x,f,ih(m04),xx,ix,ih)
      if (iredir(5) .ne. 0) call cshf(natom,x,f)
      if (iredir(7) .ne. 0) call pcshift(natom,x,f)
C
C Calculate the NMR restraint energy contributions, if requested.
C
      IF (NMROPT.GT.0) THEN
          CALL NMRCAL(X,F,IH(M04),IH(M02),IX(I02),XX(L20),ENMR,DEVDIS,
     *                DEVANG,DEVTOR,TEMP0,TAUTP,CUT,NTB,XX(LNMR01),
     *                IX(INMR02),XX(L95),31,6,RK,TK,PK,CN1,
     *                CN2,ASOL,BSOL,XX(L15),NUMBND,NUMANG,NPTRA-NIMPRP,
     *                NIMPRP,NTTYP,NPHB,NATOM,NATOM,NTYPES,NRES,RAD,
     *                WEL,RADHB,WELHB,RWELL,ISFTRP,-1,'CALC')
      END IF
      IF (NMRAFA.EQ.1) CALL AFRDIF(F,0,1,NATOM)
      IF (NMRAFA.EQ.1) CALL AFRADD(0,1,NATOM)
C
C ------------------------------------------------------------------------
C Calculate the ligand grid stuff, if requested
C
      LGCLLS = LGCLLS + 1
      IF (LGRON.EQ.1 .AND. (LGCLLS.GT.LGSKP .OR. LGSKP.LE.0)) THEN
C
C First call LGRFCL to determine the non-bonded interaction energies of
C the ligand atoms for reference (stored in RLGNRG(2,I)).
C
         IF (ABS(RLGREP).GT.1.0D-7) THEN
C           TINE LGRFCL(X         ,IAR1      ,IPAIRS    ,CG        ,
C    *       IAC       ,MARK      ,NATOM     ,TEMP0     ,FW        ,
C    *       XWIJ      ,RW        ,JPW       ,ICO       ,NTYPES    ,
C    *       CN1       ,CN2)
            CALL LGRFCL(X         ,IX(I78)   ,IPAIRS    ,XX(L15)   ,
     *       IX(I04)   ,IX(I01)   ,NATOM     ,TEMP0     ,XX(LW15)  ,
     *       XX(L95)   ,XX(LW10)  ,IH(M06)   ,IX(I06)   ,NTYPES    ,
     *       CN1       ,CN2)
         END IF

C       TINE LGCALL    (X         ,IAR1      ,IPAIRS    ,CG        ,
C    *       IAC       ,MARK      ,NATOM     ,TEMP0     ,FW        ,
C    *       XWIJ      ,RW        ,JPW       ,LGPAK     ,ISTART    ,
C    *       IEND)
         CALL LGCALL   (X         ,IX(I78)   ,IPAIRS    ,XX(L15)   ,
     *       IX(I04)   ,IX(I01)   ,NATOM     ,TEMP0     ,XX(LW15)  ,
     *       XX(L95)   ,XX(LW10)  ,IH(M06)   ,IH(M08)   ,0         ,
     *       0)
      END IF
C ------------------------------------------------------------------------
#ifdef COSMO
c
c     calculate COSMO solvation energy, result into epolar (for now)
c       current version requires Michel Sanners "msms" program
c       to be in the path.
c
      if (cosmod.gt.0.0)
     *     call cosmo(natom,x,xx(L15),ih(m04),cosmod,epolar,f)
#endif
C
C     ----- CALCULATE TOTAL ENERGY AND GROUP THE COMPONENTS -----
C
C
C If secondary cutoff used, add the secondary non-bonded energies
C into the total energy accumulators now. Also include the virial
C due to the secondary region.
C
      IF (CUT2ND.GT.0.01D0) THEN
         ENE(2) = ENE(2) + ENB2ND
         ENE(3) = ENE(3) + EEL2ND
         ENE(4) = ENE(4) + EHB2ND
         VIR(1) = VIR(1) + VIR2ND(1)
         VIR(2) = VIR(2) + VIR2ND(2)
         VIR(3) = VIR(3) + VIR2ND(3)
      END IF
c
      DO 152 M = 2,NREP
  152 ENE(1) = ENE(1)+ENE(M)
      ene(1) = ene(1) + epolar + e3bod
C
      ENE(5) = ENE(6)+ENE(7)+ENE(17)
      ENE(6) = ENE(8)+ENE(9)-ENE(18)
      ENE(7) = ENE(10)+ENE(13)-ENE(19)
      ENE(8) = ENE(11)+ENE(14)
      ENE(9) = ENE(12)+ENE(15)
c   dac change 1/90: add in nmr energies to constraints:
c     ENE(10) = ENE(17)+ENE(18)+ENE(19)+ENE(20)
      ENE(10) = ENE(17)+ENE(18)+ENE(19)+ENE(20)+eshf+enoe+
     *           epcshf+enmr(1)+enmr(2)+enmr(3)                
      ENE(1) = ENE(1)+ENE(10)-ENE(18)-ENE(19)
C
C     ----- TRANSFER THE ENERGIES TO THE ARRAY ENER -----
C
c  Vectorizing this loop in Unicos version 8.0.2 results in a crash.
c  On the other hand, de-vectorizing this loop in 8.0.3 results in a
c  crash, so the "scalar" directive is commented out.
c  Unicos 9.0 is fine either way.
c  - Kottalam, CRI 4/6/95 - Ross 5/3/95
c cdir$ next scalar
      DO 155 M = 1,10
        ENER(M) = ENE(M)
  155 CONTINUE
      ener(11) = epolar
      ener(12) = aveper
      ener(13) = aveind
      ener(14) = avetot
      if(instep.gt.0) then
         ener(15) = emcur/instep
      else
         ener(15) = 0.0d0
      endif
      ener(16) = e3bod
C
C If a secondary cutoff is being used (CUT2ND>0), add the secondary
C range forces into the primary forces here. XX(LVM01) is the start
C of the F2ND array. Also add these forces
C into the FAF forces used for averaged forces (routine STRAVF).
C
      IF (CUT2ND.GT.0.01D0) THEN
#ifdef MPI
        istart = iparpt(mytaskid) + 1
        iend = iparpt(mytaskid+1)
        istart3 = 3*istart -2
        iend3 = 3*iend
        DO 185 I = istart3,iend3
#else
        DO 185 I = 1,3*NATOM
#endif
          F(I) = F(I) + XX(LVM01+I-1)
  185   CONTINUE
        CALL STRAVF(XX(LVM01),NATOM)
      END IF
C
C     ----- IF BELLY IS ON THEN SET THE BELLY ATOM FORCES TO ZERO -----
C
      IF (BELLY) CALL BELLYF(NATOM,IX(I62),F)
  165 CONTINUE
      RETURN
  478 FORMAT(t2,'NB-update: NPAIRS =',I8,'  HBPAIR =',I8)
  479 FORMAT(t2,'Parallel : NPAIRS =',I8,'  HBPAIR =',I8)
  449 FORMAT(/,'2ndary cut: NPAIRS =',I8,'  HBPAIR =',I8,
     *         '  ENERGY=',F12.3)
      END
#endif /* SHARED_MEMORY */

#ifdef MEM_ALLOC
      subroutine pairr (npair,nhb,iar1,ipptr)
#else
      subroutine pairr (npair,nhb,iar1,ipairs)
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
#ifdef MEM_ALLOC
      pointer(ipptr, ipairs)
      dimension ipairs(1)
      integer crealloc
      external crealloc
#else
      dimension ipairs(*)
#endif
      dimension iar1(*)
c
#include "memory.h"
#include "sizes.h"
#include "iewald.h"
c
#ifdef PAIRFORM
      call amopen(56, 'PRLIST', 'O', 'F', 'R')
      read(56,3) latom,npair,nhb
#else
      call amopen(56, 'PRLIST', 'O', 'U', 'R')
      read(56) latom,npair,nhb
#endif
      if (latom.ne.natom) then
          write(6,*) ' Reading nonbon list header: atoms ',latom,
     +               ' pairs ',npair,' nhb ',nhb
          write(6,*) ' -- Atoms .ne. parm ',natom
          call mexit(6, 1)
      endif
      
      if (npair.gt.maxnb) then
#ifdef MEM_ALLOC
           maxnb = npair
           if ( iewald.eq.0 ) then
              idum = maxnb / nwdvar
           else
              idum = maxnb
           endif
c
c          -- if there is packing, an extra natoms space is needed
c
           if (ipack.eq.1) idum = idum + natom
c
c          -- convert to bytes
c
           ibsize = 6
           call i1mach(ibsize)
           nbytes = ibsize * idum
C
C          -- Reallocate the pairlist array
C
           ipptr = crealloc(ipptr, nbytes)
           write(6,'(a,i14,a)')
     +                 '| * force: pairlist limit has grown to',
     +                 maxnb,' pairs'
           if (ipptr.eq.0) then
             write(6,*) ' ** failed to allocate, bytes: ',nbytes
             call mexit(6,1)
           endif
#else
           write(6,*) ' Pairlist from file (',npair,
     +                ') greater than array size (',maxnb,
     +                ') redimension MAXPR in sizes.h and recompile'
           call mexit(6,1)
#endif
      endif
      write(6,1) npair
    1 format(' Reading ',i6,' pairs from file PRLIST')
#ifdef PAIRFORM
        read(56,3) (iar1(i),i=1,natom)
        read(56,3) (ipairs(i),i=1,npair)
    3 format (12(x,i5))
#else
        read(56) (iar1(i),i=1,natom)
        read(56) (ipairs(i),i=1,npair)
#endif
      close (56)
c
      return
      end
c------------------------------------------------------------------------
      subroutine pairw (natom,npair,nhb,iar1,ipairs)
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
      dimension iar1(*),ipairs(*)
c
#ifdef PAIRFORM
        call amopen(56, 'PRLIST', 'U', 'F', 'W')
        write(56,2) natom,npair,nhb
        write(56,2) (iar1(i),i=1,natom)
        write(56,2) (ipairs(i),i=1,npair)
        close (56)
    2 format (12(x,i5))
#else
        call amopen(56, 'PRLIST', 'U', 'U', 'W')
        write(56) natom,npair,nhb
        write(56) (iar1(i),i=1,natom)
        write(56) (ipairs(i),i=1,npair)
        close (56)
#endif
c
      write(6,1) npair
    1 format('0WRITING ',i6,' PAIRS TO FILE PRLIST (unit 56)')
c
      return
      end
C
      SUBROUTINE DIPLE(X,CG,IPRES,NRES,PTOT,EMSQ,EMTOT,EMCUR,AVEPER,
     *                 AVETOT)
C
C Routine calculates some elementary dipole quantities for systems
C where no polarization calculation being carried out... -jcc/dap
C
#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
      DIMENSION X(*),CG(*),IPRES(*)
C
#ifdef MPI
# include "parallel.h"
# include "mpif.h"
      dimension temp1(4), temp2(4)
#endif
      AVEPER = 0.0D0
      AVEIND = 0.0D0
      AVETOT = 0.0D0
      PPER = 0.0D0
      PTOT = 0.0D0
      EMX = 0.0D0
      EMY = 0.0D0
      EMZ = 0.0D0
      TERM1 = 4.8D0/18.2223D0
      ICNT = 1
C
cforcevector
#ifdef CRAY_MP
cmic$ do all
cmic$*       shared (CG, EMX, EMY, EMZ, ICNT, IPRES, NRES, PPER)
cmic$*       shared (X)
cmic$*       private (I, IULIM, J, LLIM, PX, PY, PZ)
#endif
      DO 10 J = 1,NRES
         PX = 0.0D0
         PY = 0.0D0
         PZ = 0.0D0
         LLIM =  IPRES(J)
         IULIM =  IPRES(J+1)-1
#ifdef CRAY_MP
         if (iulim-llim .gt. 3) then
#endif
         DO 20 I =LLIM,IULIM
            PX = PX + CG(I)*X(ICNT  )
            PY = PY + CG(I)*X(ICNT+1)
            PZ = PZ + CG(I)*X(ICNT+2)
            ICNT = ICNT + 3
   20    CONTINUE
#ifdef CRAY_MP
         else
cdir$ nextscalar
         DO 21 I =LLIM,IULIM
            PX = PX + CG(I)*X(ICNT  )
            PY = PY + CG(I)*X(ICNT+1)
            PZ = PZ + CG(I)*X(ICNT+2)
            ICNT = ICNT + 3
   21    CONTINUE
         endif
#endif
         PPER = PPER + SQRT(PX*PX + PY*PY + PZ*PZ)
         EMX = EMX + PX
         EMY = EMY + PY
         EMZ = EMZ + PZ
   10 continue
c
c
      PTOT = PPER
      EMSQ = ((EMX*TERM1)**2 + (EMY*TERM1)**2 + (EMZ*TERM1)**2)
      EMTOT = EMTOT + EMSQ/(NRES)
      EMCUR = EMTOT
      AVEPER = (PPER / NRES) * TERM1
      AVETOT = (PTOT / NRES) * TERM1
C
      RETURN
      END
