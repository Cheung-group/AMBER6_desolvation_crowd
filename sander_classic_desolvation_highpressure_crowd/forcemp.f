#ifdef SGI_MP
# define SHARED_MEMORY
#endif
#ifdef CRAY_MP
# define SHARED_MEMORY
#endif
#ifdef CNX_MP
# define SHARED_MEMORY
#endif
#include "vector.h"
c-----------------------------------------------------------------------
      SUBROUTINE FORCE(xx,ix,ih,ipairs,X,F,ENER,VIR)
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
c     Mods for Rev A:
c     Core arrays moved to call args
c     Hollerith mapping switched to /memlc/
c     scratch array partitioning moved out of if block
c     scratch array management completely revamped for better efficiency
c     resnba and nnbond args changed
c     pairlist overwrite check moved to resnba, where it will actually work.
c     #ifdef implicit doubles statement
c     NOCRST hardwired to .false.  this means we always image bonds, ang, dih.
c     Added support for residue based periodic imaging
C
C     Mods for NMR (version 4):
C     Added calls to NMRCAL
C     Added appropriate accumulations of nmr-restraint related quantities -- dap
c     Fixed bug in calculation of L15 partition. Might have caused check on
c       sufficient real storage in LOCMEM to fail to flag problem...-- dap
C
C     Mods for Version 4.1:
C       Add code to allow use of a secondary cutoff range -- dap (12/92)
C       Added and verified Tom Dardens EWALD code -- tec3 (10/94)
c     
c         EWALD: bypass RESNBA and NONBON for special pairlist and 
c         nonbond force and energy routines...  NOTE: this is somewhat 
c         tricky in this routine, forcemp.f, due to the need to prevent 
c         parallelism of the nonbon call during shared memory...
c
c     bugfixes
c       31
c
#ifdef DPREC 
      implicit double precision (a-h,o-z)
      parameter  (nwdbyt = 8)
#else
      parameter  (nwdbyt = 4)
#endif
#ifdef MPI
# include "parallel.h"
#endif
#include "extra.h"
      parameter (NT_LIM=5000)
#ifdef ISTAR2
      integer*2 itrp
#endif
      LOGICAL BELLY,NOCRST,FSTCLL
      DATA FSTCLL/.TRUE./
      SAVE FSTCLL
#ifdef DPREC
      double precision atvir(6), molvir(6), subvir(6)
#else
      real atvir(6), molvir(6), subvir(6)
#endif
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
C
      dimension xx(*)
      integer   ix(*), ih(*), ipairs(*)
c
      dimension itrp(3,NT_LIM)
      DIMENSION X(*),F(*),ENE(30),ENER(*),VIR(*),VIR2ND(4)
      save itrp

c     SGI_MP: Setup and allocation of scratch memory
c     The following is the space needed:
c       ene:     30
c       vir:     3
c       force:   3*natom
c       xwij:    9*natom (since qiktip needs this much)
c       rw:      1*natom
c       fw:      3*natom
c       rw2:     1*natom (qiktip)
c       rw3:     1*natom (qiktip)
c       jpw:     1*natom (integer)
c
#if defined SGI_MP || defined CRAY_MP
      pointer (ptr_f,p_f)
      dimension p_f(1)
      pointer (ptr_i,ip_t)
      dimension ip_t(1)
#endif
#ifdef CNX_MP
      integer*4 ptr_f, ptr_i
      allocatable(p_f(:), ip_t(:))
#endif
      common /SGI_COM/ptr_f,ptr_i,NumProc,NT_FSIZ,NT_ISIZ
c     total real is 30 + 3 + (3+9+1+3+1+1)*NATOM
      NB_FSIZ = 33 + 18*NATOM
c     total integer is 3*natom
      NB_ISIZ = 3*NATOM

#ifdef CNX_MP
      NCPU = NumProc
      allocate(p_f(NB_FSIZ*(NCPU-1)), ip_t(NB_ISIZ*(NCPU-1)))
#endif
#ifdef SGI_MP     
      CALL SGI_GETMEM( 0 , 0 )
99    NCPU = NumProc
      CALL SGI_GETMEM( NB_FSIZ*(NCPU-1) , NB_ISIZ*(NCPU-1) )
      IF (NCPU.NE.NumProc) GOTO 99     
#endif
c     SGI_MP: end setup and allocation of scratch memory
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
c      L95       3*Natom  Real  XWIJ(3,*)       XWIJ(3,*)
c      LW10        Natom  Real  RW(*)           RRW(*)
c      LW11        Natom  Real  RW2(*)
c      LW12        Natom  Real  RW3(*)
c      LW15      3*Natom  Real  FW(3,*)
c      m06         Natom  Int   JPW(*)          IWH(*)
c      m08         Natom  Int                   IWA(*)
c      m10         Natom  Int                   IEXW(*)
c      iw10        Natom  Int                   IRP(*)
c      
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
     +                    IX(I02), IX(I66), IX(I04), IX(I06),  IX(I08),
     +                    IX(I10), IX(I78), ipairs,  X, 
     +                    IX(I62), IH(M06), IH(M08), ADUM,     IH(IW10), 
     +                    NTYPES,  XX(L95), XX(LW10),IH(M10),  CUT2ND,   
     +                    NWDVAR)
              CALL SECOND(TIME1)
              TIMSTS(1) = TIMSTS(1) + (TIME1-TIME0)
C
              CALL ZEROUT(3*NATOM,XX(LVM01))
              CALL SECOND(TIME0)

#ifdef SGI_MP
C$doacross 
C$&  LOCAL(iproc,iepert,if_off,ii_off,istart,iend,npp,
C$&        DUM,time1,time0),
C$&  SHARE(p_f,ip_t,NumProc,NB_FSIZ,NB_ISIZ,
C$&        natom,ix,x,f,cn1,cn2,asol,bsol,xx,ih,vir2nd,ntypes,iptatm,
C$&        imgslt,ntf,ene,nocrst,ipairs,nsolw)
#endif
#ifdef CRAY_MP
        ncpu = numproc
        call hpalloc(ptr_f,ncpu*NB_FSIZ,ierr,1)
        call hpalloc(ptr_i,ncpu*NB_ISIZ,ierr,1)

cmic$ parallel autoscope
cmic$*private(iproc,iepert,if_off,ii_off,istart,iend,npp)
cmic$*private(DUM,time1,time0)
cmic$*shared(p_f,ip_t,NumProc,NB_FSIZ,NB_ISIZ)
cmic$*shared(natom,ix,x,f,cn1,cn2,asol,bsol,xx,ih,vir2nd,ntypes,iptatm)
cmic$*shared(imgslt,ntf,ene,nocrst)
cmic$*shared(ipairs)
cmic$*shared(nsolw)
cmic$ do parallel
#endif
#ifdef CNX_MP
C$DIR NO_PEEL
C$DIR LOOP_PARALLEL(IVAR=iproc)
C$DIR LOOP_PRIVATE(iepert,if_off,ii_off)
#endif

        do iproc = 0, numproc-1

           if_off = (IProc-1)*NB_FSIZ
           ii_off = (IProc-1)*NB_ISIZ
           if (iproc .ne. 0) then
              do iepert = 1, 33 + 3*NATOM
                 p_f(iepert+if_off) = 0.0
              enddo
           end if

           if (iproc .eq. 0) then
             CALL NONBON(NATOM,   IX(I78),  ipairs,  IX(I04),  IX(I06),
     +                   X,       XX(LVM01),cn1,     cn2,      asol,
     +                   bsol,    XX(L15),  ENE(2),  ENE(4),   ENE(3),
     +                   XX(L95), XX(LW10), XX(L45), XX(LW15), IH(M06),
     +                   VIR2ND,  NTYPES,
     +                   IX(I01), NSOLW,    XX(LW11),XX(LW12), IX(I02),
     .                   iproc,   numproc)
           else
             CALL NONBON(NATOM,   IX(I78),  ipairs,  IX(I04),  IX(I06),
     +                   X,       p_f(34+if_off),
     .                                      cn1,     cn2,      asol,
     +                   bsol,    XX(L15),  
     .                   p_f(2+if_off), p_f(4+if_off), p_f(3+if_off),
     +                   p_f(34+if_off+3*natom),
     .                   p_f(34+if_off+12*natom),
     .                                      XX(L45),
     .                   p_f(34+if_off+13*natom),
     .                   ip_t(1+ii_off),
     +                   p_f(31+if_off),    NTYPES,
     +                   IX(I01), NSOLW,    
     .                   p_f(34+if_off+16*natom),
     .                   p_f(34+if_off+17*natom),              IX(I02),
     .                   iproc,   numproc)
           endif
c       SGI_MP: end of multiprocessing loop
        enddo
#ifdef CRAY_MP
cmic$ enddo
cmic$ end parallel
#endif

        if (NumProc .gt. 1) then
           do IProc=1,NumProc-1
              if_off = (IProc-1)*NB_FSIZ
              ene(2)  = ene(2)  + p_f(2+if_off)
              ene(3)  = ene(3)  + p_f(3+if_off)
              ene(4)  = ene(4)  + p_f(4+if_off)
              vir2nd(1)  = vir2nd(1)  + p_f(31+if_off)
              vir2nd(2)  = vir2nd(2)  + p_f(32+if_off)
              vir2nd(3)  = vir2nd(3)  + p_f(33+if_off)
           enddo
#ifdef SGI_MP
c$doacross local(i,iproc,if_off),share(natom,f,p_f,NB_FSIZ)
           do i = 0 , (3*natom)-1
              do IProc=1,NumProc-1
#endif
#ifdef CRAY_MP
cmic$ parallel autoscope
cmic$*private(i,iproc,if_off)
cmic$*shared(natom,f,p_f,NB_FSIZ)
cmic$ do parallel
           do IProc=1,NumProc-1
              do i = 0 , (3*natom)-1
#endif
#ifdef CNX_MP
C$DIR PREFER_PARALLEL 
           do IProc=1,NumProc-1
              do i = 0 , (3*natom)-1
#endif
                 if_off = (IProc-1)*NB_FSIZ
                 XX(LVM01+i) = XX(LVM01+i) + p_f(34+i+if_off)
              enddo
           enddo
#ifdef CRAY_MP
cmic$ enddo
cmic$ end parallel
#endif
        endif

#ifdef MPI
c     =========================== AMBER/MPI ===========================
c
c     JV add force, ene, vir, npair, nhb copies from all nodes
c
       call fdist(XX(LVM01),XX(LW15),ene,vir2nd,npair,nhb)
c
c     ========================= END AMBER/MPI =========================
#endif
#ifdef CRAY_MP
      call hpdeallc(ptr_f,ierr,1)
      call hpdeallc(ptr_i,ierr,1)
#endif
              CALL SECOND(TIME1)
              TIMSTS(2) = TIMSTS(2) + (TIME1-TIME0)
              ENB2ND = ENE(2)
              EEL2ND = ENE(3)
              EHB2ND = ENE(4)
              ESECND = ENE(2) + ENE(3) + ENE(4)
c              if (master)WRITE(6,449) NPAIR,NHB,ESECND
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
            CALL PAIRR(NPAIR,NHB,IX(I78),ipairs)
            IPRR = 0
         ELSE IF (IEWALD .EQ. 1) THEN
            CALL SECOND(TIME0)
            call ewald_list(x,ix(i04),ix(i06),ix(i08),ix(i10),
     .           ntypes,natom,xx,ix,ipairs,xx(l15))
         ELSE
            CALL SECOND(TIME0)
            CALL RESNBA(MAXNB,   NATOM,   NRES,    NPAIR,    NHB,
     +                  IX(I02), IX(I66), IX(I04), IX(I06),  IX(I08),
     +                  IX(I10), IX(I78), ipairs,  X, 
     +                  IX(I62), IH(M06), IH(M08), ADUM,     IH(IW10), 
     +                  NTYPES,  XX(L95), XX(LW10),IH(M10),  ZERO,   
     +                  NWDVAR)
#  ifndef ARCHIVE
#  ifndef MPI
c            WRITE(6,478) NPAIR,NHB
#  endif
#  endif
            CALL SECOND(TIME1)
            TIMSTS(1) = TIMSTS(1)+(TIME1-TIME0)
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
      if (ipol.gt.0) then
          call politr(natom,x,f,xx(l15),xx(l25),ix(i78),ipairs,
     $                xx(l05),xx(l10),xx(l65),epolar,xx(l95),xx(l100),
     $                ix(i02),aveper,aveind,avetot,emtot,instep,nres
     $               ,ih(m12),ih(m14),ih(m16),ix(i01)
     $               ,ih(m06),xx(lw15))
          ene(21) = epolar
          call second(time0)
c
c     subroutine polder(natom,x,f,p,    q,      iara,   iarb,
c         xrc,    xij,    r2,     fw,      vt, n14,    ni14,
c         iarx,   mark,   jpw)
c
          call polder(natom,x,f,xx(l10),xx(l15),ix(i78),ipairs,
     $      xx(l45),xx(l85),xx(l95),xx(l100),vt,ih(m12),ih(m14),
     $      ih(m16),ix(i01),ih(m06))
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
c     EWALD: things get a little hairy here, at least
c     in terms of how the code flow goes.  With EWALD, we call special 
c     routines to do the nonbond energy and force evaluation.  This 
c     requires not calling NONBON below.  
c     
c     ...now we need to call the special EWALD routines when EWALD
c     is active (this must be done prior to entering the multiprocessor
c     loop below!)
c
      IF(NTF.LT.8 .and. iewald .eq. 1) THEN
         CALL SECOND(TIME0)
         call ewald_force(X,natom,ix(i04),ix(i06),ntypes,
     $        XX(l15),cn1,cn2,asol,bsol,
     $        eelt,evdw,ehb,F,XX,IX,ipairs,
     $        atvir,molvir,subvir,XX(l45),virvsene)
         ene(2) = evdw
         ene(3) = eelt
         ene(4) = ehb
c MEA CULPA kluge to slip in runtime error estimate
         ener(21) = virvsene
c for now use molecular virial. For atomic need to put in
c corrections for intramolecular forces  i.e. add fi.ri for intramolecular
c forces fi on atom i. The molecular virial is calculated using the
c coords wrt center of molecule in XX(L45).  
         vir(1) = 0.5d0*molvir(1)
         vir(2) = 0.5d0*molvir(4)
         vir(3) = 0.5d0*molvir(6)
         CALL SECOND(TIME1)
         TIMSTS(2) = TIMSTS(2)+(TIME1-TIME0)
      ENDIF
c
#ifdef SGI_MP
c     ========================= SHARED MEMORY =========================
c
c     SGI_MP: major parallel loop
c     if_off -- offset into allocated scratch float array for iproc
c     ii_off -- offset into allocated scratch integer array for iproc
c     istart -- start of chunk for iproc
c     iend   -- end of chunk for iproc
c     npp    -- actual number per processor
c     chunk  -- size of suggested work per iproc
c     time1, time0 for the timer
c     p_f    -- scratch real array (size NB_FSIZ)
c     p_i    -- scratch integer array (size NB_ISIZ)
c
C$doacross 
C$&  LOCAL(IProc,iepert,if_off,ii_off,istart,iend,npp,
C$&        DUM,chunk,time1,time0),
C$&  SHARE(p_f,ip_t,NumProc,NB_FSIZ,NB_ISIZ,
C$&        natom,ix,x,f,cn1,cn2,asol,bsol,xx,ih,vir,ntypes,iptatm,
C$&        idiel,imgslt,ntf,nbonh,mbona,nbona,ene,nocrst,kbond,
C$&        ntheth,ntheta,mtheta,nphih,nphia,mphia,iewald,ipairs,nsolw)
#endif
#ifdef CRAY_MP
      ncpu=numproc
      call hpalloc(ptr_f,ncpu*NB_FSIZ,ierr,1)
      call hpalloc(ptr_i,ncpu*NB_ISIZ,ierr,1)

cmic$ parallel autoscope
cmic$*private(IProc,iepert,if_off,ii_off,istart,iend,npp)
cmic$*private(DUM,chunk,time1,time0)
cmic$*shared(p_f,ip_t,NumProc,NB_FSIZ,NB_ISIZ)
cmic$*shared(natom,ix,x,f,cn1,cn2,asol,bsol,xx,ih,vir,ntypes,iptatm)
cmic$*shared(idiel,imgslt,ntf,nbonh,mbona,nbona,ene,nocrst,kbond)
cmic$*shared(ntheth,ntheta,mtheta,nphih,nphia,mphia,iewald)
cmic$*shared(ipairs)
cmic$*shared(nsolw)
cmic$ do parallel
#endif
#ifdef CNX_MP
C$DIR NO_PEEL
C$DIR LOOP_PARALLEL(IVAR=iproc)
C$DIR LOOP_PRIVATE(iepert,if_off,ii_off,istart,iend,npp,DUM,chunk)
#endif
      do 1000 IProc=0,NumProc-1

         if_off = (IProc-1)*NB_FSIZ
         ii_off = (IProc-1)*NB_ISIZ
         if (IProc .ne. 0) then
           do iepert = 1, 33 + 3*NATOM
              p_f(iepert+if_off) = 0.0
           enddo
         end if

      IF(NTF.LT.8 .and. iewald .eq. 0) THEN
c       SGI_MP: this is a hack throughout this routine such that only
c       one process calls the timer so this may not report actual time!
        if (iproc .eq. 0) then
           CALL SECOND(TIME0)
        endif
C
        if (iproc .eq. 0) then
           CALL AFRCMK(F,iproc,1,NATOM)
           CALL Nonbon(NATOM,   IX(I78), ipairs,  IX(I04),  IX(I06),
     +                 X,       F,       cn1,     cn2,      asol,
     +                 bsol,    XX(L15), ENE(2),  ENE(4),   ENE(3),
     +                 XX(L95), XX(LW10),XX(L45), XX(LW15), IH(M06),
     +                 VIR,     NTYPES,
     +                 IX(I01), NSOLW,   XX(LW11),XX(LW12), IX(I02),
     .                 iproc,   numproc)
           CALL AFRDIF(F,iproc,1,NATOM)
        else
           CALL AFRCMK(p_f(34+if_off),iproc,1,NATOM)
           CALL Nonbon(NATOM,   IX(I78), ipairs,  IX(I04),  IX(I06),
     +                 X,       p_f(34+if_off),
     .                                   cn1,     cn2,      asol,
     +                 bsol,    XX(L15), 
     .                 p_f(2+if_off), p_f(4+if_off), p_f(3+if_off),
     .                 p_f(34+if_off+3*natom),
     .                 p_f(34+if_off+12*natom),
     +                                   XX(L45),
     .                 p_f(34+if_off+13*natom),
     .                 ip_t(1+ii_off),
     +                 p_f(31+if_off),     NTYPES,
     +                 IX(I01), NSOLW,   
     .                 p_f(34+if_off+16*natom),  
     .                 p_f(34+if_off+17*natom),             IX(I02),
     .                 iproc,   numproc)
           CALL AFRDIF(p_f(34+if_off),iproc,1,NATOM)
        endif
        if (iproc .eq. 0) then
           CALL SECOND(TIME1)
           TIMSTS(2) = TIMSTS(2)+(TIME1-TIME0)
        endif
      END IF
C
C ----------------------------------------------------------------
C Calculate the other contributions
C ----------------------------------------------------------------
C
      if (iproc .eq. 0) then
      CALL SECOND(TIME0)
      endif
c
c     SGI_MP: moved 3 body contribution outside parallel loop for now
c
c     SGI_MP: delay NTF=8 until later
      GOTO (41,42,43,44,45,46,150,150),NTF
C
C     ----- BOND ENERGY CONTRIBUTION -----
C
   41 IF(NBONH.GT.0) THEN
c     SGI_MP: bonds with hydrogen
         chunk = (nbonh+0.0)/(numproc+0.0)
         istart = iproc*chunk
         iend = (iproc+1)*chunk
         if ((iproc+1).eq.numproc) iend = nbonh
         npp = iend - istart
         if (npp .gt. 0) then
            if (iproc .eq. 0) then
      CALL BOND(npp,IX(I12),IX(I14),IX(I16),X,F,ENE(6),NOCRST)
            else
      CALL BOND(npp,IX(I12+istart),IX(I14+istart),IX(I16+istart),X,
     .          p_f(34+if_off),p_f(6+if_off),NOCRST)
            endif
         endif
      END IF
   42 IF(MBONA.GT.0) THEN
c     SGI_MP: bonds without hydrogen + constraints
         chunk = (mbona+0.0)/(numproc+0.0)
         istart = iproc*chunk
         iend = (iproc+1)*chunk
         if ((iproc+1).eq.numproc) iend = mbona
         npp = iend - istart
         if (npp .gt. 0) then
            if (iproc .eq. 0) then
      CALL BOND(npp,IX(I18),IX(I20),IX(I22),X,F,ENE(7),NOCRST)
            else
      CALL BOND(npp,IX(I18+istart),IX(I20+istart),IX(I22+istart),X,
     .          p_f(34+if_off),p_f(7+if_off),NOCRST)
            endif
         endif
      END IF
      if (iproc .eq. 0) then
      CALL SECOND(TIME1)
      TIMSTS(3) = TIMSTS(3)+(TIME1-TIME0)
      TIME0 = TIME1
      endif
C
C     ----- CALCULATE THE CONSTRAINT BOND CONTRIBUTION SINCE
C           CONSTRAINED BONDS ARE NOT INCLUDED FOR SHAKE -----
C
   43 CONTINUE
      IF(KBOND.GT.0) THEN
c     SGI_MP: constraint bond contribution, not in parallel!
      if (iproc .eq. 0) then
      CALL CBOND(KBOND,IX(I18+MBONA),IX(I20+MBONA),
     +           IX(I22+MBONA),X,F,ENE(17),NOCRST)
      endif
      END IF
C
C     ----- ANGLE ENERGY CONTRIBUTION -----
C
      IF(NTHETH.GT.0) THEN
c     SGI_MP: angles with hydrogen
         chunk = (ntheth+0.0)/(numproc+0.0)
         istart = iproc*chunk
         iend = (iproc+1)*chunk
         if ((iproc+1).eq.numproc) iend = ntheth
         npp = iend - istart
         if (npp .gt. 0) then
            if (iproc .eq. 0) then
      CALL ANGL(npp,IX(I24),IX(I26),IX(I28),IX(I30),X,F,
     +          ENE(8),npp,DUM,NOCRST)
            else
      CALL ANGL(npp,IX(I24+istart),IX(I26+istart),IX(I28+istart),
     .          IX(I30+istart),X,p_f(34+if_off),
     +          p_f(8+if_off),npp,DUM,NOCRST)
            endif
         endif
      END IF
   44 IF(NTHETA.GT.0) THEN
c     SGI_MP: angles without hydrogen + constraints
         chunk = (ntheta+0.0)/(numproc+0.0)
         istart = iproc*chunk
         iend = (iproc+1)*chunk
         if ((iproc+1).eq.numproc) iend = ntheta
         npp = iend - istart
         if (npp .gt. 0) then
            if (iproc .eq. 0) then
      CALL ANGL(npp,IX(I32),IX(I34),IX(I36),IX(I38),X,F,
     +          ENE(9),max0(0,mtheta-istart),ENE(18),NOCRST)
            else
      CALL ANGL(npp,IX(I32+istart),IX(I34+istart),IX(I36+istart),
     .          IX(I38+istart),X,p_f(34+if_off), p_f(9+if_off),
     .          max0(0,mtheta-istart),p_f(18+if_off),NOCRST)
            endif
         endif
      END IF
      if (iproc .eq. 0) then
      CALL SECOND(TIME1)
      TIMSTS(4) = TIMSTS(4)+(TIME1-TIME0)
      TIME0 = TIME1
      endif
C
C     ----- DIHEDRAL ENERGY CONTRIBUTION -----
C
   45 IF(NPHIH.GT.0) THEN
c     SGI_MP: dihedrals with hydrogen
         chunk = (nphih+0.0)/(numproc+0.0)
         istart = iproc*chunk
         iend = (iproc+1)*chunk
         if ((iproc+1).eq.numproc) iend = nphih
         npp = iend - istart
         if (npp .gt. 0) then
            if (iproc .eq. 0) then
      CALL EPHI(npp,IX(I40),IX(I42),IX(I44),IX(I46),IX(I48),
     +          XX(L15),IX(I04),X,F,
     +     ENE(10),ENE(11),ENE(12),npp,DUM,NOCRST)
            else
      CALL EPHI(npp,IX(I40+istart),IX(I42+istart),IX(I44+istart),
     .          IX(I46+istart),IX(I48+istart),
     +          XX(L15),IX(I04),X,p_f(34+if_off),
     +     p_f(10+if_off),p_f(11+if_off),p_f(12+if_off),npp,DUM,NOCRST)
            endif
         endif
      END IF
   46 IF(NPHIA.GT.0) THEN
c     SGI_MP: dihedrals without hydrogen + constraints
         chunk = (nphia+0.0)/(numproc+0.0)
         istart = iproc*chunk
         iend = (iproc+1)*chunk
         if ((iproc+1).eq.numproc) iend = nphia
         npp = iend - istart
         if (npp .gt. 0) then
            if (iproc .eq. 0) then
      CALL EPHI(npp,IX(I50),IX(I52),IX(I54),IX(I56),IX(I58),
     +          XX(L15),IX(I04),X,F,
     +     ENE(13),ENE(14),ENE(15),max0(0,mphia-istart),ENE(19),NOCRST)
            else
      CALL EPHI(npp,IX(I50+istart),IX(I52+istart),IX(I54+istart),
     .          IX(I56+istart),IX(I58+istart),
     +          XX(L15),IX(I04),X,p_f(34+if_off),
     .       p_f(13+if_off),p_f(14+if_off),p_f(15+if_off),
     .       max0(0,mphia-istart),p_f(19+if_off),NOCRST)
            endif
         endif
      END IF
      if (iproc .eq. 0) then
      CALL SECOND(TIME1)
      TIMSTS(5) = TIMSTS(5)+(TIME1-TIME0)
      endif
  150 CONTINUE
c     SGI_MP: end of major parallel loop
 1000 continue
#ifdef CRAY_MP
cmic$ end do
cmic$ end parallel
#endif
c     SGI_MP: now add up all the partial results

C
C add the forces into the averaged force arrays, FAF
C
        do Iproc=0,NumProc-1
           CALL AFRADD(iproc,1,natom)
        end do

      if (NumProc .gt. 1) then
         do IProc=1,NumProc-1
            if_off = (IProc-1)*NB_FSIZ
            do i = 1, 19
               ene(i)  = ene(i)  + p_f(i+if_off)
            enddo
            vir(1)  = vir(1)  + p_f(31+if_off)
            vir(2)  = vir(2)  + p_f(32+if_off)
            vir(3)  = vir(3)  + p_f(33+if_off)
         enddo
#ifdef SGI_MP
C$doacross local(i,IProc,if_off),share(natom,f,p_f,NB_FSIZ)
#endif
#ifdef CRAY_MP
cmic$ parallel autoscope
cmic$*private(i,IProc,if_off)
cmic$*shared(natom,f,p_f,NB_FSIZ)
cmic$ do parallel
#endif
#ifdef CNX_MP
C$DIR PREFER_PARALLEL
#endif
         do i = 0 , (3*natom)-1
            do IProc=1,NumProc-1
               if_off = (IProc-1)*NB_FSIZ
               f(1+i) = f(1+i) + p_f(34+i+if_off)
            enddo
         enddo
#ifdef CRAY_MP
cmic$ end do
cmic$ end parallel
#endif
      end if
      vir(1) = vir(1) + vt(1)
      vir(2) = vir(2) + vt(2)
      vir(3) = vir(3) + vt(3)
      vir(4) = vir(4) + vt(4)
c
c     SGI_MP: this was moved from above outside the parallel loop.
c     This should be parallelized and put back
      e3bod = 0.0d0
      if(n3b.gt.0) then
           CALL SECOND(TIME0)
           call threeb(x,f,itrp,e3bod,nt_lim)
           CALL SECOND(TIME1)
             TIMSTS(10) = TIMSTS(10)+(TIME1-TIME0)
      endif
C     SGI_MP: fix up the fact that NTF.eq.8 was skipped
      if (NTF.eq.8) goto 165
C
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
c     =========================== AMBER/MPI ===========================
c
c     JV add force, ene, vir, npair, nhb copies from all nodes
c
      call second(time0)
      call fdist(F,XX(LW15),ene,vir,npair,nhb)
      CALL SECOND(TIME1)
      TIMSTS(22) = TIMSTS(22) + (TIME1-TIME0)
c
c     ========================= END AMBER/MPI =========================
#endif
C
C If doing averaged forces and NMRAFA = 1, then we add the NMR
C forces into the averaged forces. The forces before calling the NMR
C routines are stored in AFTMP and the differences after calling
C the nmr routines are added into FAF afterwards. 
C
      CALL AFRCMK(F,0,1,NATOM)
c
c   dac change 1/90: add support for nmr-based refinement:
c
      eshf = 0.0
      enoe = 0.0
      epcshf = 0.0
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
C
#ifdef SGI_MP
C$doacross 
C$&  LOCAL(IProc,iepert,if_off,ii_off,istart,iend,npp,
C$&        DUM,chunk,time1,time0),
C$&  SHARE(p_f,ip_t,NumProc,NB_FSIZ,NB_ISIZ,
C$&        X,IX,XX,IH,ipairs,natom,temp0 )
#endif
        do 4200 iproc = 0,numproc-1
           if_off = (IProc-1)*(8*NATOM)
           ii_off = (IProc-1)*(3*NATOM)
           chunk = lgpts/numproc
           istart = iproc*chunk+1
           iend = (iproc+1)*chunk
           if (iproc.eq.numproc-1) iend = lgpts

           if (iproc.eq.0) then
C       TINE LGCALL    (X         ,IAR1      ,IPAIRS    ,CG        ,
C    *       IAC       ,MARK      ,NATOM     ,TEMP0     ,FW        ,
C    *       XWIJ      ,RW        ,JPW       ,LGPAK     ,ISTART    ,
C    *       IEND)
         CALL LGCALL   (X         ,IX(I78)   ,IPAIRS    ,XX(L15)   ,
     *       IX(I04)   ,IX(I01)   ,NATOM     ,TEMP0     ,XX(LW15)  ,
     *       XX(L95)   ,XX(LW10)  ,IH(M06)   ,IH(M08)   ,ISTART    ,
     *       IEND)
           else
         CALL LGCALL   (X         ,IX(I78)   ,IPAIRS    ,XX(L15)   ,
     *       IX(I04)   ,IX(I01)   ,NATOM     ,TEMP0     ,
     *       p_f(1+if_off),
     *       p_f(1+if_off+3*natom),
     *       p_f(1+if_off+6*natom),
     *       ip_t(1+ii_off),
     *       ip_t(1+ii_off+natom),
     *                                                   ISTART    ,
     *       IEND)
           end if
 4200   continue
      END IF
C ------------------------------------------------------------------------
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
     *           enmr(1)+enmr(2)+enmr(3)                
      ENE(1) = ENE(1)+ENE(10)-ENE(18)-ENE(19)
C
C     ----- TRANSFER THE ENERGIES TO THE ARRAY ENER -----
C
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
         DO 185 I = 1,3*NATOM
            F(I) = F(I) + XX(LVM01+I-1)
  185    CONTINUE
         CALL STRAVF(XX(LVM01),NATOM)
      END IF
C
C     ----- IF BELLY IS ON THEN SET THE BELLY ATOM FORCES TO ZERO -----
C
      IF(BELLY) CALL BELLYF(NATOM,IX(I62),F)
  165 CONTINUE
#ifdef CRAY_MP
      call hpdeallc(ptr_f,ierr,1)
      call hpdeallc(ptr_i,ierr,1)
#endif
#ifdef CNX_MP
      deallocate (p_f, ip_t)
#endif
      RETURN
  478 FORMAT(t2,'NB-update: NPAIRS =',I8,'  HBPAIR =',I8)
  479 FORMAT(t2,'Parallel : NPAIRS =',I8,'  HBPAIR =',I8)
  449 FORMAT(/,'2ndary cut: NPAIRS =',I8,'  HBPAIR =',I8,
     *         '  ENERGY=',F12.3)
      END
#ifdef SGI_MP
C
C     SGI_MP: Routine to allocate memory into pointers stored
C     in common variables
C
      SUBROUTINE SGI_GETMEM(MEM_FT,MEM_IT)
#     ifndef SGI_MP_MAX_PROCESSORS
#     define SGI_MP_MAX_PROCESSORS 36
#     endif
      pointer (ptr_f,p_f)
      dimension p_f(1)
      pointer (ptr_i,ip_t)
      dimension ip_t(1)
      common /SGI_COM/ptr_f,ptr_i,NumProc,NT_FSIZ,NT_ISIZ
      SAVE SGI_NEW,MAX_TOT
      DATA SGI_NEW,MAX_TOT/0,0/
# ifdef DPREC 
      ibyte = 8
# else
      ibyte = 4
# endif
      IF (SGI_NEW.EQ.0) THEN
        NT_FSIZ = 0
        NT_ISIZ = 0
        NumProc = min(mp_numthreads(),SGI_MP_MAX_PROCESSORS)
        SGI_NEW = 1
      ENDIF

      MEM_TOT = MEM_FT*IBYTE + MEM_IT*4

      IF (MEM_TOT.EQ.0) GOTO 9999

      IF ( MEM_TOT .GT. MAX_TOT ) THEN
        IF (MAX_TOT.GT.0) CALL FREE(PTR_F)
        PTR_F = MALLOC(MEM_TOT)
        IF (PTR_F.NE.0) THEN
          NT_FSIZ = MEM_FT
          NT_ISIZ = MEM_IT
          MAX_TOT = MEM_TOT
        ELSE
          write (6, 9991) MEM_TOT
          write (6, 9992) MEM_FT, MEM_IT
          NumProc = NumProc - 1
          NT_FSIZ = 0
          NT_ISIZ = 0
          MAX_TOT = 0
          write (6, 9993) NumProc
          GOTO 9999
        ENDIF
      END IF

      PTR_I = PTR_F + IBYTE*MEM_FT

9999  CALL FLUSH(6)
9991  format (/,'SGI_GETMEM: Allocation rejected.  Trying to get ',i5)
9992  format ('bytes.  Floats: ', i10, 'Integers: ', i10)
9993  format ('Try to reduce the number of processors to ', i4)

      RETURN
      END
#endif

