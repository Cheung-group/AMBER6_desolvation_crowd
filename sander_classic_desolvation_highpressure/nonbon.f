#ifdef SGI_MP
# define SHARED_MEMORY
#endif
#ifdef CRAY_MP
# define SHARED_MEMORY
#endif
#ifdef CNX_MP
# define SHARED_MEMORY
#endif
#ifdef SHARED_MEMORY
c     ========================= SHARED MEMORY =========================
c
c     NOTE: the only changes for the shared memory implementation are
c     additional parameters for the subroutine NONBON, representing
c     the current processor number and total number of processors,
c     and also SGI_MP load balancing for the main loop over residues.
c
c     ======================= END SHARED MEMORY =======================
#endif
#include "vector.h"
#ifdef FPS264
c fps264 uses cray cpp switch to trigger pairlist packing
# define CRAYFISH
#endif
      SUBROUTINE NONBON(NATOM,  IAR1,    ipairs,  IAC,      ICO,
     *                  X,      F,       CN1,     CN2,      ASOL,
     *                  BSOL,   CG,      ENB,     EHB,      EEL,
     *                  XWIJ,   RW,      XRC,     FW,       JPW,
     *                  VIR,    NTYPES,  
     *                  MARK,   NSOLW,   RW2,     RW3,      IPRES
#ifndef SHARED_MEMORY
# ifndef MPI
     * )
# else
     * ,natomnum)
# endif
#else
     .                 ,Iproc,  NumProc)
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
c
c     --- Evaluates energy and gradient due to nonbonded interactions ---
c
c     Authors:  This routine is based on the 3.0 nonbon by U.C. Singh.
c     Energy routines optimized by Rad Olson and Bill Swope, modified
c     by George Seibel.   Vectorized residue-based imaging by George Seibel.
c     Current revision of this routine is 3.0 Rev A.
c
c     (8/90):
c     For Version 4.0 (NMR), code to allow "soft-repulsion" non-bonded
c     terms added by David Pearlman.
C
C     (7/91)
C     Modified soft-repulsion h-bond loops to use r*(mixed) from the
C     h-bond ASOL & BSOL coefficients, not from the 6-12 r* values. -- D. Pearlman
C
C     (3/93)
C     Added code to handle the special case of TIP3P-TIP3P water
C     interactions. This code is much faster than the standard loops.
c
c     Machine dependencies are handled by cpp directives.  This includes
c     selection of single/double precision, Cray/Convex compiler directives,
c     integer*2 pairlist, and Cray/FPS264/generic packing.
c
c     The following cpp defines are used, and may be mixed where appropriate:
c     /lib/cpp -DCRAYFISH to make packed cray source with Cray directives.
c     /lib/cpp -DFPS264 for packed fps264 source
c     /lib/cpp -DDPREC for implicit double precision.
c     /lib/cpp -DISTAR2 for Integer*2 pairlist.
c     /lib/cpp -DFPS500 for fps500 compiler directives
c     Default (no defines) is unpacked with Convex/Stellar compiler
c     directives, Ansi integer pairlist, single precision.
c
c     11-Jun-89 Changed shape of iar1() to iar1(2,*) to reduce paging.
c
#ifdef DPREC
      implicit double precision (a-h,o-z)
#endif
#ifdef MPI
#include "parallel.h"
#include "mpif.h"
       dimension natomnum(*)
#endif
#ifdef ISTAR2 
      integer*2 ipairs 
#endif

#ifdef FPS264
      real ipairs
#endif
c
      dimension ipairs(*)
c        ... pairlist
      logical dield
c        ... true if distance dependent dielectric
      LOGICAL VIRIAL
c        ... constant pressure if true
      logical sltimg
#ifdef TFP
      logical ltmp1
#endif
#include "box.h"
#ifdef LES
#  include "les.h"
#endif
C Double well potential. margaret
#include "douwell.h"

C
#ifdef CRAY_MP
      DIMENSION IAR1(3,*),IAC(*),X(3,*),F(3,*),CG(*)
#else
      DIMENSION IAR1(2,*),IAC(*),X(3,*),F(3,*),CG(*)
#endif
      DIMENSION CN1(*),CN2(*),ASOL(*),BSOL(*)
      DIMENSION XWIJ(3,*),RW(*),ICO(*),XRC(3,*),VIR(3),RW2(*),RW3(*)
      dimension jpw(*)
      DIMENSION FW(3,*)
      DIMENSION VIRT(3),IPRES(*)
      integer mark(natom)

      integer ibctype
c
      data ZERO,HALF,ONE,FOUR/0.0d0,0.5d0,1.0d0,4.0d0/
      data SIX,TEN,TWELVE/6.0d0,10.0d0,12.0d0/
C
      dield  = (idiel .le. 0)
      sltimg = (imgslt .le. 0)
      VIRIAL = (IABS(NTB) .GE. 2)
      LPAIR  = 0
      lpack  = 1
      ENBT   = ZERO
      EELT   = ZERO
      EHBT   = ZERO
      VIRT(1) = ZERO
      VIRT(2) = ZERO
      VIRT(3) = ZERO

c     Used for bc tests. Note that cut is really the cutoff distance
c     squared. rcut=rcut*rcut in subroutine MDREAD. Adding on 4.0
c     just to be safe
      radius = sqrt(cut) + 4.0

C
C IAR1D1 is the number of dimensions for the IAR1 array.
C
#ifdef CRAY_MP
      IAR1D1 = 3
#else
      IAR1D1 = 2
#endif
C
C ---------------------------------------------------------------------------
C     *-- MAIN LOOP OVER THE ATOMS --*
C
C NSOLW is a pointer to the first atom of the solvent if 
C
C   1) the solvent is 3-point (TIP3P) water;
C   2) the solvent is packed contiguously from NSOLW to the end of atom list;
C   3) each water molecule is packed in the order O-H-H. 
C   4) JFASTW.NE.2 or 3 in input file.
C
C In this case, we can use much faster loops to calculate the water 
C interactions. Otherwise (NSOLW=0), use the standard loops for all 
C non-bonded interactions. NSOLW is set in the calling routine (FORCE).
C ---------------------------------------------------------------------------
C
#ifdef MPI
      DO 270 ITMP = 1,natomcnt
         I = natomnum(ITMP)
         if (i.gt.nsolw .and. doqik) go to 271
#else
      IF (NSOLW.EQ.0) THEN
        ILAST = NATOM - 1
      ELSE
        ILAST = NSOLW
      end if
C
      DO 270 I = 1,ILAST
#endif
#ifdef LES
        lestmp=nlesty*(lestyp(i)-1)
#endif
C
        LPR = IAR1(1,I)
        NPR = LPR+IAR1(2,I)
#ifdef CRAYFISH
        npack = 1
#endif
C
C       --- IF NPR.EQ.0 THEN SKIP TO NEXT ATOM ---
C
        IF(NPR.EQ.0) GO TO 260
C
#ifdef SHARED_MEMORY
c     ========================= SHARED MEMORY =========================
c
c       --- Load balancing for shared memory, BEGIN ---
c
        if ( mod(i, NumProc) .eq. Iproc ) then
c
c     ======================= END SHARED MEMORY =======================
#endif /* SGI_MP */
#ifdef CRAYFISH
C
C       --- unpack the pairlist for atom i ---
C
        npack = NPR/NUMPK+1
#  ifdef FPS264
        CALL VUSI32(ipairs(LPACK),1,JPW,1,npack)
#  else
#    ifdef CRAY_MP
        CALL UNPACK(ipairs(IAR1(3,i)),NBIT,JPW,NUMPK*npack)
#    else
        CALL UNPACK(ipairs(LPACK),NBIT,JPW,NUMPK*npack)
#    endif
#  endif
#endif
C
C
        IACI = NTYPES*(IAC(I)-1)
        if (dield) then
             CGI = -(CG(I)+CG(I))
        else
             CGI = -CG(I)
        endif
        IF (NTB.ne.0) then
C
C         --- PERIODIC IMAGING ---
C
          DO 10 JN = 1,NPR
#ifdef CRAYFISH
             j = jpw(jn)
#else
             J = ipairs(JN+LPAIR)
#endif
             XWIJ(1,JN) = X(1,I)-X(1,J)
             XWIJ(2,JN) = X(2,I)-X(2,J)
             XWIJ(3,JN) = X(3,I)-X(3,J)
c
c            --- fw is used as scratch for marker displacements ---
c                later it is used as scratch for forces.
c
             fw(1,jn) = X(1,mark(i))-X(1,mark(j))
             fw(2,jn) = X(2,mark(i))-X(2,mark(j))
             fw(3,jn) = X(3,mark(i))-X(3,mark(j))
   10     CONTINUE
c
c         --- iptatm is 0 unless iftres=0 (do all solute nonbonds), so
c             in that case only, no solute-solute imaging is ever done,
c             because solute might be .gt. box/2.  This loop is likely
c             to be slow and can be avoided if no solute-solvent 
c             imaging is done (imgslt=1) but then solute *must* stay in 
c             center of box ---
c
          if (i .le. iptatm) then
c
c             --- IPTATM.gt.0 => IFTRES=0 ---
c
c             --- i<iptatm => ATOM I IS SOLUTE ---
c
              if (sltimg) then
c
c                 --- IMGSLT=0 => imaging solute with solvent ---
c
#ifdef CRAYFISH
                  call bound1a(npr,fw,jpw,xwij)
#else
                  call bound1a(npr,fw,ipairs(lpair+1),xwij)
#endif
              else
c
c                 --- no imaging of solute => 
c                     solute better be centered ---
c

              endif
          else
c
c             --- ALL ATOMS IMAGED ---
c
             ibctype = 0
             if(abs(x(1,mark(i))-boxh(1)).gt.boxh(1)-radius) 
     *            ibctype=ibctype+4
             if(abs(x(2,mark(i))-boxh(2)).gt.boxh(2)-radius) 
     *            ibctype=ibctype+2
             if(abs(x(3,mark(i))-boxh(3)).gt.boxh(3)-radius) 
     *            ibctype=ibctype+1

             call bound2a(npr,fw,xwij,ibctype)
c              call bound2a(npr,fw,xwij,7)
          endif

c     Calculate inverse distances (constant dielectric) or inverse squared
c     distances (distance dependent dielectric). Note that calculation of
c     rw has been pulled out of bound1a() and bound2a().

c     In loops for the constant dielectric case, RW is used in place of
c     sqrt(rw). After the operations involving the inverse distance, RW is
c     immediately squared for use in later operations

          if(dield) then
             DO JN = 1,NPR
                RW(JN) = ONE/
     *               (XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2) 
             ENDDO
          else

#ifdef IBM
             do jn=1,npr
                rw(jn) = XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2
             enddo
             call vrsqrt(rw,rw,npr)
#else
             DO JN = 1,NPR
                RW(JN) = ONE/
     *               sqrt(XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2)
             ENDDO
#endif

          endif
          

C
C         --- CALCULATE THE ENERGY CONTRIBUTIONS FOR PERIODIC CASE
C             (NTB.NE.0)
C
C         --    "NORMAL" NONBONDS
C

          if(virial) then !! constant pressure simulation
             
             dumx = ZERO 
             dumy = ZERO
             dumz = ZERO

             if (dield) then

                if (isftrp.le.0) then

                   DO 70 JN = 1,LPR !! 6-12,DDIEL,VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      EELT = EELT+DF2
                      IC = ICO(IACI+IAC(J))
                      R6 = RW(JN)**3
#ifdef LES
                      F2 = CN2(IC)*R6*lfac
                      F1 = CN1(IC)*(R6*R6)*lfac
#else 
                      F2 = CN2(IC)*R6
                      F1 = CN1(IC)*(R6*R6)
#endif
                      ENBT = ENBT + (F2-F1)
                      DF = (DF2+SIX*((F2-F1)-F1))*RW(JN)
                      FW(1,JN) = XWIJ(1,JN)*DF
                      FW(2,JN) = XWIJ(2,JN)*DF
                      FW(3,JN) = XWIJ(3,JN)*DF
                      dumx = dumx + fw(1,jn)
                      dumy = dumy + fw(2,jn)
                      dumz = dumz + fw(3,jn)
                      F(1,J) = F(1,J)+FW(1,JN)
                      F(2,J) = F(2,J)+FW(2,JN)
                      F(3,J) = F(3,J)+FW(3,JN)
 70                CONTINUE
                else

                   IC = IAC(I)

                   DO 80 JN = 1,LPR !! SOFT,DDIEL,VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
                      JC = IAC(J)
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      RK = RWELL*lfac
#else 
                      RK = RWELL
#endif
                      RST = RAD(IC) + RAD(JC)
                      RR = ONE/RW(JN)
                      ARG = MAX(ZERO,(RST*RST-RR))
                      F1 = RK * ARG*ARG
                      ENBT = ENBT - F1
                      DF1 = -FOUR * RK * ARG

#ifdef LES
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      EELT = EELT+DF2
                      DF = DF1 + DF2*RW(JN)

                      FW(1,JN) = XWIJ(1,JN)*DF
                      FW(2,JN) = XWIJ(2,JN)*DF
                      FW(3,JN) = XWIJ(3,JN)*DF
                      dumx = dumx + fw(1,jn)
                      dumy = dumy + fw(2,jn)
                      dumz = dumz + fw(3,jn)
                      F(1,J) = F(1,J)+FW(1,JN)
                      F(2,J) = F(2,J)+FW(2,JN)
                      F(3,J) = F(3,J)+FW(3,JN)
 80                CONTINUE
                endif
             else

                if (isftrp.le.0) then

                   DO 90 JN = 1,LPR !! 6-12,CDIEL,VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      rw(jn) = rw(jn)*rw(jn)
                      
                      EELT = EELT+DF2
                      IC = ICO(IACI+IAC(J))
                      R6 = RW(JN)**3
#ifdef LES
                      F2 = CN2(IC)*R6*lfac
                      F1 = CN1(IC)*(R6*R6)*lfac
#else 
                      F2 = CN2(IC)*R6
                      F1 = CN1(IC)*(R6*R6)
#endif
                      ENBT = ENBT + (F2-F1)
                      DF = (DF2+SIX*((F2-F1)-F1))*RW(JN)
                      FW(1,JN) = XWIJ(1,JN)*DF
                      FW(2,JN) = XWIJ(2,JN)*DF
                      FW(3,JN) = XWIJ(3,JN)*DF
                      dumx = dumx + fw(1,jn)
                      dumy = dumy + fw(2,jn)
                      dumz = dumz + fw(3,jn)
                      F(1,J) = F(1,J)+FW(1,JN)
                      F(2,J) = F(2,J)+FW(2,JN)
                      F(3,J) = F(3,J)+FW(3,JN)
 90                CONTINUE
                else

                   IC = IAC(I)

                   DO 100 JN = 1,LPR !! SOFT,CDIEL,VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
                      JC = IAC(J)
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      RK = RWELL*lfac
#else 
                      RK = RWELL
#endif

#ifdef LES
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      rw(jn) = rw(jn)*rw(jn)
                      
                      
                      RST = RAD(IC) + RAD(JC)
                      RR = ONE/RW(JN)
                      ARG = MAX(ZERO,(RST*RST-RR))
                      F1 = RK * ARG*ARG
                      ENBT = ENBT - F1
                      DF1 = -FOUR * RK * ARG
                      
                      EELT = EELT+DF2
                      DF = DF1 + DF2*RW(JN)

                      FW(1,JN) = XWIJ(1,JN)*DF
                      FW(2,JN) = XWIJ(2,JN)*DF
                      FW(3,JN) = XWIJ(3,JN)*DF
                      dumx = dumx + fw(1,jn)
                      dumy = dumy + fw(2,jn)
                      dumz = dumz + fw(3,jn)
                      F(1,J) = F(1,J)+FW(1,JN)
                      F(2,J) = F(2,J)+FW(2,JN)
                      F(3,J) = F(3,J)+FW(3,JN)
 100               CONTINUE
                endif
             endif

             F(1,I) = F(1,I)-dumx
             F(2,I) = F(2,I)-dumy
             F(3,I) = F(3,I)-dumz

             dumx = ZERO
             dumy = ZERO
             dumz = ZERO
             if (dield) then

                if (isftrp.le.1) then

                   DO 110 JN = LPR+1,NPR !! 10-12 HB,DDIEL,VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      EELT = EELT+DF2
                      IC = -ICO(IACI+IAC(J))
                      R10 = RW(JN)**5
#ifdef LES
                      F1 = ASOL(IC)*R10*RW(JN)*lfac
                      F2 = BSOL(IC)*R10*lfac
#else 
                      F1 = ASOL(IC)*R10*RW(JN)
                      F2 = BSOL(IC)*R10
#endif
                      DF = (DF2-TWELVE*F1+TEN*F2)*RW(JN)
                      EHBT = EHBT+F1-F2
                      FW(1,JN) = XWIJ(1,JN)*DF
                      FW(2,JN) = XWIJ(2,JN)*DF
                      FW(3,JN) = XWIJ(3,JN)*DF
                      dumx = dumx + fw(1,jn)
                      dumy = dumy + fw(2,jn)
                      dumz = dumz + fw(3,jn)
                      F(1,J) = F(1,J)+FW(1,JN)
                      F(2,J) = F(2,J)+FW(2,JN)
                      F(3,J) = F(3,J)+FW(3,JN)
 110               CONTINUE
                else

                   DO 120 JN = LPR+1,NPR !! SOFT HB,DDIEL,VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      RK = RWELL*lfac
#else 
                      RK = RWELL
#endif
                      
                      IC = -ICO(IACI+IAC(J))
                      RST = RADHB(IC)
                      RR = ONE/RW(JN)
                      ARG = MAX(ZERO,(RST*RST - RR))
                      F1 = RK * ARG*ARG
                      EHBT = EHBT + F1
                      DF1 = -FOUR * RK * ARG

#ifdef LES
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      EELT = EELT+DF2

                      DF = DF1 + DF2*RW(JN)
                      FW(1,JN) = XWIJ(1,JN)*DF
                      FW(2,JN) = XWIJ(2,JN)*DF
                      FW(3,JN) = XWIJ(3,JN)*DF
                      dumx = dumx + fw(1,jn)
                      dumy = dumy + fw(2,jn)
                      dumz = dumz + fw(3,jn)
                      F(1,J) = F(1,J)+FW(1,JN)
                      F(2,J) = F(2,J)+FW(2,JN)
                      F(3,J) = F(3,J)+FW(3,JN)
 120               CONTINUE
                endif
             else

                if (isftrp.le.1) then

                   DO 130 JN = LPR+1,NPR !! 10-12 HB,CDIEL,VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      rw(jn) = rw(jn)*rw(jn)
                      
                      EELT = EELT+DF2
                      IC = -ICO(IACI+IAC(J))
                      R10 = RW(JN)**5
#ifdef LES
                      F1 = ASOL(IC)*R10*RW(JN)*lfac
                      F2 = BSOL(IC)*R10*lfac
#else 
                      F1 = ASOL(IC)*R10*RW(JN)
                      F2 = BSOL(IC)*R10
#endif
                      DF = (DF2-TWELVE*F1+TEN*F2)*RW(JN)
                      EHBT = EHBT+F1-F2
                      FW(1,JN) = XWIJ(1,JN)*DF
                      FW(2,JN) = XWIJ(2,JN)*DF
                      FW(3,JN) = XWIJ(3,JN)*DF
                      dumx = dumx + fw(1,jn)
                      dumy = dumy + fw(2,jn)
                      dumz = dumz + fw(3,jn)
                      F(1,J) = F(1,J)+FW(1,JN)
                      F(2,J) = F(2,J)+FW(2,JN)
                      F(3,J) = F(3,J)+FW(3,JN)
 130               CONTINUE
                else

                   DO 140 JN = LPR+1,NPR !! SOFT HB,CDIEL,VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      RK = RWELL*lfac
#else 
                      RK = RWELL
#endif
                      

#ifdef LES
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      rw(jn) = rw(jn)*rw(jn)
                      
                      IC = -ICO(IACI+IAC(J))
                      RST = RADHB(IC)
                      RR = ONE/RW(JN)
                      ARG = MAX(ZERO,(RST*RST - RR))
                      F1 = RK * ARG*ARG
                      EHBT = EHBT + F1
                      DF1 = -FOUR * RK * ARG
                      EELT = EELT+DF2

                      DF = DF1 + DF2*RW(JN)
                      FW(1,JN) = XWIJ(1,JN)*DF
                      FW(2,JN) = XWIJ(2,JN)*DF
                      FW(3,JN) = XWIJ(3,JN)*DF
                      dumx = dumx + fw(1,jn)
                      dumy = dumy + fw(2,jn)
                      dumz = dumz + fw(3,jn)
                      F(1,J) = F(1,J)+FW(1,JN)
                      F(2,J) = F(2,J)+FW(2,JN)
                      F(3,J) = F(3,J)+FW(3,JN)
 140               CONTINUE
                endif
             endif

             F(1,I) = F(1,I)-dumx
             F(2,I) = F(2,I)-dumy
             F(3,I) = F(3,I)-dumz
             
          else !! ----- NON-VIRIAL (NON-CONSTANT PRESSURE) CASE -----
             
             dumx = ZERO       
             dumy = ZERO
             dumz = ZERO

             if (dield) then

                if (isftrp.le.0) then

                   DO 71 JN = 1,LPR !! 6-12,DDIEL,NON-VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      EELT = EELT+DF2
                      IC = ICO(IACI+IAC(J))
                      R6 = RW(JN)**3
#ifdef LES
                      F2 = CN2(IC)*R6*lfac
                      F1 = CN1(IC)*(R6*R6)*lfac
#else 
                      F2 = CN2(IC)*R6
                      F1 = CN1(IC)*(R6*R6)
#endif
                      ENBT = ENBT + (F2-F1)
                      DF = (DF2+SIX*((F2-F1)-F1))*RW(JN)
                      FW1 = XWIJ(1,JN)*DF
                      FW2 = XWIJ(2,JN)*DF
                      FW3 = XWIJ(3,JN)*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
 71                CONTINUE
                else

                   IC = IAC(I)

                   DO 81 JN = 1,LPR !! SOFT,DDIEL,NON-VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
                      JC = IAC(J)
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      RK = RWELL*lfac
#else 
                      RK = RWELL
#endif
                      RST = RAD(IC) + RAD(JC)
                      RR = ONE/RW(JN)
                      ARG = MAX(ZERO,(RST*RST-RR))
                      F1 = RK * ARG*ARG
                      ENBT = ENBT - F1
                      DF1 = -FOUR * RK * ARG

#ifdef LES
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      EELT = EELT+DF2
                      DF = DF1 + DF2*RW(JN)

                      FW1 = XWIJ(1,JN)*DF
                      FW2 = XWIJ(2,JN)*DF
                      FW3 = XWIJ(3,JN)*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
 81                CONTINUE
                endif
             else

                if (isftrp.le.0) then

                   DO 91 JN = 1,LPR !! 6-12,CDIEL,NON-VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      rw(jn) = rw(jn)*rw(jn)
                      
                      EELT = EELT+DF2
                      IC = ICO(IACI+IAC(J))
                      R6 = RW(JN)**3
#ifdef LES
                      F2 = CN2(IC)*R6*lfac
                      F1 = CN1(IC)*(R6*R6)*lfac
#else 
                      F2 = CN2(IC)*R6
                      F1 = CN1(IC)*(R6*R6)
#endif
                      ENBT = ENBT + (F2-F1)
                      DF = (DF2+SIX*((F2-F1)-F1))*RW(JN)
                      FW1 = XWIJ(1,JN)*DF
                      FW2 = XWIJ(2,JN)*DF
                      FW3 = XWIJ(3,JN)*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
 91                CONTINUE
                else

                   IC = IAC(I)

                   DO 101 JN = 1,LPR !! SOFT,CDIEL,NON-VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
                      JC = IAC(J)
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      RK = RWELL*lfac
#else 
                      RK = RWELL
#endif


#ifdef LES
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      rw(jn) = rw(jn)*rw(jn)
                      
                      
                      RST = RAD(IC) + RAD(JC)
                      RR = ONE/RW(JN)
                      ARG = MAX(ZERO,(RST*RST-RR))
                      F1 = RK * ARG*ARG
                      ENBT = ENBT - F1
                      DF1 = -FOUR * RK * ARG
                      
                      EELT = EELT+DF2
                      DF = DF1 + DF2*RW(JN)

                      FW1 = XWIJ(1,JN)*DF
                      FW2 = XWIJ(2,JN)*DF
                      FW3 = XWIJ(3,JN)*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
 101               CONTINUE
                endif
             endif

             F(1,I) = F(1,I)-dumx
             F(2,I) = F(2,I)-dumy
             F(3,I) = F(3,I)-dumz

             dumx = ZERO
             dumy = ZERO
             dumz = ZERO
             if (dield) then

                if (isftrp.le.1) then

                   DO 111 JN = LPR+1,NPR !! 10-12 HB,DDIEL,NON-VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      EELT = EELT+DF2
                      IC = -ICO(IACI+IAC(J))
                      R10 = RW(JN)**5
#ifdef LES
                      F1 = ASOL(IC)*R10*RW(JN)*lfac
                      F2 = BSOL(IC)*R10*lfac
#else 
                      F1 = ASOL(IC)*R10*RW(JN)
                      F2 = BSOL(IC)*R10
#endif
                      DF = (DF2-TWELVE*F1+TEN*F2)*RW(JN)
                      EHBT = EHBT+F1-F2
                      FW1 = XWIJ(1,JN)*DF
                      FW2 = XWIJ(2,JN)*DF
                      FW3 = XWIJ(3,JN)*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
 111               CONTINUE
                else

                   DO 121 JN = LPR+1,NPR !! SOFT HB,DDIEL,NON-VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      RK = RWELL*lfac
#else 
                      RK = RWELL
#endif
                      
                      IC = -ICO(IACI+IAC(J))
                      RST = RADHB(IC)
                      RR = ONE/RW(JN)
                      ARG = MAX(ZERO,(RST*RST - RR))
                      F1 = RK * ARG*ARG
                      EHBT = EHBT + F1
                      DF1 = -FOUR * RK * ARG

#ifdef LES
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      EELT = EELT+DF2

                      DF = DF1 + DF2*RW(JN)
                      FW1 = XWIJ(1,JN)*DF
                      FW2 = XWIJ(2,JN)*DF
                      FW3 = XWIJ(3,JN)*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
 121               CONTINUE
                endif
             else

                if (isftrp.le.1) then

                   DO 131 JN = LPR+1,NPR !! 10-12 HB,CDIEL,NON-VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      rw(jn) = rw(jn)*rw(jn)
                      
                      EELT = EELT+DF2
                      IC = -ICO(IACI+IAC(J))
                      R10 = RW(JN)**5
#ifdef LES
                      F1 = ASOL(IC)*R10*RW(JN)*lfac
                      F2 = BSOL(IC)*R10*lfac
#else 
                      F1 = ASOL(IC)*R10*RW(JN)
                      F2 = BSOL(IC)*R10
#endif
                      DF = (DF2-TWELVE*F1+TEN*F2)*RW(JN)
                      EHBT = EHBT+F1-F2
                      FW1 = XWIJ(1,JN)*DF
                      FW2 = XWIJ(2,JN)*DF
                      FW3 = XWIJ(3,JN)*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
 131               CONTINUE
                else

                   DO 141 JN = LPR+1,NPR !! SOFT HB,CDIEL,NON-VIRIAL,PERIODIC
#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      RK = RWELL*lfac
#else 
                      RK = RWELL
#endif
                      

#ifdef LES
                      DF2 = CGI*CG(J)*RW(JN)*lfac
#else 
                      DF2 = CGI*CG(J)*RW(JN)
#endif
                      rw(jn) = rw(jn)*rw(jn)
                      
                      IC = -ICO(IACI+IAC(J))
                      RST = RADHB(IC)
                      RR = ONE/RW(JN)
                      ARG = MAX(ZERO,(RST*RST - RR))
                      F1 = RK * ARG*ARG
                      EHBT = EHBT + F1
                      DF1 = -FOUR * RK * ARG
                      EELT = EELT+DF2

                      DF = DF1 + DF2*RW(JN)
                      FW1 = XWIJ(1,JN)*DF
                      FW2 = XWIJ(2,JN)*DF
                      FW3 = XWIJ(3,JN)*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
 141               CONTINUE
                endif
             endif

             F(1,I) = F(1,I)-dumx
             F(2,I) = F(2,I)-dumy
             F(3,I) = F(3,I)-dumz
             
             
          endif

C
C         --- CALCULATE THE VIRIAL IF constant pressure simulation ---
C
          IF (VIRIAL) then
            DO 150 JN = 1,NPR
#ifdef CRAYFISH
                j = jpw(jn)
#else
                J = ipairs(JN+LPAIR)
#endif
                VIRT(1)=VIRT(1)+FW(1,JN)*(XWIJ(1,JN)-XRC(1,I)+XRC(1,J))
                VIRT(2)=VIRT(2)+FW(2,JN)*(XWIJ(2,JN)-XRC(2,I)+XRC(2,J))
                VIRT(3)=VIRT(3)+FW(3,JN)*(XWIJ(3,JN)-XRC(3,I)+XRC(3,J))
  150       continue
          endif
c         -- end of periodic case, ntb.ne.0 (used to be goto 250)

        else
C
C         --- CALCULATE THE ENERGY CONTRIBUTIONS FOR NO PERIODICITY ---
C
C         -- "NORMAL" NONBONDS
C
          dumx = ZERO
          dumy = ZERO
          dumz = ZERO
          if (dield) then
c
c             --- distance dependent dielectric
c
              if (isftrp.le.0) then
c                 --- 6-12 form, distance dependent dielectric, 
c                                non-periodic case ---
cforcevector
                  DO 170 JN = 1,LPR
#ifdef CRAYFISH
                      j = jpw(jn)
#else
                      J = ipairs(JN+LPAIR)
#endif
                      XW1 = X(1,I)-X(1,J)
                      XW2 = X(2,I)-X(2,J)
                      XW3 = X(3,I)-X(3,J)
                      R2INV = ONE/(XW1**2+XW2**2+XW3**2)
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      DF2 = CGI*CG(J)*R2INV*lfac
#else
                      DF2 = CGI*CG(J)*R2INV
#endif
                      EELT = EELT+DF2
                      IC = ICO(IACI+IAC(J))
                      R6 = R2INV**3
#ifdef LES
                      F2 = CN2(IC)*R6*lfac
                      F1 = CN1(IC)*(R6*R6)*lfac
#else
                      F2 = CN2(IC)*R6
                      F1 = CN1(IC)*(R6*R6)
#endif
                      ENBT = ENBT + (F2-F1)
                      DF = (DF2+SIX*((F2-F1)-F1))*R2INV
                      FW1 = XW1*DF
                      FW2 = XW2*DF
                      FW3 = XW3*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
  170             CONTINUE
              else
c                 --- soft-repulsion form, distance dep dielectric, 
c                                non-periodic case ---
                  IC = IAC(I)
cforcevector
                  DO 180 JN = 1,LPR
#ifdef CRAYFISH
                      j = jpw(jn)
#else
                      J = ipairs(JN+LPAIR)
#endif
                      XW1 = X(1,I)-X(1,J)
                      XW2 = X(2,I)-X(2,J)
                      XW3 = X(3,I)-X(3,J)
                      RR = (XW1**2+XW2**2+XW3**2)
                      R2INV = ONE/RR
c
                      JC = IAC(J)
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      RK = RWELL*lfac
#else
                      RK = RWELL
#endif
                      RST = RAD(IC) + RAD(JC)
                      ARG = MAX(ZERO,(RST*RST-RR))
                      F1 = RK * ARG*ARG
                      ENBT = ENBT - F1
                      DF1 = -FOUR * RK * ARG
c
#ifdef LES
                      DF2 = CGI*CG(J)*R2INV*lfac
#else
                      DF2 = CGI*CG(J)*R2INV
#endif
                      EELT = EELT+DF2
c
                      DF = DF1 + DF2*R2INV
                      FW1 = XW1*DF
                      FW2 = XW2*DF
                      FW3 = XW3*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
  180             CONTINUE
              endif
          else
c
c             --- constant dielectric
c
c             idiel=1 isftrp=0:default margaret

              if (isftrp.le.0) then
c                 --- 6-12 form, constant dielectric, 
c                                non-periodic case ---
cforcevector
                  DO 190 JN = 1,LPR
#ifdef CRAYFISH
                      j = jpw(jn)
#else
                      J = ipairs(JN+LPAIR)
#endif
                      XW1 = X(1,I)-X(1,J)
                      XW2 = X(2,I)-X(2,J)
                      XW3 = X(3,I)-X(3,J)
                      R2INV = ONE/(XW1**2+XW2**2+XW3**2)

c RD:margaret
		      RD = SQRT(XW1**2+XW2**2+XW3**2)
c end of margaret

#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      DF2 = CGI*CG(J)*SQRT(R2INV)*lfac
#else
                      DF2 = CGI*CG(J)*SQRT(R2INV)
#endif
                      EELT = EELT+DF2
                      IC = ICO(IACI+IAC(J))
                      R6 = R2INV**3
#ifdef LES
                      F2 = CN2(IC)*R6*lfac
                      F1 = CN1(IC)*(R6*R6)*lfac
#else
                      F2 = CN2(IC)*R6
                      F1 = CN1(IC)*(R6*R6)
#endif



c change force here. margaret. 8.15.00 .
c double well potential
c of the negative sign. 
c  This file is made to conjugate three potential for desolvation
c from r=0 to r0
c   U(r)=D0*Z(r)*(Z(r)-2) with Z(r)=pow(r0/r,k)
c from r0 to r1
c   U(r)=C*pow(Y(r),n)*(pow(Y(r),n)/2-pow(r1-r0,2n))/(2*n)+D1
c from r1 to infty
c   U(r)=-B(Y(r)-h1)/(pow(Y(r),m)+h2) where Y(r)=pow(r-r1,2)
c by Guo Chin Lin
 
c check enbt:-enb (energy)  
cccccccccccccc ENERGY part
          if(mapgo(i,j).eq.1) then

             F2=0
             F1=0

             ki=kpower
             ni=npower
             t=ni
             mi=mpower
             u=(r1(i,j)-r0(i,j))**(2*ni)
                                         

             if(rd.lt.r0(i,j)) then       

                s=r0(i,j)/rd
                s=s**ki
c Sep. 15, 2014, Jianfa Chen:
c                fr=DB0(i,j)*s*(s-2.0)
                if(DB0(i,j).gt.0) then
                       fr=DB0(i,j)*s*(s-2.0)
                else
                       fr=0-DB0(i,j)*s*(s-2.0)-2*DB0(i,j)
                endif
cccEnd Jianfa Chen

                ENBT=ENBT-fr          

 
             else  

                y=rd-r1(i,j)
                y=y*y    
 
                if(rd.lt.r1(i,j)) then
 
                 s=y**ni
                 fr=CB1(i,j)*0.5*s*(s*0.5-u)/t+DB1(i,j)
                 ENBT=ENBT-fr   
                else
                   s=y**mi
                   fr=-B2(i,j)*(y-h1(i,j))/(s+h2(i,j))
                 ENBT=ENBT-fr
 
                endif
                          
             endif      

c                   if(mapgo(j,i).eq.1) then
c                       write(79,*)j,i,RD,fr
c                      endif
c                    if(j.eq.57 .and. i.eq.2) then
c                      write(79,*)rd,fr
c                     endif          

c original            ENBT = ENBT + (F2-F1)

cccccccccccccccccccccccccccccccc FORCE    


c origin              DF = (DF2+SIX*((F2-F1)-F1))*R2INV

c margaret new DF

                   if(rd.lt.r0(i,j)) then
 
                      ri=r0(i,j)/rd
                      DF=1-ri**ki      
                      DF=DF*2*ki*(r0(i,j)**ki)
                      DF=DF/(rd**(ki+1))
c Sep.15, 2014 Jianfa Chen:
c                      DF=DB0(i,j)*DF/rd
                      if(DB0(i,j).gt.0) then
                                DF=DB0(i,j)*DF/rd
                      else
                                DF=0-DB0(i,j)*DF/rd
                      endif
cccEnd Jianfa Chen
 
                      else                                               

                        ri=rd-r1(i,j) 

                         if(rd.lt.r1(i,j)) then 

                          rj=r1(i,j)-r0(i,j)
                          temp1=(rj**(2*ni))*(ri**(2*ni-1))
                          temp1=ri**(4*ni-1)-temp1
                          DF=CB1(i,j)*temp1/rd
                                                  
                         else    

                          mi2=mi*2.0D0
                          temp1=(1-mi)*(ri**mi2)+h2(i,j)
                          temp2=h1(i,j)*mi*(ri**(mi2-2))
                          temp1=temp1+temp2
                          temp2=(ri**mi2)+h2(i,j)
                          temp2=temp2*temp2
                          DF=-B2(i,j)*2*ri*temp1/temp2
                          DF=DF/rd

                        endif       
                      endif  

c                    if(j.eq.57 .and. i.eq.2) then
c                      write(80,*)rd, dff
c                     endif
c check the sign  of the force: df goes with dU/dr    


	  else  
c else the regular L-J hard core.

          ENBT = ENBT + (F2-F1)    	 
	  DF = (DF2+SIX*((F2-F1)-F1))*R2INV   
         
          endif

c end of margaret:default
              
		


                      FW1 = XW1*DF
                      FW2 = XW2*DF
                      FW3 = XW3*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
  190             CONTINUE
              else
c                 --- soft-repulsion form, constant dielectric, 
c                                non-periodic case ---
                  IC = IAC(I)
cforcevector
                  DO 200 JN = 1,LPR
#ifdef CRAYFISH
                      j = jpw(jn)
#else
                      J = ipairs(JN+LPAIR)
#endif
                      XW1 = X(1,I)-X(1,J)
                      XW2 = X(2,I)-X(2,J)
                      XW3 = X(3,I)-X(3,J)
                      RR = (XW1**2+XW2**2+XW3**2)
                      R2INV = ONE/RR
c
                      JC = IAC(J)
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      RK = RWELL*lfac
#else
                      RK = RWELL
#endif
                      RST = RAD(IC) + RAD(JC)
                      ARG = MAX(ZERO,(RST*RST-RR))
                      F1 = RK * ARG*ARG
                      ENBT = ENBT - F1
                      DF1 = -FOUR * RK * ARG
c
#ifdef LES
                      DF2 = CGI*CG(J)*SQRT(R2INV)*lfac
#else
                      DF2 = CGI*CG(J)*SQRT(R2INV)
#endif
                      EELT = EELT+DF2
c
                      DF = DF1 + DF2*R2INV
                      FW1 = XW1*DF
                      FW2 = XW2*DF
                      FW3 = XW3*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
  200             CONTINUE
              endif
          endif
c         -- end of dielectric for normal pairs
          F(1,I) = F(1,I)-dumx
          F(2,I) = F(2,I)-dumy
          F(3,I) = F(3,I)-dumz
C
C         --- H-BOND PAIRS 10-12 POTENTIAL ---
C
          dumx = ZERO
          dumy = ZERO
          dumz = ZERO
          if (dield) then
c
c             --- distance dependent dielectric
c
              if (isftrp.le.1) then
c                 --- 10-12 form, distance dependent dielectric, 
c                                 non-periodic case ---
cforcevector
                  DO 210 JN = LPR+1,NPR
#ifdef CRAYFISH
                      j = jpw(jn)
#else
                      J = ipairs(JN+LPAIR)
#endif
                      XW1 = X(1,I)-X(1,J)
                      XW2 = X(2,I)-X(2,J)
                      XW3 = X(3,I)-X(3,J)
                      R2INV = ONE/(XW1**2+XW2**2+XW3**2)
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      DF2 = CGI*CG(J)*R2INV*lfac
#else
                      DF2 = CGI*CG(J)*R2INV
#endif
                      EELT = EELT+DF2
                      IC = -ICO(IACI+IAC(J))
                      R10 = R2INV**5
#ifdef LES
                      F1 = ASOL(IC)*R10*R2INV*lfac
                      F2 = BSOL(IC)*R10*lfac
#else
                      F1 = ASOL(IC)*R10*R2INV
                      F2 = BSOL(IC)*R10
#endif
                      DF = (DF2-TWELVE*F1+TEN*F2)*R2INV
                      EHBT = EHBT+F1-F2
                      FW1 = XW1*DF
                      FW2 = XW2*DF
                      FW3 = XW3*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
  210             CONTINUE
              else
c                 --- soft-repulsion form, distance dependent dielectric, 
c                                          non-periodic case ---
cforcevector
                  DO 220 JN = LPR+1,NPR
#ifdef CRAYFISH
                      j = jpw(jn)
#else
                      J = ipairs(JN+LPAIR)
#endif
                      XW1 = X(1,I)-X(1,J)
                      XW2 = X(2,I)-X(2,J)
                      XW3 = X(3,I)-X(3,J)
                      RR = (XW1**2+XW2**2+XW3**2)
                      R2INV = ONE/RR
C
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      RK = RWELL*lfac
#else
                      RK = RWELL
#endif
                      IC = -ICO(IACI+IAC(J))
                      RST = RADHB(IC)
                      ARG = MAX(ZERO,(RST*RST - RR))
                      F1 = RK * ARG*ARG
                      EHBT = EHBT + F1
                      DF1 = -FOUR * RK * ARG
C
#ifdef LES
                      DF2 = CGI*CG(J)*R2INV*lfac
#else
                      DF2 = CGI*CG(J)*R2INV
#endif
                      EELT = EELT+DF2
C
                      DF = DF1 + DF2*R2INV
                      FW1 = XW1*DF
                      FW2 = XW2*DF
                      FW3 = XW3*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
  220             CONTINUE
              endif
          else
c
c             --- constant dielectric
c
              if (isftrp.le.1) then
c                 --- 10-12 hbond form, constant dielectric, 
c                                       non-periodic case ---
cforcevector
                  DO 230 JN = LPR+1,NPR
#ifdef CRAYFISH
                      j = jpw(jn)
#else
                      J = ipairs(JN+LPAIR)
#endif
                      XW1 = X(1,I)-X(1,J)
                      XW2 = X(2,I)-X(2,J)
                      XW3 = X(3,I)-X(3,J)
                      R2INV = ONE/(XW1**2+XW2**2+XW3**2)
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      DF2 = CGI*CG(J)*SQRT(R2INV)*lfac
#else
                      DF2 = CGI*CG(J)*SQRT(R2INV)
#endif
                      EELT = EELT+DF2
                      IC = -ICO(IACI+IAC(J))
                      R10 = R2INV**5
#ifdef LES
                      F1 = ASOL(IC)*R10*R2INV*lfac
                      F2 = BSOL(IC)*R10*lfac
#else
                      F1 = ASOL(IC)*R10*R2INV
                      F2 = BSOL(IC)*R10
#endif
                      DF = (DF2-TWELVE*F1+TEN*F2)*R2INV
                      EHBT = EHBT+F1-F2
                      FW1 = XW1*DF
                      FW2 = XW2*DF
                      FW3 = XW3*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
  230             CONTINUE
              else
c                 --- soft-repulsion hbond form, constant dielectric, 
c                                       non-periodic case ---
cforcevector
                  DO 240 JN = LPR+1,NPR
#ifdef CRAYFISH
                      j = jpw(jn)
#else
                      J = ipairs(JN+LPAIR)
#endif
                      XW1 = X(1,I)-X(1,J)
                      XW2 = X(2,I)-X(2,J)
                      XW3 = X(3,I)-X(3,J)
                      RR = (XW1**2+XW2**2+XW3**2)
                      R2INV = ONE/RR
C
#ifdef LES
                      lfac=lesfac(lestmp+lestyp(j))
                      RK = RWELL*lfac
#else
                      RK = RWELL
#endif
                      IC = -ICO(IACI+IAC(J))
                      RST = RADHB(IC)
                      ARG = MAX(ZERO,(RST*RST - RR))
                      F1 = RK * ARG*ARG
                      EHBT = EHBT + F1
                      DF1 = -FOUR * RK * ARG
C
#ifdef LES
                      DF2 = CGI*CG(J)*SQRT(R2INV)*lfac
#else
                      DF2 = CGI*CG(J)*SQRT(R2INV)
#endif
                      EELT = EELT+DF2
C
                      DF = DF1 + DF2*R2INV
                      FW1 = XW1*DF
                      FW2 = XW2*DF
                      FW3 = XW3*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
  240             CONTINUE
              endif
          endif
c         -- end dielectric for hbonds
          F(1,I) = F(1,I)-dumx
          F(2,I) = F(2,I)-dumy
          F(3,I) = F(3,I)-dumz
        endif
c       -- end periodicity / no periodicity (used to be 250)
C
C       --- END OF PAIRS INVOLVING ATOM I ---
C
#ifdef SHARED_MEMORY
c     ========================= SHARED MEMORY =========================
c
c       --- END of load balancing shared memory ---
c
      endif
c
c     ======================= END SHARED MEMORY =======================
#endif
  260   continue
c       -- 260 is for skipping an atom w/out pairs
#ifdef CRAYFISH
        lpack = lpack+npack
#endif
        LPAIR = LPAIR+NPR
  270 CONTINUE
C
C--------------------------------------------------------------------
C If the conditions have been met to carry out the special-case
C fast loop TIP3P-TIP3P water interaction calculations, make a call
C to do these now.
C--------------------------------------------------------------------
C
#ifdef MPI
  271 IF (doqik) THEN
#else
      IF (NSOLW.NE.0) THEN
#endif
            CALL QIKTIP(    CG        ,IAC       ,ICO       ,CN1       ,
     *           CN2       ,X         ,F         ,IAR1      ,ipairs    ,
     *           JPW       ,XRC       ,FW        ,XWIJ      ,RW        ,
     *           RW2       ,RW3       ,NTYPES    ,NATOM     ,NUMPK     ,
     *           NBIT      ,LPACK     ,LPAIR     ,NTB       ,IAR1D1    ,
     *           DIELD     ,VIRIAL    ,EELT      ,ENBT      ,VIRT      ,
#ifndef SHARED_MEMORY
# ifndef MPI
     *        NSOLW     ,NATOM-1)
# else
     *        NSOLW     ,NATOMCNT     ,ITMP      ,natomnum)
# endif
#else 
     *        NSOLW     ,NATOM-1       ,Iproc     ,NumProc)
#endif
      END IF

      if (dield) then
           EEL = -EELT*HALF
      else
           EEL = -EELT
      endif
      ENB = -ENBT
      EHB = EHBT
      VIR(1) = VIRT(1)*HALF
      VIR(2) = VIRT(2)*HALF
      VIR(3) = VIRT(3)*HALF

      RETURN
      END







