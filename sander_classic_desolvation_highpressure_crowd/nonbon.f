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
Cc
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
c Andrei added 1.7.2016 (double well potential)
#include "douwell.h"
c end Andrei
c Margaret
#include "HB.h"
        REAL cosphi,cosphip,acosphi
        REAL XIJ,YIJ,ZIJ,XKJ,YKJ,ZKJ
        REAL XLK,YLK,ZLK,TX,TY,TZ,UX,UY,UZ
        REAL DTXX, DTY,DTZ,DUX,DUY,DUZ
	REAL CHI, U
#include "CHI.h"
#include "CROWD.h"
c Margaret
      integer NPP 



#ifdef LES
#  include "les.h"
#endif
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

c Margaret
      PARAMETER (PI = 3.141592653589793)
      data pt999 /0.9999d0/
      data pt000 /0.0001d0/
c      data Amp /1d0/
c      data phia /0.465999d0/
      data cosphia /0.893373d0/
      data MAXM /10d5/


c Margaret
      data MAXLIM/5.0D2/
C
      dield  = (idiel .le. 0)
      sltimg = (imgslt .le. 0)
      VIRIAL = (IABS(NTB) .GE. 2)
      LPAIR  = 0
      lpack  = 1
      ENBT   = ZERO
      EELT   = ZERO
      EHBT   = ZERO

c! Margaret add EHBA for backbone Hydrogen Bond 6-12 energies
c! add EHBV for van der waal's 6-12 energies
      EHBV = ZERO
      EHBA   = ZERO
      ECHI   = ZERO

c Margaret
c ENPP, 6-12 energy between proteins. ENPC 6-12 energy between protein
c and crowding agents
c NPP number of proteins
      ENPC   = ZERO
      ENCC   = ZERO
      NPP    = NTYPES-1

c end Margaret

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
c	        write(76,*) i,j,iptatm
c
              if (sltimg) then
c
c                 --- IMGSLT=0 => imaging solute with solvent ---
c
c	        write(76,*) i,j,iptatm

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

c	   write(76,*)i,j
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

c constrant D. constant Volume. 6-12. periodic 
c

#ifdef CRAYFISH
                      j = jpw(jn)
#else 
                      J = ipairs(JN+LPAIR)
#endif

c Jan 7, 2016, Andrei Gasic:
            RD2 = XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2
            RD = SQRT(RD2)
            R2INV = ONE/(RD2)
c end Andrei

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

c Margaret change here, periodic, 
                  if(maphb(i,j).lt.1) then
c                  EHBV=EHBV+F1-F2

c Jan 7, 2016 Andrei Gasic:
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

                   if(DB0(i,j).gt.0) then
                       fr=DB0(i,j)*s*(s-2.0)
                   else
                       fr=0-DB0(i,j)*s*(s-2.0)-2*DB0(i,j)
                   endif

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

ccc check energy for debuging purposes
c                   if(mapgo(j,i).eq.1) then
c                       write(999,*)j,i,RD,fr
c                     endif

cccccccccccccccccccccccccccccccc FORCE
               if(rd.lt.r0(i,j)) then

                 ri=r0(i,j)/rd
                 DF=1-ri**ki
                 DF=DF*2*ki*(r0(i,j)**ki)

                 if(DB0(i,j).gt.0) then
                     DF=DB0(i,j)*DF/rd
ccc check energy for debuging purposes
c                       DF=DF*rd
c                       write(998,*)j,i,RD,DF
c                       DF=DF/rd

                 else
                     DF=0-DB0(i,j)*DF/rd
ccc check energy for debuging purposes
c                       DF=DF*rd
c                       write(997,*)j,i,RD,DF
c                       DF=DF/rd

                 endif


               else

                  ri=rd-r1(i,j)

                  if(rd.lt.r1(i,j)) then

                      rj=r1(i,j)-r0(i,j)
                      temp1=(rj**(2*ni))*(ri**(2*ni-1))
                      temp1=ri**(4*ni-1)-temp1
                      DF=CB1(i,j)*temp1/rd

ccc check energy for debuging purposes
c                       DF=DF*rd
c                       write(996,*)j,i,RD,DF
c                       DF=DF/rd


                  else

                      mi2=mi*2.0D0
                      temp1=(1-mi)*(ri**mi2)+h2(i,j)
                      temp2=h1(i,j)*mi*(ri**(mi2-2))
                      temp1=temp1+temp2
                      temp2=(ri**mi2)+h2(i,j)
                      temp2=temp2*temp2
                      DF=-B2(i,j)*2*ri*temp1/temp2

ccc check energy for debuging purposes
c                       write(995,*)j,i,RD,DF


                      DF=DF/rd

                   endif
                 endif

ccc check energy for debuging purposes
c                   if(mapgo(j,i).eq.1) then
c                       DF=DF*rd
c                       write(999,*)j,i,RD,DF
c                     endif


c check the sign  of the force: df goes with dU/dr 

          else
c else the regular L-J hard core.
                      ENBT = ENBT + (F2-F1)
                      DF = (DF2+SIX*((F2-F1)-F1))*RW(JN)
         endif
ccc End Andrei Gasic

                      FW1 = XWIJ(1,JN)*DF
                      FW2 = XWIJ(2,JN)*DF
                      FW3 = XWIJ(3,JN)*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3

c crowding

                     if (i.le.NPP.and.j.gt.NPP) then
                     ENPC = ENPC - (F2-F1)
c                    write(71,*)i,j,ENPC,mark(i),mark(j)

                     elseif (i.gt.NPP.and.j.le.NPP) then
                     ENPC = ENPC - (F2-F1)
c                    write(72,*)i,j,ENPC,mark(i),mark(j)

                     elseif (i.gt.NPP.and.j.gt.NPP) then
                     ENCC = ENCC - (F2-F1)
c                    write(73,*)i,j,ENCC,mark(i),mark(j)
                     elseif (i.le.NPP.and.j.le.NPP) then
                     EHBV=EHBV+F1-F2
c                    write(73,*)i,j,ENCC,mark(i),mark(j)

                     endif





                  else


              num= maphb(i,j)
              U=F1-F2
                      DF = (DF2+SIX*((F2-F1)-F1))*RW(JN)
                      FW1 = XWIJ(1,JN)*DF
                      FW2 = XWIJ(2,JN)*DF
                      FW3 = XWIJ(3,JN)*DF
              ixi=NXI(num)
              ixj=NXI1(num)
              ixk=NXI2(num)
              ixl=NXI3(num)

c check if ixk=i,ixj=j, OK
c Margaret
                      XIJ = X(1,ixi)-X(1,ixj)
                      YIJ = X(2,ixi)-X(2,ixj)
                      ZIJ = X(3,ixi)-X(3,ixj)
                      XKJ = XWIJ(1,JN)
                      YKJ = XWIJ(2,JN)
                      ZKJ = XWIJ(3,JN)
                      XLK = X(1,ixl)-X(1,ixk)
                      YLK = X(2,ixl)-X(2,ixk)
                      ZLK = X(3,ixl)-X(3,ixk)

c imaging.this has to be itself. because xij<<box

c		write (84,*)xij,yij,zij
		  xij = xij - anint(xij*boxi(1)) *box(1)
		  yij = yij - anint(yij*boxi(2)) *box(2)
		  zij = zij - anint(zij*boxi(3)) *box(3)
c		write (84,*)xij,yij,zij
 
		  xlk = xlk - anint(xlk*boxi(1)) *box(1)
		  ylk = ylk - anint(ylk*boxi(2)) *box(2)
		  zlk = zlk - anint(zlk*boxi(3)) *box(3)


        TX=(YIJ*ZKJ-YKJ*ZIJ)
        TY=(XKJ*ZIJ-XIJ*ZKJ)
        TZ=(XIJ*YKJ-XKJ*YIJ)

        UX=(YLK*ZKJ-YKJ*ZLK)    
        UY=(XKJ*ZLK-XLK*ZKJ)
        UZ=(XLK*YKJ-XKJ*YLK)
        
        temp=TX*UX+TY*UY+TZ*UZ
        RT=TX*TX+TY*TY+TZ*TZ    
        RT=sqrt(RT)
        RU=UX*UX+UY*UY+UZ*UZ    
        RU=sqrt(RU)
        cosphi=temp/(RT*RU)


c ----- calculate structural factor
cccc energy
c the unperturbed cosphi

c           sinphi=(1-cosphi**2)
c           sinphi=sqrt(sinphi)

       tempA=1-cosphi
       tempB=cosphi+1
       tempC=cosphi/cosphia
       tempC=1-tempC

       dummy=tempA*tempB*tempC
       dummy=1+Amp*dummy*dummy
       Chi=1/dummy


c ENBT has a (-) change at the end to ENB
	       ENBT = ENBT-Chi*U


               EHBA = EHBA+Chi*U

c       if(i.eq.99.and.j.eq.115) then
c          write(76,*)Chi,U,Chi*U,cosphi,CONVERT*acos(cosphi)
c       endif


cccccc Forces

c dChi/dcosphi
       temp=-tempB*tempC+tempA*tempC-tempA*tempB/cosphia
       temp=temp*tempA*tempB*tempC
       temp=temp*2*Amp
       dChi=-temp/(dummy**2)

c partial derivatives

        RT2=RT*RT
        RU2=RU*RU
        RTU=RT*RU
        DTXX=UX/RTU-cosphi*TX/RT2
        DTY=UY/RTU-cosphi*TY/RT2
        DTZ=UZ/RTU-cosphi*TZ/RT2
        
        DUX=TX/RTU-cosphi*UX/RU2
        DUY=TY/RTU-cosphi*UY/RU2
        DUZ=TZ/RTU-cosphi*UZ/RU2


c force dcos/dxi 
        FXI=DTY*(-ZKJ)+DTZ*(YKJ)

c force dcos/dyi
        FYI=DTXX*(ZKJ)+DTZ*(-XKJ)

c force dcos/dzi
        FZI=DTXX*(-YKJ)+DTY*(XKJ)


c force dcos/dxj
        FXJ=DTY*(ZKJ-ZIJ)+DTZ*(-YKJ+YIJ)
     *     +DUY*(-ZLK) + DUZ*(YLK)

c force dcos/dyj
        FYJ=DTXX*(ZIJ-ZKJ)+DTZ*(XKJ-XIJ)
     *     +DUX*(ZLK) + DUZ*(-XLK)

c force dcos/dzj
        FZJ=DTXX*(YKJ-YIJ)+DTY*(XIJ-XKJ)
     *     +DUX*(-YLK) + DUY*(XLK)



c force dcos/dxk
        FXK=DTY*(ZIJ)+DTZ*(-YIJ)
     *     +DUY*(ZKJ+ZLK) + DUZ*(-YKJ-YLK)

c force dcos/dyk
        FYK=DTXX*(-ZIJ)+DTZ*(XIJ)
     *     +DUX*(-ZKJ-ZLK) + DUZ*(XLK+XKJ)

c force dcos/dzk
        FZK=DTXX*(YIJ)+DTY*(-XIJ)
     *     +DUX*(YLK+YKJ) + DUY*(-XLK-XKJ)


c force dcos/dxl
        FXl= DUY*(-ZKJ) + DUZ*(YKJ)

c force dcos/dyl
        FYl= DUX*(ZKJ) + DUZ*(-XKJ)

c force dcos/dzl
        FZl= DUX*(-YKJ) + DUY*(XKJ)




c force -dP/dx

        FFxi=-U*dChi*FXI
        FFyi=-U*dChi*FYI
        FFzi=-U*dChi*FZI
        FFxj=-U*dChi*FXJ
        FFyj=-U*dChi*FYJ
        FFzj=-U*dChi*FZJ
        FFxk=-U*dChi*FXK
        FFyk=-U*dChi*FYK
        FFzk=-U*dChi*FZK
        FFxl=-U*dChi*FXL
        FFyl=-U*dChi*FYL
        FFzl=-U*dChi*FZL


c -dp/dxi
                      N=NXI(num)
                      F(1,N) = F(1,N)+FFxi
                      F(2,N) = F(2,N)+FFyi
                      F(3,N) = F(3,N)+FFzi


c -dp/dxj
                      N=NXI1(num)
               F(1,N) = F(1,N)+FFxj+Chi*FW1
               F(2,N) = F(2,N)+FFyj+Chi*FW2
               F(3,N) = F(3,N)+FFzj+Chi*FW3

c-dp/dxk
                      N=NXI2(num)
               F(1,N) = F(1,N)+FFxk-Chi*FW1
               F(2,N) = F(2,N)+FFyk-Chi*FW2
               F(3,N) = F(3,N)+FFzk-Chi*FW3
        
c dp/dxl
                      N=NXI3(num)
                      F(1,N) = F(1,N)+FFxl
                      F(2,N) = F(2,N)+FFyl
                      F(3,N) = F(3,N)+FFzl


c Mark
c       if(U.gt.5.or.(U*Chi).gt.5) then
c               write(82,*)i,j,U,U*Chi,acosphi
c       endif







                  endif

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

c Andrei Gasic
                      RD = SQRT(XW1**2+XW2**2+XW3**2)
c end of Andrei 

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

c Margaret change here



c Note: HB force field U is modifed by pseudo dihedral C.
c potential term P
c P=C(cos(pseudo dihedral phi(xi,xj,xk,xl) )U(r_hb)
c since we dont care about sign of the angle, use 
c formula in J.Comp.Chem.vol13.1992.585
c consider rij=xi-xj,rkj=xk-xj,rlk=xl-xk
c XIJ=XI-XJ, YIJ=YI-YJ, ZIJ=ZI-ZJ,
c XKJ=XK-XJ, YKJ=YK-YJ, ZKJ=ZK-ZJ,
c XLK=XL-XK, YLK=YL-YK, ZLK=ZL-ZK 
c t=|rkj x rkj|= i(YIJZKJ-YKJZIJ) - j(XIJZIJ-XKJZIJ) + k(XIJYKJ-XKJYIJ)
c u=|rlk x rkj|= i(YLKZKJ-YKJZLK) - j(XLKZKJ-XKJZLK) +k(XLKYKJ-XKJYLK)
c cosphi=t*u/|t||u|
c r_hb=|rkj|
c U=U(rkj)
c A=
c phia= 0.465999 (rad) [i-1,i,i+4,i+5]
c C=1/1+A[(1-cosphi)(1+cosphi)(1-cosphi/cosphia)]**2
c Grad_i C(phi)= (dC(phi)/dcosphi)*(Grad_i cosphi)
c
c B=1+A[(1-cosphi)(1+cosphi)(1-cosphi/cosphia)]**2
c tempA=(1-cosphi),tempB=(1+cosphi),tempC=(1-cosphi/cosphia)
c dC/dcosphi = 2A(tempA*tempB*tempC))*
c {-tempB*tempC+tempA*tempC-tempA*tempB/cosphia}/B**2
c
c dcosphi/dxi = dcosphi/dty(-ZKJ) + dcosphi/dtz(YKJ) 
c dcosphi/dyi = dcosphi/dtx(ZKJ) + dcosphi/dtz(-XKJ) 
c dcosphi/dzi = dcosphi/dtx(-YKJ) + dcosphi/dtz(XKJ) 
c 
c dcosphi/dxj = dcosphi/dty(ZKJ-ZIJ) + dcosphi/dtz(-YKJ+YIJ) 
c              +dcosphi/duy(-ZLK) + dcosphi/duz(YLK)
c dcosphi/dyj = dcosphi/dtx(-ZKJ+ZIJ) + dcosphi/dtz(XKJ-XIJ) 
c              +dcosphi/dux(ZLK) + dcosphi/duz(-XLK)
c dcosphi/dzj = dcosphi/dtx(YKJ-YIJ) + dcosphi/dty(XIJ-XKJ) 
c              +dcosphi/dux(-YLK) + dcosphi/duy(XLK)
 
c dcosphi/dxk = dcosphi/dty(ZIJ) + dcosphi/dtz(-YIJ) 
c              +dcosphi/duy(ZLK+XKJ) + dcosphi/duz(-YLK-YKJ)
c dcosphi/dyk = dcosphi/dtx(-ZIJ) + dcosphi/dtz(XIJ) 
c              +dcosphi/dux(-ZLK-ZKJ) + dcosphi/duz(XLK+XKJ)
c dcosphi/dzk = dcosphi/dtx(YIJ) + dcosphi/dty(-XIJ) 
c              +dcosphi/dux(YLK+YKJ) + dcosphi/duy(-XLK-XKJ)
c
c dcosphi/dxl = dcosphi/duy(-ZKJ) + dcosphi/duz(YKJ)
c dcosphi/dyl = dcosphi/dux(ZKJ) + dcosphi/duz(-XKJ)
c dcosphi/dzl = dcosphi/dux(-YKJ) + dcosphi/duy(+XKJ)

c dU/dx=dF
c dF={-60CN1 (r**-13) + 60CN2 (r**-11)} * (r2)/|r2|
c -dP/dxi= U(rjk)*dC/dcosphi*dcosphi/dxi
c -dP/dxj=C*dF + U(rjk)*dC/dcosphi*dcosphi/dxj  (xj corresponding to J)
c -dP/dxk=-C*dF + U(rjk)*dC/dcosphi*dcosphi/dxk  (xk corresponding to I)
c -dP/dxl=U(r2)*dC/dcosphi*dcosphi/dxl
c 


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



                  if(maphb(i,j).lt.1) then
c no HB
c F_J= dU/dx=dF F1:rep, F2:attr
c                   EHBV=EHBV+F1-F2
c modified by Andrei Gasic 
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

ccc check energy for debuging purposes
c                   if(mapgo(j,i).eq.1) then
c                       write(78,*)j,i,RD,fr
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

c end Andrei 


c crowding

                     if (i.le.NPP.and.j.gt.NPP) then
                     ENPC = ENPC - (F2-F1)
c                    write(71,*)i,j,ENPC

                     elseif (i.gt.NPP.and.j.le.NPP) then
                     ENPC = ENPC - (F2-F1)
c                    write(72,*)i,j,ENPC

                     elseif (i.gt.NPP.and.j.gt.NPP) then
                     ENCC = ENCC - (F2-F1)
c                    write(73,*)i,j,ENCC

                     elseif (i.le.NPP.and.j.le.NPP) then
                     EHBV=EHBV+F1-F2
c                    write(73,*)i,j,ENCC,mark(i),mark(j)


                     endif


    
                      FW1 = XW1*DF
                      FW2 = XW2*DF
                      FW3 = XW3*DF
                      dumx = dumx + fw1
                      dumy = dumy + fw2
                      dumz = dumz + fw3
                      F(1,J) = F(1,J)+FW1
                      F(2,J) = F(2,J)+FW2
                      F(3,J) = F(3,J)+FW3
                  else
c yes HB
              num= maphb(i,j)
              U=F1-F2
                      DF = (DF2+SIX*((F2-F1)-F1))*R2INV
                      FW1 = XW1*DF
                      FW2 = XW2*DF
                      FW3 = XW3*DF
              ixi=NXI(num)
              ixj=NXI1(num)
              ixk=NXI2(num)
              ixl=NXI3(num)

c         write(72,*)i,j,ixi,ixj,ixk,ixl 
c check if ixk=i,ixj=j, OK.


                      XIJ = X(1,ixi)-X(1,ixj)
                      YIJ = X(2,ixi)-X(2,ixj)
                      ZIJ = X(3,ixi)-X(3,ixj)
                      XKJ = XW1
                      YKJ = XW2
                      ZKJ = XW3
                      XLK = X(1,ixl)-X(1,ixk)
                      YLK = X(2,ixl)-X(2,ixk)
                      ZLK = X(3,ixl)-X(3,ixk)

        TX=(YIJ*ZKJ-YKJ*ZIJ)
        TY=(XKJ*ZIJ-XIJ*ZKJ)
        TZ=(XIJ*YKJ-XKJ*YIJ)

        UX=(YLK*ZKJ-YKJ*ZLK)    
        UY=(XKJ*ZLK-XLK*ZKJ)
        UZ=(XLK*YKJ-XKJ*YLK)
        
        temp=TX*UX+TY*UY+TZ*UZ
        RT=TX*TX+TY*TY+TZ*TZ    
        RT=sqrt(RT)
        RU=UX*UX+UY*UY+UZ*UZ    
        RU=sqrt(RU)
        cosphi=temp/(RT*RU)


c ----- calculate structural factor
cccc energy
c the unperturbed cosphi

c           sinphi=(1-cosphi**2)
c           sinphi=sqrt(sinphi)

       tempA=1-cosphi
       tempB=cosphi+1
       tempC=cosphi/cosphia
       tempC=1-tempC

       dummy=tempA*tempB*tempC
       dummy=1+Amp*dummy*dummy
       Chi=1/dummy


c ENBT has a (-) change at the end to ENB
	       ENBT = ENBT-Chi*U


               EHBA = EHBA+Chi*U

c       if(i.eq.99.and.j.eq.115) then
c          write(76,*)Chi,U,Chi*U,cosphi,CONVERT*acos(cosphi)
c       endif


cccccc Forces

c dChi/dcosphi
       temp=-tempB*tempC+tempA*tempC-tempA*tempB/cosphia
       temp=temp*tempA*tempB*tempC
       temp=temp*2*Amp
       dChi=-temp/(dummy**2)


c partial derivatives

        RT2=RT*RT
        RU2=RU*RU
        RTU=RT*RU
        DTXX=UX/RTU-cosphi*TX/RT2
        DTY=UY/RTU-cosphi*TY/RT2
        DTZ=UZ/RTU-cosphi*TZ/RT2
        
        DUX=TX/RTU-cosphi*UX/RU2
        DUY=TY/RTU-cosphi*UY/RU2
        DUZ=TZ/RTU-cosphi*UZ/RU2


c force dcos/dxi 
        FXI=DTY*(-ZKJ)+DTZ*(YKJ)

c force dcos/dyi
        FYI=DTXX*(ZKJ)+DTZ*(-XKJ)

c force dcos/dzi
        FZI=DTXX*(-YKJ)+DTY*(XKJ)


c force dcos/dxj
        FXJ=DTY*(ZKJ-ZIJ)+DTZ*(-YKJ+YIJ)
     *     +DUY*(-ZLK) + DUZ*(YLK)

c force dcos/dyj
        FYJ=DTXX*(ZIJ-ZKJ)+DTZ*(XKJ-XIJ)
     *     +DUX*(ZLK) + DUZ*(-XLK)

c force dcos/dzj
        FZJ=DTXX*(YKJ-YIJ)+DTY*(XIJ-XKJ)
     *     +DUX*(-YLK) + DUY*(XLK)



c force dcos/dxk
        FXK=DTY*(ZIJ)+DTZ*(-YIJ)
     *     +DUY*(ZKJ+ZLK) + DUZ*(-YKJ-YLK)

c force dcos/dyk
        FYK=DTXX*(-ZIJ)+DTZ*(XIJ)
     *     +DUX*(-ZKJ-ZLK) + DUZ*(XLK+XKJ)

c force dcos/dzk
        FZK=DTXX*(YIJ)+DTY*(-XIJ)
     *     +DUX*(YLK+YKJ) + DUY*(-XLK-XKJ)


c force dcos/dxl
        FXl= DUY*(-ZKJ) + DUZ*(YKJ)

c force dcos/dyl
        FYl= DUX*(ZKJ) + DUZ*(-XKJ)

c force dcos/dzl
        FZl= DUX*(-YKJ) + DUY*(XKJ)




c force -dP/dx

        FFxi=-U*dChi*FXI
        FFyi=-U*dChi*FYI
        FFzi=-U*dChi*FZI
        FFxj=-U*dChi*FXJ
        FFyj=-U*dChi*FYJ
        FFzj=-U*dChi*FZJ
        FFxk=-U*dChi*FXK
        FFyk=-U*dChi*FYK
        FFzk=-U*dChi*FZK
        FFxl=-U*dChi*FXL
        FFyl=-U*dChi*FYL
        FFzl=-U*dChi*FZL


c -dp/dxi
                      N=NXI(num)
                      F(1,N) = F(1,N)+FFxi
                      F(2,N) = F(2,N)+FFyi
                      F(3,N) = F(3,N)+FFzi

c -dp/dxj
                      N=NXI1(num)
               F(1,N) = F(1,N)+FFxj+Chi*FW1
               F(2,N) = F(2,N)+FFyj+Chi*FW2
               F(3,N) = F(3,N)+FFzj+Chi*FW3
c -dp/dxk
                      N=NXI2(num)
               F(1,N) = F(1,N)+FFxk-Chi*FW1
               F(2,N) = F(2,N)+FFyk-Chi*FW2
               F(3,N) = F(3,N)+FFzk-Chi*FW3
        
c dp/dxl
                      N=NXI3(num)
                      F(1,N) = F(1,N)+FFxl
                      F(2,N) = F(2,N)+FFyl
                      F(3,N) = F(3,N)+FFzl




c Mark
c       if(U.gt.5.or.(U*Chi).gt.5) then
c               write(82,*)i,j,U,U*Chi,acosphi
c       endif



                  endif

c end of Margaret
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
c	Margaret:default
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
c Margaret change here
c	           if(maphb(i,j).lt.1) then
c		   EHBV=EHBV+F1-F2
c	 	   else
c		   EHBA=EHBA+F1-F2
c	 	   endif
c end of Margaret
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


c Margaret
c add chiral term
c U(chiral)=0.5*xkchi*(triple-triple_o)^2
c triple is the triple scaler product
c triple=(AxB)*C, where A=CB-CA,B=NC-CA, C=CC-CA
c CA=xi,CB=x(i+1),NC=x(i+2),CC=x(i+3)
c Ax=x(i+1)-xi
c Ay=y(i+1)-yi
c Az=z(i+1)-zi
c Bx=x(i+2)-xi
c By=y(i+2)-yi
c Bz=z(i+2)-zi
c Cx=x(i+3)-xi
c Cy=y(i+3)-yi
c Cz=z(i+3)-zi
c F_xi=-dU(chiral)/dxi
c     =-dU/dAx * dAx/dxi - dU/dBx* dBx/dxi - dU/dCx* dCx/dxi
c     = dU/dAx+dU/dBx+dU/dCx
c F_xi1=-dU(chiral)/dxi1
c       = -dU/dAx * dAx/dxi1
c       = -dU/dAx 
c F_xi+2=-dU(chiral)/dxi2
c       = -dU/dBx * dBx/dxi2
c       = -dU/dBx 
c F_xi+3=-dU(chiral)/dxi3
c       = -dU/dCx * dCx/dxi3
c       = -dU/dCx 

	do i=1,ICHI
	
	ixi=ICA(i)
	ixi1=ICB(i)
	ixi2=INC(i)
	ixi3=ICC(i)

	XAX = X(1,ixi1)-X(1,ixi)
        XAY = X(2,ixi1)-X(2,ixi)
        XAZ = X(3,ixi1)-X(3,ixi)

	XBX = X(1,ixi2)-X(1,ixi)
        XBY = X(2,ixi2)-X(2,ixi)
        XBZ = X(3,ixi2)-X(3,ixi)

	XCX = X(1,ixi3)-X(1,ixi)
        XCY = X(2,ixi3)-X(2,ixi)
        XCZ = X(3,ixi3)-X(3,ixi)


	xtriple=  XAX*XBY*XCZ+ XAZ*XBX*XCY
     +		+ XAY*XBZ*XCX- XAZ*XBY*XCX
     +		- XAY*XBX*XCZ- XAX*XBZ*XCY

c calculate energy

	temp=xkchi*(xtriple-XCHI(i))

c        write(85,*) 'temp ', xkchi, xtriple, XCHI(i) 


	ECHII=0.5*temp*(xtriple-XCHI(i))
	ECHI=ECHII+ECHI

c	if(i.eq.1) then
c	write(81,*)xchi(i),xtriple,ECHII,temp
c	write(81,*)ixi,ixi1,ixi2,ixi3
c	endif
c calculate force

	Fxi1= -(XBY*XCZ-XBZ*XCY)*temp
	Fxi2= -(XAZ*XCY-XAY*XCZ)*temp
	Fxi3= -(XAY*XBZ-XAZ*XBY)*temp

c	if(i.eq.1) then
c	write(82,*)fxi1,fxi2,fxi3
c	endif

	Fxi = -(Fxi1+Fxi2+Fxi3)

	Fyi1= -(XBZ*XCX-XBX*XCZ)*temp
	Fyi2= -(XAX*XCZ-XAZ*XCX)*temp
	Fyi3= -(XAZ*XBX-XAX*XBZ)*temp
	Fyi=   -(Fyi1+Fyi2+Fyi3)

c	if(i.eq.1) then
c	write(83,*)fyi1,fyi2,fyi3
c	endif

	Fzi1= -(XBX*XCY-XBY*XCX)*temp
	Fzi2= -(XAY*XCX-XAX*XCY)*temp
	Fzi3= -(XAX*XBY-XAY*XBX)*temp
	Fzi= -(Fzi1+Fzi2+Fzi3)

cc	if(i.eq.1) then
c	write(79,*)fzi1,fzi2,fzi3
c	endif

	
	if(Fxi.gt.MAXLIM) then
c 	write(84,*) 'ixi,x ', ixi,fxi,temp,xtriple
	Fxi=MAXLIM
	endif
	if(Fyi.gt.MAXLIM) then
 	write(84,*) 'ixi,y ', ixi,fyi,temp,xtriple
	Fyi=MAXLIM
	endif
	if(Fzi.gt.MAXLIM) then
	Fzi=MAXLIM
 	write(84,*) 'ixi,z ', ixi,fzi,temp,xtriple
	endif

	if(Fxi1.gt.MAXLIM) then
 	write(84,*) 'ixi1,x ', ixi1,fxi1,temp,xtriple
	Fxi1=MAXLIM
	endif
	if(Fyi1.gt.MAXLIM) then
 	write(84,*) 'ixi1,y ', ixi1,fyi1,temp,xtriple
	Fyi1=MAXLIM
	endif
	if(Fzi1.gt.MAXLIM) then
 	write(84,*) 'ixi1,z ', ixi1,fzi1,temp,xtriple
	Fzi1=MAXLIM
	endif

	if(Fxi2.gt.MAXLIM) then
 	write(84,*) 'ixi2,x ', ixi2,fxi2,temp,xtriple
	Fxi2=MAXLIM
	endif
	if(Fyi2.gt.MAXLIM) then
 	write(84,*) 'ixi2,y ', ixi2,fyi2,temp,xtriple
	Fyi2=MAXLIM
	endif
	if(Fzi2.gt.MAXLIM) then
 	write(84,*) 'ixi2,z ', ixi2,fzi2,temp,xtriple
	Fzi2=MAXLIM
	endif

	if(Fxi3.gt.MAXLIM) then
 	write(84,*) 'ixi3,x ', ixi3,fxi3,temp,xtriple
	Fxi3=MAXLIM
	endif
	if(Fyi3.gt.MAXLIM) then
 	write(84,*) 'ixi3,y ', ixi3,fyi3,temp,xtriple
	Fyi3=MAXLIM
	endif
	if(Fzi3.gt.MAXLIM) then
 	write(84,*) 'ixi3,z ', ixi3,fzi3,temp,xtriple
	Fzi3=MAXLIM
	endif



        F(1,ixi) = F(1,ixi)+Fxi
        F(2,ixi) = F(2,ixi)+Fyi
        F(3,ixi) = F(3,ixi)+Fzi
        F(1,ixi1) = F(1,ixi1)+Fxi1
        F(2,ixi1) = F(2,ixi1)+Fyi1
        F(3,ixi1) = F(3,ixi1)+Fzi1
        F(1,ixi2) = F(1,ixi2)+Fxi2
        F(2,ixi2) = F(2,ixi2)+Fyi2
        F(3,ixi2) = F(3,ixi2)+Fzi2
        F(1,ixi3) = F(1,ixi3)+Fxi3
        F(2,ixi3) = F(2,ixi3)+Fyi3
        F(3,ixi3) = F(3,ixi3)+Fzi3

	enddo	



      RETURN
      END






