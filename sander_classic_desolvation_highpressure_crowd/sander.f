c-----------------------------------------------------------------------
      PROGRAM SANDER
C
C        SANDER, version 5    
C
C The Molecular Dynamics/NMR Refinement/Modeling Module of the AMBER Package.
C
C     This is a modified version of the AMBER 3.0, Rev. A MD Module which
C     includes an extensive suite of utilites for use with NMR refinement,
C     other modeling calculations, as well as other new options within 
C     the bulk of the MD program (minimization is also included).
C
C     The NMR refinement/modeling suite, hooks thereto, and the new options
C     were written by
C
C             David A. Pearlman (UCSF)
C             David A. Case (Scripps) &
C             Ping Yip      (Scripps)
C
C     Version 4.1 also include the capability to carry out polarization
C     calculations. The polarizability code was written by Jim Caldwell
C     and Liem X. Dang (UCSF).  Truncated octahedral periodicity was
c     added by Thomas Huber of Ludwig Maximilian Universitaet Muenchen.
c
C
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995, 1997                **
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
C  Significant changes for Version 4.1:
C
C     1) Inclusion of polarization code. 
C        (J. Caldwell and L. Dang)
C     2) Use of fast analytical shake for 3 point waters.
C        (D. Pearlman and S. Miyamoto)
C     3) New much faster routines to handle TIP3P-TIP3P water interactions
C        (D. Case and D. Pearlman)
C     4) Inclusion of standard and time-averaged J-coupling restraints
C        (D. Pearlman)
C     5) New methods for NMR Intensity refinement and ring current calcs.
C        (D. Case)
C     6) Allow for dual non-bonded cutoffs.
C        (D. Pearlman)
C     7) New VLIMIT option to limit max. atomic velocity.
C        (D. Pearlman)
C     8) Incorporation of PEACS constant nrg contour conf. search capability
C        (D. Case)
C     9) Plus various minor modifications for clearer output, cleanup, 
C        to incorperate bugfixes, etc. 
C    10) New code for truncated octahedral periodic boundary conditions
C        (Thomas Huber, Ludwig Maximilian Universitaet Muenchen,
C         email: thuber@Physik.TU-Muenchen.de) and reorganization/relocation
C         of imaging routines to period.f (B. Ross).
C    11) Particle Mesh Ewald code written and contributed by Tom Darden
C        of the National Institute of Environmental Health Sciences
C        (NIEHS) based on AMBER 3.  (See the manual for citations.)
C        Ported to AMBER 4.1 by Tom Darden with the assistance of 
C        T. Cheatham.
C    12) Support for shared memory parallelization.
C        The original shared memory implementation was done by Roberto
C        Gomperts (SGI) and Michael Schlenkrich (SGI), assisted by 
C        Thomas Cheatham, under AMBER 4.0 and then ported to AMBER 4.1, 
C        cleaned-up and optimized.  CRAY_MP support was added by Mike Page
C        and Jeyapandian Kottalam of Cray Research, with the assistance 
C        of Thomas Cheatham, based on the SGI_MP code.  
C    13) A Message Passing Interface version of AMBER 4.0 and 4.1 was 
C        contributed by James Vincent and Ken Merz of the Pennsylvania 
C        State University.  (See the manual for proper citations.)
C        Optimization and enhancement of the code was performed by 
C        James Vincent, Dave Case, Thomas Cheatham, Jeyapandian 
C        Kottalam (CRI), Asiri Nanayakkara (PSC/CRI), Thomas Huber and 
C        Michael Crowley (PSC).
C
C  Changes for Version 4 (NMR):
C
C     Aside from the various hooks required to
C     integrate the NMR package, 
C     1) the nonbonded code in nonbon and ephi
C        was modified to allow a "soft repulsion" non-bonded potential in
C        place of 6-12 or 10-12 vdw interactions; 
C     2) several new temperature-coupling options were added in RUNMD.
C        These will be particularly useful when carrying out simulations
C        where the internal energy of the molecular system is changing
C        very quickly (such as in some MD/NMR refinement schemes).
C     3) The "NMR" package itself, which allows
C        a large number of simulation protocols appropriate for NMR/MD
C        refinements and general modeling work, and offers a relatively 
C        flexible and easy-to-use interface. See the SANDER refinement 
C        manual for more details.
C
C ============================================================================
C
#ifdef DPREC
      implicit double precision (a-h,o-z)
#endif
      LOGICAL SKIP, BELLY, erstop, convgd
#include "files.h"
#include "sizes.h"
#include "memory.h"
#include "nmr.h"
#include "box.h"
#include "md.h"
#include "parms.h"
#include "other.h"
#include "iewald.h"
#include "extra.h"
#include "lgcom.h"
#include "avfrc.h"
c Margaret Langevin
#include "LAN.h"
#include "HB.h"
#include "CHI.h"
c end of Margaret
c Andrei double well potential
#include "douwell.h"
c end of Andrei 
#ifdef MPI
c     =========================== AMBER/MPI ===========================
#include "parallel.h"
# include "ew_parallel.h"
#include "mpif.h"
#  ifdef MPI_BUFFER_SIZE
      integer*4 mpibuf(MPI_BUFFER_SIZE)
#  endif

      integer imes
      dimension ener(30),vir(4)
c     ========================= END AMBER/MPI =========================
#endif
C
#ifdef CNX_MP
#define SHARED_MEMORY
c     ========================= SHARED MEMORY =========================
c     ...shared memory common block, necessary since in the Convex SPP 
c     implementation, the number of processors is set here, rather than
c     in forcemp.f as in the SGI implementation.
c
      integer*4 ptr_f, ptr_i
      common /SGI_COM/ptr_f,ptr_i,NumProc,NT_FSIZ,NT_ISIZ
c
c     ======================= END SHARED MEMORY =======================
#endif
#ifdef SGI_MP
#define SHARED_MEMORY
#endif
#ifdef CRAY_MP
#define SHARED_MEMORY
c     ========================= SHARED MEMORY =========================
c     ...shared memory common block, necessary since in the cray 
c     implementation, the number of processors is set here, rather than
c     in forcemp.f as in the SGI implementation.
c
      common /SGI_COM/ptr_f,ptr_i,NumProc,NT_FSIZ,NT_ISIZ
c
c     ======================= END SHARED MEMORY =======================
#endif
c
c     -- macro compatibility checking
c
#ifdef SHARED_MEMORY
#define MALLOC_CHECK
#endif
#ifdef MPI
#  ifndef CSPP
#    define MALLOC_CHECK
#  endif
#endif
#ifdef MEM_ALLOC
# ifdef MALLOC_CHECK

 ccccccccccccccccccccccccccccccccccccccccccc

        -DMEM_ALLOC not compatible with
        parallel versions for now, This
        message is intended to produce
        a compiler error.. remove the
        -D*_MP, -DMPI or -DMEM_ALLOC flag 
        from your MACHINE file or use a non-
        parallel one & try again.

 ccccccccccccccccccccccccccccccccccccccccccc
# endif
#endif
      dimension ene(30)
      integer itim0, itim1, itim2
c
c     ---- DECLARE THE MEMORY OF THE MAIN REAL, INTEGER, HOLLERITH
c          AND NONBONDED PAIRLIST ARRAYS  ----
c          These arrays are a very old-fashioned memory-greedy 
c          device that permeates the style of amber - statically 
c          dimensioned per-type arrays are suballocated within 
c          the process (in the locmem subroutine) according to the
c          needs of the current job.
c
#ifdef MEM_ALLOC
C 
C     ---- OMUPD rkw 12/01/95 Pointer declarations for dynamic 
C          memory use
C
      pointer(xptr,  x)
      pointer(ixptr, ix)
      pointer(ihptr, ih)
      pointer(ipptr, ipairs)
      dimension x(1)
      integer   ix(1), ih(1), ipairs(1), cmalloc, p2p_shmalloc
c     common/ptrs/xptr,ixptr,ihptr,ipptr
C OMUPD end
#else
c
c     ---- The purpose of putting the per-type arrays in a
c          common block (which is not used anywhere else) is 
c          to guarantee static and aligned allocation, thus
c          "anchoring" the location of these arrays in all
c          processes when running in parallel since data
c          sharing relies on the shared memory being in the
c          same place in all processes. If the statement is 
c          _not_ there (among the highest-level declarations), 
c          the T3D (MPI) parallel version will crash and SGI 
c          shared memory parallel version will give wrong
c          results; and if the real array X is not first, the 
c          T3D still crashes.
c          
      common/root/x,ix,ih,ipairs
c
      dimension X(MAXREA)
      integer ix(MAXINT)
      integer ih(MAXHOL)
      integer ipairs(MAXPR)
#endif
c
c     ---- HERE BEGIN THE EXECUTABLE STATEMENTS ----
c
#ifdef CTSS
C  Following is CTSS-only call. Not needed Unicos. Collectors item now.
C  CTSS=Cray Time Sharing System; "dropfile" was a restartable image of 
C  the program. 
      CALL DROPFILE(0)
#endif
C
C Initialize the cpu timer. Needed for machines where returned cpu times 
C are relative.
C
#ifdef CSPP
      t0 = my_walltime( 0.0 )
#endif
      CALL TIMIT(3,SKIP,6)
      call wall(itim0)
C
#if defined HP || defined SPP
c
c     --- set up certain underflow operations
c
      on double precision underflow call trapud
      on real underflow call trapu
#endif
c
c     --- Code Configuration Section ---
c
c     Revision A code is designed to port easily to machines of
c     varying wordsize.   We have changed the memory mapping scheme
c     of version 3.0 to improve portabilty.  Instead of one large
c     array there are now three separate arrays.  The following points
c     are notable:
c     1)  Three arrays X, IH, and IX are passed as arguments instead
c         of COMMON.  These arrays are for Reals, Hollerith ints, and
c         Integers, respectively.  The Reals may be 32 or 64 bit, and
c         Integers may be 32 or 64 bit as well.   The Hollerith int
c         and numerical int arrays share the same starting address,
c         by virtue of the equivalence in main.  This is *not* non-portable,
c         do not get excited.  All of the structural arrays in the program
c         are mapped into the three large arrays, and are referenced
c         by the offsets passed in commons MEMLA for reals, MEMLB and
c         MEMLD for ints, and MEMLC for Hollerith data.
c     2)      Two variables are used to describe the packing of the
c         nonbonded pairlist.   IPACK is set to 1 if explicit word-
c         packing routines are called, as on Cray or FPS.  Otherwise
c         IPACK = 0.  NWDVAR is the number of NB pair pointers that is
c         held in a default integer word.   This will usually be 4 for
c         the Cray, 2 or 4 for the FPS, and on a 32 bit machine will
c         be 1 if the entire word is used, or 2 if an integer*2 declaration
c         is used for the pairlist.
c         NATIVE is the number of bits in a default integer on the target
c         machine.  It is usually 32 or 64.
c     3)      Memory Requirements:  The three parameters MAXREA, MAXINT,
c         and MAXDUP are used to control memory use.
c
c         Array         Use       Parameter      Typical Value
c           X       floating pt    MAXREA         ~ 23 * Natom
c          IX       Integers       MAXINT        ~ 150 * Natom
c         (various) dihedral dup   MAXDUP      0 - 1000, data dependent
c
c         The typical values given are only rough estimates.  The Integer
c         memory requirement consists of a "static" requirement, which
c         is topology dependent and does not vary throughout the run,
c         and a variable amount for the nonbonded pairlist pointers.
c         The value of ~150*Natom includes both the static requirement
c         and the pairlist requirement for a "typical" system assuming
c         a full word is used to store a pairlist pointer.   To determine
c         the actual Integer memory requirement, add the static requirement
c         reported at the start of a run to the pairlist requirement.  The
c         pairlist requirement is the total number of nonbonded pairs (this
c         is geometry and cutoff dependent) divided by NWDVAR.  If IPACK
c         is not 0, you should add Natom.  Since the pairlist can grow
c         during a run (and often does) it is a good idea to increase the
c         room for it by ~10%.   The maximum number of nonbonded pairs
c         for the value of MAXINT used will be reported at the top of the
c         output.
c
c         IPACK=0:   MAXINT = Static Int + (NPAIRS/NWDVAR)*1.1
c         IPACK=1:   MAXINT = Static Int + (NPAIRS/NWDVAR)*1.1+NATOM
c
c         All cases: MAXREA = reported static output.
c
c         Dihedrals that have more than one fourier term will have their
c         pointers duplicated for the vectorized dihedral routine.  This
c         is done twice; once for heavy atom dihedrals and once for diheds
c         involving H-atoms.   MAXDUP must be at least as large as the
c         larger number reported in the output.  The actual amount of space
c         allocated is 10*MAXDUP, so it should not be set too large in a
c         tight memory environment.
c
#ifdef MPI
c     =========================== AMBER/MPI ===========================
c
c     Set up parallel execution 
c
       call mpi_init(ierr)
       call mpi_comm_rank(MPI_COMM_WORLD,mytaskid,ierr)
       call mpi_comm_size(MPI_COMM_WORLD,numtasks,ierr)
       world_comm=MPI_COMM_WORLD
       call mpi_comm_group(MPI_COMM_WORLD,world_group,ierr)

#  ifdef MPI_BUFFER_SIZE
       call mpi_buffer_attach(mpibuf, MPI_BUFFER_SIZE*4, ierr)
#  endif
c
c     Make PE 0 the master
c
      master = mytaskid.EQ.0
      if(master)print *,"All processors started"
c
c
c
c     ========================= END AMBER/MPI =========================
#else
c
c     in the single-threaded version, the 1 process is master
c
      master = .true.
#endif
      NRU = 0
      erstop = .false.
c     --- generic packing scheme ---
      nwdvar = 1
      native = 32
#ifdef CRAYFISH
c     --- Cray packing scheme ---
      nwdvar = 4
      native = 64
#endif
#ifdef FPS264
c     --- FPS264 packing scheme ---
      nwdvar = 2
      native = 64
#endif
#ifdef ISTAR2
c     --- Int*2 packing scheme ---
      nwdvar = 2
#endif
      numpk = nwdvar
      nbit = native/numpk
      ONE = 1.0D0
      SMALL = 1.0D-4
c
c     ----- Only the master node (only node when single-process)
c           performs the initial setup and reading/writing -----
c
      if (master) then
c
c        --- get file names ---
c
         call mdfil
C
C        ----- READ THE NECESSARY DATA TO INITIATE THE RUN -----
C
#ifdef MEM_ALLOC
         CALL MDREAD(xptr,ixptr,ihptr,ipptr)
#else
         CALL MDREAD(x,ix,ih,ipairs)
#endif
C
C        ----- EVALUATE SOME CONSTANTS FROM MDREAD SETTINGS -----
C
         NRPT = NRP
         NR = NRPT+NSM*NRAM
         NR3 = 3*NR
         BELLY = IBELLY.GT.0
c
c        --- seed the random number generator ---
c
         call amrset(ig)
c
         if (nbit .lt. 32 .and. nr .gt. 32767) then
              write(6, *) '  Too many atoms for 16 bit pairlist -'
              write(6, *) '    Recompile without ISTAR2'
              call mexit(6, 1)
         endif
c
c        --- this check important because of alloc of L45 & 
c            its use in runmd ---
c
         IF (NTP.GT.0.AND.IABS(NTB).NE.2) GOTO 1000
C
C        ----- READ COORDINATES AND VELOCITIES -----
C
         CALL GETCOR(NR,X(L30),X(L55),X(L40),X(L35),NTX,BOX,IREST,T)
c
c        --- Set up principal (marker) atom list for res based imaging ---
c
         if (ntb .ne. 0) 
     +        call setmrk(natom,nres,ix(i02),x(L30),ix(i01))
c
         IF (INIT.EQ.4.AND.NTX.LT.4) INIT = 3
C
C        ----- SET THE INITIAL VELOCITIES -----
C
         IF (NTX.LE.3) CALL SETVEL(NRP,NR,NTU,NTX,X(L30),
     +             X(L40),X(L35),X(L20),TEMPI,HEAT,DT,INIT,IG,
     +             iscale,scalm)
C
         IF (BELLY) CALL BELLYF(NATOM,IX(I62),X(L40))
c
c        --- remove com motion if required ---
c
         if (imin.eq.0 .and. ntcm.ne.0 .and. ntb.ne.0) then
           vx = 0.d0
           vy = 0.d0
           vz = 0.d0
  
           j = 1
           do i = 1, 3*natom,3
              aamass = x(l70+j-1)
              vx = vx + aamass * x(l40+i-1)
              vy = vy + aamass * x(l40+i)
              vz = vz + aamass * x(l40+i+1)
              j = j + 1
           enddo
 
           vx = vx * tmassinv
           vy = vy * tmassinv
           vz = vz * tmassinv
 
           vel = sqrt(vx*vx + vy*vy + vz*vz)
           if ( vel.gt.0.01d0 ) then
             write (6,'(a,f9.2)')
     .           '  removed COM velocity: ', vel
           endif

           do i = 1, 3*natom, 3
              x(l40+i-1) = x(l40+i-1) - vx
              x(l40+i)   = x(l40+i)   - vy
              x(l40+i+1) = x(l40+i+1) - vz
           enddo
         endif
C
C        --- If we are reading NMR restraints/weight changes, 
C            read them now:
C
         IF (NMROPT.ge.1) THEN
            CALL NMRCAL(X(L30),F,IH(M04),IH(M02),IX(I02),X(L20),ENMR,
     *               DEVDIS,DEVANG,DEVTOR,TEMP0,TAUTP,CUT,NTB,X(LNMR01),
     *               IX(INMR02),X(L95),5,6,RK,TK,PK,CN1,CN2,
     *               AG,BG,CG,NUMBND,NUMANG,NUMPHI,NIMPRP,
     *               NTTYP,NHB,NATOM,NATOM,NTYPES,NRES,RAD,WEL,RADHB,
     *               WELHB,RWELL,ISFTRP,-1,'READ')
C
C           --- Determine how many of the torsional parameters 
C               are impropers
C
            CALL IMPNUM(IX(I46),IX(I56),IX(I48),IX(I58),NPHIH,NPHIA,
     *               0,NPTRA,NIMPRP)
         ENDIF
C
C        --- Set the van der Waals radii (for soft-repulsion form):
C
         CALL NMRRAD(RAD,WEL,CN1,CN2,NTYPES,0,0.0D0)
         CALL DECNVH(ASOL,BSOL,NPHB,RADHB,WELHB)
C
C        --- Call FASTWT, which will tag those bonds which are part 
C            of 3-point water molecules. Constraints will be effected 
C            for these waters using a fast analytic routine -- dap.
C
C        SUBROUTINE FASTWT(IGRAPH    ,NRES      ,IPRES     ,LBRES     ,
C    2          NBONH     ,NBONA     ,NBPER     ,IB        ,JB        ,
C    3          NBRST     ,NARST     ,NPRST     ,NATRST    ,ITORTY    ,
C    4          IBELLY    ,IGRP      ,NCORC     ,ICMPDR    ,IAPER     ,
C    5          IWTNM     ,IOWTNM    ,IHWTNM    ,JFASTW    ,IFSTWT    ,
C    6          IFSTWR    ,IBGWAT    ,IENWAT    ,IORWAT    ,IWATPR    ,
C    7          IOUT)
C
C        SUBROUTINE GETWDS(IGRAPH    ,NRES      ,IPRES     ,LBRES     ,
C    2          NBONH     ,NBONA     ,NBPER     ,IB        ,JB        ,
C    3          IWTNM     ,IOWTNM    ,IHWTNM    ,JFASTW    ,ICB       ,
C    4          REQ       ,WINV      ,RBTARG    ,IBELLY    ,IGRP      ,
C    5          IOUT)
C
C        SUBROUTINE FSWIND(IBGWAT    ,IENWAT    ,JFASTW    ,IORWAT    ,
C    2          IWATPR    ,NRES      ,NTB       ,NTM       ,ISFTRP    ,
C    3          NTYPES    ,IAC       ,ICO       ,IPRES     ,NSOLW     ,
C    4          IBELLY)
C
         CALL FASTWT(IH(M04)   ,NRES      ,IX(I02)   ,IH(M02)   ,
     2       NBONH     ,NBONA     ,0         ,IX(I12)   ,IX(I14)   ,
     3       0         ,0         ,0         ,DUMM      ,DUMM      ,
     4       IBELLY    ,IX(I62)   ,0         ,0         ,DUMM      ,
     5       IWTNM     ,IOWTNM    ,IHWTNM    ,JFASTW    ,IX(IVM01) ,
     6       IX(IVM02) ,IBGWAT    ,IENWAT    ,IORWAT    ,IDUMM     ,
     7       6)
c
         CALL GETWDS(IH(M04)   ,NRES      ,IX(I02)   ,IH(M02)   ,
     2       NBONH     ,NBONA     ,0         ,IX(I12)   ,IX(I14)   ,
     3       IWTNM     ,IOWTNM    ,IHWTNM    ,JFASTW    ,IX(I16)   ,
     4       REQ       ,X(L20)    ,RBTARG    ,IBELLY    ,IX(I62)   ,
     5       6)
C
         CALL FSWIND(IBGWAT    ,IENWAT    ,JFASTW    ,IORWAT    ,
     *       IWATPR    ,NRES      ,NTB       ,0         ,ISFTRP    ,
     *       NTYPES    ,IX(I04)   ,IX(I06)   ,IX(I02)   ,NSOLW     ,
     *       IBELLY)
C
C        --- Set up the solute/solvent pointers:
C
         CALL SOLPNT(NSOLW      ,IBGWAT    ,IPTRES    ,IPTSOL    ,
     *       NATRCM    ,IPTATM     ,IFTRES    ,ISOLVP    ,NATOM     ,
     *       NRES      ,NSOLUT     ,NTT       ,NSPSOL    ,NSPSTR    ,
     *       IX(I02)    ,6)
C
C If we're calculating a ligand grid, call LGSET to do the setup
C
C           TINE LGSET (X         ,IB        ,JB        ,NBOND     ,
C    *       NAMES     ,IOUT      ,NTYPES    ,CN1       ,CN2       ,
C    *       IBLO      ,INB)
C
         IF (LGRON.EQ.1) THEN
            NBOND = NBONH + NBONA
            CALL LGSET (X(L30)    ,IX(I12)   ,IX(I14)   ,NBOND     ,
     *       IH(M04)   ,6         ,NTYPES    ,CN1       ,CN2       ,
     *       IX(I08)   ,IX(I10)   ,66)

            IF (LGREST.GT.0) 
     *      CALL LGRST (66        ,1         ,NTYPES    ,CN1       ,
     *       CN2)

         END IF
C
C --------------------------------------------------------------------------
C If RLMAVF > 0 (using averaged guiding forces), call MKSBST to initialize
C the substructures list, and set all GS0 values initially to 0.0
C Note that currently we don't include any added restraints in the
C topology for determining substructures.
C --------------------------------------------------------------------------
C
      IF (RLMAVF.GT.1.0D-6) THEN
C        TINE MKSBST   (NBONH     ,NBON      ,NBPER     ,NBRST     ,
C    *       IB        ,JB        ,NATRST    ,IDEEP     ,RIMASS    ,
C    *       NATOM)
         CALL MKSBST   (NBONH     ,NBONA     ,0         ,0         ,
     *       IX(I12)   ,IX(I14)   ,IDUMM     ,ISGDEP    ,X(L20)    ,
     *       NATOM)
      END IF
C
C        --- OPEN THE DATA DUMPING FILES AND POSITION IT DEPENDING
C            ON THE TYPE OF RUN -----
C
         CALL OUTOPN
C
C        --- end of master process setup ---
      endif
C
  120 CONTINUE
C
C ----------------------------------------------------------------------
C Main MD loop over NRUN
C ----------------------------------------------------------------------
C
      if(master) then

C -------------------------------------------------------------------------
C ------------------ READ INPUT COORDINATE ARCHIVE ----------------------
C If IRDARC.NE.0, we don't use the standard input coordinates, but rather
C coordinates read from an input coordinate archive. Read these now,
C before we start the MIN/MD loops.
C
      IF (IRDARC.NE.0) THEN
         IACSET = IACSET + ABS(IRDARC)
         WRITE(6,1070) IACSET
C        TINE    REDARC(     X         ,XSCR      ,ISET      ,NATOM    ,
C    *            NSTART    ,NEND      ,IFORM     ,IUNIT     ,IOUT     ,
C    *            ISTAT)
         CALL    REDARC(     X(L30)    ,X(L95)    ,IACSET    ,NATOM    ,
     *            IRD1ST    ,IRDLST    ,IRDARC    ,IURDAR    ,6        ,
     *            ISTATC)
C
C ISTAT=0 if read was successful. If read was unsuccessful (attempt
C to read a set that's not there), we're done reading from archive.
C Jump to full exit.
C
         IF (ISTATC.EQ.1) THEN
            WRITE(6,1071)
            GO TO 250
         END IF
C
C Reset the step and average accumulators for the NMR routines
C
         CALL RESETN
C
C Assign a new set of random velocities for the new coordinate set if
C we're doing MD (IMIN=0).
C
         INIT = 3
         NTX = 1
         IF (IMIN.EQ.0)
     *   CALL SETVEL(NPM,NRP,NR,NTU,NTX,X(L30),X(L40),X(L35),X(L20),
     *               TEMPI,HEAT,DT,INIT,IG,iscale,scalm)
         IF (BELLY) CALL BELLYF(NATOM,IX(I62),X(L40))
      END IF
C ---------------- END READ INPUT COORDINATE ARCHIVE ----------------------
C -------------------------------------------------------------------------
      endif
C
  220 CONTINUE
      SKIP = .FALSE.
      NRU = NRU+1
C
c start Margaret for langevin.

 
        JSEED(1)=mod(IG,4096)
        JSEED(2)=mod(IG*3,4096)
        JSEED(3)=mod(IG*7,4096)
        IRAN = mod(IG*5,4096)
        IF(mod(IRAN,2).EQ.0) THEN
          JSEED(4)=IRAN+1
        ELSE
          JSEED(4)=IRAN
        ENDIF      

        PI =  3.141592653589793
c       water visco = 0.1439 kcal*ps/mol*(AA**3) * 20.455 
        WVISCO = 0.1439 * 20.455 
        BOLTZ = 1.9872d-3        

c high friction limit
c	RRADIUS = 3
c	eta= 6*PI*RRADIUS*WVISCO   
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c low friction limit
c read friction
      open(77,file='rsize.inp')
      read(77,*)nrsize

      if(nrsize.lt.natom.or.nrsize.gt.natom) then
      write(756,*) 'error: nrsize not equal natom, check' 
      endif

        do i=1,nrsize
	read(77,*)RRADIUS		
	xeta(i)= 6*PI*RRADIUS*WVISCO   
        write(80,*)xeta(i),RRADIUS   
	enddo

c read triple.inp

      open(78,file='triple.inp')
      read(78,*)nbeta,xxkchi

        XKCHI=xxkchi
        ICHI=nbeta
      write(80,*)ICHI,XKCHI

      do j=1,nbeta
        read(78,*)IICA,IICB,IINC,IICC,xxchi
        write(80,*)IICA,IICB,IINC,IICC,xxchi
       ICA(j) = iica 
       ICB(j) = iicb 
       INC(j) = iinc
       ICC(j) = iicc
       XCHI(j) = xxchi 
      enddo

c read hb pairs

      ZERO = 0.0d0
      open(79,file='pseudo.inp')
      read(79,*)nhb_pair,Amp
      write(80,*)nhb_pair,Amp
      if(natom.gt.10000) then
      write(755,*) 'error:natom gt 10000,check HB.h'
      endif

      if(nhb_pair.gt.95000) then
      write(755,*) 'error:nhb_pair gt 95000,check HB.h'
      endif

      do j=1,natom
       do k=1,natom
       maphb(j,k) = ZERO
       enddo
      enddo



      do j=1,nhb_pair
       NXI(j) = ZERO
       NXI1(j) = ZERO
       NXI2(j) = ZERO
       NXI3(j) = ZERO
      enddo



c ixi1<ixi2
      do j=1,nhb_pair
        read(79,*)ixi,ixi1,ixi2,ixi3
        write(80,*)ixi,ixi1,ixi2,ixi3

       maphb(ixi1,ixi2) = j
       maphb(ixi2,ixi1) = j

       NXI(j) = ixi
       NXI1(j) = ixi1
       NXI2(j) = ixi2
       NXI3(j) = ixi3

      enddo

!      if(ntb.eq.1.and.iftres.eq.1) then
!      write(755,*) 'error:ntb ne 1, iftres to be 0'
!      endif


c        write(78,*)eta    
         close(77)
         close(78)
         close(79)
         close(80)
c end of Margaret for langevin

c Andrei for double well.

        ZERO = 0.0d0
        open(77,file='douwellwide.inp')
        read(77,*)mpower,kpower,npower
        read(77,*)ngo_pair

        write(75,*)mpower,kpower,npower

        do j=1,ntypes-1
         do k=1,ntypes-1 
         DB0(j,k) = ZERO
         DB1(j,k) = ZERO
         CB1(j,k) = ZERO
         B2(j,k) = ZERO
         h1(j,k) = ZERO
         h2(j,k) = ZERO
         r0(j,k) = ZERO
         r1(j,k) = ZERO 
         r2(j,k) = ZERO
         mapgo(j,k) = ZERO
c         write(1000,*)j, k, mapgo(j,k)
c         write(1001,*)j, k, DB0(j,k)
c         write(1002,*)j, k, DB1(j,k)
         enddo
        enddo

        do j=1,ngo_pair
         read(77,*)jgo,kgo,dbb0,dbb1,cbb1,dbb2,hh1,hh2
         mapgo(jgo,kgo) = 1
         mapgo(kgo,jgo) = 1
         DB0(jgo,kgo) = dbb0
         DB1(jgo,kgo) = dbb1
         CB1(jgo,kgo) = cbb1
         B2(jgo,kgo) = dbb2
         h1(jgo,kgo) = hh1
         h2(jgo,kgo) = hh2
c         write(1003,*)j, jgo, kgo, mapgo(jgo,kgo)
c         write(1004,*)j, jgo, kgo, DB0(jgo,kgo)
c         write(1005,*)j, jgo, kgo, DB1(jgo,kgo)
        enddo

        close(77)

        open(77,file='r0r1r2.inp')
        read(77,*)ngo_pair
        write(75,*) ngo_pair

        do j=1,ngo_pair
         read(77,*) jgo,kgo,rr0,rr1,rr2
         r0(jgo,kgo) = rr0
         r1(jgo,kgo) = rr1
         r2(jgo,kgo) = rr2
         write(75,*)rr0,rr1,rr2
        enddo

        close(77)
        close(75)
c        close(1000)
c        close(1001)
c        close(1002)
c        close(1003)
c        close(1004)
c        close(1005)
c end of Andrei 
      CALL MDBOX

#ifdef MPI
c     =========================== AMBER/MPI ===========================
c
c     NOTE: in the current AMBER/MPI implementation, two means of 
c     running in parallel within sander are supported.  
c     The original implementation had each of the non-master nodes
c     entering a never ending loop which repeatedly called force().
c     This implementation is general and supports most of the input 
c     options within sander.  Later, a more optimal approach was
c     implemented which involves calling runmd() in parallel.
c     Calling runmd() in parallel allows not only the parallelization 
c     of the energy and force evaluation, but also the dynamics 
c     integration steps.
c
c     Previously, which mode used was determined at compile time.
c     In the current code, a the value of mpi_orig (set below and
c     hardcoded) determines which approach is used.
c     Currently this is turned on when minimization (imin .ne. 0)
c     or when using NMR restraints.  Eventually runmin() will be 
c     adapted to be called in parallel, the NMR code updated for
c     parallel operation, and the mpi_orig .eq. true case probably 
c     disabled.  In the meantime, the use of the mpi_orig variable
c     obfuscates the source code in this routine, force.f, qiktip.f,
c     runmd.f, and perhaps others.
c
c     [When running the mpi_orig case, a variable notdone is now
c     set by the master and determines when to exit the force()
c     loop.  When the master has finished calling force, the
c     master changes notdone to 0 and broadcasts the data one more
c     time to signal end of the loop.  force() is modified so that
c     in the mpi_orig case, an initial broadcast is done to receive
c     the value from the master to decide whether to do the work or
c     simply exit (tec3)]
c
c     ...set up initial data and send all needed data to other nodes, 
c     now that the master has it
c
      NRPT = NRP
      NR = NRPT+NSM*NRAM
      NR3 = 3*NR
      BELLY = IBELLY.GT.0
c
      call startup_groups(ierr)
c
      if (NRU .eq. 1) then
#ifdef MEM_ALLOC
        call mpi_bcast(nrealb,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if(.not.master) then
c
c     -- allocate!
c
# ifdef CSPP
           xptr = p2p_shmalloc(%val(nrealb))
c          xptr = cmalloc(nrealb)
# else
           xptr = cmalloc(nrealb)
# endif
           if (xptr.eq.0) then
             write(6,*) ' ** failed to allocate real array, bytes: ',
     .          nrealb
             call mexit(6,1)
           endif
           ixptr = cmalloc(nintb)
           if (ixptr.eq.0) then
             write(6,*) ' ** failed to allocate integer array, bytes: ',
     .           nintb
             call mexit(6,1)
           endif
           ihptr = cmalloc(nholb)
           if (ihptr.eq.0) then
             write(6,*) ' ** failed to allocate holl array, bytes: ',
     .            nholb
             call mexit(6,1)
           endif
           ipptr = cmalloc(npairb)
           if (ipptr.eq.0) then
             write(6,*) ' ** failed to allocate pair array, bytes: ',
     .            npairb
             call mexit(6,1)
           endif
         endif
c
#endif
         call startup(x,ix,ih)

         if (NMROPT.ge.1) then
           CALL NMRCAL(X(L30),F,IH(M04),IH(M02),IX(I02),X(L20),ENMR,
     *               DEVDIS,DEVANG,DEVTOR,TEMP0,TAUTP,CUT,NTB,X(LNMR01),
     *               IX(INMR02),X(L95),5,6,RK,TK,PK,CN1,CN2,
     *               AG,BG,CG,NUMBND,NUMANG,NUMPHI,NIMPRP,
     *               NTTYP,NHB,NATOM,NATOM,NTYPES,NRES,RAD,WEL,RADHB,
     *               WELHB,RWELL,ISFTRP,-1,'MPI ')
         endif

c        if (.not. master) call amrset(ig)
      endif

      if (imin.ne.0 .or. plevel.eq.0) then
        mpi_orig = .true.
        notdone = 1
      else
        mpi_orig = .false.
      end if
      if (mpi_orig .and. .NOT.master) then
c       
c        ...all nodes only do the force calculations (JV)
c
         do while( notdone .eq. 1 )
            call FORCE(x,ix,ih,ipairs,X(L30),X(L35),ENER,VIR)
         enddo
#ifndef TERRA
         call profil2
#endif
         call mexit(0,0)
      endif
c
      if (master) write(6, '(a,i4,a,/)') 
     .      '|  Running AMBER/MPI version on ',numtasks, ' nodes'
c
c     ========================= END AMBER/MPI =========================
#endif
#ifdef SHARED_MEMORY
c     ========================= SHARED MEMORY =========================
c
c     This code is a sanity check to make sure that the desired number
c     of processors does not exceed the number specified in the MACHINE
c     script.  SHARED_MEMORY_MAX_PROCESSORS is a CPP #define used in 
c     a few places in the SHARED_MEMORY code to allocate small scratch 
c     arrays.
c
# ifndef SHARED_MEMORY_MAX_PROCESSORS
#   ifdef SGI_MP
#     define SHARED_MEMORY_MAX_PROCESSORS 36
#   else
#     ifdef CRAY_MP
#       define SHARED_MEMORY_MAX_PROCESSORS 16
#     else
#       ifdef CNX_MP
#         define SHARED_MEMORY_MAX_PROCESSORS 128 
#       endif
#     endif
#   endif
# endif

      numproc = 1
# ifdef SGI_MP
c$    numproc = mp_numthreads()
      if (numproc .gt. SHARED_MEMORY_MAX_PROCESSORS) then
         write (6, 21991) numproc, SHARED_MEMORY_MAX_PROCESSORS
         write (6, 21992)
         write (6, 21993)
c$       call mp_set_numthreads(SHARED_MEMORY_MAX_PROCESSORS)
      endif
# endif
# ifdef CRAY_MP
      numproc = ncpu()
      if (numproc .gt. SHARED_MEMORY_MAX_PROCESSORS) then
         write (6, 21991) numproc, SHARED_MEMORY_MAX_PROCESSORS
         write (6, 21992)
         numproc = SHARED_MEMORY_MAX_PROCESSORS
      endif
# endif
# ifdef CNX_MP
      numproc = num_threads()
      if (numproc .gt. SHARED_MEMORY_MAX_PROCESSORS) then
         write (6, 21991) numproc, SHARED_MEMORY_MAX_PROCESSORS
         write (6, 21992)
         numproc = SHARED_MEMORY_MAX_PROCESSORS
      endif
      write(6, 21994) numproc
# endif


21991 format('ERROR: cannot run with ', i4, ' threads. ',
     .     'Compiled maximum is ', i4)
21992 format('WARNING: resetting number of processors to ',
     .     'maximum number allowed...')
21993 format('NOTE: to run with more threads, modify the ',
     .     'MACHINE script and recompile.')
21994 format('Convex SPP running with ', i4, ' threads. ')
c
      write(6, '(a,i4,a,/)') 
     .   '|  Running shared memory version on ',numproc,' processors'
c

c     ======================= END SHARED MEMORY =======================
#endif
      call second(wtim0)
      call wall(itim1)
C
C ----------------------------------------------------------------------
C Now do the dynamics or minimization.
C ----------------------------------------------------------------------
C
      IF (IMIN.EQ.0) THEN
C
C        --- Dynamics:
C
         CALL RUNMD(x,ix,ih,
#ifdef MEM_ALLOC
     +           ipptr,
#else
     +           ipairs,
#endif
     +           X(L30),X(L20),X(L70),X(L35),X(L40),X(L45),X(L55),
     +           X(L50),X(L95),IX(I70),X(L75),erstop)
c
         if (master) call amflsh(6)
c
         IF (ERSTOP) THEN
            if (master) then
               WRITE(6, *) 'FATAL ERROR'
               CALL MEXIT(6,1)
            else
               CALL MEXIT(0,0)
            endif
         END IF
C
C        --- At the end of each NSTLIM steps, if NTT=0, adjust the temp.
C
         INIT = 4
         IF (NSCM.LE.NSTLIM) NTCM = 1
         IF (NSNB.LE.NSTLIM) NTNB = 1
         IF (NTCM.NE.0) INIT = 3
C
         IF (NTT.EQ.0 .AND.
     *       ABS(TEMP-TEMP0).GT.DTEMP .AND. 
     *       ABS(TEMP).GT.ONE .AND.
     *       ABS(DTEMP).GT.SMALL) THEN
           HEAT = SQRT(TEMP0/TEMP)
           DO 260 I = 1,NR3
                X(I+L40-1) = X(I+L40-1)*HEAT
  260      CONTINUE
           INIT = 3
           if (master)WRITE(6,1030)
         ENDIF
C
C        --- Check time remaining. Exit if max. time exceeded.
C
         CALL TIMIT(1,SKIP,6)
         IF (SKIP) WRITE(6,1020)
C
C        --- If time not exceeded, and specified NRUN runs 
C            not completed, go do another run:
C
         IF (.NOT. SKIP .AND. NRU .LT.NRUN) GO TO 220
C
      ELSE
C
C        --- Minimization:
C
c        SUBROUTINE RUNMIN(xx,ix,ih,ipairs,
c    .                     X,FG,W,IGRAPH,LBRES,IPRES,IB,JB,CONP,
c    .                     WINV,IGRP,SKIPS,NSP,ERSTOP,CONVGD,ENE)
         CALL RUNMIN(x,ix,ih,
#ifdef MEM_ALLOC
     +               ipptr,
#else
     +               ipairs,
#endif
     +               X(L30),X(L35),X(L40),ih(m04),ih(m02),
     +               IX(I02),IX(I12),IX(I14),X(L50),X(L20),IX(I62),
     +               X(L95),IX(I70),ERSTOP,CONVGD,ene)
c
         if (master) then
C           --- Write the restart file:
C
            CALL MINRIT(X(L30))
            CALL TIMIT(1,SKIP,6)
         endif
C
      END IF
C
C At this point, we're essentially done in the standard (one coordinate
C set input case). However, if we're reading multiple sets of coordinates
C from an input archive file (IRDARC.NE.0), we need to read another
C and start again. Return to step 120, where we'll read another coordinate
C set and then do the MIN/MD stuff again.
C
      IF (IRDARC.NE.0) GO TO 120
c
c     -- calc time spent running vs setup
c
  250 call second(wtim1)
      call wall(itim2)
      timsts(NUMSTS) = wtim1 - wtim0

#ifdef MPI
c     =========================== AMBER/MPI ===========================
c
c     Set and broadcast notdone in mpi_orig case to inform 
c     other nodes that we are finished calling force(). (tec3)
c
      if ( mpi_orig ) then
         notdone = 0
         call MPI_BCAST(notdone,1,MPI_INTEGER,0,
     .        MPI_COMM_WORLD,ierr)
      endif

#ifndef TERRA
      call profil2 
#endif
c
c     JV Need to close to get last output on SP1 (maybe due to AFS)
c     JV This could be an error in regular code!!!
c
      if(master) then 
        if(NTWX.GT.0) close(12)
        if(NTWV.GT.0) close(13)
        if(NTWE.GT.0) close(15)
      else
        call mexit(0,0)
      endif
c     ========================= END AMBER/MPI =========================
#endif
C
C When run is over, call profil to write timings:
C
      CALL PROFIL
      write(6, '(/''|'',5x,''Setup wallclock   '',i9,'' seconds'')') 
     .        itim1-itim0
      write(6, '(''|'',5x,''Nonsetup wallclock'',i9,'' seconds'')') 
     .        itim2-itim1
      CALL MEXIT(6, 0)
c
 1000 WRITE(6,1010)
 1010 FORMAT(/ /,'INPUT NTP/NTB INCONSISTENT')
 1020 FORMAT(/ /5X,'CPU TIME LIMIT EXCEEDED')
 1030 FORMAT(/5X,'VELOCITIES HAVE BEEN RESCALED',/)
c1050 FORMAT(/10X,'ENTERING INTO DYNAMICS',/)
 1070 FORMAT(/,78('#'),/,T15
     *       'Reading set ',I5,' from input coordinate archive',
     *       /,78('#'),/)
 1071 FORMAT(T23,'End of coordinate archive detected',
     *       /,78('#'),/)
      CALL MEXIT(6, 1)
      END
