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
c=======================================================================
c
c     This is the main "driver" routine for carrying out molecular
c     dynamics.
c
c=======================================================================
#ifdef MEM_ALLOC
      SUBROUTINE RUNMD(xx,ix,ih,ipptr,X,WINV,amass,F,
     $                 V,XR,XC,CONP,SKIP,NSP,TMA,erstop)
# ifdef DPREC
      implicit double precision (a-h,o-z)
# endif
      pointer(ipptr, ipairs)
      integer   ipairs(1)
#else
      SUBROUTINE RUNMD(xx,ix,ih,ipairs,X,WINV,amass,F,
     $                 V,XR,XC,CONP,SKIP,NSP,TMA,erstop)
# ifdef DPREC
      implicit double precision (a-h,o-z)
# endif
      integer   ipairs(*)
#endif
      dimension xx(*)
      integer   ix(*), ih(*)
      double precision td
      LOGICAL SKIP,BELLY,LOUT,LOUTFM,erstop,vlim,peac,RECENT
#ifdef MPI
c     =========================== AMBER/MPI ===========================
c
c     NOTE: this routine contains MPI functionality to update the
c     positions and velocities for all the atoms on a given node.
c     After the positions are updated, communication (where necessary)
c     is performed to update the coordinates and velocities on each
c     processor.  Also note that wrappers to all the I/O routines are
c     present to prevent all nodes except the master from executing
c     I/O.  In addition, extra timing variables and calls are defined 
c     to allow profiling of this routine.
c
#  include "parallel.h"
#  include "mpif.h"
      dimension enewt(2),ekcmtt(3)
#  ifdef T3D
      common/ekcmtt/enewt,ekcmtt
#  endif
c     ========================= END AMBER/MPI =========================
#endif
c
c     ...for EWALD the only real changes in this routine
c     are that subroutine wrap_molecules() is called instead of SHIAG(),
c     volume is calculated using a special routine get_current_volume(),
c     during pressure coupling, and a set of routines is called
c     instead of MDBOX and PSCALE to update the box/unit cell and 
c     perform pressure scaling.
c
c     Special notes: 
c
c     (1) NREN is updated to 43 intead of 42 to store an extra energy
c     value in the EWALD case which represents the accuracy...
c     (2) the MDBOX/PSCALE equivalents in EWALD currently may not be
c     correct is NTP > 1 (i.e. non-isotropic scaling) with non-
c     orthogonal boxes!
c
#include "files.h"
#include "iewald.h"
#include "md.h"
#include "box.h"
#include "nmr.h"
#include "memory.h"
#include "other.h"
#include "extra.h"
#include "parms.h"
#include "lgcom.h"
#include "avfrc.h"
C Margaret:Langevin
#include "LAN.h"
#include "HB.h"
#include "CHI.h"
#include "CROWD.h"

      DOUBLE PRECISION NOISE(25001) 
c
      DIMENSION X(*),WINV(*),amass(*),F(*),V(*),XR(*),XC(*),CONP(*)
      DIMENSION ENERT(51),ENERT2(51),ENER(51),VIR(4),EKCMT(4)
      DIMENSION PRES(4),RMU(3),FAC(3),XCM(3),VCM(3),ACM(3),OCM(3)
      DIMENSION SKIP(*),TMA(*),NSP(*),ENEW(2)
      DIMENSION IDUMAR(4)
C
      EQUIVALENCE (SCALP,ENER(5)),(SCALS,ENER(6)),(VOL,ENER(10))
      EQUIVALENCE (PRES(1),ENER(11)),(EKCMT(1),ENER(15))
      EQUIVALENCE (VIR(1),ENER(19))
C
      SAVE SMALL,EKKEEP
      DATA SMALL/1.0D-7/
c Margaret
      REAL XXX,xkinetic

c  Runmd operates in kcal/mol units for energy, amu for masses,
c     and angstoms for distances.  To convert the input time parameters 
c     from picoseconds to internal units, multiply by 20.455 
c     (which is 10.0*sqrt(4.184)).
c   
C
C     ----- INITIALIZE SOME VARIABLES -----
C
      vlim = vlimit.gt.small
      peac = tauv0.gt.small .and. tauv.gt.small
      tauv = 20.455*tauv
      td = t
      intcm = max(nstlim/10,ntpr)
      if (intcm .lt. 100) intcm = 100
      NTCMT = 0
      IZERO = 0
      IONE  = 1
      ITWO  = 2
      IFOUR = 4
      BELLY = IBELLY.GT.0
      LOUT = .TRUE.
      LOUTFM = IOUTFM.LE.0
      NRPT = NRP
      NRPT3 = 3*NRPT
      NR = NRPT+NSM*NRAM
      NR3 = 3*NR
#ifdef MPI
c     =========================== AMBER/MPI ===========================
c
c   ---- divide atoms up among the processors, always splitting on 
c        residue boundaries:
c
      call setpar(ix(I02), nsp, nspm)
      if ( mpi_orig ) then
         istart = 1
         iend = natom
      else
         istart = iparpt(mytaskid) + 1
         iend = iparpt(mytaskid+1)
      endif
      istart3 = 3*istart -2
      iend3 = 3*iend
c
c     ========================= END AMBER/MPI =========================
#endif
c
c If NTWPRT.NE.0, only print the solute atoms in the coordinate/vel. archives
C If NTWPR0 > 0, print only atoms >= NTWPR0.
c
      NRX  = NR3
      NRX0 = 1
      IF (NTWPRT.LT.0 .AND. NSPSTR.GT.0) NRX = NSPSTR*3
      IF (NTWPRT.LT.0 .AND. NSPSTR.LE.0) NRX = NATOM*3
      IF (NTWPRT.GT.0) NRX = NTWPRT*3
      IF (NTWPR0.GT.0) NRX0 = 3*(NTWPR0-1)+1
C
C     ----- CLEANUP THE VELOCITY IF BELLY RUN -----
C
      IF(BELLY) CALL BELLYF(NR,IX(I62),V)
C
c=======================================================================
c
C Determine system degrees of freedom (for T scaling, reporting)
C
C Call DEGCNT to get the actual number of degrees of freedom for the
C solute and solvent. This call is new to AMBER4 and for the first time
C returns the correct numbers for belly simulations and simulations
C with separate solute/solvent scaling -- dap
C
C "IDUMAR" is dummy array. Used since this routine also used w/ GIBBS.
C
      CALL DEGCNT(IBELLY,NR,IX(I62),NSOLUT,NBONH,NBONA,0,
     *            IX(I12),IX(I14),IX(I18),IX(I20),IDUMAR,
     *            IDUMAR,NTC,IDUMAR,0,0,0,
     *            IDUMAR,IBELSV,RNDFP,RNDFS)
C
C RNDFP = # degrees of freedom for solute
C RNDFS = # degrees of freedom for solvent
C RNDF = total number of degrees of freedom.
C    modify RNDFP to reflect NDFMIN (set in user-input)
C
      RNDFP = RNDFP - FLOAT(NDFMIN)
      RNDF = RNDFP+RNDFS
C
      NSOLV = NRP - NSOLUT
      IF (IBELLY.GT.0) NSOLV = IBELSV
C
C End of degrees of freedom stuff
C
c=======================================================================
c
      ONET = 1.d0/3.d0
C
C     (PCONV IS THE FACTOR TO CONVERT THE PRESSURE KCAL/MOLE TO BAR)
C

      BOLTZ2 = 8.31441d-3 * 0.5d0
      PCONV = 1.6604345d+04
c
c  convert to "kcal" units
c
      BOLTZ2 = BOLTZ2/4.184d0
      DTX = DT*20.455d+00
      dtxinv = 1.0d0 / dtx
      DT5 = DTX * 0.5d0
      PCONV = PCONV*4.184d0
c
      FAC(1) = BOLTZ2*RNDF
      FAC(2) = BOLTZ2*RNDFP
      IF(RNDFP.LT.0.1d0) FAC(2) = 1.d-6
      FAC(3) = BOLTZ2*RNDFS
      IF(RNDFS.LT.0.1d0) FAC(3) = 1.d-6
      FACTT = RNDF/(RNDF+FLOAT(NDFMIN))
      ekin0  = fac(1)*temp0
      EKINP0 = FAC(2)*TEMP0
      EKINS0 = FAC(3)*TEMP0
      DTMP0  = FAC(1)*DTEMP
C
C Set tauts = tautp, if not already set.
C
      if (tautp.lt.1.0d-7) tautp = 1.0d-7
      if (tauts.lt.1.0d-7) tauts = tautp
C
      IF (NTT.NE.0) DTTP = DT/TAUTP
      IF (NTT.NE.0) DTTS = DT/TAUTS
      IF (NTP.GT.0) DTCP = COMP * 1.0d-06 * DT / TAUP
C
      NREN = 42
c
c     EWALD: NREN=43 is to put use the first unused space in
c     the ener array to hold the accuracy (or quick error 
c     estimate) of the current ewald run.
c     NOTE: if nren changes, or you want to change the place this
c     energy is stored, you will need to change the ene(21) reference
c     in force.f, and the ener(43) everywhere...
c
      NREN = 43
c
      NREK = 4
      NREP = 15
C
      NSTEP = 0
      ISTEP = 0
      FIT = 0.d0
      FITI = 0.d0
      FIT2 = 0.d0
C
      DO 30 I = 1,NREN
        ENER(I) = 0.0d0
        ENERT(I) = 0.0d0
        ENERT2(I) = 0.0d0
   30 CONTINUE
C
      ENER(5) = 1.d0
      ENER(6) = 1.d0
      DO 40 M = 1,3
        ENER(M+6) = BOX(M)
   40 CONTINUE
C
      NITP = 0
      NITS = 0
      if (iscale.gt.0) then
        do 50 im=1,iscale
          x(nr3+im) = bscale(im)
   50   continue
      end if
c
c=======================================================================
C
C     ----- MAKE A FIRST DYNAMICS STEP -----
C
c=======================================================================
c
      GOTO (60,60,110,200),INIT
c
c--------------
c   init = 1,2:
c--------------
c
   60 CONTINUE
      if (master) then
        write(6,*) ' --- attempt to execute dead code: exiting'
        call mexit(6,1)
      else
        call mexit(0,0)
      endif
c
c--------------
c  init = 3:
c--------------
c
  110 IF (NTCM.ne.0) then
        CALL CENMAS(NRP,IZERO,X,V,tmass,tmassinv,IZERO,amass,EKCM,
     +             XCM,VCM,ACM,EKROT,OCM,IFOUR,LOUT)
        CALL STOPCM(NR,X,V,XCM,VCM,OCM)
        NTCMT = NTCM
        NTCM = 0
      end if
      CALL CENMAS(NRP,IZERO,X,V,tmass,tmassinv,IZERO,amass,EKCM,
     +             XCM,VCM,ACM,EKROT,OCM,IFOUR,LOUT)
C
C     ----- EXTRA ENERGY EVALUATION IF CENTER OF MASS MOTION IS
C           REMOVED -----
C
  120 CONTINUE
      if (NTB.ne.0) then
        if (iewald .eq. 1) then
           call wrap_molecules(nspm,nsp,x)
        else
C
C          -- Recenter box on this SHIAG call unless ITRSLU = 2
C
           RECENT = .false.
           IF (NTR.EQ.0 .AND. IBELLY.EQ.0 .AND. ITRSLU.NE.2) 
     *        RECENT = .true.
           CALL SHIAG(NSPM,NSPSTR,NSP,X,NR,RECENT)
        endif
      else
C
C routine POPBCK makes sure a simulation run in vacuuo and without a belly
C will be recentered on {0,0,0} if it gets too close to the boundaries
C that would cause overflows in coordinate writes
C
         IF (NTR.EQ.0 .AND. IBELLY.EQ.0 .AND. ITRSLU.NE.2) THEN
            CALL POPBCK(X,F,NATOM,100.0D0,6)
         END IF
C
      endif
      if (NTP.gt.0) then
        EKCMT(1) = 0.d0
        EKCMT(2) = 0.d0
        EKCMT(3) = 0.d0
        DO 130 I3 = 1,NR3
          XR(I3) = X(I3)
  130   CONTINUE
C
C     ----- CALCULATE THE CENTER OF MASS ENERGY AND THE COORDINATES
C           OF THE SUB-MOLECULES WITH RESPECT TO ITS OWN CENTER OF
C           MASS -----
C
#ifdef MPI
        CALL EKCMR(NSPM,NSP,TMA,EKCMT,XR,V,amass,1,natom)
#else
        CALL EKCMR(NSPM,NSP,TMA,EKCMT,XR,V,amass)
#endif
      end if
C
C     ----- CALCULATE THE FORCE -----
C
      instep = nstep
#ifdef MEM_ALLOC
      CALL FORCE(xx,ix,ih,ipptr,X,F,ENER(23),VIR)
#else
      CALL FORCE(xx,ix,ih,ipairs,X,F,ENER(23),VIR)
#endif

c Margaret
        ENER(23)=ENER(23)+ECHI
C
C This FORCE call does not count as a "step". CALL NMRDCP to decrement
C local NMR step counter
C
      CALL NMRDCP
C
C Reset quantities depending on TEMP0 and TAUTP (which may have been 
C changed by MODWT during FORCE call).
C
      EKINP0 = FAC(2) * TEMP0
      EKINS0 = FAC(3) * TEMP0
      EKIN0 = FAC(1) * TEMP0
      IF (NTT.NE.0) DTTP = DT / TAUTP
      DTMP0 = FAC(1) * DTEMP
C
      if (NTP.gt.0) then
        VOL = BOX(1) * BOX(2) * BOX(3)
        if (ifbox.eq.2) VOL = VOL * 0.5d0
        if (iewald .eq. 1) call get_current_volume(vol)
        ENER(42) = tmass / (0.602204d0*VOL)
        EKCMT(4) = 0.d0
        VIR(4) = 0.d0
        PRES(4) = 0.d0
        DO 140 M = 1,3
          EKCMT(M) = EKCMT(M) * 0.5d0
          EKCMT(4) = EKCMT(4) + EKCMT(M)
          VIR(4) = VIR(4) + VIR(M)
          PRES(M) = (PCONV+PCONV) * (EKCMT(M)-VIR(M)) / VOL
          PRES(4) = PRES(4) + PRES(M)
  140   CONTINUE
        PRES(4) = PRES(4) / 3.d0
      end if
C
      NTNB = 0
      I3 = 0
      TEMPSU = 0.0d0
      TEMPSV = 0.0d0
      DO 160 J = 1,NRP
        WINF = WINV(J) * DT5
        aamass = amass(j)
        DO 150 M = 1,3
          I3 = I3+1
          RTERM = V(I3)*V(I3) * aamass
          IF (J.LE.NSOLUT) THEN
            TEMPSU = TEMPSU + RTERM
          ELSE
            TEMPSV = TEMPSV + RTERM
          END IF
          V(I3) = V(I3) - F(I3) * WINF
          if (vlim) v(i3) = sign(min(abs(v(i3)),vlimit),v(i3))
  150   CONTINUE
  160 CONTINUE
      if (iscale.gt.0) then
        do 180 im=1,iscale
            v(nr3+im) = v(nr3+im) - f(nr3+im) * dt5 / scalm
            tempsu = tempsu + scalm * v(nr3+im)*v(nr3+im)
  180   continue
      end if
      ENER(3) = TEMPSU * 0.5d0
      ENER(4) = TEMPSV * 0.5d0
      ENER(2) = ENER(3)+ENER(4)
      ENER(1) = ENER(2)+ENER(23)
      NSHAK = 1
      IF(NTT.EQ.0) GOTO 200
      EKGP = ENER(3)
      EKGS = ENER(4)
      IF(EKGP.LT.1.d-4) EKGP = EKINP0
      IF(EKGS.LT.1.d-4) EKGS = EKINS0
C
c------------
c init = 4:
c------------
c
  200 EOLD3 = 0.0d0
      EOLD4 = 0.0d0
      I3 = 0
      DO 220 J = 1,NRP
        aamass = amass(j)
        DO 210 M = 1,3
          I3 = I3+1
          RTERM = V(I3)*V(I3) * aamass
          IF (J.LE.NSOLUT) THEN
            EOLD3 = EOLD3 + RTERM
          ELSE
            EOLD4 = EOLD4 + RTERM
          END IF
  210   CONTINUE
  220 CONTINUE
      if (iscale.gt.0) then
          do 240 im=1,iscale
            eold3 = eold3 + scalm*v(nr3+im)*v(nr3+im)
  240     continue
      end if
      EOLD3 = EOLD3 * 0.5d0
      EOLD4 = EOLD4 * 0.5d0
      EKKEEP = EOLD3 + EOLD4
      if (INIT.eq.4) then
        if (NTT.ne.0) then
          EKGP = EOLD3
          EKGS = EOLD4
        end if
        NSHAK = 0
      else
C
c-------------------------------------------------------------------
C           PRINT THE INITIAL ENERGIES AND TEMPERATURES 
c-------------------------------------------------------------------
C
        IF (NSTEP.LE.0 .and. master) THEN
          call opinfo(7)
          CALL PRNTMD(NSTEP,NITP,NITS,T,ENER,FAC,bscale,iscale,7)
          IF (NMROPT.GT.0) THEN
            CALL NMRPTX(6)
            CALL NMRPTX(7)
          END IF
          call closc(7)
        END IF
        IF (NSTLIM.EQ.0) RETURN
      end if
c
c=======================================================================
c
C     ----- MAIN LOOP FOR PERFORMING THE DYNAMICS STEP -----
c
c=======================================================================
C
C Note that every 500 steps we recenter the box on the solute.
C
C (If ITRSLU < 0, then translate waters back into box and
C recenter box every ABS(ITRSLU) steps).
C
  260 ICKSTT = 500
      IF (ITRSLU.LT.0) ICKSTT = ABS(ITRSLU)
C
      IF (NTB.NE.0) then
        if (iewald .eq. 1) then
           call wrap_molecules(nspm,nsp,x)
        else
C
C Recenter box on this SHIAG call every 500 steps unless ITRSLU = 2
C
           RECENT = .false.
           IF (MOD(NSTEP,ICKSTT).EQ.0 .AND. NTR.EQ.0 .AND. 
     *         IBELLY.EQ.0 .AND.
     *         ITRSLU.NE.2) RECENT = .true.
           CALL SHIAG(NSPM,NSPSTR,NSP,X,NR,RECENT)
C
        endif
      ELSE
C
C routine POPBCK makes sure a simulation run in vacuuo and without a belly
C will be recentered on {0,0,0} if it gets too close to the boundaries
C that would cause overflows in coordinate writes
C
         IF (MOD(NSTEP,ICKSTT).EQ.0 .AND. NTR.EQ.0 .AND. IBELLY.EQ.0 
     *   .AND. ITRSLU.NE.2) THEN
            CALL POPBCK(X,F,NATOM,100.0D0,6)
         END IF

      END IF
C
C Update XR(I) (which is the same as XRC() or XX(L45)). XR is set
C to the distance from the COM in routine EKCMR.
C
      if (NTP.gt.0) then
        EKCMT(1) = 0.d0
        EKCMT(2) = 0.d0
        EKCMT(3) = 0.d0
        DO 270 I3 = 1,NR3
          XR(I3) = X(I3)
  270   CONTINUE
C
C     ----- CALCULATE THE CENTER OF MASS ENERGY AND THE COORDINATES
C           OF THE SUB-MOLECULES WITH RESPECT TO ITS OWN CENTER OF
C           MASS -----
C
#ifdef MPI
        if (plevel .eq. 1) then
           CALL EKCMR(NSPM,NSP,TMA,EKCMT,XR,V,amass,1,natom)
c           write(6,*) 'plevel is 1, ekcmt: ',mytaskid,ekcmt(1)
        else
           CALL EKCMR(NSPM,NSP,TMA,EKCMT,XR,V,amass,istart,iend)
           call mpi_allreduce(ekcmt,ekcmtt,3,
     +          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
           ekcmt(1) = ekcmtt(1)
           ekcmt(2) = ekcmtt(2)
           ekcmt(3) = ekcmtt(3)
c           write(6,*) 'plevel is 2, ekcmt: ',mytaskid,ekcmt(1)
        endif
c     ========================= END AMBER/MPI =========================
#else
        CALL EKCMR(NSPM,NSP,TMA,EKCMT,XR,V,amass)
#endif
      end if
C
C     ----- CALCULATE THE FORCE -----
C
      print = mod(nstep,ntpr).eq.0
      instep = nstep
#ifdef MEM_ALLOC
      CALL  FORCE(xx,ix,ih,ipptr,X,F,ENER(23),VIR)
#else
      CALL FORCE(xx,ix,ih,ipairs,X,F,ENER(23),VIR)
#endif


c Margaret
        ENER(23)=ENER(23)+ECHI


#ifdef MPI
      call second(timeg0)
#endif

C If averaged forces are being calculated and added in, do this
C now (since total normal forces were just determined in FORCE)
C Routine GUIDFC determines guiding forces, GS.
C Then calculate the factor ChiE, which is used to scale the velocities.
C Finally, AVEFUD adds the guiding forces into the regular forces.
C
      CHIE = 1.0D0
      IF (ABS(RLMAVF).GT.SMALL) THEN
         NATUSE = NATOM
         IF (ILTAVF.GT.0) NATUSE = ILTAVF
         CALL GUIDFC(FAF,DT,XX(L20),NATUSE)
         CALL CHISCL(CHIE,EKKEEP,F,V,DT,GS,RLMAVF,XX(L20),NATUSE)
         CALL AVEFUD(F,GS,RLMAVF,NATUSE)
      END IF

C
C Reset quantities depending on TEMP0 and TAUTP (which may have been 
C changed by MODWT during FORCE call).
C
      EKINP0 = FAC(2)*TEMP0
      EKINS0 = FAC(3)*TEMP0
      EKIN0 = FAC(1)*TEMP0
      IF (NTT.NE.0) DTTP = DT/TAUTP
      DTMP0 = FAC(1)*DTEMP
C
c --set the PEACS parameters at initial values:
c
      if (peac) then
        if (nstep.le.1) then
          vmin = ener(23)
          if (vzero.eq.0.0d0) vzero = ener(23)
          if (master)
     +       write(6,*) 'PEACS: initial vzero,vmin: ',vzero,vmin
        end if
C
c -- update the PEACS parameters
c
        vmin = min(vmin,ener(23))
        vzero = vzero*exp(-dt/tauv) + vmin*(1. - exp(-dt/tauv))
c           ********need a broadcast here???************
      end if
c
c  Pressure coupling:
c
      if (NTP.gt.0) then
        VOL = BOX(1)*BOX(2)*BOX(3)
        IF (ifbox.eq.2) VOL = VOL*0.5d0
        if (iewald .eq. 1) call get_current_volume(vol)
        ENER(42) = tmass / (0.602204d0*VOL)
        EKCMT(4) = 0.d0
        VIR(4) = 0.d0
        PRES(4) = 0.d0
        DO 280 M = 1,3
          EKCMT(M) = EKCMT(M)*0.5d0
          EKCMT(4) = EKCMT(4)+EKCMT(M)
          VIR(4) = VIR(4)+VIR(M)
          PRES(M) = (PCONV+PCONV)*(EKCMT(M)-VIR(M))/VOL
          PRES(4) = PRES(4)+PRES(M)
  280   CONTINUE
        PRES(4) = PRES(4)/3.d0
      end if
C
c=======================================================================
c
C Temperature coupling:
c
C Several new schemes have been added. HOWEVER, the user should be
C aware that some of these are rather ad-hoc and really only made
C available for md runs where one is not concerned with generating
C a thermodynamically relevant ensemble (e.g. NMR refinement). If 
C one is concerned with the thermodynamics of the system, one is
C best-advised to stick with the tried-and-true(?) options.
C                                             D. Pearlman
C
c-------------------------------------------------------------------
C NTT < 0:
c-------------------------------------------------------------------
C Assign new random velocities, if necessary (NTT < 0):
C Only reassign if outside of DTEMP range of target or every
C ABS(NTT) steps, if NTT.NE.-1. All velocities are reassigned.
C
      IF (NTT.LT.0) THEN
          NSEL = 0
          IF (ABS(EKIN0-(EKGP+EKGS)).GT.DTMP0 .OR.
     *        (NTT.NE.-1 .AND. MOD(NSTEP,ABS(NTT)).EQ.0)) THEN
             if (master) then
                CALL NEWVEL(NTCMT,NRP,IZERO,X,V,tmass,tmassinv,
     *                  amass,WINV,IX(I62),
     *                  IBELLY,EKCM,XCM,VCM,ACM,EKROT,OCM,IFOUR,NR,
     *                  DT5,ENER,F,EKGP,EKGS,EKINP0,EKINS0,EOLD3,EOLD4,
     *                  NSHAK,NTU,TEMP0*FACTT,IG,NSOLUT,NTT,IH(M10),
     *                  NSEL)
#ifdef MPI
c
c               ...this is a hack to fix up velocity reassignment
c               in the absence of a parallelized newvel.  This 
c               results in communication of all the velocities
c               from the master.  There is probably a better way to
c               do this, but I have not found it yet...  [The 
c               communication will be more efficient if it is done
c               as a scatter, however the T3D version does not have
c               this builtin yet...]
c
                 call mpi_bcast(v, 3*natom, MPI_DOUBLE_PRECISION,
     .                          0, MPI_COMM_WORLD, ierr)
#endif
              endif
          ENDIF
          SCALP = 1.0d0
          SCALS = SCALP
C
c-------------------------------------------------------------------
C Constant Energy (NTT=0):
c-------------------------------------------------------------------
C
      ELSE IF (NTT.EQ.0) THEN
          GO TO 290
C
c-------------------------------------------------------------------
C Constant Temperature; Berendsen Algorithm (NTT=1):
C Single scaling factor:
c-------------------------------------------------------------------
C
      ELSE IF (NTT.EQ.1) THEN
          SCALP =  SQRT(1.0d0+DTTP*(EKIN0/(EKGP+EKGS)-1.0d0))
          SCALS = SCALP
C
c-------------------------------------------------------------------
C Constant Temp; Berendsen Algorithm. Only consider solute energy (NTT=2):
C Single scaling factor:
C Could result in solvent atoms having very high temp.
c-------------------------------------------------------------------
C
      ELSE IF (NTT.EQ.2) THEN
          SCALP =  SQRT(1.0d0+DTTP*(EKINP0/EKGP-1.0d0))
          SCALS = SCALP
C
c-------------------------------------------------------------------
C Constant Temperature; Berendsen Algorithm, but only scale when
C temperature falls outside of DTEMP region (NTT=3)
C Single scaling factor:
c-------------------------------------------------------------------
C
      ELSE IF (NTT.EQ.3) THEN
          IF (ABS(EKIN0-(EKGP+EKGS)).GT.DTMP0) THEN
              SCALP =  SQRT(1.0d0+DTTP*(EKIN0/(EKGP+EKGS)-1.0d0))
          ELSE
              SCALP = 1.0d0
          END IF
          SCALS = SCALP
C
c-------------------------------------------------------------------
C When temperature falls outside of DTEMP region, do one quick scale
C back to target temperature (NTT=4):
c-------------------------------------------------------------------
C
      ELSE IF (NTT.EQ.4) THEN
          IF (ABS(EKIN0-(EKGP+EKGS)).GT.DTMP0) THEN
              SCALP = SQRT(EKIN0/(EKGP+EKGS))
          ELSE
              SCALP = 1.0d0
          END IF
          SCALS = SCALP
C
c-------------------------------------------------------------------
C Berendsen algorithm (NTT=5):
C Separate solute/solvent scaling factors:
c-------------------------------------------------------------------
C
      ELSE IF (NTT.EQ.5) THEN
          SCALP = 1.0d0
          IF (EKGP.GT.SMALL) 
     .       SCALP = SQRT(1.0d0+DTTP*(EKINP0/(EKGP)-1.0d0))
          IF (NSOLV.GT.0) 
     .       SCALS = SQRT(1.0d0+DTTS*(EKINS0/(EKGS)-1.0d0))
C
c-------------------------------------------------------------------
C Unrecognized Temp. option; stop:
c-------------------------------------------------------------------
C
      ELSE
          WRITE(6,570) NTT
          call mexit(6, 1)
      END IF
C
c=======================================================================
c
c  ------ready for leap-frog trajectory step:
c
c=======================================================================

cc skip MD.
cc Margaret: Langevin dynamics starts here.
cc with large viscosity limit.
cc a. high friction region
cc x=x+dt*(F+G)/eta
cc b. low friction region
cc
cc eta=6*pi*radius*water
cc water=0.1439 kcal*ps/mol*AA**3 * 20.455 
cc If NR3>15000 need to reset the array of NOISE

	if (NR3.gt.25000) then
	open(76,file='alert')
	write(76,*) 'reset size of NOISE NR3= ',NR3
	end if 

 
        CALL DLARNV(3,JSEED,25000,NOISE)   

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc 	HIGH FRICTION REGION
cc     fix the first atom	
cc origin   DO I3  = 1,NR3+iscale
cc        DO I3  = 4,NR3+iscale
cc     end of fix
cc   nofix
c          DO I3  = 1,NR3+iscale
c	F(I3) = F(I3)*dtx/eta
c        tmpp= 2*boltz*temp0*dtx
cc        write(75,*)JSEED(1),JSEED(2),JSEED(3),JSEED(4)
cc        CALL DLARNV(3,JSEED,1,NOISE)   
cc        CALL DLARUV(JSEED,1,NOISE)   
c	gammaa = NOISE(I3)*sqrt(tmpp/eta)
c        x(I3)=x(I3)+F(I3)+gammaa 
c       enddo                       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c  	low friction region
c 	The verlet-langevin is taken from Chinlin Guo's JCP paper
c	xtemp=eta*dtx/2
c       x(t+dt)-x(t)=(1-xtemp/m)/(1+xtemp/m) {x(t)-x(t-dt)} + 
c	  		dt^2/(m+xtemp)*{F+G}
c	F=-dU(r)/dr
c	G= random*sqrt(2eta*KT/dtx}	

        xkinetic = 0.0

	DO I = 1, NR

	eta=xeta(I)
        tmpp= 2*eta*boltz*temp0/dtx
	xmass=amass(I)
	xtemp=eta*dtx/2

	DO J = 1, 3

	I3=(I-1)*3+J
	gammaa = NOISE(I3)*sqrt(tmpp)

	if(nstep.gt.1) then
	XXX=XDEL(I3)
	else
c 	initial first step
	XXX=0.0
	endif

	temp=((1-xtemp/xmass)*XXX)/(1+xtemp/xmass)
	temp=temp+((dtx**2)/(xmass+xtemp))*(F(I3)+gammaa)
	XM=X(I3)
c 	update coordinate and velocity
	X(I3)=X(I3)+temp
	V(I3)=(X(I3)-XM+XXX)/(2*dtx)
c	store the old value to the next iteration	
	XDEL(I3)=X(I3)-XM


c Kinetic energy (ener(2)) should be hacked here.
        xkinetic = xkinetic + V(I3)*V(I3) * aamass


c              if (mod(nstep,ntpr) .eq. 0 ) then  
c                 write(75,*)nstep,I3,x(i3),F(I3),gammaa
c                 write(76,*)nstep,I3,xmass,v(i3),xdel(i3)
c             endif

       enddo                       
       enddo                       


c Kinetic energy (ener(2)) should be hacked here.
      ENER(2)=0.5*xkinetic

       goto 245

cc END of Langevin go to printout
cc skip the rest
c
c-----------------------------------------------------
c  ---Step 1: update the velocities using the forces:
c-----------------------------------------------------
C
  290 continue
c
c     --- note: looping separately over solute,
c         solvent, having peacs/nonpeacs versions
c         and checking vlimit in its own loop speeds 
c         up plain non-peacs, non-vlimit by ~2% for
c         6K atoms w/ 9A cut. -bross
c
      if (peac) then
c
c       --- peacs version, solute
c
#ifdef MPI
        i3 = 3*(istart-1)
        ibiggest = min( iend, nsolut )
        do j=istart,ibiggest
#else
        I3 = 0
        DO J = 1,nsolut
#endif
c
c         --- PEACS correction to wfac, Eq. 12 of Schaik, et al.
c
          fsq = f(i3+1)*f(i3+1) + f(i3+2)*f(i3+2) + f(i3+3)*f(i3+3)
          wfac = WINV(J)*DTX - (vzero-ener(23)) / (natom*tauv*fsq)
c
c         --- update the atom x,y,z
c
          DO M = 1,3
            I3 = I3+1
            V(I3) = (V(I3) + F(I3)*wfac) * SCALP
          enddo
        enddo
c
c       --- solvent
c
#ifdef MPI
        ismallest = max( nsolut+1, istart )
        do j=ismallest,iend
#else
        do j=nsolut+1,nrp
#endif
c
c         --- PEACS correction to wfac, Eq. 12 of Schaik, et al.
c
          fsq = f(i3+1)*f(i3+1) + f(i3+2)*f(i3+2) + f(i3+3)*f(i3+3)
          wfac = WINV(J) * DTX  -  (vzero-ener(23)) / (natom*tauv*fsq)
c
c         --- update the atom x,y,z
c
          do m=1,3
            I3 = I3+1
            V(I3) = (V(I3) + F(I3)*wfac) * SCALS
          enddo
        enddo
      else
c
c       --- non-peacs; solute
c
#ifdef MPI
        i3 = 3*(istart-1)
        ibiggest = min( iend, nsolut )
        do j=istart,ibiggest
#else
        I3 = 0
        DO J = 1,nsolut
#endif
          wfac = WINV(J) * DTX
c
c         --- update the atom x,y,z
c
          DO M = 1,3
            I3 = I3+1
            V(I3) = (V(I3) + F(I3)*wfac) * SCALP
          enddo
        enddo
c
c       --- solvent
c
#ifdef MPI
        ismallest = max( nsolut+1, istart )
        do j=ismallest,iend
#else
        do j=nsolut+1,nrp
#endif
          wfac = WINV(J) * DTX
c
c         --- update the atom x,y,z
c
          do m=1,3
            I3 = I3+1
            V(I3) = (V(I3) + F(I3)*wfac) * SCALS
          enddo
        enddo
      endif
c
c apply CHIE
c
      if (rlmavf.le.small) then
         iltavc = 0
      else if (iltavf.le.0) then
         iltavc = 3*(natom)
      else
         iltavc = 3*iltavf
      end if
#ifdef MPI
      do i=3*istart,min(i3,iltavc)
#else
      do i=1,min(i3,iltavc)
#endif
         v(i) = v(i) * chie
      end do
c
c     --- consider vlimit
c
      if (vlim) then
        vmax = 0.0d0
#ifdef MPI
        do i=3*istart,i3
#else
        do i=1,i3
#endif
          vmax = max(vmax,abs(v(i)))
          v(i) = sign(min(abs(v(i)),vlimit),v(i))
        enddo
        if (vmax.gt.vlimit .and. master) write(6,*) 
     .      'vlimit exceeded for step ',nstep,'; vmax = ',vmax
      endif
c
c     --- iscale
c
      if (iscale.gt.0) then
        do im=1,iscale
          v(nr3+im) = (v(nr3+im) + f(nr3+im)*dtx/scalm)*scalp
        enddo
      end if
c
c-------------------------------------------------------------------
c  Step 2: Calculate quantities needed to get kinetic energy:
c          (this section if SHAKE is not being carried out)
c-------------------------------------------------------------------
C
      ENEW(1) = 0.0d0
      ENEW(2) = 0.0d0
C
      if (NTC.eq.1) then
c
c       --- solute
c
#ifdef MPI
        i3 = 3*(istart-1)
        ibiggest = min( iend, nsolut )
        do j=istart,ibiggest
#else
        I3 = 0
        DO J=1,nsolut
#endif
          aamass = amass(j)
          DO M = 1,3
            I3 = I3+1
            ENEW(1) = ENEW(1) + V(I3)*V(I3) * aamass
          enddo
        enddo
c
c       --- solvent
c
#ifdef MPI
        ismallest = max( nsolut+1, istart )
        do j=ismallest,iend
#else
        DO j=nsolut+1,NRP
#endif
          aamass = amass(j)
          DO M = 1,3
            I3 = I3+1
            ENEW(2) = ENEW(2) + V(I3)*V(I3) * aamass
          enddo
        enddo
#ifdef MPI
c     =========================== AMBER/MPI ===========================
c
c  --- sum up the partial kinetic energies:
c
        if ( .not. mpi_orig .and. numtasks .gt. 1 ) then
           call mpi_allreduce(enew,enewt,2,
     +          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
           enew(1) = enewt(1)
           enew(2) = enewt(2)
        endif
c
c     ========================= END AMBER/MPI =========================
#endif
        if (iscale.gt.0) then
          do im=1,iscale
            ENEW(1) = ENEW(1) + scalm*v(nr3+im)*v(nr3+im)
          enddo
        end if
        ENEW(1) = ENEW(1) * 0.5d0
        ENEW(2) = ENEW(2) * 0.5d0
      end if
c
c-------------------------------------------------------------------
c   Step 3: update the positions, putting the "old" positions into F:
c-------------------------------------------------------------------
c
#ifdef MPI
      DO I3 = istart3,iend3
#else
      DO I3 = 1,NR3+iscale
#endif
        F(I3) = X(I3)
        X(I3) = X(I3)+V(I3)*DTX
      enddo
#ifdef MPI
      call second(timeg1)
      timsts(26) = timsts(26) + timeg1 - timeg0 
      timeg0 = timeg1
#endif
C
c
c-------------------------------------------------------------------
c   Step 4: if shake is being used, update the new positions to fix the
c      bond lengths, and then re-estimate the velocities from
c      differences in positions:
c-------------------------------------------------------------------
c
      IF (NTC.ne.1) then
        call second(time0)
        CALL SHAKE(NRP,NBONH,MBONA,0,IX(I12),IX(I14),IX(I62),
     +           WINV,CONP,SKIP,F,X,NITP,BELLY,IX(IVM01))
        IF(NITP.EQ.0) then
          erstop = .true.
          GOTO 480
        endif
        CALL QUICK3(F,X,IX(IVM02),NATOM,NRES,IX(I02),XX(L85),XX(L90))
        call second(time1)
        timsts(6) = timsts(6) + time1 - time0
        timeg0 = time1
        ENEW(1) = 0.0d0
        ENEW(2) = 0.0d0
c
c       --- solute
c
#ifdef MPI
        i3 = 3*(istart-1)
        ibiggest = min( iend, nsolut )
        do j=istart,ibiggest
#else
        I3 = 0
        DO J = 1,nsolut
#endif
          aamass = amass(j)
          DO M = 1,3
            I3 = I3+1
            V(I3) = (X(I3)-F(I3)) * dtxinv
            if (vlim) v(i3) = sign(min(abs(v(i3)),vlimit),v(i3))
            ENEW(1) = ENEW(1) + V(I3)*V(I3) * aamass
          enddo
        enddo
c
c       --- solvent
c
#ifdef MPI
        ismallest = max( nsolut+1, istart )
        do j=ismallest,iend
#else
        DO j=nsolut+1,NRP
#endif
          aamass = amass(j)
          DO M = 1,3
            I3 = I3+1
            V(I3) = (X(I3)-F(I3)) * dtxinv
            if (vlim) v(i3) = sign(min(abs(v(i3)),vlimit),v(i3))
            ENEW(2) = ENEW(2) + V(I3)*V(I3) * aamass
          enddo
        enddo
#ifdef MPI
c     =========================== AMBER/MPI ===========================
c
c  --- sum up the partial kinetic energies:
c
        if ( .not. mpi_orig .and. numtasks .gt. 1 ) then
           call mpi_allreduce(enew,enewt,2,
     +          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
           enew(1) = enewt(1)
           enew(2) = enewt(2)
        endif
c
c     ========================= END AMBER/MPI =========================
#endif
        if (iscale.gt.0) then
          do im=1,iscale
            v(nr3+im) = (x(nr3+im) - f(nr3+im)) / dtx
            ENEW(1) = ENEW(1) + scalm*v(nr3+im)*v(nr3+im)
          enddo
        end if
        ENEW(1) = ENEW(1) * 0.5d0
        ENEW(2) = ENEW(2) * 0.5d0
      end if
#ifdef MPI
c     =========================== AMBER/MPI ===========================
c
c  --- now distribute the coordinates:
c      NOTE: this is *only* necessary when mpi_orig is false!
c
      call second(timeg1)
      timsts(26) = timsts(26) + timeg1 - timeg0 
      timeg0 = timeg1

      if ( .not. mpi_orig .and. numtasks .gt. 1 ) then
        do i = 0, numtasks-1
          rcvcnt(i) = iparpt3(i+1) - iparpt3(i)
        enddo
#  ifdef CSPP
        call MPI_ALLGATHERV(
     .     %val(0),
     .     rcvcnt(mytaskid),
     .     MPI_DOUBLE_PRECISION,
     .     x,
     .     rcvcnt,
     .     iparpt3,
     .     MPI_DOUBLE_PRECISION,
     .     MPI_COMM_WORLD, ierr)
#  else
        call MPI_ALLGATHERV(
     .     x(iparpt3(mytaskid)+1),
     .     rcvcnt(mytaskid),
     .     MPI_DOUBLE_PRECISION,
     .     xx(l95),
     .     rcvcnt,
     .     iparpt3,
     .     MPI_DOUBLE_PRECISION,
     .     MPI_COMM_WORLD, ierr)
        do i = 1, 3*natom
          x(i) = xx(l95+i-1)
        enddo
#  endif
c
c     It is necessary for every processor to have a full copy of the
c     velocities currently for constant pressure simulations, 
c     removal of center of mass motion, etc.  (This perhaps can be 
c     optimized by fixing EKCMR() to do parallel and allreducing 
c     the EKCMT values).
c     Also, the master needs velocities when coordinate
c     archiving and at the final run; this can be accomplished 
c     with a simple gather to the master node rather than the 
c     presumably more expensive allgatherv().  This is not important
c     for the T3D and other fast communication machines but may
c     provide some limited advantage on clusters.
c
c     For simplicity for now, an allgatherv is done...
c
        if ( (plevel .eq. 1 .and. NTP .gt. 0) .or.
     .     (NTWV.GT.0 .AND. NTWV.LE.NTWVM) .or.
     .     (NSTEP+1 .ge. NSTLIM) .or.
     .     (mod(nstep+1,nscm).eq.0) ) then
#  ifdef CSPP
          call MPI_ALLGATHERV(
     .        %val(0),
     .        rcvcnt(mytaskid),
     .        MPI_DOUBLE_PRECISION,
     .        v,
     .        rcvcnt,
     .        iparpt3,
     .        MPI_DOUBLE_PRECISION,
     .        MPI_COMM_WORLD, ierr)
#  else
          call MPI_ALLGATHERV(
     .        v(iparpt3(mytaskid)+1),
     .        rcvcnt(mytaskid),
     .        MPI_DOUBLE_PRECISION,
     .        xx(l100),
     .        rcvcnt,
     .        iparpt3,
     .        MPI_DOUBLE_PRECISION,
     .        MPI_COMM_WORLD, ierr)

          do i = 1, 3*natom
             v(i) = xx(l100+i-1)
          enddo
#  endif
        endif
c
c     ...end of the check if .not. mpi_orig
c
      endif

      call second(timeg1)
      timsts(22) = timsts(22) + timeg1 - timeg0
      timeg0 = timeg1
c
c     ========================= END AMBER/MPI =========================
#endif
c
c-------------------------------------------------------------------
c   Step 4a: zero COM velocity if requested; used for preventing 
c   ewald "block of ice flying thru space" phenomenon
c-------------------------------------------------------------------
c
      if (mod(nstep+1,nscm).eq.0 .and. NTB.ne.0) then
c
        vcmx = 0.d0
        vcmy = 0.d0
        vcmz = 0.d0

        j = 1
        do i = 1, 3*natom,3
           aamass = amass(j)

           vcmx = vcmx + aamass * v(i)
           vcmy = vcmy + aamass * v(i+1)
           vcmz = vcmz + aamass * v(i+2)

           j = j + 1
        enddo

        vcmx = vcmx * tmassinv
        vcmy = vcmy * tmassinv
        vcmz = vcmz * tmassinv

        vel2 = vcmx*vcmx + vcmy*vcmy + vcmz*vcmz
        atempdrop = 0.5d0 * tmass * vel2 / fac(1)
        vel = sqrt(vel2)
        if ( master .and. vel.gt.0.01d0 ) then
             write (6,'(a,f9.2,a,f9.2,a)')
     .           '  removed COM velocity ',vel,' approx temp drop ',
     .           atempdrop, ' K'
        endif

        do i = 1, 3*natom, 3
             v(i)   = v(i)   - vcmx
             v(i+1) = v(i+1) - vcmy
             v(i+2) = v(i+2) - vcmz
        enddo
      endif
c
c-------------------------------------------------------------------
c  Step 5: Zero out the non-moving forces if a belly is active:
c-------------------------------------------------------------------
C
      IF (BELLY) CALL BELLYF(NR,IX(I62),V)
c
c-------------------------------------------------------------------
c  Step 6: scale coordinates if constant pressure run:
c-------------------------------------------------------------------
c
      if (ntp.eq.1) then
        RMU(1) = (1.d0-DTCP*(PRES0-PRES(4)))**ONET
        RMU(2) = RMU(1)
        RMU(3) = RMU(1)
      else if (ntp.gt.1) then
        RMU(1) = (1.d0-DTCP*(PRES0-PRES(1)))**ONET
        RMU(2) = (1.d0-DTCP*(PRES0-PRES(2)))**ONET
        RMU(3) = (1.d0-DTCP*(PRES0-PRES(3)))**ONET
      end if
      if (ntp.gt.0) then
        do m=1,3
          BOX(M) = BOX(M)*RMU(M)
          ENER(M+6) = BOX(M)
        enddo
        if (iewald .eq. 1) then
c
c         ...a series of routines is called when IEWALD is active 
c         instead of MDBOX and PSCALE.
c
c         WARNING!!   This is not correct for non-orthogonal boxes if
c         NTP > 1 (i.e. non-isotropic scaling).  Currently general cell
c         updates which allow cell angles to change are not implemented.
c         The viral tensor computed for ewald is the general Nose Klein,
c         however the cell response needs a more general treatment.
c
           call redo_ucell(rmu)
           call ew_pscale(natom,x,amass,nspm,nsp,npscal)
           if (ntr.gt.0 .and. nrc.gt.0)
     .          call ew_pscale(natom,xc,amass,nspm,nsp,npscal)
        else
C
C         ----- RESET THE BOX PARAMETERS -----
C
          CALL MDBOX
C
C         ----- SCALE THE COORDINATES APPROPRIATELY -----
C
          CALL PSCALE(NSPM,NSP,NR,RMU,X,TMA,amass,NPSCAL)
          IF(NTR.gt.0.and.NRC.gt.0) 
     .      CALL PSCALE(NSPM,NSP,NR,RMU,XC,TMA,amass,NPSCAL)
        endif
      endif
c

      if (NSHAK.eq.0) then
        ENER(3) = (EOLD3+ENEW(1)) * 0.5d0
        ENER(4) = (EOLD4+ENEW(2)) * 0.5d0
      end if

      ENER(2) = ENER(3)+ENER(4)
      ENER(1) = ENER(2)+ENER(23)
      EOLD3 = ENEW(1)
      EOLD4 = ENEW(2)
      EKKEEP = EOLD3 + EOLD4
      if (NTT.ne.0) then
        EKGP = ENEW(1)
        EKGS = ENEW(2)
      end if

cccc Margaret Langevin
cccc printout
C
c-------------------------------------------------------------------
c  Step 7:  update the step counter and the integration time:
c-------------------------------------------------------------------
c
cc Langevin needs this. Margaret
  
  245  continue


   
 
      NSTEP = NSTEP+1
      Td = Td+DT
      t = td
      ISTEP = ISTEP+1

c Margaret folding
c Margaret sh3 (-53)
c Margaret im9 (-53)
c Margaret im7 (-52)
c Margaret ci2 (-48)
c Margaret ci2.i+4 (-48)
c Margaret ilb (-173)
c Margaret im7,i+4 (-69)
c Margaret u1a, (-76)
c 6-12LJ
c       if(ENER(24).lt.-25.5) then
c	NSTEP = NSTLIM
c 	NRU = NRUN
c      endif


c Margaret shg v58t, (-51)
c Margaret shg v44t, (-50)
cccccccccccccccccccccccccccccccccccccccccccccccc
c 10-12LJ
c      if(ENER(26).lt.-50) then
c	NSTEP = NSTLIM
c      endif

c Margaret unfolding
c Margaret sh3 (-5)
c Margaret im9 (-5)
c Margaret im7 (-5)
c Margaret ci2 (-5)
c      if(ENER(24).gt.-5) then
c	NSTEP = NSTLIM
c      endif
c end of Margaret

      FIT = FIT + ENER(1)
      FITI = FITI + ENER(1)*ISTEP
      FIT2 = FIT2 + ENER(1)*ENER(1)
      DO M = 1,NREN
        ENERT(M) = ENERT(M)+ENER(M)
        ENERT2(M) = ENERT2(M) + ENER(M)*ENER(M)
      enddo
C
      NTNB = 0
      NSHAK = 0
      if (mod(nstep,nsnb) .eq. 0) ntnb = 1
      lout = mod(nstep,ntpr) .eq. 0

cc Langevin goto printout Margaret
c  Need  to skip step 8. set NSCM==0 to skip 
C
c-------------------------------------------------------------------
C  Step 8: stop motion about center of mass (ntb.eq.0) if requested:
c-------------------------------------------------------------------
C
C     Note: no one seems to know how this works: it appears to
C     involve taking an extra step and going back to init=3 code.
C     See the ntb.ne.0 version above for simpler, hopefully
C     correct version. Modified so that it happens before writing
C     restrt & mdvel, with the "goto init=3" remaining after the
C     file writing. 11/96
C     -Bill Ross
C
      if (mod(nstep,nscm).eq.0 .and. NTB.eq.0) then
c
c       -- no box: stop translation & rotation
c
        DO I3 = 1,NR3
          F(I3) = X(I3)-V(I3)*DT5
        enddo
        CALL CENMAS(NRP,IZERO,F,V,tmass,tmassinv,IZERO,amass,EKCM,
     +             XCM,VCM,ACM,EKROT,OCM,IFOUR,LOUT)
        DO I3 = 1,NR3
          X(I3) = F(I3)
        enddo

        CALL STOPCM(NR,X,V,XCM,VCM,OCM)
        CALL CENMAS(NRP,IZERO,X,V,tmass,tmassinv,IZERO,amass,EKCM,
     +              XCM,VCM,ACM,EKROT,OCM,IFOUR,LOUT)
      endif
C
c-------------------------------------------------------------------
C  Step 9:  output from this step if required:
c-------------------------------------------------------------------
c
c     ...only the master needs to do the output
c

cc Langein output. Margaret

cc  245  continue
      if (master) then
C
C        -- restrt:
C
         iprint = 0
         if (ntwr.ne.0) then
           if (mod(nstep,ntwr) .eq. 0 ) then
             iprint = 1
           endif
         endif
         if (nstep.eq.nstlim) iprint = 1
         if (iprint.eq.1) then
            NR = NRP+NSM*NRAM
            CALL MDWRIT(nstep,NRP,NR,NRES,NTXO,NTR,NTB,X,V,
     +                  XX(L55),BOX,ih(m04),ih(m02),IX(I02),T)
         endif
C
C        -- Coordinate archive:
C
c     IF (NTWX.GT.0 .AND. NTWX.LE.NTWXM) THEN  changed 9/95 -- dap
         IF (NTWX.GT.0 .AND. nstep.ge.NTWXM) THEN
           IF (mod(nstep,ntwx) .eq. 0) then
             CALL CORPAC(X,NRX0,NRX,12,LOUTFM)
             IF(NTB.GT.0) then
               if ( iewald .eq. 1 ) 
     .           call find_ucell(box(1),box(2), box(3),
     .                           alpha,betta,gamma)
               CALL CORPAC(BOX,1,3,12,LOUTFM)
             endif
           END IF
         END IF
C
C Velocity archive:
C
c        IF (NTWV.GT.0 .AND. NTWV.LE.NTWVM) THEN  changed 9/95 -- dap
         IF (NTWV.GT.0 .and. nstep.ge.ntwvm) THEN
            IF (NTWV.LE.NTWVM .AND. MOD(NSTEP,NTWV).EQ.0)
     +           CALL CORPAC(V,NRX0,NRX,13,LOUTFM)
         ENDIF
C
C Energy archive:
C
c        IF (NTWE.GT.0 .AND. NSTEP.ge.NTWEM) THEN  changed 9/95 -- dap
         IF (NTWE.GT.0 .and. nstep.ge.ntwem) THEN
            IF (MOD(NSTEP,NTWE).EQ.0) then
                CALL MDENG(15,NSTEP,T,ENER,FAC)
c Margaret add to print out twvdw and hbvdw
c ener(30) includes scnb, 1-4 van der waals

              write(37,*)EHBA
              write(36,*)EHBV
              write(38,*)ECHI
              ener(42)=ECHI

              write(39,*)ENPC
              write(40,*)ENCC
	     endif
         END IF
C
         IF (LOUT .or. (init .eq. 4 .and. nstep .eq. 1)) THEN
           call opinfo(7)
           CALL PRNTMD(NSTEP,NITP,NITS,T,ENER,FAC,bscale,iscale,7)
           IF (NMROPT.GT.0) THEN
              CALL NMRPTX(6)
              CALL NMRPTX(7)
           END IF
           CALL CLOSC(7)
         END IF
C
C If IWRTSM.GT.0, write energy/restrained coord archive file:
C
         IF (IWRTSM.GT.0) THEN
            IF (MOD(NSTEP,IWRTSM).EQ.0) THEN
C          TINE  WRTRSL(     X         ,RIMASS    ,NTB       ,MDSTEP   ,
C    *            NITP      ,NITS      ,T         ,ENER      ,FAC      ,
C    *            BSCALE    ,ISCALE    ,NCALLS    ,RMS       ,FDMAX    ,
C    *            IATMAX    ,LABMAX    ,ISET      ,ITYP      ,IWRT     ,
C    *            WORK      ,IWORK     )
           CALL  WRTRSL(     X         ,WINV      ,NTB       ,NSTEP    ,
     *            NITP      ,NITS      ,T         ,ENER      ,FAC      ,
     *            BSCALE    ,ISCALE    ,NCALLS    ,dumm      ,dumm     ,
     *            idumm     ,idumm     ,IACSET    ,1         ,IUWRTS   ,
     *            XX(LNMR01),IX(INMR02))
            END IF
         END IF
C
C If Ligand Grid run being performed (LGRON = 1) write the ligand grid
C restart file every LGWRIT steps
C
         IF (LGRON.EQ.1 .AND. MOD(NSTEP,LGWRIT).EQ.0) THEN
C           TINE LGRST(IUN       ,IFUNC     ,NTYPES    ,CN1       ,
C    *       CN2)
            CALL LGRST(66        ,0         ,NTYPES    ,CN1       ,
     *       CN2)
         END IF
C
         IF (MOD(NSTEP,INTCM).EQ.0 .OR. NSTEP.EQ.NSTLIM) THEN
C
C           --- PERFORM A LEAST SQUARES FIT OF THE TOTAL ENERGY TO
C              A STRAIGHT LINE ONCE IN 100 STEPS -----
C
            CALL LSQFIT(ISTEP,FIT,FITI,FIT2,.true.)
            ISTEP = 0
            FIT = 0.d0
            FITI = 0.d0
            FIT2 = 0.d0
          end if
c
c         --- end masters output ---
c
      end if
c
c=======================================================================
C
c  ---major cycle back to new step unless we have reached our limit:
c
      IF (NSTEP.LT.NSTLIM) then
        if (mod(nstep,nscm).eq.0 .and. NTB.eq.0) goto 120
        GOTO 260
      endif
  480 CONTINUE
c
c=======================================================================
C
C     ----- PRINT AVERAGES -----
C
c=======================================================================
c
c
c     ...update the velocities if necessary and then return unless
c     you are the master (who needs to print the final information)
c     
      if (master) then
         TSPAN = NSTEP
         IF (NSTEP.GT.0) THEN
           DO M = 1,NREN
             ENERT(M) = ENERT(M)/TSPAN
             ENERT2(M) = ENERT2(M)/TSPAN - ENERT(M)*ENERT(M)
             IF(ENERT2(M).LT.0.d0) ENERT2(M) = 0.d0
             ENERT2(M) =  SQRT(ENERT2(M))
           enddo
C
           WRITE(6,540) NSTEP
           CALL PRNTMD(NSTEP,IZERO,IZERO,T,ENERT,FAC,bscale,iscale,0)
           IF (NMROPT.GT.0) CALL NMRPTX(6)
           WRITE(6,550)
           CALL PRNTMD(NSTEP,IZERO,IZERO,T,ENERT2,FAC,bscale,iscale,0)
           IF (NMROPT.EQ.1) THEN
              WRITE(6,500)
  500         FORMAT(/,' NMR restraints on final step:'/)
              CALL NDVPTX(X,F,IH(M04),IH(M02),IX(I02),NRES,XX(L95),
     +                   NATOM, XX(L20),NTB,XX(LNMR01),IX(INMR02),6)
           END IF
C
           DO M = 2,NREK
             ENERT(M) = ENERT(M)/FAC(M-1)
             ENERT2(M) = ENERT2(M)/FAC(M-1)
           enddo
           TEMP = ENERT(2)
         END IF
c
C If IWRTSM.NE.0, write energy/restrained coord archive file before returning
C
         IF (IWRTSM.NE.0) THEN
            IF (IWRTSM.LT.0 .OR. 
     *      (IWRTSM.GT.0 .AND. MOD(NCALLS,IWRTSM).NE.0)) THEN
            CALL WRTRSL(     X         ,WINV      ,NTB       ,NSTEP    ,
     *            NITP      ,NITS      ,T         ,ENER      ,FAC      ,
     *            BSCALE    ,ISCALE    ,NCALLS    ,dumm      ,dumm     ,
     *            idumm     ,idumm     ,IACSET    ,1         ,IUWRTS   ,
     *            XX(LNMR01),IX(INMR02))
            END IF
         END IF
C
C Before returning, write the ligand grid stuff if necessary
C
         IF (LGRON.EQ.1) THEN
            CALL LGWRT(XX(L95)   ,66        ,TEMP0)

C           TINE LGRST(IUN       ,IFUNC     ,NTYPES    ,CN1       ,
C    *       CN2)

            CALL LGRST(66        ,0         ,NTYPES    ,CN1       ,
     *       CN2)
         END IF
C
c        --- end of masters velocity stuff ---
c
      endif
C
  520 FORMAT('   SHAKE COORDINATES, NITER  = ',I4)
  530 FORMAT('   SHAKE VELOCITIES,  NITER  = ',I4)
  540 FORMAT(/5X,' A V E R A G E S   O V E R ',I5,' S T E P S',/)
  550 FORMAT(/5X,' R M S  F L U C T U A T I O N S',/)
  560 FORMAT(5E15.8)
  570 FORMAT(' ERROR: NTT=',I5,'; Unrecognized Option.')
      RETURN
      END
