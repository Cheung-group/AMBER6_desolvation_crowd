c-----------------------------------------------------------------------
      SUBROUTINE CENMAS(NRP,NPA,X,V,tmass,tmassinv,NPW,amass,
     +                  EKCM,XCM,VCM,ACM,EKROT,OCM,ICM,LOUT)
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
#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
#ifdef MPI
#include "parallel.h"
#endif     
#include "extra.h"
c Margaret
#include "HB.h"
#include "CHI.h"
      LOGICAL LOUT
C
C     ----- ROUTINE TO CALCULATE THE TRANSLATIONAL AND ROTATIONAL
C           KINETIC ENERGIES AND VELOCITIES -----
C
c     mods for Rev A by gls:
c     Simplified output.  Fixed bug in report of EKtrans and EKrot.
c        same bug applied to total mass; uncleared 1P edit descriptor
c        was making these three values 10X too big.
c        New output is same in 32/64 bits.
c     cpp-switchable double/single prec
c     changed call to 'aminv' to 'matinv'
c     removed winvs(*) reference, since it wasnt used.
c
      DIMENSION X(*),V(*),amass(*),XCM(*),VCM(*),ACM(*),OCM(*)
      DIMENSION TCM(3,3),LH(3),MH(3)
C
      DATA CRIT/1.0d-06/
C
      IF (ICM.EQ.0) return
C
      IF (NPA.NE.0.OR.NPW.NE.0) then
         write(6,9990)
 9990    format ('%CENMAS-F-error, NPA and NPW must both be zero, ',
     -     'check input') 
         call mexit(6,1)
      endif
C
      NRPT = NRP
      NR = NRPT
      I0 = 3*NPA
C
C     ----- CALCULATE THE CENTER OF MASS COORDINATES -----
C
      XCM(1) = 0.d0
      XCM(2) = 0.d0
      XCM(3) = 0.d0

      I = I0
      DO J = 1,NRP
        aamass = amass(NPW+J)
        DO M = 1,3
          I = I+1
          XCM(M) = XCM(M) + X(I)*aamass
        enddo
      enddo

      XCM(1) = XCM(1) * tmassinv
      XCM(2) = XCM(2) * tmassinv
      XCM(3) = XCM(3) * tmassinv
      
*     IF(LOUT) WRITE(6,301) (XCM(M),M = 1,3),tmass
C
      IF (IABS(ICM).LT.2) return
C
C     ----- CALCULATE THE VELOCITY AND THE TRANSLATIONAL 
C           KINETIC ENERGY OF THE CENTRE OF MASS -----
C
      EKCM = 0.d0
      VCM(1) = 0.0d0
      VCM(2) = 0.0d0
      VCM(3) = 0.0d0

      I = I0
      DO J = 1,NRP
        aamass = amass(NPW+J)
        DO M = 1,3
          I = I+1
          VCM(M) = VCM(M) + V(I)*aamass
        enddo
      enddo

      DO M = 1,3
        VCM(M) = VCM(M) * tmassinv
        EKCM = EKCM + VCM(M)*VCM(M)
      enddo

      EKCM = EKCM * tmass * 0.5d0
      comvel = sqrt(vcm(1)*vcm(1)+vcm(2)*vcm(2)+vcm(3)*vcm(3))
c     IF (ICM.ge.0 .and. lout) WRITE(6,302) (VCM(M),M = 1,3),EKCM
C
      IF(IABS(ICM).LT.3) return
c
      IF (NPA.NE.0.OR.NPW.NE.0) then
         write(6,9991)
 9991    format ('%CENMAS-F-error2, NPA and NPW must both be zero, ',
     -     'check input') 
         call mexit(6,1)
      endif
c
C     ----- CALCULATE THE ANGULAR MOMENTUM ABOUT THE 
C           CENTER OF MASS ----
C
      ACM(1) = 0.0d0
      ACM(2) = 0.0d0
      ACM(3) = 0.0d0

      I = 0
      DO J = 1,NR
        aamass = amass(J)
        ACM(1) = ACM(1) + (X(I+2)*V(I+3)-X(I+3)*V(I+2)) * aamass
        ACM(2) = ACM(2) + (X(I+3)*V(I+1)-X(I+1)*V(I+3)) * aamass
        ACM(3) = ACM(3) + (X(I+1)*V(I+2)-X(I+2)*V(I+1)) * aamass
        I = I+3
      enddo

      ACM(1) = ACM(1) - (XCM(2)*VCM(3)-XCM(3)*VCM(2)) * tmass
      ACM(2) = ACM(2) - (XCM(3)*VCM(1)-XCM(1)*VCM(3)) * tmass
      ACM(3) = ACM(3) - (XCM(1)*VCM(2)-XCM(2)*VCM(1)) * tmass
c     IF (ICM.ge.0 .and. lout) WRITE(6,303) (ACM(M),M = 1,3)

      IF (IABS(ICM).LT.4) return
C
C     ----- CALCULATE THE INERTIA TENSOR -----
C
      XX = 0.d0
      XY = 0.d0
      XZ = 0.d0
      YY = 0.d0
      YZ = 0.d0
      ZZ = 0.d0

      I = 0
      DO J = 1,NR
        X1 = X(I+1)-XCM(1)
        X2 = X(I+2)-XCM(2)
        X3 = X(I+3)-XCM(3)
        aamass = amass(J)
        XX = XX+X1*X1*aamass
        XY = XY+X1*X2*aamass
        XZ = XZ+X1*X3*aamass
        YY = YY+X2*X2*aamass
        YZ = YZ+X2*X3*aamass
        ZZ = ZZ+X3*X3*aamass
        I = I+3
      enddo
      TCM(1,1) = YY+ZZ
      TCM(2,1) = -XY
      TCM(3,1) = -XZ
      TCM(1,2) = -XY
      TCM(2,2) = XX+ZZ
      TCM(3,2) = -YZ
      TCM(1,3) = -XZ
      TCM(2,3) = -YZ
      TCM(3,3) = XX+YY
C
C     ----- INVERT THE INERTIA TENSOR -----
C
      CALL matinv(TCM,3,D,LH,MH)
      IF(ABS(D).le.CRIT) then
        WRITE(6,307)
  307   format(/5x,'%CENMAS-F-INERTIA_TENSOR,  determinant',
     -           ' is zero ... stop')
        call mexit(6,1)
      endif
C
C     ----- CALCULATE THE ANGULAR VELOCITY ABOUT THE CENTER OF 
C           MASS AND THE ROTATIONAL KINETIC ENERGY -----
C
      EKROT = 0.d0
      DO 240 M = 1,3
        OCM(M) = 0.d0
        DO 230 N = 1,3
  230     OCM(M) = OCM(M)+TCM(M,N)*ACM(N)
  240   EKROT = EKROT+OCM(M)*ACM(M)
      EKROT = EKROT * 0.5d0

      IF (ICM.lt.0) return
      if (master)
     +   write(6,'(/3x,a,f11.4,3x,a,f11.4,3x,a,f12.6)') 'KE Trans =',
     +          ekcm, 'KE Rot =', ekrot,'C.O.M. Vel =',comvel
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE corpac(R,ISTART,N,NF,LOUTFM)
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
C
C Revisions:
C    Check max. coordinate. If it is < -999.999 or > 9999.999, change to
C    F8.2 format for formatted writes to prevent overflows. -- dap, 5/93
C
C    Add ISTART to call list. Allows one to write a subrange of
C    values in the middle of the array.

#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
      PARAMETER (RMAX = 9999.99D0)
      PARAMETER (RMIN = -999.99D0)
C
      LOGICAL LOUTFM
      DIMENSION R(N)

      IF (ISTART.GT.N) RETURN
 
C Unformatted writes:

      IF (.NOT.LOUTFM) THEN
         WRITE(NF) (R(J),J=ISTART,N)
      ELSE

C Formatted writes:

         IMAX = 0
         DO 10 I = ISTART,N
           IF (R(I).GT.RMAX .OR. R(I).LT.RMIN) THEN
              IMAX = 1
              GO TO 20
           END IF
   10    CONTINUE
C
   20    IF (IMAX.EQ.0) THEN
            WRITE(NF,1000) (R(J),J=ISTART,N)
         ELSE
            WRITE(NF,1001) (R(J),J=ISTART,N)
         END IF
      END IF


 1000 FORMAT(10F8.3)
 1001 FORMAT(10F8.2)
      RETURN
      END
c-----------------------------------------------------------------------
#ifdef MPI
      SUBROUTINE EKCMR(NSPM,NSP,TMA,EKCMT,XR,V,amass,istart,iend)
#else
      SUBROUTINE EKCMR(NSPM,NSP,TMA,EKCMT,XR,V,amass)
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
#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
C
C     ----- ROUTINE TO CALCULATE THE TOTAL KINETIC ENERGY OF THE
C           CENTER OF MASS OF THE SUB-MOLECULES AND ALSO THE 
C           COORDINATES OF THE MOLECULES RELATIVE TO THE CENTER OF
C           MASS -----
C
c     Rev A mods: cleanup, DP constants
#include "box.h"
C
      DIMENSION NSP(*),TMA(*),EKCMT(*),XR(*),V(*),amass(*)
      DIMENSION XCM(3),VCM(3),XI(3),EKCML(3)
#ifdef CRAY_MP
      real ekcmlx(3,nspm)
      integer i3(nspm), iat(nspm)
#endif
      data zero, one /0.0d0, 1.0d0/
C
      I3 = 0
      IAT = 0

#ifdef CRAY_MP
      i3(1) = 0
      iat(1) = 0
      do i = 2, nspm
         i3(i)  = 3*nsp(i-1)+  i3(i-1)
         iat(i) =   nsp(i-1)+ iat(i-1)
      enddo
#endif

      EKCML(1) = zero
      EKCML(2) = zero
      EKCML(3) = zero
C
#ifdef CRAY_MP
cmic$ do all autoscope if(nspm.gt.1200)
      do i=1,nspm
        ekcmlx(1,i)=0.
        ekcmlx(2,i)=0.
        ekcmlx(3,i)=0.
      end do

cmic$ parallel autoscope
cmic$*private(nn,xi,xcm,vcm)
cmic$*shared(xr,beta)
cmic$ do parallel
#endif
c
      DO 900 N = 1,NSPM
        NN = NSP(N)
        TMN = one/TMA(N)
C
C       ----- FIRST GATHER THE ATOMS OF A MOLECULE BY APPLYING
C             THE BOUNDARY CONDITIONS SUCH THAT EACH ATOM LIES WITHIN
C             BOX/2 OR BOX*SQRT(3)/4 OF ITS PREDECESSOR. THE FIRST ATOM
C             OF THE MOLECULE IS COMPARED TO THE REFERENCE POSITION XR ---
C
        IF(NTB.ne.0) then
#ifdef CRAY_MP
          J3 = i3(n)
#else
          J3 = I3
#endif
          XI(1) = XR(J3+1)
          XI(2) = XR(J3+2)
          XI(3) = XR(J3+3)
          call bound3(nn,xi,xr,j3)
        endif
C
C       ----- NOW CALCULATE THE CENTER OF MASS AND THEN MOVE EACH
C             SUB-MOLECULE TO ITS CENTER OF MASS -----
C
#ifdef CRAY_MP
        j3 = i3(n)
#else
        J3 = I3
#endif
        XCM(1) = zero
        XCM(2) = zero
        XCM(3) = zero
        VCM(1) = zero
        VCM(2) = zero
        VCM(3) = zero
        DO 220 J = 1,NN
#ifdef CRAY_MP
          aamass = amass(IAT(n)+J)
#else
          aamass = amass(IAT+J)
#endif
          XCM(1) = XCM(1)+XR(J3+1)*aamass
          XCM(2) = XCM(2)+XR(J3+2)*aamass
          XCM(3) = XCM(3)+XR(J3+3)*aamass
#ifdef MPI
          if (iat+j.ge.istart .and. iat+j.le.iend) then
#endif
          VCM(1) = VCM(1)+V(J3+1)*aamass
          VCM(2) = VCM(2)+V(J3+2)*aamass
          VCM(3) = VCM(3)+V(J3+3)*aamass
#ifdef MPI
          endif
#endif

          J3 = J3+3
  220   CONTINUE
        XCM(1) = XCM(1)*TMN
        XCM(2) = XCM(2)*TMN
        XCM(3) = XCM(3)*TMN
C
#ifdef CRAY_MP
        J3 = i3(n)
#else
        J3 = I3
#endif
        DO 240 J = 1,NN
          XR(J3+1) = XR(J3+1)-XCM(1)
          XR(J3+2) = XR(J3+2)-XCM(2)
          XR(J3+3) = XR(J3+3)-XCM(3)
          J3 = J3+3
  240   CONTINUE
C
#ifdef CRAY_MP
        EKCMLx(1,n) = TMN*VCM(1)*VCM(1)
        EKCMLx(2,n) = TMN*VCM(2)*VCM(1)
        EKCMLx(3,n) = TMN*VCM(3)*VCM(1)
  900 CONTINUE
cmic$ end do
cmic$ end parallel

        do i=1,nspm
          ekcml(1)=ekcml(1)+ekcmlx(1,i)
          ekcml(2)=ekcml(2)+ekcmlx(2,i)
          ekcml(3)=ekcml(3)+ekcmlx(3,i)
        enddo
#else
        EKCML(1) = EKCML(1)+TMN*VCM(1)*VCM(1)
        EKCML(2) = EKCML(2)+TMN*VCM(2)*VCM(2)
        EKCML(3) = EKCML(3)+TMN*VCM(3)*VCM(3)
C
C       ----- END OF CALCULATION FOR EACH SUB-MOLECULE -----
C
        I3 = J3
        IAT = IAT+NN
  900 CONTINUE
#endif
C
      EKCMT(1) = EKCML(1)
      EKCMT(2) = EKCML(2)
      EKCMT(3) = EKCML(3)
C
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE OUTOPN
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
#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
#include "files.h"
#include "md.h"
C
C     -- ROUTINE TO OPEN THE DUMPING and restart FILES ---
C
C
c     subr amopen(lun,fname,fstat,fform,facc)
      IF(IOUTFM.le.0) then
C
C     ----- FORMATTED DUMPING -----
C
           if (ntwx .gt. 0) then
                call amopen(12,mdcrd,owrite,'F','W')
                WRITE(12,1000) ITITL
           endif
           if (ntwv .gt. 0) then
                call amopen(13,mdvel,owrite,'F','W')
                WRITE(13,1000) ITITL
           endif
           if (ntwe .gt. 0) then
                call amopen(15,mden,owrite,'F','W')
c               WRITE(15,1000) ITITL
           endif
c Margaret adds
           if (ntwe .gt. 0) then
                call amopen(37,twhb,owrite,'F','W')
                WRITE(37,1000) ITITL
           endif
           if (ntwe .gt. 0) then
                call amopen(36,twvdw,owrite,'F','W')
                WRITE(36,1000) ITITL
           endif
           if (ntwe .gt. 0) then
                call amopen(38,twchi,owrite,'F','W')
                WRITE(38,1000) ITITL
           endif
           if (ntwe .gt. 0) then
                call amopen(39,fenpc,owrite,'F','W')
                WRITE(39,1000) ITITL
           endif
           if (ntwe .gt. 0) then
                call amopen(40,fencc,owrite,'F','W')
                WRITE(40,1000) ITITL
           endif

      else
C
C     ----- UNFORMATTED DUMPING -----
C
           if (ntwx .gt. 0) then
                call amopen(12,mdcrd,owrite,'U','W') 
                WRITE(12) ITITL
           endif 
           if (ntwv .gt. 0) then
                call amopen(13,mdvel,owrite,'U','W')   
                WRITE(13) ITITL
           endif
           if (ntwe .gt. 0) then
                call amopen(15,mden,owrite,'U','W')
                WRITE(15) ITITL
           endif
      endif
C
C     ----- OPEN THE RESTART FILE -----
C
      if (ntxo.le.0) then
           call amopen(16,restrt,owrite,'U','W') 
      else
           call amopen(16,restrt,owrite,'F','W') 
      endif
 1000 FORMAT(20A4)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE OPINFO(NF)
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
#include "files.h"
c     --- opens the info file for md/gibbs ---
c     subr amopen(lun,fname,fstat,fform,facc)
      call amopen(nf,mdinfo,'U','F','W')
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE PRNTMD(NSTEP,NITP,NITS,TIME,ENER,FAC,BSCALE,ISCALE, 
     $                  iout7)
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
c     Mods for Rev A:
c     added call to flush (wrapper).
c     changed output format from e14.7 to e13.6 because last digit
c      was meaningless most of the time.
C
C 4.1 mods:
C     Density (average and rms) now calculated correctly in RUNMD and
C     passed through ENER(42). -- dap
C
C 4.2 mods:
C     IOUT7 < 0 can be used to do print ONLY to unit IOUT7 without rewinding
C     of this unit -- dap 9/95
c
c     EWALD: print out quick error estimate referenced in ener(43)
c
C
#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
#include "md.h"
#include "iewald.h"
      DIMENSION ENER(*),FAC(*)
c
c   ---dac change 2/90:  support bscale:
c
      dimension bscale(*)
C
C     ----- DEFINE VARIOUS TERMS -----
C
      ETOT   = ENER(1)
      EKTOT  = ENER(2)
      TEMP   = EKTOT/FAC(1)
      EKSOLT = ENER(3)/fac(2)
      if(ntt.eq.5) then 
          EKSOLV = ENER(4)/fac(3)
      else 
          eksolv = 0.0d0
      endif
C
      SCALP  = ENER(5)
      SCALS  = ENER(6)
C
      BOXX   = ENER(7)
      BOXY   = ENER(8)
      BOXZ   = ENER(9)
      VOLUME = ENER(10)
      DENSIT = ENER(42)
C
      PRESX  = ENER(11)
      PRESY  = ENER(12)
      PRESZ  = ENER(13)
      PRESS  = ENER(14)
C
      EKCMX  = ENER(15)
      EKCMY  = ENER(16)
      EKCMZ  = ENER(17)
      EKCMT  = ENER(18)
C
      VIRX   = ENER(19)
      VIRY   = ENER(20)
      VIRZ   = ENER(21)
      VIRT   = ENER(22)
C
      EPOT   = ENER(23)
      ENONB  = ENER(24)
      ENELE  = ENER(25)
      EHBOND = ENER(26)
      EBOND  = ENER(27)
      EANGLE = ENER(28)
      EDIHED = ENER(29)
      ENB14  = ENER(30)
      EEL14  = ENER(31)
      ECONST = ENER(32)
      EPOL   = ENER(33)
      aveper   = ENER(34)
      aveind   = ENER(35)
      avetot   = ENER(36)
      e3bod = ener(38)
      if (iewald .eq. 1) virvsene = ener(43)
C
C Write out info to standard output UNLESS IOUT7 < 0. In this case,
C ONLY output to the (absolute value of the) unit named in IOUT7
C and do not rewind IOUT7.
C
      IF (IOUT7.LT.0) GO TO 10
C
      WRITE(6,9018) NSTEP,TIME,TEMP,PRESS
      WRITE(6,9028) ETOT,EKTOT,EPOT
      WRITE(6,9038) EBOND,EANGLE,EDIHED
c Margaret
c      WRITE(6,9048) ENB14,EEL14,ENONB
      WRITE(6,9048) ENB14,ener(42),ENONB
      WRITE(6,9058) ENELE,EHBOND,ECONST
      if (econst.ne.0.0) WRITE(6,9076) epot-econst
      if (VOLUME.ne.0.0) WRITE(6,9078) EKCMT,VIRT,VOLUME
      if (ntt .eq. 5) 
     $        WRITE(6,9068) EKSOLT,EKSOLV
      if (EPOL.ne.0.0 .or. e3bod .ne.0.0) 
     $        WRITE(6,9070) EPOL,e3bod 
      if (EPOL.ne.0.0) WRITE(6,9088) aveper,aveIND,avetot
      if (iscale.ne.0) write(6,9069) (bscale(im),im=1,min(50,iscale))
      if (volume.ne.0.0) write(6,9079) DENSIT
      if (iewald .eq. 1) write(6,9188) virvsene
      write(6,8088)
c
c     --- flush i/o buffer ---
c
      call amflsh(6)
      if (iout7.eq.0) return
C
C       ----- OUTPUT THE INFO FILE if requested -----
C             Write to unit IOUT7
C
      rewind(iout7)
   10 IOUT = ABS(IOUT7)
      WRITE(iout,9018) NSTEP,TIME,TEMP,PRESS
      WRITE(iout,9028) ETOT,EKTOT,EPOT
      WRITE(iout,9038) EBOND,EANGLE,EDIHED
c      WRITE(iout,9048) ENB14,EEL14,ENONB
c Margaret
      WRITE(iout,9048) ENB14,ener(42),ENONB
      WRITE(iout,9058) ENELE,EHBOND,ECONST
      if (econst.ne.0.0.OR.IOUT7.LT.0) WRITE(iout,9076) epot-econst
      if (VOLUME.ne.0.0) WRITE(iout,9078) EKCMT,VIRT,VOLUME
      if (ntt .eq. 5) 
     $        WRITE(iout,9068) EKSOLT,EKSOLV
      if (EPOL.ne.0.0 .or. e3bod .ne.0.0) 
     $        WRITE(iout,9070) EPOL,e3bod 
      if (EPOL.ne.0.0) WRITE(iout,9088) aveper,aveind,avetot
      if (iscale.ne.0) 
     $        write(iout,9069) (bscale(im),im=1,min(50,iscale))
      if (volume.ne.0.0) write(iout,9079) DENSIT
      if (iewald .eq. 1) write(iout,9188) virvsene
C

 8088 format(t2,78('-'),/)
 9018 FORMAT(/1X, 'NSTEP =',I6,2X,'TIME(PS) =',F9.3,2X,
     .        'TEMP(K) =',F9.2,2x,'PRESS =',F10.2)
 9028 FORMAT (1X,'Etot   = ',f12.4,2X,'EKtot   = ',f12.4,2X,
     +        'EPtot      = ',f12.4)
 9038 FORMAT (1X,'BOND   = ',f12.4,2X,'ANGLE   = ',f12.4,2X,
     +        'DIHED      = ',f12.4)
 9048 FORMAT (1X,'1-4 NB = ',f12.4,2X,' Chiral = ',f12.4,2X,
     +        'VDWAALS    = ',f12.4)
c Margaret
c 9048 FORMAT (1X,'1-4 NB = ',f12.4,2X,'1-4 EEL = ',f12.4,2X,
c     +        'VDWAALS    = ',f12.4)
 9058 FORMAT (1X,'EELEC  = ',f12.4,2X,'EHBOND  = ',f12.4,2X,
     +        'CONSTRAINT = ',f12.4)
 9076 FORMAT (1X,'EAMBER (non-constraint) = ',f12.4)
 9078 FORMAT (1X,'EKCMT  = ',f12.4,2X,'VIRIAL  = ',f12.4,2X,
     +        'VOLUME     = ',f12.4)
 9079 format (48x,'Density    = ',f12.4)
c
 9068 FORMAT (1X,'T_SOLUTE =',F11.4,2X, 'T_SOLVENT =',F11.4)
 9070 FORMAT (1X,'EPOLZ  = ',F12.4,2X,'E3BODY  = ',F12.4)
c
 9069 format (1x,'SCALE  = ',f12.4,2x,'          ',f13.4,2x,
     +        '             ',f13.4)
 9088 FORMAT (1X,'DIPOLE MOMENTS/RESIDUE:',/
     $       ,1X,'PERMEN = ',F8.3,2X,'INDUCED = ',F8.3,
     $        2X,'VECTOR SUM =',F8.3)
 9188 FORMAT (1X,'Ewald error estimate: ', e12.4)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE PSCALE(NSPM,NSP,NR,RMU,X,TMA,amass,INDEX)
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
#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
C
C     ----- ROUTINE TO DO THE SCALING OF THE COORDINATES OF
C           THE SYSTEM IN THE CASE OF CONSTANT PRESSURE RUN -----
C
c     Mods for Rev A by gls:
c     - this was made a library routine, so the sense of INDEX
c       had to be made consistent between gibbs and md.  It was
c       changed to the gibbs sense where INDEX=0 means uniform coordinate
c       scaling.
c     - cpp switchable precision and DP constants
c
      DIMENSION NSP(*),RMU(*),X(3,*),amass(*),TMA(*)
      DIMENSION XCM(3),XCMS(3)
#ifdef CRAY_MP
      integer j1(nspm), jj1(nspm)
#endif
      data zero, one /0.0d0, 1.0d0/
C
C     ----- BRANCH DEPENDING ON THE THE VALUE OF INDEX -----
C
      IF(INDEX .le. 0) GO TO 200
C
C     ----- SUB MOLECULE CENTER OF MASS SCALING -----
C
      J = 0
      JJ = 0
#ifdef CRAY_MP
      j1(1) = 0
      jj1(1) = 0
      do i = 2, nspm
         j1(i)  = nsp(i-1) +  j1(i-1)
         jj1(i) = nsp(i-1) + jj1(i-1)
      enddo
cmic$ parallel autoscope
cmic$*private(xcm,xcms)
cmic$*shared(x)
cmic$ do parallel
#endif
c
      DO 100 I = 1,NSPM
        NN = NSP(I)
        TMN = one/TMA(I)
        XCM(1) = zero
        XCM(2) = zero
        XCM(3) = zero
        DO K = 1,NN
#ifdef CRAY_MP
          j = j1(i) + k
#else
          J = J+1
#endif
          aamass = amass(J)
          XCM(1) = XCM(1) + X(1,J) * aamass
          XCM(2) = XCM(2) + X(2,J) * aamass
          XCM(3) = XCM(3) + X(3,J) * aamass
        enddo
        XCM(1) = XCM(1) * TMN
        XCM(2) = XCM(2) * TMN
        XCM(3) = XCM(3) * TMN
        XCMS(1) = XCM(1)*(RMU(1)-one)
        XCMS(2) = XCM(2)*(RMU(2)-one)
        XCMS(3) = XCM(3)*(RMU(3)-one)
C
C       ----- MOVE THE MOLECULE TO THE NEW CENTER OF MASS -----
C
        DO K = 1,NN
#ifdef CRAY_MP
          jj = jj1(i) + k
#else
          JJ = JJ+1
#endif
          X(1,JJ) = X(1,JJ)+XCMS(1)
          X(2,JJ) = X(2,JJ)+XCMS(2)
          X(3,JJ) = X(3,JJ)+XCMS(3)
        enddo
  100 continue
#ifdef CRAY_MP
cmic$ end do
cmic$ end parallel
#endif
      RETURN
C
C     ----- UNIFORM COORDINATES SCALING -----
C
  200 CONTINUE
      DO 220 I = 1,NR
        X(1,I) = X(1,I)*RMU(1)
        X(2,I) = X(2,I)*RMU(2)
        X(3,I) = X(3,I)*RMU(3)
  220 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE SETVEL(NRP,NR,NTU,NTX,X,V,F,WINV,
     +                  TEMPI,HEAT,DT,INIT,IG,iscale,scalm)
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
c     sets up initial velocities.
c     Mods for Rev A by gls:
c
c     - cpp selectable single/double prec.
c     - changed to DP constants
#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
      DIMENSION X(*),V(*),WINV(*),F(*)
C
      NR3 = 3*NR
      IF(NTU.NE.2 ) GOTO 131
      CRX = 0.1d0
      DO 110 I = 1,NR3+iscale
  110  X(I) = X(I)*CRX
      IF(NTX-3) 141,121,111
  111 CRV =  SQRT(4.184d0)
      DO 120 I = 1,NR3+iscale
  120   V(I) = V(I)*CRV
      GOTO 141
  121 CONTINUE
      DO 130 I = 1,NR3+iscale
  130   F(I) = F(I)*CRX
C
C     ----- CALCULATE THE VELOCITIES FROM THE POSITONS,
C           WHEN REQUIRED -----
C
  131 IF(NTX.NE.3) GOTO 141
      dtinv = 1.0d0 / dt
      DO 140 I = 1,NR3+iscale
  140   V(I) = (F(I)-X(I)) * dtinv
C
C     ----- TAKE THE VELOCITIES FROM A MAXWELLIAN,
C           WHEN REQUIRED -----
C
  141 IF(TEMPI.LT.1.d-6) GOTO 171
      BOLTZ = 8.31441d-3
      IF(NTU.EQ.1) BOLTZ = BOLTZ/4.184d0
      BOLTZ = BOLTZ*TEMPI
      I = 0
c
c   dac change, 1.18.91: apply tempi to all atoms, not just
c     the solute, making this more compatible with previous
c     versions of amber:
c
      do j=1,nr
        SD =  SQRT(BOLTZ*WINV(J))
        DO M = 1,3
          CALL GAUSS(0.d0,SD,VG)
          I = I+1
          V(I) = VG
        enddo
      enddo
      if (iscale.gt.0) then
        SD =  SQRT(BOLTZ/SCALM)
        do 164 j=1,iscale
          CALL GAUSS(0.d0,SD,VG)
          I = I+1
          V(I) = VG
  164   CONTINUE
      end if
      IF(INIT.EQ.4) INIT = 3
C
C     ----- SCALE VELOCITIES, WHEN REQUIRED -----
C
  171 IF( ABS(HEAT).LT.1.d-6) GOTO 181
      DO 175 I = 1,NR3
  175   V(I) = V(I)*HEAT
      IF(INIT.EQ.4) INIT = 3
  181 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE STOPCM(NR,X,V,XCM,VCM,OCM)
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
C
C     ----- ROUTINE TO STOP THE TRANSLATIONAL MOTION -----
C
#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
#ifdef MPI
# include "parallel.h"
#endif
#include "extra.h"
      DIMENSION X(*),V(*),XCM(*),VCM(*),OCM(*)
C
C     ----- STOP THE CENTER OF MASS TRANSLATION -----
C
      I = 0
      DO J = 1,NR
        DO M = 1,3
          I = I+1
          V(I) = V(I)-VCM(M)
        enddo
      enddo
C
C     ----- STOP THE ROTATION ABOUT THE CENTER OF MASS -----
C
      I = 0
      DO J = 1,NR
        X1 = X(I+1)-XCM(1)
        X2 = X(I+2)-XCM(2)
        X3 = X(I+3)-XCM(3)
        V(I+1) = V(I+1)-OCM(2)*X3+OCM(3)*X2
        V(I+2) = V(I+2)-OCM(3)*X1+OCM(1)*X3
        V(I+3) = V(I+3)-OCM(1)*X2+OCM(2)*X1
        I = I+3
      enddo
      if (master) WRITE(6,9008)
 9008 FORMAT(/5X,'TRANSLATIONAL AND ROTATIONAL MOTION REMOVED')
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE MDENG(nf,NSTEP,TIME,ENER,FAC)
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
#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
#include "box.h"
      DIMENSION ENER(*),FAC(*)
      logical first
      character*16 labs(40)
      save first
      save labs
      data first/.true./
      data labs/'Nsteps  ','time(ps)  ','Etot  ','EKinetic  ',
     .          'Temp  ','T_solute ','T_solv  ','Pres_scal_solu ',
     .          'Pres_scal_solv ','BoxX  ','BoxY  ','BoxZ  ',
     .          'volume  ','pres_X  ','pres_Y  ','pres_Z  ',
     .          'Pressure ','EKCoM_x ','EKCoM_y ','EKCoM_z',
     .          'EKComTot ','VIRIAL_x ','VIRIAL_y ','VIRIAL_z ',
     .          'VIRIAL_tot ','E_pot  ','E_vdw+Echi  ','E_el  ',
     .          'E_hbon  ','E_bon  ','E_angle  ','E_dih  ',
     .          'E_14vdw  ','E_14el  ','E_const  ','E_pol  ',
     .          'AV_permMoment ','AV_indMoment ','AV_totMoment ',
c Margaret
c     .                                            'Density'/
     .                                            'Chiral'/
C
C     ----- DEFINE VARIOUS TERMS -----
C
      if (first) then
c       -- up to Ekinetic:
        write(nf,1) 'L0 ', (labs(i),i=1,4)
c       -- up to Pres_scal_solu:
        write(nf,1) 'L1 ', (labs(i),i=5,8)
c       -- up to boxZ:
        write(nf,1) 'L2 ', (labs(i),i=9,12)
c       -- up to pres_Z:
        write(nf,1) 'L3 ', (labs(i),i=13,16)
c       -- up to EKCoM_z:
        write(nf,1) 'L4 ', (labs(i),i=17,20)
c       -- up to VIRIAL_z:
        write(nf,1) 'L5 ', (labs(i),i=21,24)
c       -- up to E_el:
        write(nf,1) 'L6 ', (labs(i),i=25,28)
c       -- up to E_dih:
        write(nf,1) 'L7 ', (labs(i),i=29,32)
c       -- up to E_pol:
        write(nf,1) 'L8 ', (labs(i),i=33,36)
c       -- up to Density:
        write(nf,1) 'L9 ', (labs(i),i=37,40)
    1   format(a,10(1x,a))
        first = .false.
      endif
c
c     ----- write values for this step -----
c
c     -- up to Ekinetic:
      write(nf, 2) 'L0 ', nstep, time, ener(1), ener(2)
c     -- up to Pres_scal_solu:
      write(nf, 3) 'L1 ', ENER(2)/FAC(1), ENER(3)/fac(2), 
     .                   ENER(4)/fac(3), ener(5)
c     -- up to boxZ:
      write(nf, 3) 'L2 ', ener(6), box(1), box(2), box(3)
c     -- up to pres_Z:
      write(nf, 3) 'L3 ', (ener(i), i=10,13)
c     -- up to EKCoM_z:
      write(nf, 3) 'L4 ', (ener(i), i=14,17)
c     -- up to VIRIAL_z:
      write(nf, 3) 'L5 ', (ener(i), i=18,21)
c     -- up to E_el:
      write(nf, 3) 'L6 ', (ener(i), i=22,25)
c     -- up to E_dih:
      write(nf, 3) 'L7 ', (ener(i), i=26,29)
c      -- up to E_pol:
      write(nf, 3) 'L8 ', (ener(i), i=30,33)
c     -- up to Density:
c Margaret
      write(nf, 3) 'L9 ', (ener(i), i=34,36), ener(42)
c
    2 format(a, i8, 20(2x,e16.10))
    3 format(a, 20(e16.10,2x))
      RETURN
      END


      SUBROUTINE MDEN2(nf,NSTEP,TIME,ENER,FAC)
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
#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
      DIMENSION ENER(*),FAC(*)
C
C     ----- DEFINE VARIOUS TERMS -----
C
      write(nf,'(''Nsteps = '',i10)')nstep
      write(nf,'(''Time               = '',f20.10,'' p.s. '')')time
      write(nf,'(''Etotal             = '',f20.10)')ENER(1)
      write(nf,'(''EKinetic           = '',f20.10)')ENER(2)
      write(nf,'(''Temperature        = '',f20.10)')ENER(2)/FAC(1)
      write(nf,'(''Temp_solute        = '',f20.10)')ENER(3)/fac(2)
      write(nf,'(''Temp_solvent       = '',f20.10)')ENER(4)/fac(3)
      write(nf,'(''Press_SCALE_solute = '',f20.10)')ENER(5)
      write(nf,'(''Press_SCALE_solvent= '',f20.10)')ENER(6)
      write(nf,'(''BOX_x              = '',f20.10)')ENER(7)
      write(nf,'(''BOX_y              = '',f20.10)')ENER(8)
      write(nf,'(''BOX_z              = '',f20.10)')ENER(9)
      write(nf,'(''VOLUME             = '',f20.10)')ENER(10)
      write(nf,'(''PRES_x             = '',f20.10)')ENER(11)
      write(nf,'(''PRES_y             = '',f20.10)')ENER(12)
      write(nf,'(''PRES_z             = '',f20.10)')ENER(13)
      write(nf,'(''PRESSURE           = '',f20.10)')ENER(14)
      write(nf,'(''EKCoM_x            = '',f20.10)')ENER(15)
      write(nf,'(''EKCoM_y            = '',f20.10)')ENER(16)
      write(nf,'(''EKCoM_z            = '',f20.10)')ENER(17)
      write(nf,'(''EKCoMTotal         = '',f20.10)')ENER(18)
      write(nf,'(''VIRIAL_x           = '',f20.10)')ENER(19)
      write(nf,'(''VIRIAL_y           = '',f20.10)')ENER(20)
      write(nf,'(''VIRIAL_z           = '',f20.10)')ENER(21)
      write(nf,'(''VIRIAL_Total       = '',f20.10)')ENER(22)
      write(nf,'(''Epotential         = '',f20.10)')ENER(23)
      write(nf,'(''vanderwaals        = '',f20.10)')ENER(24)
      write(nf,'(''electrostatic      = '',f20.10)')ENER(25)
      write(nf,'(''h-bond             = '',f20.10)')ENER(26)
      write(nf,'(''bond               = '',f20.10)')ENER(27)
      write(nf,'(''angle              = '',f20.10)')ENER(28)
      write(nf,'(''dihedral           = '',f20.10)')ENER(29)
      write(nf,'(''1-4v.d.w.          = '',f20.10)')ENER(30)
      write(nf,'(''1-4electrostatic   = '',f20.10)')ENER(31)
      write(nf,'(''constraint         = '',f20.10)')ENER(32)
      write(nf,'(''polarization       = '',f20.10)')ENER(33)
      write(nf,'(''ave perm moment    = '',f20.10)')ENER(34)
      write(nf,'(''ave ind  moment    = '',f20.10)')ENER(35)
      write(nf,'(''ave total moment   = '',f20.10)')ENER(36)
      write(nf,'(''density            = '',f20.10)')ENER(42)
      write(nf,'(''-----------------------------------------'')')
c
      RETURN
      END
