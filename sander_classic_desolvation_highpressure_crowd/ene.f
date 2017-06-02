c-----------------------------------------------------------------------
      SUBROUTINE BOND(NB,IB,JB,ICB,X,F,EB,NOCRST)
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
#ifdef MPI
#  include "parallel.h"
#endif
#include "nmr.h"
      LOGICAL SKIP,NOCRST
C
C     ----- ROUTINE TO GET BOND ENERGY AND FORCES FOR THE POTENTIAL
C           OF CB*(B-B0)**2
C
#include "box.h"
c Margaret
#include "HB.h"
#include "LAN.h"
#include "parms.h"
      DIMENSION XIJ(190),YIJ(190),ZIJ(190),XKJ(190),YKJ(190),
     + ZKJ(190),CST(190),SNT(190),EAW(190),RIJ(190),RKJ(190),RIK(190),
     + DFW(190),ANT(190),EPW(190)
      DIMENSION IB(*),JB(*),ICB(*),X(*),F(*)
#ifdef MPI
      integer piece,start,end,nbtmp,newnb

      piece = nb/numtasks
      start = mytaskid*piece+1
      end   = mytaskid*piece+piece
      if(mytaskid .EQ. (numtasks-1)) end = nb
      nbtmp = nb
c     JV Use actual count this PE will do, reset at end of routine
      nb = end
      IST = start -1
#else
      IST = 0
#endif
C
      EBL = 0.0D+00
C
C     ----- GRAND LOOP FOR THE bond STUFF -----
C
 4200 CONTINUE
        MAXLEN = 190
        SKIP = (IST+MAXLEN).GT.NB
        IF (SKIP) MAXLEN = NB-IST
        IF (MAXLEN.gt.0) then
C
          DO 100 JN = 1,MAXLEN
            I3 = IB(JN+IST)
            J3 = JB(JN+IST)
C
C           ----- CALCULATION OF THE bond vector -----
C
            XIJ(JN) = X(I3+1)-X(J3+1)
            YIJ(JN) = X(I3+2)-X(J3+2)
            ZIJ(JN) = X(I3+3)-X(J3+3)
  100     CONTINUE
          IF(NTB.EQ.0 .OR. NOCRST) GOTO 120
C
C           ----- APPLY PERIODIC BOUNDARY CONDITION -----
C
            CALL PIMAG(MAXLEN,XIJ,YIJ,ZIJ)
            IF(NTB.GT.0) GO TO 120
            CALL IMAGT(MAXLEN,XIJ,YIJ,ZIJ)
  120     CONTINUE
C
          DO 160 JN = 1,MAXLEN
            RIJ0 = XIJ(JN)*XIJ(JN)+YIJ(JN)*YIJ(JN)+ZIJ(JN)*ZIJ(JN)
            RIJ(JN) = SQRT(RIJ0)
  160     CONTINUE
C
C         ----- CALCULATION OF THE ENERGY AND DER -----
C
          DO 180 JN = 1,MAXLEN
            IC = ICB(JN+IST)
            RIJ0 = RIJ(JN)
            DA = RIJ0-REQ(IC)
c                                 for rms deviation from ideal bonds:
            ebdev = ebdev + da*da
            DF = RK(IC)*DA
            EAW(JN) = DF*DA
            DFW(JN) = (DF+DF)/RIJ0
  180     CONTINUE
C
C         ----- CALCULATION OF THE FORCE -----
C
          DO 200 JN = 1,MAXLEN
            I3 = IB(JN+IST)
            J3 = JB(JN+IST)
            DF = DFW(JN)
            XA = DF*XIJ(JN)
            YA = DF*YIJ(JN)
            ZA = DF*ZIJ(JN)
            F(I3+1) = F(I3+1)-XA
            F(I3+2) = F(I3+2)-YA
            F(I3+3) = F(I3+3)-ZA
            F(J3+1) = F(J3+1)+XA
            F(J3+2) = F(J3+2)+YA
            F(J3+3) = F(J3+3)+ZA
  200     CONTINUE
          do 201 ksum = 1, maxlen
            ebl = ebl + eaw(ksum)
  201     continue
*         EBL = EBL+SSUM(MAXLEN,EAW(1),1)
          IST = IST+MAXLEN
        endif
      IF(.NOT.SKIP) GO TO 4200
C
C     ----- ALL DONE -----
C
      EB = EBL
#ifdef MPI
c     JV Set NB back to correct value
      nb = nbtmp
#endif
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE ANGL(NBA,IT,JT,KT,ICT,X,F,EBA,MBA,ECN,NOCRST)
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
#ifdef MPI
#include "parallel.h"
#endif
      LOGICAL SKIP,NOCRST
C
C     ----- ROUTINE TO GET THE BOND ENERGIES AND FORCES FOR THE
C           POTENTIAL OF THE TYPE CT*(T-T0)**2
C
#include "box.h"
#include "parms.h"
#include "nmr.h"
#ifdef CHARMM
      dimension ebw(190)
#endif
      DIMENSION XIJ(190),YIJ(190),ZIJ(190),XKJ(190),YKJ(190),
     + ZKJ(190),CST(190),SNT(190),EAW(190),RIJ(190),RKJ(190),RIK(190),
     + DFW(190),ANT(190),EPW(190)
      DIMENSION IT(*),JT(*),KT(*),ICT(*),X(*),F(*)
#ifdef MPI
      integer piece,start,end,nbtmp,newnb
#endif
      data pt999 /0.9990d0/
c
#ifdef MPI
      piece = nba/numtasks
      start = mytaskid*piece+1
      end   = mytaskid*piece+piece
      if(mytaskid .EQ. (numtasks-1)) end = nba
      nbtmp = nba
c     JV Use actual count this PE will do, reset at end of routine
      nba = end
      IST = start -1
#else
      IST = 0
#endif
C
      EBAL = 0.0
      ECNL = 0.0
      EUB  = 0.0
C
C     ----- GRAND LOOP FOR THE angle STUFF -----
C
      ISTC = MBA
 4200 CONTINUE
        MAXLEN = 190
        SKIP = (IST+MAXLEN).GT.NBA
        IF (SKIP) MAXLEN = NBA-IST
        IF (MAXLEN.LE.0) GO TO 220
C
          DO 100 JN = 1,MAXLEN
            I3 = IT(JN+IST)
            J3 = JT(JN+IST)
            K3 = KT(JN+IST)
C
C           ----- CALCULATION OF THE angle -----
C
            XIJ(JN) = X(I3+1)-X(J3+1)
            YIJ(JN) = X(I3+2)-X(J3+2)
            ZIJ(JN) = X(I3+3)-X(J3+3)
            XKJ(JN) = X(K3+1)-X(J3+1)
            YKJ(JN) = X(K3+2)-X(J3+2)
            ZKJ(JN) = X(K3+3)-X(J3+3)
  100     CONTINUE
          IF(NTB.NE.0 .AND. .NOT.NOCRST) THEN
C
C           ----- APPLY PERIODIC BOUNDARY CONDITION -----
C
              CALL PIMAG(MAXLEN,XIJ,YIJ,ZIJ)
              CALL PIMAG(MAXLEN,XKJ,YKJ,ZKJ)
C
              IF(NTB.LE.0) THEN
                CALL IMAGT(MAXLEN,XIJ,YIJ,ZIJ)
                CALL IMAGT(MAXLEN,XKJ,YKJ,ZKJ)
              END IF
          END IF
C
          DO 160 JN = 1,MAXLEN
            RIJ0 = XIJ(JN)*XIJ(JN)+YIJ(JN)*YIJ(JN)+ZIJ(JN)*ZIJ(JN)
            RKJ0 = XKJ(JN)*XKJ(JN)+YKJ(JN)*YKJ(JN)+ZKJ(JN)*ZKJ(JN)
            RIK0 = SQRT(RIJ0*RKJ0)
            CT0 = (XIJ(JN)*XKJ(JN)+YIJ(JN)*YKJ(JN)+ZIJ(JN)*ZKJ(JN))/RIK0
            CT1 = MAX(-pt999,CT0)
            CT2 = MIN(pt999,CT1)
            CST(JN) = CT2
            ANT(JN) = ACOS(CT2)
            RIJ(JN) = RIJ0
            RKJ(JN) = RKJ0
            RIK(JN) = RIK0
  160     CONTINUE
C
C         ----- CALCULATION OF THE ENERGY AND DER -----
C
          DO 180 JN = 1,MAXLEN
            IC = ICT(JN+IST)
            ANT0 = ANT(JN)
            DA = ANT0-TEQ(IC)
c                                 for rms deviation from ideal angles:
            eadev = eadev + da*da
            DF = TK(IC)*DA
            EAW(JN) = DF*DA
            DFW(JN) = -(DF+DF)/SIN(ANT0)
  180     CONTINUE
C
C         ----- CALCULATION OF THE FORCE -----
C
          DO 200 JN = 1,MAXLEN
            I3 = IT(JN+IST)
            J3 = JT(JN+IST)
            K3 = KT(JN+IST)
C
            ST = DFW(JN)
            STH = ST*CST(JN)
            CIK = ST/RIK(JN)
            CII = STH/RIJ(JN)
            CKK = STH/RKJ(JN)
            DT1 = CIK*XKJ(JN)-CII*XIJ(JN)
            DT2 = CIK*YKJ(JN)-CII*YIJ(JN)
            DT3 = CIK*ZKJ(JN)-CII*ZIJ(JN)
            DT7 = CIK*XIJ(JN)-CKK*XKJ(JN)
            DT8 = CIK*YIJ(JN)-CKK*YKJ(JN)
            DT9 = CIK*ZIJ(JN)-CKK*ZKJ(JN)
            DT4 = -DT1-DT7
            DT5 = -DT2-DT8
            DT6 = -DT3-DT9
C
            F(I3+1) = F(I3+1)-DT1
            F(I3+2) = F(I3+2)-DT2
            F(I3+3) = F(I3+3)-DT3
            F(J3+1) = F(J3+1)-DT4
            F(J3+2) = F(J3+2)-DT5
            F(J3+3) = F(J3+3)-DT6
            F(K3+1) = F(K3+1)-DT7
            F(K3+2) = F(K3+2)-DT8
            F(K3+3) = F(K3+3)-DT9
  200     CONTINUE
#ifdef CHARMM
c
c         --- include Urey-Bradley terms:
c
          DO 300 JN = 1,MAXLEN
            I3 = IT(JN+IST)
            J3 = KT(JN+IST)
C
C           ----- CALCULATION OF THE bond vector -----
C
            XIJ(JN) = X(I3+1)-X(J3+1)
            YIJ(JN) = X(I3+2)-X(J3+2)
            ZIJ(JN) = X(I3+3)-X(J3+3)
  300     CONTINUE
          IF(NTB.NE.0 .and. .not.NOCRST) then
C
C           ----- APPLY PERIODIC BOUNDARY CONDITION -----
C
            CALL PIMAG(MAXLEN,XIJ,YIJ,ZIJ)
            IF(NTB.LE.0) CALL IMAGT(MAXLEN,XIJ,YIJ,ZIJ)
          END IF
C
          DO 360 JN = 1,MAXLEN
            RIJ0 = XIJ(JN)*XIJ(JN)+YIJ(JN)*YIJ(JN)+ZIJ(JN)*ZIJ(JN)
            RIJ(JN) = SQRT(RIJ0)
  360     CONTINUE
C
C         ----- CALCULATION OF THE ENERGY AND DER -----
C
          DO 380 JN = 1,MAXLEN
            IC = ICT(JN+IST)
            RIJ0 = RIJ(JN)
            DA = RIJ0-RUB(IC)
            DF = RKUB(IC)*DA
            EBW(JN) = DF*DA
            DFW(JN) = (DF+DF)/RIJ0
  380     CONTINUE
C
C         ----- CALCULATION OF THE FORCE -----
C
          DO 400 JN = 1,MAXLEN
            I3 = IT(JN+IST)
            J3 = KT(JN+IST)
            DF = DFW(JN)
            XA = DF*XIJ(JN)
            YA = DF*YIJ(JN)
            ZA = DF*ZIJ(JN)
            F(I3+1) = F(I3+1)-XA
            F(I3+2) = F(I3+2)-YA
            F(I3+3) = F(I3+3)-ZA
            F(J3+1) = F(J3+1)+XA
            F(J3+2) = F(J3+2)+YA
            F(J3+3) = F(J3+3)+ZA
  400     CONTINUE
          do 401 ksum = 1, maxlen
           eub = eub + ebw(ksum)
  401     continue
*         EUB  = EUB +SSUM(MAXLEN,EBW(1),1)
#endif
          do 201 ksum=1, maxlen
           ebal = ebal + eaw(ksum)
  201     continue
*         EBAL = EBAL+SSUM(MAXLEN,EAW(1),1)
          IST = IST+MAXLEN
          IF(IST.LE.ISTC) GO TO 220
          LENC = IST-ISTC
          do 202 ksum = maxlen-lenc+1, maxlen
           ecnl = ecnl + eaw(ksum)
  202     continue
*         ECNL = ECNL+SSUM(LENC,EAW(MAXLEN-LENC+1),1)
          ISTC = ISTC+LENC
  220   CONTINUE
      IF(.NOT.SKIP) GO TO 4200
      EBA = EBAL
#ifdef CHARMM
      write(6,*) 'Urey-Bradley energy: ',eub
      eba = eba + eub
#endif
      ECN = ECNL
#ifdef MPI
      nba = nbtmp
#endif
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE EPHI(NPHI,IP,JP,KP,LP,ICP,CG,IAC,X,F,EP,ENBP,EELP,
     +                MPHI,ECN,NOCRST)
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
C  Mods for 4.1
C  (tec3) Deleted all the equivalence statements that were inhibiting 
c  parallelism in the shared memory implementation
C
#ifdef DPREC
      implicit double precision (a-h,o-z)
#endif
#ifdef MPI
#include "parallel.h"
#endif
#ifdef LES
#  include "les.h"
#endif

c     Mods for Rev A by GLS:
c     removed "getlen" nonsense
c     cpp selectable single/double prec, DP constants
c     3-Dec-88: made dual cray/generic source.  
c     Use cpp -DCRAYFISH for Cray source.
c     removed ssum calls, newen.
      LOGICAL DIELD,SKIP,NOCRST
C
#include "parms.h"
#include "box.h"
c Margaret
#include "HB.h"
#include "md.h"
      DIMENSION XIJ(190),YIJ(190),ZIJ(190),XKJ(190),YKJ(190),
     + ZKJ(190),XKL(190),YKL(190),ZKL(190),DX(190),DY(190),DZ(190),
     + GX(190),GY(190),GZ(190),CT(190),CPHI(190),SPHI(190),Z1(190),
     + Z2(190),FXI(190),FYI(190),FZI(190),FXJ(190),FYJ(190),FZJ(190),
     + FXK(190),FYK(190),FZK(190),FXL(190),FYL(190),FZL(190),DF(190)
C
      DIMENSION EPW(190)
      DIMENSION GMUL(10)
C
      DIMENSION IP(*),JP(*),KP(*),LP(*),ICP(*)
      DIMENSION CG(*),IAC(*),X(*),F(*)
C
      DATA GMUL/0.0d+00,2.0d+00,0.0d+00,4.0d+00,0.0d+00,6.0d+00,
     +          0.0d+00,8.0d+00,0.0d+00,10.0d+00/
      DATA TM24,TM06,tenm3/1.0d-18,1.0d-06,1.0d-03/
      DATA PI/3.141592653589793d+00/
      data zero,one,four,six,twelve/0.0d0,1.0d0,4.0d0,6.0d0,12.0d0/
C
C     ---- ARRAYS GAMC = PK*COS(PHASE) AND GAMS = PK*SIN(PHASE) ----
C
#ifdef MPI
      integer piece,start,end,nbtmp,newnb

      piece = NPHI/numtasks
      start = mytaskid*piece+1
      end   = mytaskid*piece+piece
      if(mytaskid .EQ. (numtasks-1)) end = NPHI
      nbtmp = NPHI
c     JV Use actual count this PE will do, reset at end of routine
      NPHI = end
      IST = start -1
#else
      IST = 0
#endif
      EPL = zero
      ECNL = zero
      ENBPL = zero
      EELPL = zero
      SCNB0 = one/SCNB
      SCEE0 = one/SCEE
      DIELD = IDIEL.LE.0
C
C     ----- GRAND LOOP FOR THE DIHEDRAL STUFF -----
C
      ISTC = MPHI
 4200 CONTINUE
        MAXLEN = 190
        SKIP = (IST+MAXLEN).GT.NPHI
        IF(SKIP) MAXLEN = NPHI-IST
        IF(MAXLEN.LE.0) GO TO 820 
C
          DO 100 JN = 1,MAXLEN
            I3 = IP(JN+IST)
            J3 = JP(JN+IST)
            K3T = KP(JN+IST)
            L3T = LP(JN+IST)
            K3 = IABS(K3T)
            L3 = IABS(L3T)
C
C           ----- CALCULATION OF ij, kj, kl VECTORS -----
C
            XIJ(JN) = X(I3+1)-X(J3+1)
            YIJ(JN) = X(I3+2)-X(J3+2)
            ZIJ(JN) = X(I3+3)-X(J3+3)
            XKJ(JN) = X(K3+1)-X(J3+1)
            YKJ(JN) = X(K3+2)-X(J3+2)
            ZKJ(JN) = X(K3+3)-X(J3+3)
            XKL(JN) = X(K3+1)-X(L3+1)
            YKL(JN) = X(K3+2)-X(L3+2)
            ZKL(JN) = X(K3+3)-X(L3+3)                                   
  100     CONTINUE
          IF(NTB.EQ.0 .OR. NOCRST) GOTO 120
C
C           ----- APPLY PERIODIC BOUNDARY CONDITION -----
C
            CALL PIMAG(MAXLEN,XIJ,YIJ,ZIJ)
            CALL PIMAG(MAXLEN,XKJ,YKJ,ZKJ)
            CALL PIMAG(MAXLEN,XKL,YKL,ZKL)
C
            IF(NTB.GT.0) GO TO 120
              CALL IMAGT(MAXLEN,XIJ,YIJ,ZIJ)
              CALL IMAGT(MAXLEN,XKJ,YKJ,ZKJ)
              CALL IMAGT(MAXLEN,XKL,YKL,ZKL)
  120     CONTINUE
C
C         ----- GET THE NORMAL VECTOR -----
C
          DO 140 JN = 1,MAXLEN
            DX(JN) = YIJ(JN)*ZKJ(JN)-ZIJ(JN)*YKJ(JN)
            DY(JN) = ZIJ(JN)*XKJ(JN)-XIJ(JN)*ZKJ(JN)
            DZ(JN) = XIJ(JN)*YKJ(JN)-YIJ(JN)*XKJ(JN)
            GX(JN) = ZKJ(JN)*YKL(JN)-YKJ(JN)*ZKL(JN)
            GY(JN) = XKJ(JN)*ZKL(JN)-ZKJ(JN)*XKL(JN)
            GZ(JN) = YKJ(JN)*XKL(JN)-XKJ(JN)*YKL(JN)
  140     CONTINUE
C
          DO 160 JN = 1,MAXLEN
            FXI(JN) = SQRT(DX(JN)*DX(JN)
     .                    +DY(JN)*DY(JN)
     .                    +DZ(JN)*DZ(JN)+TM24)
            FYI(JN) = SQRT(GX(JN)*GX(JN)
     .                    +GY(JN)*GY(JN)
     .                    +GZ(JN)*GZ(JN)+TM24)
            CT(JN) = DX(JN)*GX(JN)+DY(JN)*GY(JN)+DZ(JN)*GZ(JN)
  160     CONTINUE
C
C         ----- BRANCH IF LINEAR DIHEDRAL -----
C                             
          DO 180 JN = 1,MAXLEN
#ifdef CRAYFISH
            BIT = one/FXI(JN)
            BIK = one/FYI(JN)
            Z10 = CVMGT(zero,BIT,tenm3.GT.FXI(JN))
            Z20 = CVMGT(zero,BIK,tenm3.GT.FYI(JN))
#else
            z10 = one/FXI(jn)
            z20 = one/FYI(jn)
            if (tenm3 .gt. FXI(jn)) z10 = zero
            if (tenm3 .gt. FYI(jn)) z20 = zero
#endif
            Z12 = Z10*Z20
            Z1(JN) = Z10
            Z2(JN) = Z20
#ifdef CRAYFISH
            FTEM = CVMGZ(zero,one,Z12)
#else
            ftem = zero
            if (z12 .ne. zero) ftem = one
#endif
            FZI(JN) = FTEM
cd          CT0 = AMIN1(one,CT(JN)*Z12)
            CT0 = MIN(one,CT(JN)*Z12)
cd          CT1 = AMAX1(-one,CT0)
            CT1 = MAX(-one,CT0)
            S = XKJ(JN)*(DZ(JN)*GY(JN)-DY(JN)*GZ(JN))+
     +          YKJ(JN)*(DX(JN)*GZ(JN)-DZ(JN)*GX(JN))+
     +          ZKJ(JN)*(DY(JN)*GX(JN)-DX(JN)*GY(JN))
            AP0 = ACOS(CT1)
            AP1 = PI-SIGN(AP0,S)
            CT(JN) = AP1
            CPHI(JN) = COS(AP1)
            SPHI(JN) = SIN(AP1)
  180     CONTINUE
C
C         ----- CALCULATE THE ENERGY AND THE DERIVATIVES WITH RESPECT TO
C               COSPHI -----
C
          DO 200 JN = 1,MAXLEN
            IC = ICP(JN+IST)
            INC = IPN(IC)
            CT0 = PN(IC)*CT(JN)
            COSNP = COS(CT0)
            SINNP = SIN(CT0)
            EPW(JN) = (PK(IC)+COSNP*GAMC(IC)+SINNP*GAMS(IC))*FZI(JN)
            DF0 = PN(IC)*(GAMC(IC)*SINNP-GAMS(IC)*COSNP)
            DUMS = SPHI(JN)+SIGN(TM24,SPHI(JN))
            DFLIM = GAMC(IC)*(PN(IC)-GMUL(INC)+GMUL(INC)*CPHI(JN))
#ifdef CRAYFISH
            DF1 = CVMGT(DFLIM,DF0/DUMS,TM06.GT.ABS(DUMS))
#else
            df1 = df0/dums
            if(tm06.gt.abs(dums)) df1 = dflim
#endif
            DF(JN) = DF1*FZI(JN)
  200     CONTINUE          
C                                     
C         ----- NOW DO TORSIONAL FIRST DERIVATIVES -----
C
          DO 220 JN = 1,MAXLEN
C
C           ----- NOW, SET UP ARRAY DC = FIRST DER. OF COSPHI W/RESPECT
C                 TO THE CARTESIAN DIFFERENCES T -----
C
            Z11 = Z1(JN)*Z1(JN)
            Z12 = Z1(JN)*Z2(JN)
            Z22 = Z2(JN)*Z2(JN)
            DC1 = -GX(JN)*Z12-CPHI(JN)*DX(JN)*Z11
            DC2 = -GY(JN)*Z12-CPHI(JN)*DY(JN)*Z11
            DC3 = -GZ(JN)*Z12-CPHI(JN)*DZ(JN)*Z11
            DC4 =  DX(JN)*Z12+CPHI(JN)*GX(JN)*Z22
            DC5 =  DY(JN)*Z12+CPHI(JN)*GY(JN)*Z22
            DC6 =  DZ(JN)*Z12+CPHI(JN)*GZ(JN)*Z22
C
C           ----- UPDATE THE FIRST DERIVATIVE ARRAY -----
C
            DR1 = DF(JN)*(DC3*YKJ(JN) - DC2*ZKJ(JN))
            DR2 = DF(JN)*(DC1*ZKJ(JN) - DC3*XKJ(JN))
            DR3 = DF(JN)*(DC2*XKJ(JN) - DC1*YKJ(JN))
            DR4 = DF(JN)*(DC6*YKJ(JN) - DC5*ZKJ(JN))
            DR5 = DF(JN)*(DC4*ZKJ(JN) - DC6*XKJ(JN))
            DR6 = DF(JN)*(DC5*XKJ(JN) - DC4*YKJ(JN))
            DRX = DF(JN)*(-DC2*ZIJ(JN) + DC3*YIJ(JN) +
     +               DC5*ZKL(JN) - DC6*YKL(JN))
            DRY = DF(JN)*( DC1*ZIJ(JN) - DC3*XIJ(JN) -
     +               DC4*ZKL(JN) + DC6*XKL(JN))
            DRZ = DF(JN)*(-DC1*YIJ(JN) + DC2*XIJ(JN) +
     +               DC4*YKL(JN) - DC5*XKL(JN))
            FXI(JN) = - DR1
            FYI(JN) = - DR2
            FZI(JN) = - DR3
            FXJ(JN) = - DRX + DR1
            FYJ(JN) = - DRY + DR2
            FZJ(JN) = - DRZ + DR3
            FXK(JN) = + DRX + DR4
            FYK(JN) = + DRY + DR5
            FZK(JN) = + DRZ + DR6
            FXL(JN) = - DR4
            FYL(JN) = - DR5
            FZL(JN) = - DR6
  220     CONTINUE
C
C         ----- END OF A STRIP OF DIHEDRALS AND START OF 1-4 NB -----
C
          DO 300 JN = 1,MAXLEN
            I3 = IP(JN+IST)
            L3T = LP(JN+IST)
            L3 = IABS(L3T)
            XIJ(JN) = X(I3+1)-X(L3+1)
            YIJ(JN) = X(I3+2)-X(L3+2)
            ZIJ(JN) = X(I3+3)-X(L3+3)
  300     CONTINUE
          IF(NTB.EQ.0 .OR. NOCRST) GO TO 320
C
C           ----- APPLY PERIODIC BOUNDARY CONDITION -----
C
            CALL PIMAG(MAXLEN,XIJ,YIJ,ZIJ)
C
            IF(NTB.GT.0) GO TO 320
              CALL IMAGT(MAXLEN,XIJ,YIJ,ZIJ)
C
C             ----- NOW LOOP OVER ALL THE DIHEDRALS (CONSTANT DIEL) -----
C
  320     CONTINUE
          DO 360 JN = 1,MAXLEN
            CT(JN) = XIJ(JN)*XIJ(JN)+YIJ(JN)*YIJ(JN)+ZIJ(JN)*ZIJ(JN)
  360     CONTINUE
C
          IF (DIELD) GO TO 720
            DO 700 JN = 1,MAXLEN
              I3 = IP(JN+IST)
              K3T = KP(JN+IST)
              L3T = LP(JN+IST)
              IC0 = ICP(JN+IST)
              IDUMI = ISIGN(1,K3T)
              IDUML = ISIGN(1,L3T)
              KDIV = (2+IDUMI+IDUML)/4
              L3 = IABS(L3T)
              FMULN = FLOAT(KDIV)*FMN(IC0)
C
              II = (I3+3)/3
              JJ = (L3+3)/3
              IA1 = IAC(II)
              IA2 = IAC(JJ)
              IBIG = MAX0(IA1,IA2)
              ISML = MIN0(IA1,IA2)   
              IC = IBIG*(IBIG-1)/2+ISML
C
C             ----- CALCULATE THE 14-EEL ENERGY -----
C
              R2 = FMULN/CT(JN)
              R1 = SQRT(R2)
#ifdef LES
              lfac=lesfac(nlesty*(lestyp(ii)-1)+lestyp(jj))
#else
              lfac = 1.d0
#endif
              G = CG(II)*CG(JJ)*R1*lfac
              SPHI(JN) = G
              if (isftrp.le.0) then
C
C               --- normal 6-12:
c
                R6 = R2*R2*R2
                R12 = R6*R6
#ifdef CHARMM
                F1 = CN114(IC)*R12*lfac
                F2 = CN214(IC)*R6*lfac
#else
                F1 = CN1(IC)*R12*lfac
                F2 = CN2(IC)*R6*lfac
#endif
                CPHI(JN) = F1-F2

c Margaret. constant dielec,non, or periodic both here.
c Margaret change here
                  if(maphb(II,JJ).lt.1) then
                  EHBV=EHBV+F1-F2
                  else
                  EHBA=EHBA+F1-F2
                  endif
c end of Margaret




                DFN =((-twelve*F1+six*F2)*SCNB0-G*SCEE0)*R2
              else
c
C               --- soft-repulsion nonbond:
c
                rrk = rwell*lfac
                rst = rad(ia1) + rad(ia2)
                arg = fmuln * max(zero,rst*rst-ct(jn))
                cphi(jn) = rrk * arg*arg
                dfn = -four*rrk*arg*scnb0 - g*scee0*r2
              end if
c
              XA = XIJ(JN)*DFN
              YA = YIJ(JN)*DFN 
              ZA = ZIJ(JN)*DFN
              FXI(JN) = FXI(JN)-XA
              FYI(JN) = FYI(JN)-YA
              FZI(JN) = FZI(JN)-ZA
              FXL(JN) = FXL(JN)+XA
              FYL(JN) = FYL(JN)+YA
              FZL(JN) = FZL(JN)+ZA
  700       CONTINUE
            GO TO 760
  720     CONTINUE
C
C           ----- DISTANCE DEPENDENT DIELEECTRIC -----
C
            DO 740 JN = 1,MAXLEN
              I3 = IP(JN+IST)
              K3T = KP(JN+IST)
              L3T = LP(JN+IST)
              IC0 = ICP(JN+IST)
              IDUMI = ISIGN(1,K3T)
              IDUML = ISIGN(1,L3T)
              KDIV = (2+IDUMI+IDUML)/4
              L3 = IABS(L3T)
              FMULN = FLOAT(KDIV)*FMN(IC0)
C
              II = (I3+3)/3  
              JJ = (L3+3)/3
              IA1 = IAC(II)
              IA2 = IAC(JJ)
              IBIG = MAX0(IA1,IA2)
              ISML = MIN0(IA1,IA2)
              IC = IBIG*(IBIG-1)/2+ISML
C
C             ----- CALCULATE THE 14-EEL ENERGY -----
C
              R2 = FMULN/CT(JN)
#ifdef LES
              lfac=lesfac(nlesty*(lestyp(ii)-1)+lestyp(jj))
#else
              lfac = 1.d0
#endif
              G = CG(II)*CG(JJ)*R2*lfac
              SPHI(JN) = G
              if (isftrp.le.0) then
C
C               --- normal 6-12:
                R6 = R2*R2*R2
                R12 = R6*R6
#ifdef CHARMM
                F1 = CN114(IC)*R12*lfac
                F2 = CN214(IC)*R6*lfac
#else
                F1 = CN1(IC)*R12*lfac
                F2 = CN2(IC)*R6*lfac
#endif
                CPHI(JN) = F1-F2
                DFN =((-twelve*F1+six*F2)*SCNB0-(G+G)*SCEE0)*R2
              else
c               -- soft-repulsion nonbond
                rrk = rwell*lfac
                rst = rad(ia1) + rad(ia2)
                arg = fmuln * max(zero,rst*rst-ct(jn))
                cphi(jn) = rrk * arg*arg
                dfn = -four*rrk*arg*scnb0 - (g+g)*scee0*r2
              end if
c
              XA = XIJ(JN)*DFN
              YA = YIJ(JN)*DFN
              ZA = ZIJ(JN)*DFN
              FXI(JN) = FXI(JN)-XA
              FYI(JN) = FYI(JN)-YA
              FZI(JN) = FZI(JN)-ZA
              FXL(JN) = FXL(JN)+XA
              FYL(JN) = FYL(JN)+YA
              FZL(JN) = FZL(JN)+ZA
  740       CONTINUE
  760     CONTINUE
C
C         ----- SUMUP ALL THE GRADIENTS -----
C
          DO 780 JN = 1,MAXLEN
            I3 = IP(JN+IST)
            J3 = JP(JN+IST)
            K3 = IABS(KP(JN+IST))
            L3 = IABS(LP(JN+IST))
C      
            F(I3+1) = F(I3+1) + FXI(JN)
            F(I3+2) = F(I3+2) + FYI(JN)
            F(I3+3) = F(I3+3) + FZI(JN)
            F(J3+1) = F(J3+1) + FXJ(JN) 
            F(J3+2) = F(J3+2) + FYJ(JN)
            F(J3+3) = F(J3+3) + FZJ(JN)
            F(K3+1) = F(K3+1) + FXK(JN)
            F(K3+2) = F(K3+2) + FYK(JN)
            F(K3+3) = F(K3+3) + FZK(JN)
            F(L3+1) = F(L3+1) + FXL(JN)
            F(L3+2) = F(L3+2) + FYL(JN)
            F(L3+3) = F(L3+3) + FZL(JN)
  780     CONTINUE
          do 440 ksum = 1,maxlen
            enbpl = enbpl+cphi(ksum)
            eelpl = eelpl+sphi(ksum)
            epl   = epl  +epw(ksum)
  440     continue
*         ENBPL = ENBPL+SSUM(MAXLEN,CPHI(1),1)
*         EELPL = EELPL+SSUM(MAXLEN,SPHI(1),1)
*         EPL   = EPL  +SSUM(MAXLEN,EPW(1),1)
C         
          IST = IST+MAXLEN
          IF(IST.LE.ISTC) GO TO 820
            LENC = IST-ISTC
            do 450 ksum = maxlen-lenc+1,maxlen
              ecnl = ecnl+epw(ksum)
  450       continue
*           ECNL = ECNL+SSUM(LENC,EPW(MAXLEN-LENC+1),1)
            ISTC = ISTC+LENC
  820   CONTINUE
      IF(.NOT.SKIP) GO TO 4200
C
C     ---- ALL DONE -----
C
      ENBP = ENBPL*SCNB0
      EELP = EELPL*SCEE0
      EP   = EPL
      ECN = ECNL
#ifdef MPI
      NPHI = nbtmp
#endif
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE CAPWAT(NAT,X,F)
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
c     Mods for Rev A by GLS:
c     cpp switchable single/double prec.
c
#ifdef DPREC
      implicit double precision (a-h,o-z)
#endif
#ifdef MPI
#include "parallel.h"
#endif
C
C     ----- ROUTINE TO CALCULATE THE CAP FORCE -----
C
#include "box.h"
      DIMENSION X(3,*),F(3,*)
      DATA TM34,zero/1.0D-34,0.0d0/
C
#ifdef MPI
      DO 100 I = NATCAP+1+mytaskid,NAT,numtasks
#else
      DO 100 I = NATCAP+1,NAT
#endif
      XA = XCAP-X(1,I)
      YA = YCAP-X(2,I)
      ZA = ZCAP-X(3,I)
      DA = SQRT(XA*XA+YA*YA+ZA*ZA+TM34)
cd    DF = FCAP*AMAX1(zero,DA-CUTCAP)/DA
      DF = FCAP*MAX(zero,DA-CUTCAP)/DA
      F(1,I) = F(1,I)+DF*XA
      F(2,I) = F(2,I)+DF*YA
      F(3,I) = F(3,I)+DF*ZA
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE CBOND(NB,IB,JB,ICB,X,F,EB,NOCRST)
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
c     cpp-switchable single/double prec.
#ifdef DPREC
      implicit double precision (a-h,o-z)
#endif
#ifdef MPI
#include "parallel.h"
#endif
      LOGICAL NOCRST
C
C     ----- ROUTINE TO GET BOND ENERGY AND FORCES FOR THE POTENTIAL
C           OF CB*(B-B0)**2
C
#include "box.h"
#include "parms.h"
C
      DIMENSION IB(*),JB(*),ICB(*),X(*),F(*)
      DIMENSION XIJ(3)
      data zero,two/0.0d0,2.0d0/
C
      EB = zero
C
#ifdef MPI
      DO 30 N = mytaskid+1,NB,numtasks
#else
      DO 30 N = 1,NB
#endif
      I3 = IB(N)
      J3 = JB(N)
      IPC = ICB(N)
      RIJ2 = zero
      DO 32 M = 1,3
      XIJ(M) = X(I3+M)-X(J3+M)
   32 RIJ2 = RIJ2+XIJ(M)*XIJ(M)
      IF(NTB.EQ.0) GOTO 41
C
C     ----- CALL THE PERIODIC CONDITION -----
C
      CALL PERCON(RIJ2,XIJ)
   41 CONTINUE
      RIJ = SQRT(RIJ2)
      DB = RIJ-REQ(IPC)
C
C     ----- SKIP THE ENERGY IF DB IS WITHIN BOUND -----
C
c     --- uncomment following line to make half harmonic ---
C     IF(DB.LE.zero) GO TO 30
C
      DF = RK(IPC)*DB
      EBH = DF*DB
      EB = EB+EBH
C
C     ----- UPDATE THE FORCE ARRAY -----
C
      DF = two*DF/RIJ
      DO 90 M = 1,3
      XH = XIJ(M)*DF
      F(I3+M) = F(I3+M)-XH
      F(J3+M) = F(J3+M)+XH
   90 CONTINUE
   30 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE XCONST(NATC,ECON,IGROUP,X,F,XC,WEIT)
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
#ifdef MPI
#include "parallel.h"
#endif
c     Mods for Rev A by GLS:
c     cpp selectable single/double prec.
C
C     ----- ROUTINE TO PUT HARMONIC CONSTRAINTS FOR POSITION -----
C
      DIMENSION IGROUP(*),X(*),F(*),XC(*),WEIT(*)
C
      ECON = 0.0D+00
#ifdef MPI
      DO 100 II = 1+mytaskid,NATC,numtasks
#else
      DO 100 II = 1,NATC
#endif
      I = IGROUP(II)
      WT = WEIT(II)
      I3 = 3*I-3
      AX = X(I3+1)-XC(I3+1)
      AY = X(I3+2)-XC(I3+2)
      AZ = X(I3+3)-XC(I3+3)
      WX = WT*AX
      WY = WT*AY
      WZ = WT*AZ
      EADD = WX*AX+WY*AY+WZ*AZ
      ECON = ECON+EADD
      F(I3+1) = F(I3+1)-(WX+WX)
      F(I3+2) = F(I3+2)-(WY+WY)
      F(I3+3) = F(I3+3)-(WZ+WZ)
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE BELLYF(NAT,IGRP,F)
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
      DIMENSION IGRP(*),F(*)
      data zero/0.0d0/
      DO 100 I = 1,NAT
      IF(IGRP(I).GT.0) GO TO 100
      I3 = 3*I-3
      F(I3+1) = zero
      F(I3+2) = zero
      F(I3+3) = zero
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE IMAGT(MAX,X,Y,Z)
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
c     Mods for Rev A by GLS:
c     cpp selectable single/double prec.
#include "box.h"
      DIMENSION X(*),Y(*),Z(*)
      data zero,one /0.0d0,1.0d0/
C
      DO 140 JN = 1,MAX
           DUMI = BOXOQ-(ABS(X(JN))+ABS(Y(JN))+ABS(Z(JN)))
           DUMJ = one
           IF(DUMI.GT.zero) DUMJ = zero
C          DUMJ = CVMGP(zero,one,DUMI)
           X(JN) = X(JN)-SIGN(BOXOH,X(JN))*DUMJ
           Y(JN) = Y(JN)-SIGN(BOXOH,Y(JN))*DUMJ
           Z(JN) = Z(JN)-SIGN(BOXOH,Z(JN))*DUMJ 
  140 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE IMAGT2(MAX,X)
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
c     Mods for Rev A by GLS:
c     cpp selectable single/double prec.
#include "box.h"
      DIMENSION X(3,*)
      data zero,one /0.0d0,1.0d0/
C
      DO 140 JN = 1,MAX
           DUMI = BOXOQ-(ABS(X(1,JN))+ABS(X(2,JN))+
     +              ABS(X(3,JN)))
           DUMJ = one
           IF(DUMI.GT.zero) DUMJ = zero
C          DUMJ = CVMGP(zero,one,DUMI)
           X(1,JN) = X(1,JN)-SIGN(BOXOH,X(1,JN))*DUMJ
           X(2,JN) = X(2,JN)-SIGN(BOXOH,X(2,JN))*DUMJ
           X(3,JN) = X(3,JN)-SIGN(BOXOH,X(3,JN))*DUMJ
  140 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE TRACO(NR,NPA,X,IT)
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
c     Mods for Rev A by GLS:
c     cpp selectable single/double prec.
      LOGICAL NEW
C
C     ----- ROUTINE TO MAKE TRANSFORMATION FROM CARTESIAN TO OBLIQUE
C           COORDINATES AND VICE VERSA -----
C
      DIMENSION X(*)
      DATA NEW/.FALSE./
C
C     ----- INITIALIZE -----
C
      IF(NEW) GO TO 8
      NEW = .TRUE.
      ONE = 1.0D0
      PYE = 4.0D0* ATAN(ONE)
      CONV = 1.8D2/PYE
      BETAR = 90.0D0/CONV
      COSB = COS(BETAR)
      SINB = SIN(BETAR)
    8 CONTINUE
C
C     ----- PERFORM THE TRANSFORMATION -----
C
      IF(IT) 41,51,31
C
   31 I3 = 3*NPA+1
      DO 40 J = 1,NR
      XH = X(I3+2)/SINB
      X(I3+2) = XH
      X(I3) = X(I3)-COSB*XH
      I3 = I3+3
   40 CONTINUE
      GO TO 51
   41 I3 = 3*NPA+1
      DO 50 J = 1,NR
      X(I3) = X(I3)+COSB*X(I3+2)
      X(I3+2) = X(I3+2)*SINB
      I3 = I3+3
   50 CONTINUE
   51 RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE TRACOL(NR,NPA,X,Y,Z,IT)
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
c     Mods for Rev A by GLS:
c     cpp selectable single/double prec.
      LOGICAL NEW
C
C     ----- ROUTINE TO MAKE TRANSFORMATION FROM CARTESIAN TO OBLIQUE
C           COORDIANTES AND VICE VERSA -----
C
      DIMENSION X(*),Y(*),Z(*)
      DATA NEW/.FALSE./
C
C     ----- INITIALIZE -----
C
      IF(NEW) GO TO 100
      NEW = .TRUE.
      ONE = 1.0D0
      PYE = 4.0D0* ATAN(ONE)
      CONV = 1.8D2/PYE
      BETAR = 90.0D0/CONV
      COSB = COS(BETAR)
      SINB = SIN(BETAR)
  100 CONTINUE
C
C     ----- PERFORM THE TRANSFORMATION -----
C
      IF(IT) 140,150,130
C
  130 CONTINUE
      I = NPA+1
      DO 160 J = I,NR
      XH = Z(J)/SINB
      Z(J) = XH
      X(J) = X(J)-COSB*XH
  160 CONTINUE
      GO TO 150
  140 CONTINUE
      I = NPA+1
      DO 180 J = I,NR
      X(J) = X(J)+COSB*Z(J)
      Z(J) = Z(J)*SINB
  180 CONTINUE
  150 RETURN
      END
c-----------------------------------------------------------------------
