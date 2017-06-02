#include "vector.h"
c TODO - put loops inside #Ifdef USENINTso that cforcevector can 
c         be used on them
c        consider ifdef DOUBLE define XNINT DNINT else define XNINT DNINT
c        check input to make sure not obliq when ifbox=2
c     ---- period.f: permutations of periodic imaging
c          over rectangular boxes, truncated octahedral boxes,
c          coords in 2D or 1D arrays, etc. These routines were
c          gathered here from various places when truncated
c          octahedral boxes were added, because oct stuff
c          increases the length of the code significantly.
c
c          The boundX() routines were originally inline code.
c
c          routines            where from originally      
c
c          mdbox()             set.f
c          pimg(), pimg2()     ene.f
c          percon()            misc.f
c          shiag()             misc.f
c          bound1()            efield.f, indip.f, polder.f, nonbon.f
c          bound2()            efield.f, indip.f, polder.f, nonbon.f
c          bound3()            dynlib.f
c          bound6()            qiktip.f
c          bound7()            resnba.f
c          bound9()            shake.f
c          bound10()           threeb.f
c
#ifdef USENINT
C
C         The sections of code within USENINT #ifdefs provide an 
c         alternate imaging scheme, in which periodic imaging is 
c         performed not by determining if xi-xj > 1/2 boxlength (the 
c         "standard" method) but by dividing xi-xj by the box length 
c         and taking the nearest integer. It has been found that this 
c         scheme is faster on some machines (e.g. Cray YMP, IBM RS/6000).
C
#endif
      SUBROUTINE MDBOX
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
      LOGICAL initd
C
C     ----- ROUTINE TO SET THE BOX RELATED VARIABLES -----
C
#include "box.h"
      data ONE  /1.0d0/
      DATA initd/.FALSE./
c     save initd
c
c     -- set arrays for box & box inverted 
c
      if (ntb.ne.0) then
          BOXI(1) = ONE/BOX(1)
          BOXI(2) = ONE/BOX(2)
          BOXI(3) = ONE/BOX(3)
      endif
c
      IF (initd) then
c
c         -- after the first time, only calc if constant pressure
c
          if (NTB.ne.2) return
c
c         -- calc size of max vector that is guaranteed to fit
c            in box w/out imaging
c
          BOXHM = BOX(1)*0.5d0
          DO 5 M = 1,3
              BOXH(M) = BOX(M)*0.5d0
    5         IF(BOXH(M).LT.BOXHM) BOXHM = BOXH(M)
          BOXHM2 = BOXHM*BOXHM
      else if (ntb.gt.0) then
c
c         -- first time and periodic conditions are involved;
c            set the box-related stuff and 
c         -- calc size of max vector that is guaranteed to fit
c            in box w/out imaging
c
          initd = .TRUE.
          BOXOQ = 1.0d+20
          NTBB = NTB
          PYE = 4.d0 * ATAN(ONE)
          CONV = 1.8d2/PYE
          NTM = 0
          BETAR = 90.0D0/CONV
          COSB = COS(BETAR)
          COSB2 = COSB+COSB
          IF( ABS(COSB).GE.1.d-4) NTM = 1
          BOXHM = BOX(1)*0.5d0
          DO 6 M = 1,3
              BOXH(M) = BOX(M)*0.5d0
    6         IF(BOXH(M).LT.BOXHM) BOXHM = BOXH(M)
          BOXHM2 = BOXHM*BOXHM
c
c         -- set convenient "four thirds" var in case octa bound
c
          fothi = 4.0d0/3.0d0
      else
c
c         -- no periodic situation: scram
c
          return
      endif

c
c     -- periodic conditions, {first time} or {later + constant press}
c        go the extra mile if truncated octahedral
c
      if (ifbox.eq.2) then
c
c         -- set "boxahedron" array
c
          BOCT(1) = box(3)/box(1)
          BOCT(2) = box(3)/box(2)
          BOCT(3) = box(3)*0.75
c
c         Truncated octahedron with different xyz box dimensions:
c         calculate the max length of a vector in the box which is 
c         never transformed by the periodicity transformations.
c
c         We view the box as centered at 0,0,0 for purposes of this 
c         one calculation (in the simulation space the _corner_ is at
c         0,0,0).  We consider the x>0,y>0,z>0 (positive) octant of 
c         this box, which extends to halfboxX,halfboxY,halfboxZ, and
c         has the corner chopped off by a plane that intersects the
c         axes at 1.5 times the (halfbox) extents: .75 * 2 * halfboxX, 
c         .75 * 2 * halfboxY, .75 * 2 * halfboxZ, or .75 * each of the
c         full box dimensions:
c
c                 .75  ( box(1)|      0|      0)  = P
c                 .75  (      0| box(2)|      0)  = Q
c                 .75  (      0|      0| box(3))  = R
c
c         so a normal vector to this plane is the cross-product of
c         the vectors PQxQR, which (using X,Y,Z for box(1),box(2),box(3)) 
c         comes to .75^2 (YZ,XZ,XY), and the equation for the plane is 
c         YZx + XZy + XYz = .75 XYZ, and the distance of this plane 
c         from the origin is .75 XYZ / sqrt((YZ)^2 + (XZ)^2 + (XY)^2).
c
          yz = box(2) * box(3)
          xz = box(1) * box(3)
          xy = box(1) * box(2)
          dum = (0.75 * box(1) * box(2) * box(3)) /
     +          sqrt( yz*yz + xz*xz + xy*xy )
          write(6, '(1x,a,f8.3,f8.3)') 
     .       'oct box min vector; ', boxhm,dum
          if (dum.lt.boxhm) then
              boxhm = dum
              boxhm2 = boxhm*boxhm
          endif
      endif
      RETURN
      END
c-----------------------------------------------------------------------

      SUBROUTINE PIMAG(MAX,X,Y,Z)
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
c     --- PIMAG() - 1D arrays for X,Y,Z
c
#ifdef DPREC
      implicit double precision (a-h,o-z)
#endif
c     Mods for Rev A by GLS:
c     cpp selectable single/double prec.
      DIMENSION X(*),Y(*),Z(*)
#include "box.h"
      data zero, one /0.0d0, 1.0d0/
C
      if (ifbox.eq.1) then
c
c         -- rectilinear box
c
cforcevector
          DO 140 JN = 1,MAX
              IF(ABS(X(JN)).gt.BOXH(1)) 
     +            X(JN) = X(JN)-SIGN(BOX(1),X(JN))
              IF(ABS(Y(JN)).gt.BOXH(2)) 
     +            Y(JN) = Y(JN)-SIGN(BOX(2),Y(JN))
              IF(ABS(Z(JN)).gt.BOXH(3)) 
     +            Z(JN) = Z(JN)-SIGN(BOX(3),Z(JN))
C             X(JN) = CVMGP(X(JN)-
C    .                 SIGN(BOX(1),X(JN)),X(JN),ABS(X(JN))-BOXH(1))
C             Y(JN) = CVMGP(Y(JN)-
C    .                 SIGN(BOX(2),Y(JN)),Y(JN),ABS(Y(JN))-BOXH(2))
C             Z(JN) = CVMGP(Z(JN)-
C    .                 SIGN(BOX(3),Z(JN)),Z(JN),ABS(Z(JN))-BOXH(3))
  140     CONTINUE
      else
c
c         -- truncated octahedral box
c
#ifdef USENINT
cforcevector
          DO 141 JN = 1,MAX
              dux = X(JN)*boxi(1)
              duy = Y(JN)*boxi(2)
              duz = Z(JN)*boxi(3)
# ifdef DPREC
              dux = dux - dnint( dux )
              duy = duy - dnint( duy )
              duz = duz - dnint( duz )
# else
              dux = dux - anint( dux )
              duy = duy - anint( duy )
              duz = duz - anint( duz )
# endif
              corr= 0.5 * aint( fothi* (abs( dux ) +
     +                                 abs( duy ) +
     +                                 abs( duz )))
              dux = dux - sign( corr,dux )
              duy = duy - sign( corr,duy )
              duz = duz - sign( corr,duz )
              X(JN) = dux * box(1)
              Y(JN) = duy * box(2)
              Z(JN) = duz * box(3)
  141     CONTINUE
#else
cforcevector
          DO 141 JN = 1,MAX
              dux = abs(X(JN))
              duy = abs(Y(JN))
              duz = abs(Z(JN))
              if( dux.gt.boxh(1) ) then
                  X(JN) = X(JN) - sign( box(1),X(JN) )
                  dux  = abs( X(JN) )
              endif
              if( duy.gt.boxh(2) ) then
                  Y(JN) = Y(JN) - sign( box(2),Y(JN) )
                  duy  = abs( Y(JN) )
              endif
              if( duz.gt.boxh(3) ) then
                  Z(JN) = Z(JN) - sign( box(3),Z(JN) )
                  duz  = abs( Z(JN) )
              endif
              if( (dux*boct(1)+duy*boct(2)+duz).gt.boct(3) ) then
                  X(JN) =
     +                    X(JN) - sign( boxh(1),X(JN) )
                  Y(JN) =
     +                    Y(JN) - sign( boxh(2),Y(JN) )
                  Z(JN) =
     +                    Z(JN) - sign( boxh(3),Z(JN) )
              endif
  141     CONTINUE
#endif
      endif
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE PIMAG2(MAX,X)
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
c     --- PIMAG2() - 2D array for X,Y,Z
c
#ifdef DPREC
      implicit double precision (a-h,o-z)
#endif
c     Mods for Rev A by GLS:
c     cpp selectable single/double prec.
      DIMENSION X(3,*)
#include "box.h"
      data zero, one /0.0d0, 1.0d0/
C
      if (ifbox.eq.1) then
c
c         -- rectilinear box
c
cforcevector
          DO 140 JN = 1,MAX
              IF(ABS(X(1,JN)).gt.BOXH(1)) 
     +            X(1,JN) = X(1,JN)-SIGN(BOX(1),X(1,JN))
              IF(ABS(X(2,JN)).gt.BOXH(2)) 
     +            X(2,JN) = X(2,JN)-SIGN(BOX(2),X(2,JN))
              IF(ABS(X(3,JN)).gt.BOXH(3)) 
     +            X(3,JN) = X(3,JN)-SIGN(BOX(3),X(3,JN))
C             X(1,JN) = CVMGP(X(1,JN)-SIGN(BOX(1),X(1,JN)),X(1,JN),
C    +                   ABS(X(1,JN))-BOXH(1))
C             X(2,JN) = CVMGP(X(2,JN)-SIGN(BOX(2),X(2,JN)),X(2,JN),
C    +                   ABS(X(2,JN))-BOXH(2))
C             X(3,JN) = CVMGP(X(3,JN)-SIGN(BOX(3),X(3,JN)),X(3,JN),
C    +                   ABS(X(3,JN))-BOXH(3))
  140     CONTINUE
      else
c
c         -- truncated octahedral box
c
#ifdef USENINT
cforcevector
          DO 141 JN = 1,MAX
              dux = X(1,JN)*boxi(1)
              duy = X(2,JN)*boxi(2)
              duz = X(3,JN)*boxi(3)
#  ifdef DPREC
              dux = dux - dnint( dux )
              duy = duy - dnint( duy )
              duz = duz - dnint( duz )
#  else
              dux = dux - anint( dux )
              duy = duy - anint( duy )
              duz = duz - anint( duz )
#  endif
              corr= 0.5 * aint( fothi* (abs( dux ) +
     +                                 abs( duy ) +
     +                                 abs( duz )))
              dux = dux - sign( corr,dux )
              duy = duy - sign( corr,duy )
              duz = duz - sign( corr,duz )
              X(1,JN) = dux * box(1)
              X(2,JN) = duy * box(2)
              X(3,JN) = duz * box(3)
  141     CONTINUE
#else
cforcevector
          DO 141 JN = 1,MAX
              dux = abs(X(1,JN))
              duy = abs(X(2,JN))
              duz = abs(X(3,JN))
              if( dux.gt.boxh(1) ) then
                  X(1,JN) = X(1,JN) - sign( box(1),X(1,JN) )
                  dux  = abs( X(1,JN) )
              endif
              if( duy.gt.boxh(2) ) then
                  X(2,JN) = X(2,JN) - sign( box(2),X(2,JN) )
                  duy  = abs( X(2,JN) )
              endif
              if( duz.gt.boxh(3) ) then
                  X(3,JN) = X(3,JN) - sign( box(3),X(3,JN) )
                  duz  = abs( X(3,JN) )
              endif
              if( (dux*boct(1)+duy*boct(2)+duz).gt.boct(3) ) then
                  X(1,JN) =
     +                  X(1,JN) - sign( boxh(1),X(1,JN) )
                  X(2,JN) =
     +                  X(2,JN) - sign( boxh(2),X(2,JN) )
                  X(3,JN) =
     +                  X(3,JN) - sign( boxh(3),X(3,JN) )
              endif
  141     CONTINUE
#endif
      endif
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE PERCON(RIJQ,XIJ)
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
#include "box.h"
      DIMENSION XIJ(3)
      data zero /0.0d0/
C
C     ----- PERIODIC BOUNDARY CONDITION THROUGH MINIMUM IMAGE
C           CONVENTION -----
C
      RIJ2 = RIJQ
      IF(RIJ2.lt.BOXHM2) then
          return
      endif
      if (ifbox.eq.1) then
c
c         --- rectangular boundaries
c
          DO 180 M = 1,3
              IF(XIJ(M).ge.BOXH(M)) then
                  XIJ(M) = XIJ(M)-BOX(M)
              else if (XIJ(M).lt.-BOXH(M)) then
                  XIJ(M) = XIJ(M)+BOX(M)
              endif
  180     CONTINUE
      else
c
c         --- truncated octahedral boundaries
c
#ifdef USENINT
          dux = xij(1)*boxi(1)
          duy = xij(2)*boxi(2)
          duz = xij(3)*boxi(3)
#  ifdef DPREC
          dux = dux - dnint( dux )
          duy = duy - dnint( duy )
          duz = duz - dnint( duz )
#  else
          dux = dux - anint( dux )
          duy = duy - anint( duy )
          duz = duz - anint( duz )
#  endif
          corr= 0.5 * aint( fothi* (abs( dux ) +
     +                             abs( duy ) +
     +                             abs( duz )))
          dux = dux - sign( corr,dux )
          duy = duy - sign( corr,duy )
          duz = duz - sign( corr,duz )
          xij(1) = dux * box(1)
          xij(2) = duy * box(2)
          xij(3) = duz * box(3)
#else
          dux = abs(xij(1))
          duy = abs(xij(2))
          duz = abs(xij(3))
          if ( dux.gt.boxh(1) ) then
              xij(1) = xij(1) - sign( box(1),xij(1) )
              dux  = abs( xij(1) )
          endif
          if ( duy.gt.boxh(2) ) then
              xij(2) = xij(2) - sign( box(2),xij(2) )
              duy  = abs( xij(2) )
          endif
          if( duz.gt.boxh(3) ) then
              xij(3) = xij(3) - sign( box(3),xij(3) )
              duz  = abs( xij(3) )
          endif
          if ( (dux*boct(1)+duy*boct(2)+duz).gt.boct(3) ) then
              xij(1) = xij(1) - sign( boxh(1),xij(1) )
              xij(2) = xij(2) - sign( boxh(2),xij(2) )
              xij(3) = xij(3) - sign( boxh(3),xij(3) )
          endif

#endif
      endif
      RIJ2 = zero
      DO 260 M = 1,3
  260     RIJ2 = RIJ2+XIJ(M)*XIJ(M)
      RIJQ = RIJ2
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE SHIAG(NSPM,NPA,NSP,X,NR,SHIFT)
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
C
C     ----- SHIAG TRANSLATES COUNTER ION  ATOMS AND SOLVENT (** now ALL **)
C           MOLECULES, APPLYING PERIODIC BOUNDARY CONDITIONS, SUCH THAT
C           THE COUNTER ION OR THE FIRST ATOM OF A SOLVENT MOLECULE
C           LIES WITHIN THE SPECIFIED PERIODIC BOX. THIS BOX LIES IN
C           THE POSITIVE QUADRANT. IT CAN BE A TRUNCATED OCTAHEDRON,
C           RECTANGULAR(BETA = 90 DEGREES) OR MONOCLINIC (BETA.NE.90)
C           THE ATOMS OF A SOLVENT MOLECULE MUST LIE WITHIN BOX/2
C           OR BOX*SQRT(3)/4 OF EACH OTHER -----
C
C     Mods for version 4.1:
C     - If "shift" flag is true, translate the entire system so that
C       the solute is centered in central "box" -- dac
C     - If IFTRES = 0 (no imaging applied to solute-solute interactions), 
C       then never translate solute atoms back into the box -- dap
C
c     Mods for Rev A by gls:
c     - hardwired to translate ALL molecules, as is required by
c       the Amber imaging method.  This removes need for rcmsol
c       hack, which was incorrect.
c     - added nr to args because it is needed for call to traco.
c     - cpp-selectable single/double prec.
c     - added check on nspsol to cover reported bug in some versions.
c       (should not be necessary but sanity checks are a good thing)
c
#include "box.h"
C
      DIMENSION NSP(*),X(*)
      DIMENSION XI(3)
      logical shift
C
      IF (NSPM.LE.0) return

C
C     ----- TRANSLATE THE SOLVENT MOLECULES -----
C
c
      if (nspsol .le. 0 .or. nr .le. 0) then
           write(6,'(/5x,a)') 'shiag: error in args'
           call mexit(6, 1)
      endif
C
C If the shift flag is set, translate the system so that the geometic
C center of the solute is at the center of the box. Only use solute
C molecules of > 1 atom (do not include ions, therefore), unless
C all molecules of the solute are 1 atom = ions.
c
      i3 = 0
      if (shift) then
          xi(1) = 0.0d0
          xi(2) = 0.0d0
          xi(3) = 0.0d0
          nute = 0
          nmax = 0
          do 6  jn=1,nspsol-1
              if (nsp(jn).gt.nmax) nmax = nsp(jn)
    6     continue

          do 12 jn = 1,nspsol-1
              if(jn.eq.1 .or. nsp(jn).gt.1 .or. nmax.eq.1) then
                  do 10 i=1,nsp(jn)
                      xi(1) = xi(1) + x(i3+1)
                      xi(2) = xi(2) + x(i3+2)
                      xi(3) = xi(3) + x(i3+3)
                      i3 = i3 + 3
                      nute = nute + 1
   10              continue
              else
                  i3 = i3 + 3
               endif
   12     continue

          if (nute.gt.0) then
              xi(1) = xi(1)/nute - box(1)/2.d0
              xi(2) = xi(2)/nute - box(2)/2.d0
              xi(3) = xi(3)/nute - box(3)/2.d0
              i3 = 0
              do 20 i=1,nr
                  x(i3+1) = x(i3+1) - xi(1)
                  x(i3+2) = x(i3+2) - xi(2)
                  x(i3+3) = x(i3+3) - xi(3)
                  i3 = i3 + 3
   20         continue
          endif
      endif
C
C     Now translate molecules back into central box. If IFTRES=1,
C     do translation over all molecules. If IFTRES=0, then no imaging
C     is being applied to solute-solute interactions. In this case,
C     do translation only over molecules of the solvent.--dap
C
      IF (IFTRES.EQ.0) THEN
          I3 = 3*NPA
          IBEG = NSPSOL
      ELSE
          I3 = 0
          IBEG = 1
      END IF
      DO 90 NN = IBEG,NSPM
          NS = 0
          NRAM = NSP(NN)
          NRAM3 = 3*NRAM-3
C
c         --- use first atom of this molecule to figure out if it needs
c             to be translated back to box center or not:
c
          if (ifbox.eq.1) then
c
c             --- cubic boundary
c
              DO 74 M = 1,3
                  I3 = I3+1
                  XI(M) = X(I3)
                  IF (X(I3).gt.BOX(M)) then 
                      dist = x(i3) - box(m)
c                      if (dist.gt.box(m)) then
c                          write(6,'(t2,a,i5,a)')
c     .                      '%SANDER-W-ESCAPE, atom # ',int(i3/3),
c     .                      ' has escaped'
c                      endif
                      X(I3) = X(I3)-BOX(M)
                      NS = 1
                  else IF(X(I3).lt.0.d0) then
c                      if (abs(x(i3)).gt.box(m)) then
c                          write(6,'(t2,a,i5,a)')
c     .                      '%SANDER-W-ESCAPE, atom # ',int(i3/3),
c     .                      ' has escaped'
c                      endif
                      X(I3) = X(I3)+BOX(M)
                      NS = 1
                  endif
   74         CONTINUE
          else
c
c             --- truncated octahedral boundary
c
              XI(1) = X(I3+1)
              XI(2) = X(I3+2)
              XI(3) = X(I3+3)
              dux = abs(x(i3+1)-boxh(1))
              duy = abs(x(i3+2)-boxh(2))
              duz = abs(x(i3+3)-boxh(3))
c
              if (dux.gt.boxh(1)) then
                  x(i3+1) = x(i3+1) - sign( box(1), x(i3+1)-boxh(1) )
                  dux = abs( x(i3+1)-boxh(1) )
                  if ( dux.gt.boxh(1) ) then
                      write(6,'(t2,a,i5,a)')
     .                      '%SANDER-W-ESCAPE, atom # ',int(i3/3+1),
     .                      ' has escaped'
                  endif
                  ns = 1
              endif
              if (duy.gt.boxh(2)) then
                  x(i3+2) = x(i3+2) - sign( box(2), x(i3+2)-boxh(2) )
                  duy = abs( x(i3+2)-boxh(2) )
                  if ( duy.gt.boxh(2) ) then
                      write(6,'(t2,a,i5,a)')
     .                      '%SANDER-W-ESCAPE, atom # ',int(i3/3+1),
     .                      ' has escaped'
                  endif
                  ns = 1
              endif
              if (duz.gt.boxh(3)) then
                  x(i3+3) = x(i3+3) - sign( box(3), x(i3+3)-boxh(3) )
                  duz = abs( x(i3+3)-boxh(3) )
                  if ( duz.gt.boxh(3) ) then
                      write(6,'(t2,a,i5,a)')
     .                      '%SANDER-W-ESCAPE, atom # ',int(i3/3+1),
     .                      ' has escaped'
                  endif
                  ns = 1
              endif

              if ( (dux*boct(1)+duy*boct(2)+duz).gt.boct(3) ) then
c
c                 since all atoms are now scaled into the cubic box,
c                 no escape from the inner octahedral box and its
c                 primary image is possible
c
                  x(i3+1) = x(i3+1) - sign( boxh(1),x(i3+1)-boxh(1) )
                  x(i3+2) = x(i3+2) - sign( boxh(2),x(i3+2)-boxh(2) )
                  x(i3+3) = x(i3+3) - sign( boxh(3),x(i3+3)-boxh(3) )
                  ns=1
              endif             
              I3 = I3+3
          endif
C
c         --- now translate the rest of the atoms of this molecule:
c
          IF (NS.ne.0.and.NRAM.gt.1) then
              I3 = I3-3
              DO 84 M = 1,3
                  I3 = I3+1
   84             XI(M) = X(I3)-XI(M)
              DO 86 J = 2,NRAM
                  DO 85 M = 1,3
                      I3 = I3+1
   85                 X(I3) = X(I3)+XI(M)
   86         CONTINUE
          else
              I3 = I3+NRAM3
          endif
C
   90 CONTINUE
C
C     ----- PERFORM THE COORDINATE TRANSFORMATION, WHEN REQUIRED ----
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine bound1a(npr,fw,itst,xwij)
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
c
c     bound1a() - like bound2a(), but tests an array element in the loop
c
#ifdef DPREC
      IMPLICIT double precision (A-H,O-Z)
#endif
#ifdef ISTAR2
      integer*2 itst
#endif
      DIMENSION XWIJ(3,*),fw(3,*)
c      dimension r2(*)
      dimension itst(*)
#include "box.h"
      data zero, one /0.0d0, 1.0d0/
C
      if (ifbox.eq.1) then
c
c         -- rectilinear box
c
cforcevector
          DO 200 JN = 1,NPR
              IF(itst(JN).gt.IPTATM) then
#ifdef USENINT
#  ifdef DPREC
                  xwij(1,jn) = xwij(1,jn)
     +                        - dnint(fw(1,jn)*boxi(1))*box(1)
                  xwij(2,jn) = xwij(2,jn)
     +                        - dnint(fw(2,jn)*boxi(2))*box(2)
                  xwij(3,jn) = xwij(3,jn)
     +                        - dnint(fw(3,jn)*boxi(3))*box(3)
#  else
                  xwij(1,jn) = xwij(1,jn)
     +                        - anint(fw(1,jn)*boxi(1))*box(1)
                  xwij(2,jn) = xwij(2,jn)
     +                        - anint(fw(2,jn)*boxi(2))*box(2)
                  xwij(3,jn) = xwij(3,jn)
     +                        - anint(fw(3,jn)*boxi(3))*box(3)
#  endif
#else
                  IF (ABS(fw(1,JN)).GT.BOXH(1))
     +                XWIJ(1,JN) = XWIJ(1,JN)-SIGN(BOX(1),fw(1,JN))
                  IF (ABS(fw(2,JN)).GT.BOXH(2))
     +                XWIJ(2,JN) = XWIJ(2,JN)-SIGN(BOX(2),fw(2,JN))
                  IF (ABS(fw(3,JN)).GT.BOXH(3))
     +                XWIJ(3,JN) = XWIJ(3,JN)-SIGN(BOX(3),fw(3,JN))
#endif
              endif
c     Commented out - calculate in subroutine nonbon
c              r2(JN) = ONE/
c     *             (XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2)
  200     CONTINUE
      else
c
c         -- truncated octahedral box
c
cforcevector
          DO 201 JN = 1,NPR
              IF(itst(JN).gt.IPTATM) then
#ifdef USENINT
                  dwx = fw(1,jn)*boxi(1)
                  dwy = fw(2,jn)*boxi(2)
                  dwz = fw(3,jn)*boxi(3)
                  dux = xwij(1,jn)*boxi(1)
                  duy = xwij(2,jn)*boxi(2)
                  duz = xwij(3,jn)*boxi(3)
# ifdef DPREC
                  dux = dux - dnint( dwx )
                  duy = duy - dnint( dwy )
                  duz = duz - dnint( dwz )
                  dwx = dwx - dnint( dwx )
                  dwy = dwy - dnint( dwy )
                  dwz = dwz - dnint( dwz )
# else
                  dux = dux - anint( dwx )
                  duy = duy - anint( dwy )
                  duz = duz - anint( dwz )
                  dwx = dwx - anint( dwx )
                  dwy = dwy - anint( dwy )
                  dwz = dwz - anint( dwz )
# endif
                  corr= 0.5 * aint( fothi* (abs( dwx ) +
     .                                     abs( dwy ) +
     .                                     abs( dwz )))
                  dux = dux - sign( corr,dwx )
                  duy = duy - sign( corr,dwy )
                  duz = duz - sign( corr,dwz )
                  xwij(1,jn) = dux * box(1)
                  xwij(2,jn) = duy * box(2)
                  xwij(3,jn) = duz * box(3)
#else
                  dwx = abs(fw(1,jn))
                  dwy = abs(fw(2,jn))
                  dwz = abs(fw(3,jn))
                  if ( dwx.gt.boxh(1) ) then
                      corr = - sign( box(1),fw(1,jn) )
                      xwij(1,jn) = xwij(1,jn) + corr
                      fw(1,jn) = fw(1,jn) + corr
                      dwx  = abs( fw(1,jn) )
                  endif
                  if( dwy.gt.boxh(2) ) then
                      corr = - sign( box(2),fw(2,jn) )
                      xwij(2,jn) = xwij(2,jn) + corr
                      fw(2,jn) = fw(2,jn) + corr
                      dwy  = abs( fw(2,jn) )
                  endif
                  if( dwz.gt.boxh(3) ) then
                      corr = - sign( box(3),fw(3,jn) )
                      xwij(3,jn) = xwij(3,jn) + corr
                      fw(3,jn) = fw(3,jn) + corr
                      dwz  = abs( fw(3,jn) )
                  endif
                  if( (dwx*boct(1)+dwy*boct(2)+dwz).gt.boct(3) ) then
                      xwij(1,jn) =
     +                      xwij(1,jn) - sign( boxh(1),fw(1,jn) )
                      xwij(2,jn) =
     +                      xwij(2,jn) - sign( boxh(2),fw(2,jn) )
                      xwij(3,jn) =
     +                      xwij(3,jn) - sign( boxh(3),fw(3,jn) )
                  endif
#endif
              endif
c     Commented out - calculate in subroutine nonbon
c              r2(JN) = ONE / 
c     *             (XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2)
  201     CONTINUE
      endif
      return
      end

c-----------------------------------------------------------------------
      subroutine bound2a(npr,fw,xwij,ibctype)
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
c
c     bound2a() - like bound1a(), but does not test an array element 
c                in the loop
c
#ifdef DPREC
      IMPLICIT double precision (A-H,O-Z)
#endif
      DIMENSION fw(3,*),XWIJ(3,*)
c      dimension R2(*)
      integer ibctype
#include "box.h"
      data zero, one /0.0d0, 1.0d0/
C
      if (ifbox.eq.1) then
c
c         -- rectilinear box
c
#ifdef USENINT
cforcevector
          DO 210 JN = 1,NPR
#  ifdef DPREC
              xwij(1,jn) = xwij(1,jn) - dnint(fw(1,jn)*boxi(1))*box(1)
              xwij(2,jn) = xwij(2,jn) - dnint(fw(2,jn)*boxi(2))*box(2)
              xwij(3,jn) = xwij(3,jn) - dnint(fw(3,jn)*boxi(3))*box(3)
#  else
              xwij(1,jn) = xwij(1,jn) - anint(fw(1,jn)*boxi(1))*box(1)
              xwij(2,jn) = xwij(2,jn) - anint(fw(2,jn)*boxi(2))*box(2)
              xwij(3,jn) = xwij(3,jn) - anint(fw(3,jn)*boxi(3))*box(3)
#  endif

c     Commented out - calculate in subroutine nonbon
c              r2(JN) =
c     +               ONE/(XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2)
  210     CONTINUE
#else
cforcevector
c         select case (ibctype)
c         case(1)
          if( ibctype .eq. 1 ) then
             DO JN = 1,NPR
                IF (ABS(fw(3,JN)).GT.BOXH(3))
     +               XWIJ(3,JN) = XWIJ(3,JN)-SIGN(BOX(3),fw(3,JN))
             enddo
c         case(2)
          else if( ibctype .eq. 2 ) then
             DO JN = 1,NPR
                IF (ABS(fw(2,JN)).GT.BOXH(2))
     +               XWIJ(2,JN) = XWIJ(2,JN)-SIGN(BOX(2),fw(2,JN))
             enddo
c         case(3)
          else if( ibctype .eq. 3 ) then
             DO JN = 1,NPR
                IF (ABS(fw(2,JN)).GT.BOXH(2))
     +               XWIJ(2,JN) = XWIJ(2,JN)-SIGN(BOX(2),fw(2,JN))
                IF (ABS(fw(3,JN)).GT.BOXH(3))
     +               XWIJ(3,JN) = XWIJ(3,JN)-SIGN(BOX(3),fw(3,JN))
             enddo
c         case(4)
          else if( ibctype .eq. 4 ) then
             DO JN = 1,NPR
                IF (ABS(fw(1,JN)).GT.BOXH(1))
     +               XWIJ(1,JN) = XWIJ(1,JN)-SIGN(BOX(1),fw(1,JN))
             enddo
c         case(5)
          else if( ibctype .eq. 5 ) then
             DO JN = 1,NPR
                IF (ABS(fw(1,JN)).GT.BOXH(1))
     +               XWIJ(1,JN) = XWIJ(1,JN)-SIGN(BOX(1),fw(1,JN))
                IF (ABS(fw(3,JN)).GT.BOXH(3))
     +               XWIJ(3,JN) = XWIJ(3,JN)-SIGN(BOX(3),fw(3,JN))
             enddo
c         case(6)
          else if( ibctype .eq. 6 ) then
             DO JN = 1,NPR
                IF (ABS(fw(1,JN)).GT.BOXH(1))
     +               XWIJ(1,JN) = XWIJ(1,JN)-SIGN(BOX(1),fw(1,JN))
                IF (ABS(fw(2,JN)).GT.BOXH(2))
     +               XWIJ(2,JN) = XWIJ(2,JN)-SIGN(BOX(2),fw(2,JN))
             enddo
c         case(7)
          else if( ibctype .eq. 7 ) then
             DO JN = 1,NPR
                IF (ABS(fw(1,JN)).GT.BOXH(1))
     +               XWIJ(1,JN) = XWIJ(1,JN)-SIGN(BOX(1),fw(1,JN))
                IF (ABS(fw(2,JN)).GT.BOXH(2))
     +               XWIJ(2,JN) = XWIJ(2,JN)-SIGN(BOX(2),fw(2,JN))
                IF (ABS(fw(3,JN)).GT.BOXH(3))
     +               XWIJ(3,JN) = XWIJ(3,JN)-SIGN(BOX(3),fw(3,JN))
             enddo
          end if
c         end select
#endif
       else
c
c         -- truncated octahedral box
c
cforcevector
          DO 200 JN = 1,NPR
#ifdef USENINT
              dwx = fw(1,jn)*boxi(1)
              dwy = fw(2,jn)*boxi(2)
              dwz = fw(3,jn)*boxi(3)
              dux = xwij(1,jn)*boxi(1)
              duy = xwij(2,jn)*boxi(2)
              duz = xwij(3,jn)*boxi(3)
# ifdef DPREC
              dux = dux - dnint( dwx )
              duy = duy - dnint( dwy )
              duz = duz - dnint( dwz )
              dwx = dwx - dnint( dwx )
              dwy = dwy - dnint( dwy )
              dwz = dwz - dnint( dwz )
# else
              dux = dux - anint( dwx )
              duy = duy - anint( dwy )
              duz = duz - anint( dwz )
              dwx = dwx - anint( dwx )
              dwy = dwy - anint( dwy )
              dwz = dwz - anint( dwz )
# endif
              corr= 0.5 * aint( fothi* (abs( dwx ) +
     +                                 abs( dwy ) +
     +                                 abs( dwz )))
              dux = dux - sign( corr,dwx )
              duy = duy - sign( corr,dwy )
              duz = duz - sign( corr,dwz )
              xwij(1,jn) = dux * box(1)
              xwij(2,jn) = duy * box(2)
              xwij(3,jn) = duz * box(3)
#else
              dwx = abs(fw(1,jn))
              dwy = abs(fw(2,jn))
              dwz = abs(fw(3,jn))
              if( dwx.gt.boxh(1) ) then
                  corr = - sign( box(1),fw(1,jn) )
                  xwij(1,jn) = xwij(1,jn) + corr
                  fw(1,jn) = fw(1,jn) + corr
                  dwx  = abs( fw(1,jn) )
              endif
              if( dwy.gt.boxh(2) ) then
                  corr = - sign( box(2),fw(2,jn) )
                  xwij(2,jn) = xwij(2,jn) + corr
                  fw(2,jn) = fw(2,jn) + corr
                  dwy  = abs( fw(2,jn) )
              endif
              if( dwz.gt.boxh(3) ) then
                  corr = - sign( box(3),fw(3,jn) )
                  xwij(3,jn) = xwij(3,jn) + corr
                  fw(3,jn) = fw(3,jn) + corr
                  dwz  = abs( fw(3,jn) )
              endif
              if( (dwx*boct(1)+dwy*boct(2)+dwz).gt.boct(3) ) then
                  xwij(1,jn) =
     +                      xwij(1,jn) - sign( boxh(1),fw(1,jn) )
                  xwij(2,jn) =
     +                      xwij(2,jn) - sign( boxh(2),fw(2,jn) )
                  xwij(3,jn) =
     +                      xwij(3,jn) - sign( boxh(3),fw(3,jn) )
              endif
#endif

c     Commented out - calculate in subroutine nonbon
c              r2(JN) = ONE/
c     *             (XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2)
  200     CONTINUE
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine bound1(npr,fw,itst,xwij,r2)
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
c
c     bound1() - like bound2(), but tests an array element in the loop
c
#ifdef DPREC
      IMPLICIT double precision (A-H,O-Z)
#endif
#ifdef ISTAR2
      integer*2 itst
#endif
      DIMENSION XWIJ(3,*),fw(3,*)
      dimension r2(*)
      dimension itst(*)
#include "box.h"
      data zero, one /0.0d0, 1.0d0/
C
      if (ifbox.eq.1) then
c
c         -- rectilinear box
c
cforcevector
          DO 200 JN = 1,NPR
              IF(itst(JN).gt.IPTATM) then
#ifdef USENINT
#  ifdef DPREC
                  xwij(1,jn) = xwij(1,jn)
     +                        - dnint(fw(1,jn)*boxi(1))*box(1)
                  xwij(2,jn) = xwij(2,jn)
     +                        - dnint(fw(2,jn)*boxi(2))*box(2)
                  xwij(3,jn) = xwij(3,jn)
     +                        - dnint(fw(3,jn)*boxi(3))*box(3)
#  else
                  xwij(1,jn) = xwij(1,jn)
     +                        - anint(fw(1,jn)*boxi(1))*box(1)
                  xwij(2,jn) = xwij(2,jn)
     +                        - anint(fw(2,jn)*boxi(2))*box(2)
                  xwij(3,jn) = xwij(3,jn)
     +                        - anint(fw(3,jn)*boxi(3))*box(3)
#  endif
#else
                  IF (ABS(fw(1,JN)).GT.BOXH(1))
     +                XWIJ(1,JN) = XWIJ(1,JN)-SIGN(BOX(1),fw(1,JN))
                  IF (ABS(fw(2,JN)).GT.BOXH(2))
     +                XWIJ(2,JN) = XWIJ(2,JN)-SIGN(BOX(2),fw(2,JN))
                  IF (ABS(fw(3,JN)).GT.BOXH(3))
     +                XWIJ(3,JN) = XWIJ(3,JN)-SIGN(BOX(3),fw(3,JN))
#endif
              endif
              r2(JN) = ONE/
     *             (XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2)
  200     CONTINUE
      else
c
c         -- truncated octahedral box
c
cforcevector
          DO 201 JN = 1,NPR
              IF(itst(JN).gt.IPTATM) then
#ifdef USENINT
                  dwx = fw(1,jn)*boxi(1)
                  dwy = fw(2,jn)*boxi(2)
                  dwz = fw(3,jn)*boxi(3)
                  dux = xwij(1,jn)*boxi(1)
                  duy = xwij(2,jn)*boxi(2)
                  duz = xwij(3,jn)*boxi(3)
# ifdef DPREC
                  dux = dux - dnint( dwx )
                  duy = duy - dnint( dwy )
                  duz = duz - dnint( dwz )
                  dwx = dwx - dnint( dwx )
                  dwy = dwy - dnint( dwy )
                  dwz = dwz - dnint( dwz )
# else
                  dux = dux - anint( dwx )
                  duy = duy - anint( dwy )
                  duz = duz - anint( dwz )
                  dwx = dwx - anint( dwx )
                  dwy = dwy - anint( dwy )
                  dwz = dwz - anint( dwz )
# endif
                  corr= 0.5 * aint( fothi* (abs( dwx ) +
     .                                     abs( dwy ) +
     .                                     abs( dwz )))
                  dux = dux - sign( corr,dwx )
                  duy = duy - sign( corr,dwy )
                  duz = duz - sign( corr,dwz )
                  xwij(1,jn) = dux * box(1)
                  xwij(2,jn) = duy * box(2)
                  xwij(3,jn) = duz * box(3)
#else
                  dwx = abs(fw(1,jn))
                  dwy = abs(fw(2,jn))
                  dwz = abs(fw(3,jn))
                  if ( dwx.gt.boxh(1) ) then
                      corr = - sign( box(1),fw(1,jn) )
                      xwij(1,jn) = xwij(1,jn) + corr
                      fw(1,jn) = fw(1,jn) + corr
                      dwx  = abs( fw(1,jn) )
                  endif
                  if( dwy.gt.boxh(2) ) then
                      corr = - sign( box(2),fw(2,jn) )
                      xwij(2,jn) = xwij(2,jn) + corr
                      fw(2,jn) = fw(2,jn) + corr
                      dwy  = abs( fw(2,jn) )
                  endif
                  if( dwz.gt.boxh(3) ) then
                      corr = - sign( box(3),fw(3,jn) )
                      xwij(3,jn) = xwij(3,jn) + corr
                      fw(3,jn) = fw(3,jn) + corr
                      dwz  = abs( fw(3,jn) )
                  endif
                  if( (dwx*boct(1)+dwy*boct(2)+dwz).gt.boct(3) ) then
                      xwij(1,jn) =
     +                      xwij(1,jn) - sign( boxh(1),fw(1,jn) )
                      xwij(2,jn) =
     +                      xwij(2,jn) - sign( boxh(2),fw(2,jn) )
                      xwij(3,jn) =
     +                      xwij(3,jn) - sign( boxh(3),fw(3,jn) )
                  endif
#endif
              endif
              r2(JN) = ONE / 
     *             (XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2)
  201     CONTINUE
      endif
      return
      end

c-----------------------------------------------------------------------
      subroutine bound2(npr,fw,xwij,r2)
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
c
c     bound2() - like bound1(), but does not test an array element 
c                in the loop
c
#ifdef DPREC
      IMPLICIT double precision (A-H,O-Z)
#endif
      DIMENSION fw(3,*),XWIJ(3,*)
      dimension R2(*)
#include "box.h"
      data zero, one /0.0d0, 1.0d0/
C
      if (ifbox.eq.1) then
c
c         -- rectilinear box
c
#ifdef USENINT
cforcevector
          DO 210 JN = 1,NPR
#  ifdef DPREC
              xwij(1,jn) = xwij(1,jn) - dnint(fw(1,jn)*boxi(1))*box(1)
              xwij(2,jn) = xwij(2,jn) - dnint(fw(2,jn)*boxi(2))*box(2)
              xwij(3,jn) = xwij(3,jn) - dnint(fw(3,jn)*boxi(3))*box(3)
#  else
              xwij(1,jn) = xwij(1,jn) - anint(fw(1,jn)*boxi(1))*box(1)
              xwij(2,jn) = xwij(2,jn) - anint(fw(2,jn)*boxi(2))*box(2)
              xwij(3,jn) = xwij(3,jn) - anint(fw(3,jn)*boxi(3))*box(3)
#  endif

              r2(JN) =
     +               ONE/(XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2)
  210     CONTINUE
#else
cforcevector
          DO 210 JN = 1,NPR
              IF (ABS(fw(1,JN)).GT.BOXH(1))
     +                 XWIJ(1,JN) = XWIJ(1,JN)-SIGN(BOX(1),fw(1,JN))
              IF (ABS(fw(2,JN)).GT.BOXH(2))
     +                 XWIJ(2,JN) = XWIJ(2,JN)-SIGN(BOX(2),fw(2,JN))
              IF (ABS(fw(3,JN)).GT.BOXH(3))
     +                 XWIJ(3,JN) = XWIJ(3,JN)-SIGN(BOX(3),fw(3,JN))

              r2(JN) =
     +               ONE/(XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2)
  210     CONTINUE
#endif
      else
c
c         -- truncated octahedral box
c
cforcevector
          DO 200 JN = 1,NPR
#ifdef USENINT
              dwx = fw(1,jn)*boxi(1)
              dwy = fw(2,jn)*boxi(2)
              dwz = fw(3,jn)*boxi(3)
              dux = xwij(1,jn)*boxi(1)
              duy = xwij(2,jn)*boxi(2)
              duz = xwij(3,jn)*boxi(3)
# ifdef DPREC
              dux = dux - dnint( dwx )
              duy = duy - dnint( dwy )
              duz = duz - dnint( dwz )
              dwx = dwx - dnint( dwx )
              dwy = dwy - dnint( dwy )
              dwz = dwz - dnint( dwz )
# else
              dux = dux - anint( dwx )
              duy = duy - anint( dwy )
              duz = duz - anint( dwz )
              dwx = dwx - anint( dwx )
              dwy = dwy - anint( dwy )
              dwz = dwz - anint( dwz )
# endif
              corr= 0.5 * aint( fothi* (abs( dwx ) +
     +                                 abs( dwy ) +
     +                                 abs( dwz )))
              dux = dux - sign( corr,dwx )
              duy = duy - sign( corr,dwy )
              duz = duz - sign( corr,dwz )
              xwij(1,jn) = dux * box(1)
              xwij(2,jn) = duy * box(2)
              xwij(3,jn) = duz * box(3)
#else
              dwx = abs(fw(1,jn))
              dwy = abs(fw(2,jn))
              dwz = abs(fw(3,jn))
              if( dwx.gt.boxh(1) ) then
                  corr = - sign( box(1),fw(1,jn) )
                  xwij(1,jn) = xwij(1,jn) + corr
                  fw(1,jn) = fw(1,jn) + corr
                  dwx  = abs( fw(1,jn) )
              endif
              if( dwy.gt.boxh(2) ) then
                  corr = - sign( box(2),fw(2,jn) )
                  xwij(2,jn) = xwij(2,jn) + corr
                  fw(2,jn) = fw(2,jn) + corr
                  dwy  = abs( fw(2,jn) )
              endif
              if( dwz.gt.boxh(3) ) then
                  corr = - sign( box(3),fw(3,jn) )
                  xwij(3,jn) = xwij(3,jn) + corr
                  fw(3,jn) = fw(3,jn) + corr
                  dwz  = abs( fw(3,jn) )
              endif
              if( (dwx*boct(1)+dwy*boct(2)+dwz).gt.boct(3) ) then
                  xwij(1,jn) =
     +                      xwij(1,jn) - sign( boxh(1),fw(1,jn) )
                  xwij(2,jn) =
     +                      xwij(2,jn) - sign( boxh(2),fw(2,jn) )
                  xwij(3,jn) =
     +                      xwij(3,jn) - sign( boxh(3),fw(3,jn) )
              endif
#endif

              r2(JN) = ONE/
     *             (XWIJ(1,JN)**2+XWIJ(2,JN)**2+XWIJ(3,JN)**2)
  200     CONTINUE
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine bound3(nn, xi, xr, j3 )
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
c
#ifdef DPREC 
      implicit double precision (a-h,o-z)
#endif
      DIMENSION XR(*)
      DIMENSION XI(3),XIJ(3)
#include "box.h"
      data zero, one /0.0d0, 1.0d0/
C
      if (ifbox.eq.1) then
c
c         -- rectilinear box
c
          xi1 = xi(1)
          xi2 = xi(2)
          xi3 = xi(3)

          DO 100 J = 1,NN
C
             xij1 = xr(J3+1) - xi1
             xij2 = xr(J3+2) - xi2
             xij3 = xr(J3+3) - xi3

             if ( xij1.ge. boxh(1) ) xij1 = xij1 - boxh(1)
             if ( xij1.lt.-boxh(1) ) xij1 = xij1 + boxh(1)

             if ( xij2.ge. boxh(2) ) xij2 = xij2 - boxh(2)
             if ( xij2.lt.-boxh(2) ) xij2 = xij2 + boxh(2)

             if ( xij3.ge. boxh(3) ) xij3 = xij3 - boxh(3)
             if ( xij3.lt.-boxh(3) ) xij3 = xij3 + boxh(3)

             xi1 = xi1 + xij1
             xi2 = xi2 + xij2
             xi3 = xi3 + xij3

             xr(J3+1) = xi1
             xr(J3+2) = xi2
             xr(J3+3) = xi3

             J3 = J3 + 3

  100     CONTINUE

          xi(1) = xi1
          xi(2) = xi2
          xi(3) = xi3

      else
c
c         -- truncated octahedral box
c
cforcevector
          do 22 j=1,nn
              XIJ(1) = XR(J3+1)-XI(1)
              XIJ(2) = XR(J3+2)-XI(2)
              XIJ(3) = XR(J3+3)-XI(3)
#ifdef USENINT
              dux = xij(1)*boxi(1)
              duy = xij(2)*boxi(2)
              duz = xij(3)*boxi(3)
# ifdef DPREC
              dux = dux - dnint( dux )
              duy = duy - dnint( duy )
              duz = duz - dnint( duz )
#  else
              dux = dux - anint( dux )
              duy = duy - anint( duy )
              duz = duz - anint( duz )
#  endif
              corr= 0.5 * aint( fothi* (abs( dux ) +
     .                                 abs( duy ) + 
     .                                 abs( duz )))
              dux = dux - sign( corr,dux )
              duy = duy - sign( corr,duy )
              duz = duz - sign( corr,duz )
              xij(1) = dux * box(1)
              xij(2) = duy * box(2)
              xij(3) = duz * box(3)
#else
              dux = abs(xij(1))
              duy = abs(xij(2))
              duz = abs(xij(3))
              if( dux.gt.boxh(1) ) then
                  xij(1) = xij(1) - sign( box(1),xij(1) )
                  dux  = abs( xij(1) )
              endif
              if( duy.gt.boxh(2) ) then
                  xij(2) = xij(2) - sign( box(2),xij(2) )
                  duy  = abs( xij(2) )
              endif
              if( duz.gt.boxh(3) ) then
                  xij(3) = xij(3) - sign( box(3),xij(3) )
                  duz  = abs( xij(3) )
              endif
              if( (dux*boct(1)+duy*boct(2)+duz).gt.boct(3) ) then
                  xij(1) = xij(1) - sign( boxh(1),xij(1) )
                  xij(2) = xij(2) - sign( boxh(2),xij(2) )
                  xij(3) = xij(3) - sign( boxh(3),xij(3) )
              endif
#endif
              DO 24 M = 1,3
                  J3 = J3+1
                  XH = XI(M)+XIJ(M)
                  XI(M) = XH
                  XR(J3) = XH
   24         CONTINUE
   22     continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine bound6(xwij, fw, lpr, npr)
c
#ifdef DPREC
      implicit double precision (a-h,o-z)
#endif
c
#include "box.h"
      dimension xwij(9,*), fw(3,*)
      data zero,one,two,six/0.0d0,1.0d0,2.0d0,6.0d0/

      if (ifbox.eq.1) then

         LN = LPR + 1

c
c         -- orthogonal boundaries
c
#ifdef USENINT
cforcevector
c     ...the following Cray directive is necessary to overcome a
c     compiler optimization bug...
cfpp$ skip
         DO 290 JN = 1,LPR 
            
#  ifdef DPREC
            FW1 = DNINT(XWIJ(1,JN)*BOXi(1))*BOX(1)
            FW2 = DNINT(XWIJ(2,JN)*BOXi(2))*BOX(2)
            FW3 = DNINT(XWIJ(3,JN)*BOXi(3))*BOX(3)
#  else
            FW1 = ANINT(XWIJ(1,JN)*BOXi(1))*BOX(1)
            FW2 = ANINT(XWIJ(2,JN)*BOXi(2))*BOX(2)
            FW3 = ANINT(XWIJ(3,JN)*BOXi(3))*BOX(3)
#  endif

               XWIJ(1,JN) = XWIJ(1,JN) - FW1
               XWIJ(2,JN) = XWIJ(2,JN) - FW2
               XWIJ(3,JN) = XWIJ(3,JN) - FW3
               XWIJ(4,JN) = XWIJ(4,JN) - FW1
               XWIJ(5,JN) = XWIJ(5,JN) - FW2
               XWIJ(6,JN) = XWIJ(6,JN) - FW3
               XWIJ(7,JN) = XWIJ(7,JN) - FW1
               XWIJ(8,JN) = XWIJ(8,JN) - FW2
               XWIJ(9,JN) = XWIJ(9,JN) - FW3

               XWIJ(1,2*JN+LPR) = XWIJ(1,2*JN+LPR) - FW1
               XWIJ(2,2*JN+LPR) = XWIJ(2,2*JN+LPR) - FW2
               XWIJ(3,2*JN+LPR) = XWIJ(3,2*JN+LPR) - FW3
               XWIJ(4,2*JN+LPR) = XWIJ(4,2*JN+LPR) - FW1
               XWIJ(5,2*JN+LPR) = XWIJ(5,2*JN+LPR) - FW2
               XWIJ(6,2*JN+LPR) = XWIJ(6,2*JN+LPR) - FW3
               XWIJ(7,2*JN+LPR) = XWIJ(7,2*JN+LPR) - FW1
               XWIJ(8,2*JN+LPR) = XWIJ(8,2*JN+LPR) - FW2
               XWIJ(9,2*JN+LPR) = XWIJ(9,2*JN+LPR) - FW3

               XWIJ(1,2*JN+LPR-1) = XWIJ(1,2*JN+LPR-1) - FW1
               XWIJ(2,2*JN+LPR-1) = XWIJ(2,2*JN+LPR-1) - FW2
               XWIJ(3,2*JN+LPR-1) = XWIJ(3,2*JN+LPR-1) - FW3
               XWIJ(4,2*JN+LPR-1) = XWIJ(4,2*JN+LPR-1) - FW1
               XWIJ(5,2*JN+LPR-1) = XWIJ(5,2*JN+LPR-1) - FW2
               XWIJ(6,2*JN+LPR-1) = XWIJ(6,2*JN+LPR-1) - FW3
               XWIJ(7,2*JN+LPR-1) = XWIJ(7,2*JN+LPR-1) - FW1
               XWIJ(8,2*JN+LPR-1) = XWIJ(8,2*JN+LPR-1) - FW2
               XWIJ(9,2*JN+LPR-1) = XWIJ(9,2*JN+LPR-1) - FW3

  290     CONTINUE
#else
C
C         This ifdef block performs periodic imaging using the 
c         "standard" method, where xi-xj is compared to 1/2 boxlength.
C
C         -- x-translations:
C
cforcevector
c     ...the following Cray directive is necessary to overcome a
c     compiler optimization bug...
cfpp$ skip
          DO 320 JN = 1,LPR

               IF (ABS(XWIJ(1,JN)).GT.BOXH(1)) THEN
                  SHIFT = SIGN(BOX(1),XWIJ(1,JN))
                  XWIJ(1,JN) = XWIJ(1,JN)-SHIFT
                  XWIJ(4,JN) = XWIJ(4,JN)-SHIFT
                  XWIJ(7,JN) = XWIJ(7,JN)-SHIFT
                  XWIJ(1,2*JN+LPR) = XWIJ(1,2*JN+LPR)-SHIFT
                  XWIJ(4,2*JN+LPR) = XWIJ(4,2*JN+LPR)-SHIFT
                  XWIJ(7,2*JN+LPR) = XWIJ(7,2*JN+LPR)-SHIFT
                  XWIJ(1,2*JN+LPR-1) = XWIJ(1,2*JN+LPR-1)-SHIFT
                  XWIJ(4,2*JN+LPR-1) = XWIJ(4,2*JN+LPR-1)-SHIFT
                  XWIJ(7,2*JN+LPR-1) = XWIJ(7,2*JN+LPR-1)-SHIFT
               END IF

               IF (ABS(XWIJ(2,JN)).GT.BOXH(2)) THEN
                  SHIFT = SIGN(BOX(2),XWIJ(2,JN))
                  XWIJ(2,JN) = XWIJ(2,JN)-SHIFT
                  XWIJ(5,JN) = XWIJ(5,JN)-SHIFT
                  XWIJ(8,JN) = XWIJ(8,JN)-SHIFT
                  XWIJ(2,2*JN+LPR) = XWIJ(2,2*JN+LPR)-SHIFT
                  XWIJ(5,2*JN+LPR) = XWIJ(5,2*JN+LPR)-SHIFT
                  XWIJ(8,2*JN+LPR) = XWIJ(8,2*JN+LPR)-SHIFT
                  XWIJ(2,2*JN+LPR-1) = XWIJ(2,2*JN+LPR-1)-SHIFT
                  XWIJ(5,2*JN+LPR-1) = XWIJ(5,2*JN+LPR-1)-SHIFT
                  XWIJ(8,2*JN+LPR-1) = XWIJ(8,2*JN+LPR-1)-SHIFT
               END IF

               IF (ABS(XWIJ(3,JN)).GT.BOXH(3)) THEN
                  SHIFT = SIGN(BOX(3),XWIJ(3,JN))
                  XWIJ(3,JN) = XWIJ(3,JN)-SHIFT
                  XWIJ(6,JN) = XWIJ(6,JN)-SHIFT
                  XWIJ(9,JN) = XWIJ(9,JN)-SHIFT
                  XWIJ(3,2*JN+LPR) = XWIJ(3,2*JN+LPR)-SHIFT
                  XWIJ(6,2*JN+LPR) = XWIJ(6,2*JN+LPR)-SHIFT
                  XWIJ(9,2*JN+LPR) = XWIJ(9,2*JN+LPR)-SHIFT
                  XWIJ(3,2*JN+LPR-1) = XWIJ(3,2*JN+LPR-1)-SHIFT
                  XWIJ(6,2*JN+LPR-1) = XWIJ(6,2*JN+LPR-1)-SHIFT
                  XWIJ(9,2*JN+LPR-1) = XWIJ(9,2*JN+LPR-1)-SHIFT
               END IF

  320     CONTINUE
#endif

      else

c
c         -- truncated octahedral boundaries
c
#ifdef USENINT
cforcevector
c     ...the following Cray directive is necessary to overcome a
c     compiler optimization bug...
cfpp$ skip
          DO 291 JN = 1,LPR
               dwx=XWIJ(1,JN)*BOXi(1)
               dwy=XWIJ(2,JN)*BOXi(2)
               dwz=XWIJ(3,JN)*BOXi(3)
#  ifdef DPREC
               FW(1,JN) = DNINT(dwx)
               FW(2,JN) = DNINT(dwy)
               FW(3,JN) = DNINT(dwz)
#  else
               FW(1,JN) = ANINT(dwx)
               FW(2,JN) = ANINT(dwy)
               FW(3,JN) = ANINT(dwz)
#  endif
               dwx=dwx-fw(1,jn)
               dwy=dwy-fw(2,jn)
               dwz=dwz-fw(3,jn)
               corr= 0.5 * aint( fothi* (abs( dwx ) +
     +                                  abs( dwy ) +
     +                                  abs( dwz )))
               fw(1,jn)=box(1)*(fw(1,jn)-sign(corr,dwx))
               fw(2,jn)=box(2)*(fw(2,jn)-sign(corr,dwy))
               fw(3,jn)=box(3)*(fw(3,jn)-sign(corr,dwz))
  291     CONTINUE
cforcevector
c     ...the following Cray directive is necessary to overcome a
c     compiler optimization bug...
cfpp$ skip
          DO 301 JN = 1,LPR
               XWIJ(1,JN) = XWIJ(1,JN) - FW(1,JN)
               XWIJ(4,JN) = XWIJ(4,JN) - FW(1,JN)
               XWIJ(7,JN) = XWIJ(7,JN) - FW(1,JN)
               XWIJ(2,JN) = XWIJ(2,JN) - FW(2,JN)
               XWIJ(5,JN) = XWIJ(5,JN) - FW(2,JN)
               XWIJ(8,JN) = XWIJ(8,JN) - FW(2,JN)
               XWIJ(3,JN) = XWIJ(3,JN) - FW(3,JN)
               XWIJ(6,JN) = XWIJ(6,JN) - FW(3,JN)
               XWIJ(9,JN) = XWIJ(9,JN) - FW(3,JN)
  301     CONTINUE
          KN = 0
cforcevector 
c     ...the following Cray directive is necessary to overcome a
c     compiler optimization bug...
cfpp$ skip
          DO 311 JN = LPR+1,NPR,2
               KN = KN +1
               XWIJ(1,JN) = XWIJ(1,JN)-FW(1,KN)
               XWIJ(4,JN) = XWIJ(4,JN)-FW(1,KN)
               XWIJ(7,JN) = XWIJ(7,JN)-FW(1,KN)
               XWIJ(2,JN) = XWIJ(2,JN)-FW(2,KN)
               XWIJ(5,JN) = XWIJ(5,JN)-FW(2,KN)
               XWIJ(8,JN) = XWIJ(8,JN)-FW(2,KN)
               XWIJ(3,JN) = XWIJ(3,JN)-FW(3,KN)
               XWIJ(6,JN) = XWIJ(6,JN)-FW(3,KN)
               XWIJ(9,JN) = XWIJ(9,JN)-FW(3,KN)
               XWIJ(1,JN+1) = XWIJ(1,JN+1)-FW(1,KN)
               XWIJ(4,JN+1) = XWIJ(4,JN+1)-FW(1,KN)
               XWIJ(7,JN+1) = XWIJ(7,JN+1)-FW(1,KN)
               XWIJ(2,JN+1) = XWIJ(2,JN+1)-FW(2,KN)
               XWIJ(5,JN+1) = XWIJ(5,JN+1)-FW(2,KN)
               XWIJ(8,JN+1) = XWIJ(8,JN+1)-FW(2,KN)
               XWIJ(3,JN+1) = XWIJ(3,JN+1)-FW(3,KN)
               XWIJ(6,JN+1) = XWIJ(6,JN+1)-FW(3,KN)
               XWIJ(9,JN+1) = XWIJ(9,JN+1)-FW(3,KN)
  311     CONTINUE
#else
c
c         "standard" periodic imaging
c
cforcevector
c     ...the following Cray directive is necessary to overcome a
c     compiler optimization bug...
cfpp$ skip
          DO 321 JN = 1,LPR
              dwx=abs(xwij(1,jn))
              dwy=abs(xwij(2,jn))
              dwz=abs(xwij(3,jn))
c
c             -- x-translations:
c
              IF (dwx.GT.BOXH(1)) THEN
                  SHIFT = SIGN(BOX(1),XWIJ(1,JN))
                  XWIJ(1,JN) = XWIJ(1,JN)-SHIFT
                  XWIJ(4,JN) = XWIJ(4,JN)-SHIFT
                  XWIJ(7,JN) = XWIJ(7,JN)-SHIFT
                  XWIJ(1,2*JN+LPR) = XWIJ(1,2*JN+LPR)-SHIFT
                  XWIJ(4,2*JN+LPR) = XWIJ(4,2*JN+LPR)-SHIFT
                  XWIJ(7,2*JN+LPR) = XWIJ(7,2*JN+LPR)-SHIFT
                  XWIJ(1,2*JN+LPR-1) = XWIJ(1,2*JN+LPR-1)-SHIFT
                  XWIJ(4,2*JN+LPR-1) = XWIJ(4,2*JN+LPR-1)-SHIFT
                  XWIJ(7,2*JN+LPR-1) = XWIJ(7,2*JN+LPR-1)-SHIFT
                  dwx=abs(xwij(1,jn))
               END IF
c
c              -- y-translations:
c
               IF (dwy.GT.BOXH(2)) THEN
                  SHIFT = SIGN(BOX(2),XWIJ(2,JN))
                  XWIJ(2,JN) = XWIJ(2,JN)-SHIFT
                  XWIJ(5,JN) = XWIJ(5,JN)-SHIFT
                  XWIJ(8,JN) = XWIJ(8,JN)-SHIFT
                  XWIJ(2,2*JN+LPR) = XWIJ(2,2*JN+LPR)-SHIFT
                  XWIJ(5,2*JN+LPR) = XWIJ(5,2*JN+LPR)-SHIFT
                  XWIJ(8,2*JN+LPR) = XWIJ(8,2*JN+LPR)-SHIFT
                  XWIJ(2,2*JN+LPR-1) = XWIJ(2,2*JN+LPR-1)-SHIFT
                  XWIJ(5,2*JN+LPR-1) = XWIJ(5,2*JN+LPR-1)-SHIFT
                  XWIJ(8,2*JN+LPR-1) = XWIJ(8,2*JN+LPR-1)-SHIFT
                  dwy=abs(xwij(2,jn))
               END IF
c
c              -- z-translations:
c
               IF (dwz.GT.BOXH(3)) THEN
                  SHIFT = SIGN(BOX(3),XWIJ(3,JN))
                  XWIJ(3,JN) = XWIJ(3,JN)-SHIFT
                  XWIJ(6,JN) = XWIJ(6,JN)-SHIFT
                  XWIJ(9,JN) = XWIJ(9,JN)-SHIFT
                  XWIJ(3,2*JN+LPR) = XWIJ(3,2*JN+LPR)-SHIFT
                  XWIJ(6,2*JN+LPR) = XWIJ(6,2*JN+LPR)-SHIFT
                  XWIJ(9,2*JN+LPR) = XWIJ(9,2*JN+LPR)-SHIFT
                  XWIJ(3,2*JN+LPR-1) = XWIJ(3,2*JN+LPR-1)-SHIFT
                  XWIJ(6,2*JN+LPR-1) = XWIJ(6,2*JN+LPR-1)-SHIFT
                  XWIJ(9,2*JN+LPR-1) = XWIJ(9,2*JN+LPR-1)-SHIFT
                  dwz=abs(xwij(3,jn))
               END IF
c
c              -- corrections
c
               if((dwx*boct(1)+dwy*boct(2)+dwz).gt.boct(3))then
c
c                 -- x correction
c
                  SHIFT = SIGN(BOXH(1),XWIJ(1,JN))
                  XWIJ(1,JN) = XWIJ(1,JN)-SHIFT
                  XWIJ(4,JN) = XWIJ(4,JN)-SHIFT
                  XWIJ(7,JN) = XWIJ(7,JN)-SHIFT
                  XWIJ(1,2*JN+LPR) = XWIJ(1,2*JN+LPR)-SHIFT
                  XWIJ(4,2*JN+LPR) = XWIJ(4,2*JN+LPR)-SHIFT
                  XWIJ(7,2*JN+LPR) = XWIJ(7,2*JN+LPR)-SHIFT
                  XWIJ(1,2*JN+LPR-1) = XWIJ(1,2*JN+LPR-1)-SHIFT
                  XWIJ(4,2*JN+LPR-1) = XWIJ(4,2*JN+LPR-1)-SHIFT
                  XWIJ(7,2*JN+LPR-1) = XWIJ(7,2*JN+LPR-1)-SHIFT
c
c                 -- y correction
c
                  SHIFT = SIGN(BOXH(2),XWIJ(2,JN))
                  XWIJ(2,JN) = XWIJ(2,JN)-SHIFT
                  XWIJ(5,JN) = XWIJ(5,JN)-SHIFT
                  XWIJ(8,JN) = XWIJ(8,JN)-SHIFT
                  XWIJ(2,2*JN+LPR) = XWIJ(2,2*JN+LPR)-SHIFT
                  XWIJ(5,2*JN+LPR) = XWIJ(5,2*JN+LPR)-SHIFT
                  XWIJ(8,2*JN+LPR) = XWIJ(8,2*JN+LPR)-SHIFT
                  XWIJ(2,2*JN+LPR-1) = XWIJ(2,2*JN+LPR-1)-SHIFT
                  XWIJ(5,2*JN+LPR-1) = XWIJ(5,2*JN+LPR-1)-SHIFT
                  XWIJ(8,2*JN+LPR-1) = XWIJ(8,2*JN+LPR-1)-SHIFT
c
c                 -- z correction
c
                  SHIFT = SIGN(BOXH(3),XWIJ(3,JN))
                  XWIJ(3,JN) = XWIJ(3,JN)-SHIFT
                  XWIJ(6,JN) = XWIJ(6,JN)-SHIFT
                  XWIJ(9,JN) = XWIJ(9,JN)-SHIFT
                  XWIJ(3,2*JN+LPR) = XWIJ(3,2*JN+LPR)-SHIFT
                  XWIJ(6,2*JN+LPR) = XWIJ(6,2*JN+LPR)-SHIFT
                  XWIJ(9,2*JN+LPR) = XWIJ(9,2*JN+LPR)-SHIFT
                  XWIJ(3,2*JN+LPR-1) = XWIJ(3,2*JN+LPR-1)-SHIFT
                  XWIJ(6,2*JN+LPR-1) = XWIJ(6,2*JN+LPR-1)-SHIFT
                  XWIJ(9,2*JN+LPR-1) = XWIJ(9,2*JN+LPR-1)-SHIFT
               endif
  321     continue
#endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine bound7(ipres,xci,yci,zci,x,rrw,i1,i2,ibctype)
#ifdef DPREC
      implicit double precision (a-h,o-z)
#endif
c
#include "box.h"
      dimension ipres(*), x(3,*), rrw(*)
      integer ibctype

      if (ifbox.eq.1) then
c
c         -- rectilinear boundary
c
cforcevector
c          DO 120 JJ = i1, i2
c              J = IPRES(JJ)
c              xw = xci - x(1,j)
c              yw = yci - x(2,j)
c              zw = zci - x(3,j)
c              IF (ABS(xw) .gt. BOXH(1)) xw = xw-SIGN(BOX(1),xw)
c              IF (ABS(yw) .gt. BOXH(2)) yw = yw-SIGN(BOX(2),yw)
c              IF (ABS(zw) .gt. BOXH(3)) zw = zw-SIGN(BOX(3),zw)
c              rrw(jj) = xw**2 + yw**2 + zw**2
c  120     CONTINUE

c         select case (ibctype)
c         case(0)
          if( ibctype .eq. 0 ) then
             DO JJ = I1, I2
                J = IPRES(JJ)
                xw = (xci - x(1,j))
                yw = (yci - x(2,j))
                zw = (zci - x(3,j))
                rrw(jj) = xw**2 + yw**2 + zw**2
             enddo
c         case(1)
          else if( ibctype .eq. 1 ) then
             DO JJ = I1, I2
                J = IPRES(JJ)
                xw = (xci - x(1,j))
                yw = (yci - x(2,j))
                zw = abs(zci - x(3,j))
                IF (zw .gt. BOXH(3)) zw = zw-BOX(3)
                rrw(jj) = xw**2 + yw**2 + zw**2
             enddo
c         case(2)
          else if( ibctype .eq. 2 ) then
             DO JJ = I1, I2
                J = IPRES(JJ)
                xw = (xci - x(1,j))
                yw = abs(yci - x(2,j))
                zw = (zci - x(3,j))
                IF (yw .gt. BOXH(2)) yw = yw-BOX(2)
                rrw(jj) = xw**2 + yw**2 + zw**2
             enddo
c         case(3)
          else if( ibctype .eq. 3 ) then
             DO JJ = I1, I2
                J = IPRES(JJ)
                xw = (xci - x(1,j))
                yw = abs(yci - x(2,j))
                zw = abs(zci - x(3,j))
                IF (yw .gt. BOXH(2)) yw = yw-BOX(2)
                IF (zw .gt. BOXH(3)) zw = zw-BOX(3)
                rrw(jj) = xw**2 + yw**2 + zw**2
             enddo
c         case(4)
          else if( ibctype .eq. 4 ) then
             DO JJ = I1, I2
                J = IPRES(JJ)
                xw = abs(xci - x(1,j))
                yw = (yci - x(2,j))
                zw = (zci - x(3,j))
                IF (xw .gt. BOXH(1)) xw = xw-BOX(1)
                rrw(jj) = xw**2 + yw**2 + zw**2
             enddo
c         case(5)
          else if( ibctype .eq. 5 ) then
             DO JJ = I1, I2
                J = IPRES(JJ)
                xw = abs(xci - x(1,j))
                yw = (yci - x(2,j))
                zw = abs(zci - x(3,j))
                IF (xw .gt. BOXH(1)) xw = xw-BOX(1)
                IF (zw .gt. BOXH(3)) zw = zw-BOX(3)
                rrw(jj) = xw**2 + yw**2 + zw**2
             enddo
c         case(6)
          else if( ibctype .eq. 6 ) then
             DO JJ = I1, I2
                J = IPRES(JJ)
                xw = abs(xci - x(1,j))
                yw = abs(yci - x(2,j))
                zw = (zci - x(3,j))
                IF (xw .gt. BOXH(1)) xw = xw-BOX(1)
                IF (yw .gt. BOXH(2)) yw = yw-BOX(2)
                rrw(jj) = xw**2 + yw**2 + zw**2
             enddo
c         case(7)
          else if( ibctype .eq. 7 ) then
             DO JJ = I1, I2
                J = IPRES(JJ)
                xw = abs(xci - x(1,j))
                yw = abs(yci - x(2,j))
                zw = abs(zci - x(3,j))
                IF (xw .gt. BOXH(1)) xw = xw-BOX(1)
                IF (yw .gt. BOXH(2)) yw = yw-BOX(2)
                IF (zw .gt. BOXH(3)) zw = zw-BOX(3)
                rrw(jj) = xw**2 + yw**2 + zw**2
             enddo
          end if
c         end select

      else


c
c         -- truncated octahedral boundary
c
cforcevector
          DO 121 JJ = i1, i2
              J = IPRES(JJ)
              xw = xci - x(1,j)
              yw = yci - x(2,j)
              zw = zci - x(3,j)
#ifdef USENINT
              dux = xw*boxi(1)
              duy = yw*boxi(2)
              duz = zw*boxi(3)
#  ifdef DPREC
              dux = dux - dnint( dux )
              duy = duy - dnint( duy )
              duz = duz - dnint( duz )
#  else
              dux = dux - anint( dux )
              duy = duy - anint( duy )
              duz = duz - anint( duz )
#  endif
              corr= 0.5 * aint( fothi* (abs( dux ) +
     +                                 abs( duy ) +
     +                                 abs( duz )))
              dux = dux - sign( corr,dux )
              duy = duy - sign( corr,duy )
              duz = duz - sign( corr,duz )
              xw = dux * box(1)
              yw = duy * box(2)
              zw = duz * box(3)
#else
              dux = abs(xw)
              duy = abs(yw)
              duz = abs(zw)
              if( dux.gt.boxh(1) ) then
                  xw = xw - sign( box(1),xw )
                  dux  = abs( xw )
              endif
              if( duy.gt.boxh(2) ) then
                  yw = yw - sign( box(2),yw )
                  duy  = abs( yw )
              endif
              if( duz.gt.boxh(3) ) then
                  zw = zw - sign( box(3),zw )
                  duz  = abs( zw )
              endif
              if( (dux*boct(1)+duy*boct(2)+duz).gt.boct(3) ) then
                  xw = xw - sign( boxh(1),xw )
                  yw = yw - sign( boxh(2),yw )
                  zw = zw - sign( boxh(3),zw )
              endif
#endif
              rrw(jj) = xw**2 + yw**2 + zw**2
  121     continue
      endif
      return
      end

c-----------------------------------------------------------------------
      subroutine bound9(xpij)
#ifdef DPREC
      implicit double precision (a-h,o-z)
#endif
c
#include "box.h"
      dimension xpij(3)

c
      if (ifbox.eq.1) then
          DO 76 M = 1,3
              IF(XPIJ(M).LT.BOXH(M)) GOTO 75
                  XPIJ(M) = XPIJ(M)-BOX(M)
                  GOTO 76
   75         IF(XPIJ(M).GE.-BOXH(M)) GOTO 76
                  XPIJ(M) = XPIJ(M)+BOX(M)
   76     CONTINUE
      else
c
c         -- truncated octahedral boundaries
c
#  ifdef USENINT
          dux = xpij(1)*boxi(1)
          duy = xpij(2)*boxi(2)
          duz = xpij(3)*boxi(3)
#  ifdef DPREC
          dux = dux - dnint( dux )
          duy = duy - dnint( duy )
          duz = duz - dnint( duz )
#  else
          dux = dux - anint( dux )
          duy = duy - anint( duy )
          duz = duz - anint( duz )
#  endif
          corr= 0.5 * aint( fothi* (abs( dux ) +
     +                             abs( duy ) +
     +                             abs( duz )))
          dux = dux - sign( corr,dux )
          duy = duy - sign( corr,duy )
          duz = duz - sign( corr,duz )
          xpij(1) = dux * box(1)
          xpij(2) = duy * box(2)
          xpij(3) = duz * box(3)
#else
          dux = abs(xpij(1))
          duy = abs(xpij(2))
          duz = abs(xpij(3))
          if ( dux.gt.boxh(1) ) then
              xpij(1) = xpij(1) - sign( box(1),xpij(1) )
              dux  = abs( xpij(1) )
          endif
          if ( duy.gt.boxh(2) ) then
              xpij(2) = xpij(2) - sign( box(2),xpij(2) )
              duy  = abs( xpij(2) )
          endif
          if ( duz.gt.boxh(3) ) then
              xpij(3) = xpij(3) - sign( box(3),xpij(3) )
              duz  = abs( xpij(3) )
          endif
          if ( (dux*boct(1)+duy*boct(2)+duz).gt.boct(3) ) then
              xpij(1) = xpij(1) - sign( boxh(1),xpij(1) )
              xpij(2) = xpij(2) - sign( boxh(2),xpij(2) )
              xpij(3) = xpij(3) - sign( boxh(3),xpij(3) )
          endif
#endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine bound10(X, Y, Z)
#ifdef DPREC
      implicit double precision (a-h,o-z)
#endif
c
#include "box.h"

      if (ifbox.eq.1) then
c
c         -- rectilinear boundary
c
          IF(ABS(X).GT.BOXH(1)) X = X-SIGN(BOX(1),X)
          IF(ABS(Y).GT.BOXH(2)) Y = Y-SIGN(BOX(2),Y)
          IF(ABS(Z).GT.BOXH(3)) Z = Z-SIGN(BOX(3),Z)
      else
c
c         -- truncated octahedral boundary
c
#ifdef USENINT
          dux = x*boxi(1)
          duy = y*boxi(2)
          duz = z*boxi(3)
#  ifdef DPREC
          dux = dux - dnint( dux )
          duy = duy - dnint( duy )
          duz = duz - dnint( duz )
#  else
          dux = dux - anint( dux )
          duy = duy - anint( duy )
          duz = duz - anint( duz )
#  endif
          corr= 0.5 * aint( fothi* (abs( dux ) +
     +                                        abs( duy ) +
     +                                        abs( duz )))
          dux = dux - sign( corr,dux )
          duy = duy - sign( corr,duy )
          duz = duz - sign( corr,duz )
          x = dux * box(1)
          y = duy * box(2)
          z = duz * box(3)
#else
          dux = abs(x)
          duy = abs(y)
          duz = abs(z)
          if( dux.gt.boxh(1) ) then
              x = x - sign( box(1),x )
              dux  = abs( x )
          endif
          if( duy.gt.boxh(2) ) then
              y = y - sign( box(2),y )
              duy  = abs( y )
          endif
          if( duz.gt.boxh(3) ) then
              z = z - sign( box(3),z )
              duz  = abs( z )
          endif
          if( (dux*boct(1)+duy*boct(2)+duz).gt.boct(3) ) then
              x = x - sign( boxh(1),x )
              y = y - sign( boxh(2),y )
              z = z - sign( boxh(3),z )
          endif
#endif
      endif
      return
      end
