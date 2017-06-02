c************************************************************************
c                              AMBER                                   **
c                                                                      **
c                  Copyright (c) 1986, 1991, 1995                      **
C             Regents of the University of California                  **
c                       All Rights Reserved.                           ** 
c                                                                      **
c  This software provided pursuant to a license agreement containing   **
c  restrictions on its disclosure, duplication, and use. This software **
c  contains confidential and proprietary information, and may not be   **
c  extracted or distributed, in whole or in part, for any purpose      **
c  whatsoever, without the express written permission of the authors.  **
c  This notice, and the associated author list, must be attached to    **
c  all copies, or extracts, of this software. Any additional           **
c  restrictions set forth in the license agreement also apply to this  **
c  software.                                                           **
c************************************************************************

c
c     -----------------------------------------------------------------
c     All of the particle mesh Ewald code was written and contributed 
c     by Tom Darden from the National Institute of Environmental Health 
c     Sciences division of the NIH.  Originally written with a modified 
c     version of AMBER 3A, the code was updated during the summer of 1994
c     to be compatible with AMBER 4.1.
c     -----------------------------------------------------------------
c
c     4.1 Revisions: tec3 (added comments, modified source to fit it with
c     default Makefile, Compile, MACHINE source handling scheme, added
c     CPP selectable precision, added CPP control of SGI multiprocessing,
c     changed output formats to rid write(6,*), streamlined source,
c     merged routines, etc...)
c
c     4.2 Revisions: tec3 (broke out source into separate files,
c     performed more streamlining, commented source, added Cray T3D
c     parallelization written by Micheal Crowley of the Pittsburgh
c     Supercomputing Center)
c
c     NOTE: Some of the routine names within here are rather long so
c     may choke lame compilers.  This code was developed on SGI
c     computers so should surely work on these; it has also been 
c     limitedly tested on Cray and HP machines.  
c
c     The following C preprocessor code (only visible in the 
c     untransformed source) is a hack to allow single precision 
c     versions of an originally all double precision code.  
c     It is recommended however that users run with double
c     precision... (by default sander is double precision)
c
#ifdef DPREC
#define _REAL_ double precision
#else
#define _REAL_ real
#endif

c-------------------------------------------------------------------
c
c     --- ADJ_MEM_PTR ---
c
c-------------------------------------------------------------------
      subroutine adj_mem_ptr(mem_ptr,assign_ptr,size)
      implicit none
      integer mem_ptr,assign_ptr,size

      assign_ptr = mem_ptr
      mem_ptr = mem_ptr + size 
      return
      end

c-------------------------------------------------------------------
c
c     --- ASSIGN_IND ---
c
c-------------------------------------------------------------------
      subroutine assign_ind(ind0,ind)
      integer ind0,ind
      ind0 = ind
      ind = ind+1
      return
      end

c-------------------------------------------------------------------
c
c     --- CHECK_NEUTRAL ---
c
c-------------------------------------------------------------------
c     ...in general, the Ewald method is only truely applicable
c     under conditions of charge neutrality.  When the system is not
c     net neutral, the direct sum and reciprocal sums are not "beta"
c     independent.  Regardless, the Ewald method can be applied 
c     with the ficticious assumption that there is a *uniform net 
c     neutralizing plasma*.  
c
c     This routine will remove any net charge resulting from
c     conversion of the low precision parm topology charges in the
c     case that the system is supposed to be net neutral.
c
      subroutine check_neutral(charge,natom,ischrgd)
      implicit none
      _REAL_ charge(*)
      integer natom,ischrgd

      integer i
      _REAL_ sum
#include "extra.h"

      sum = 0.d0
      do 100 i = 1,natom
        sum = sum + charge(i)
100   continue

      if (master) write(6, '(/,5x,a,f12.8)')
     .     'Sum of charges from parm topology file = ', 
     .     sum / 18.2223d0
      if ( ischrgd .eq. 1 ) then
        if (master) write(6, '(5x,a)')
     .        'Assuming uniform neutralizing plasma'
        return
      endif
      if (master) write(6, '(5x,a)') 'Forcing neutrality...'
      sum = sum/natom
      do 200 i = 1,natom
        charge(i) = charge(i) - sum
200   continue
      return
      end

c-------------------------------------------------------------------
c
c     --- CUBSPL ---
c
c-------------------------------------------------------------------
c     ...this code is from netlib...
c
      subroutine cubspl ( tau, c, n, ibcbeg, ibcend )
c  from  * a practical guide to splines *  by c. de boor    
c     ************************  input  ***************************
c     n = number of data points. assumed to be .ge. 2.
c     (tau(i), c(1,i), i=1,...,n) = abscissae and ordinates of the
c        data points. tau is assumed to be strictly increasing.
c     ibcbeg, ibcend = boundary condition indicators, and
c     c(2,1), c(2,n) = boundary condition information. specifically,
c        ibcbeg = 0  means no boundary condition at tau(1) is given.
c           in this case, the not-a-knot condition is used, i.e. the
c           jump in the third derivative across tau(2) is forced to
c           zero, thus the first and the second cubic polynomial pieces
c           are made to coincide.)
c        ibcbeg = 1  means that the slope at tau(1) is made to equal
c           c(2,1), supplied by input.
c        ibcbeg = 2  means that the second derivative at tau(1) is
c           made to equal c(2,1), supplied by input.
c        ibcend = 0, 1, or 2 has analogous meaning concerning the
c           boundary condition at tau(n), with the additional infor-
c           mation taken from c(2,n).
c     ***********************  output  **************************
c     c(j,i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
c        of the cubic interpolating spline with interior knots (or
c        joints) tau(2), ..., tau(n-1). precisely, in the interval
c        (tau(i), tau(i+1)), the spline f is given by
c           f(x) = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
c        where h = x - tau(i). the function program *ppvalu* may be
c        used to evaluate f or its derivatives from tau,c, l = n-1,
c        and k=4.
      implicit none
      integer ibcbeg,ibcend,n,   i,j,l,m
      _REAL_ c(4,n),tau(n),   divdf1,divdf3,dtau,g
c****** a tridiagonal linear system for the unknown slopes s(i) of
c  f  at tau(i), i=1,...,n, is generated and then solved by gauss elim-
c  ination, with s(i) ending up in c(2,i), all i.
c     c(3,.) and c(4,.) are used initially for temporary storage.
      l = n-1
compute first differences of tau sequence and store in c(3,.). also,
compute first divided difference of data and store in c(4,.).
      do 10 m=2,n
         c(3,m) = tau(m) - tau(m-1)
   10    c(4,m) = (c(1,m) - c(1,m-1))/c(3,m)
construct first equation from the boundary condition, of the form
c             c(4,1)*s(1) + c(3,1)*s(2) = c(2,1)
      if (ibcbeg-1)                     11,15,16
   11 if (n .gt. 2)                     go to 12
c     no condition at left end and n = 2.
      c(4,1) = 1.d0
      c(3,1) = 1.d0
      c(2,1) = 2.d0*c(4,2)
                                        go to 25
c     not-a-knot condition at left end and n .gt. 2.
   12 c(4,1) = c(3,3)
      c(3,1) = c(3,2) + c(3,3)
      c(2,1) =((c(3,2)+2.d0*c(3,1))*c(4,2)*c(3,3)+
     $          c(3,2)*c(3,2)*c(4,3))/c(3,1)
                                        go to 19
c     slope prescribed at left end.
   15 c(4,1) = 1.d0
      c(3,1) = 0.d0
                                        go to 18
c     second derivative prescribed at left end.
   16 c(4,1) = 2.d0
      c(3,1) = 1.d0
      c(2,1) = 3.d0*c(4,2) - c(3,2)/2.d0*c(2,1)
   18 if(n .eq. 2)                      go to 25
c  if there are interior knots, generate the corresp. equations and car-
c  ry out the forward pass of gauss elimination, after which the m-th
c  equation reads    c(4,m)*s(m) + c(3,m)*s(m+1) = c(2,m).
   19 do 20 m=2,l
         g = -c(3,m+1)/c(4,m-1)
         c(2,m) = g*c(2,m-1) + 3.d0 * 
     $           (c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
   20    c(4,m) = g*c(3,m-1) + 2.d0*(c(3,m) + c(3,m+1))
construct last equation from the second boundary condition, of the form
c           (-g*c(4,n-1))*s(n-1) + c(4,n)*s(n) = c(2,n)
c     if slope is prescribed at right end, one can go directly to back-
c     substitution, since c array happens to be set up just right for it
c     at this point.
      if (ibcend-1)                     21,30,24
   21 if (n .eq. 3 .and. ibcbeg .eq. 0) go to 22
c     not-a-knot and n .ge. 3, and either n.gt.3 or  also not-a-knot at
c     left end point.
      g = c(3,n-1) + c(3,n)
      c(2,n) = ((c(3,n)+2.d0*g)*c(4,n)*c(3,n-1)
     *            + c(3,n)*c(3,n)*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
      g = -g/c(4,n-1)
      c(4,n) = c(3,n-1)
                                        go to 29
c     either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
c     knot at left end point).
   22 c(2,n) = 2.d0*c(4,n)
      c(4,n) = 1.d0
                                        go to 28
c     second derivative prescribed at right endpoint.
   24 c(2,n) = 3.d0*c(4,n) + c(3,n)/2.d0*c(2,n)
      c(4,n) = 2.d0
                                        go to 28
   25 if (ibcend-1)                     26,30,24
   26 if (ibcbeg .gt. 0)                go to 22
c     not-a-knot at right endpoint and at left endpoint and n = 2.
      c(2,n) = c(4,n)
                                        go to 30
   28 g = -1.d0/c(4,n-1)
complete forward pass of gauss elimination.
   29 c(4,n) = g*c(3,n-1) + c(4,n)
      c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
carry out back substitution
   30 j = l 
   40    c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
         j = j - 1
         if (j .gt. 0)                  go to 40
c****** generate cubic coefficients in each interval, i.e., the deriv.s
c  at its left endpoint, from value and slope at its endpoints.
      do 50 i=2,n
         dtau = c(3,i)
         divdf1 = (c(1,i) - c(1,i-1))/dtau
         divdf3 = c(2,i-1) + c(2,i) - 2.d0*divdf1
         c(3,i-1) = 2.d0*(divdf1 - c(2,i-1) - divdf3)/dtau
   50    c(4,i-1) = (divdf3/dtau)*(6.d0/dtau)
                                        return
      end


c
c     --- ERFC_TABLE_MEM ---
c
c     ...patterned after locmem in amber code...
c
      subroutine erfc_table_mem(startreal,endreal,
     $        mxerftab,lerf_arr,ltau,ew_coeff,cutoffnb,erftbdns)

      implicit none
      integer startreal,endreal
      integer mxerftab,lerf_arr,ltau
      _REAL_ ew_coeff,cutoffnb,erftbdns

      integer mem_ptr

c erf table: assume all nonbond distances are less than 1.5 X cutoff
c between nonbond updates; i.e. no excess motion
      mxerftab = int(ew_coeff*erftbdns*cutoffnb*1.5d0)
      write(6, '(a,6X,a,i10)') 
     .     '|','Size of ERFTABLE                 = ', mxerftab

c do real array offsets
      mem_ptr = startreal

c permanent or heap REAL storage
      call adj_mem_ptr(mem_ptr,lerf_arr,4*mxerftab)

c stack REAL storage
      call adj_mem_ptr(mem_ptr,ltau,mxerftab)

      endreal =  mem_ptr

      return
      end


#ifdef HAS_FTN_ERFC
c-------------------------------------------------------------------
c
c     --- ERFCFUN ---
c
c-------------------------------------------------------------------
c     NOTE: this routine is a wrapper for the ERFC() routine
c     which is commonly contained in C libraries, yet for some
c     strange reason less often in Fortran libraries.  It seems
c     necessary to link the C version to this program via a C
c     wrapper routine (e.g. we cannot take on the whole C math 
c     library). Therefore the Fortran version, when available,
c     needs to have a wrapper with the same name, and this is it.
c     If the C route is taken, the C wrapper is placed in the
c     src/Machine/$SYSDIR/sys.a via an erfcfun.c file, and
c     -DHAS_FTN_ERFC is _NOT_ included in MACHINEs MACHINEFLAGS.
c
      subroutine erfcfun(x,y)
      real x,y

      y = erfc(x)
      return
      end
c
      subroutine derfcfun(x,y)
      double precision x,y,derfc
#ifdef CRAYFISH
      y = erfc(x)
#else
#   ifdef T3D
      y = erfc(x)
#   else
      y = derfc(x)
#   endif
#endif
      return
      end
#endif


c-------------------------------------------------------------------
c
c     --- EWALD_BOMB ---
c
c-------------------------------------------------------------------
      subroutine ewald_bomb(routine,string1,string2)
      implicit none
      character*(*) routine,string1,string2

      write(6, '(1x,3a)')
     .     'EWALD BOMB in subroutine ', routine, ' (ewald.f):'
      write(6, '(1x,a)') string1
      write(6, '(1x,a)') string2
      call mexit(6,1)
      end


c-------------------------------------------------------------------
c
c     --- EWALD_MEM ---
c
c-------------------------------------------------------------------
c     ...this routine is patterned after locmem in amber code, 
c     in order to allow seamless integration...
c
      subroutine ewald_mem(maxreal,maxint,maxpr,numatoms,numbad,
     $      startreal,endreal,startint,endint)
c     ARGUMENT LIST:
c     maxreal (INPUT) is the size of your real array X. Used to check 
c       here for overrun.
c     maxint (INPUT), the size of IX, is used to check for overrun
c     maxpr (INPUT), the size of the pairlist array
c       NOTE: currently the nonbond list is not packed into smaller words!
c       This is because the nonbond routines use preimaging which may 
c       generate larger indices in imaged atoms, or in other words 16 bits 
c       may not suffice even if numatoms < 32K!!
c     numatoms (INPUT) is the number of atoms in the system
c     startreal (INPUT) is the starting offset in X of the ewald specific 
c       real memory  e.g. in amber many other things such as coords and 
c       forces are stored in X(1) thru X(startreal-1)
c     endreal (OUTPUT) is computed as the end of that memory

      implicit none
      integer numatoms,numbad,iserr
      integer startreal,endreal,startint,endint,maxint,maxreal,maxpr
#     include "ew_pme_recip.h"
#     include "ew_adj.h"
#     include "ew_erfc_spline.h"
#     include "ew_numtasks.h"
#     include "ew_unitcell.h"
#     include "ew_localnb.h"
#     include "ew_cntrl.h"
      integer startR,startI
#ifdef MPI
#     include "parallel.h"
#endif
#     include "extra.h"

      if (master) write(6, '(/,3x,a,/)') 'EWALD MEMORY USE:'
c
c     --- FIRST GET SOME SIZES
c
c     -- pme sizes
c
      call pmesh_kspace_get_sizes(
     $     nfft1,nfft2,nfft3,numatoms,order,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack,sizscr)
c
c     -- adjust ewald code: mask size
c
      mxadjmsk = numbad
      if (master) write(6, '(a,4x,a,i10)') 
     .     '|','Adjacent nonbond minimum mask    = ', mxadjmsk
c
c     -- local nonbond code sizes
c
      call local_nb_get_sizes(enumtasks,nghb,maxnptrs,numbad,
     $     cutoffnb,reclng,dirlng,nucgrd1,nucgrd2,nucgrd3,nghb1,
     $     nghb2,nghb3,nimgrd1,nimgrd2,nimgrd3,verbose,nucgmax,
     $     nimgmax,numatoms,maximage,mxlstmsk)
c     
c     --- NOW DO POINTER OFFSETS
c
      if (master) write(6, '(3X,a)') 'EWALD LOCMEM POINTER OFFSETS'
      startR = startreal
      startI = startint
c
c     -- reals
c
      call pme_mem(numatoms,startR,endreal,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizscr,
     $     nfft1,nfft2,nfft3,order,
     $     lfftable,lbspmod1,lbspmod2,lbspmod3,
     $     lQarray,lffwork,ltheta1,ltheta2,ltheta3,ldtheta1,
     $     ldtheta2,ldtheta3,lfr1,lfr2,lfr3,lvscr)
      if (master) write(6, '(a,6X,a,i10)')
     .     '|','Real memory needed by PME        = ', endreal-startR
      startR = endreal

      call erfc_table_mem(startR,endreal,
     $     mxerftab,lerf_arr,ltau,ew_coeff,cutoffnb,erftbdns)
      if (master) write(6, '(a,6X,a,i10)')
     .     '|','Real memory needed by ERFTABLE   = ', endreal-startR
      startR = endreal
c
c     -- integers
c
      endint = startI
      call adj_mem_ptr(endint,imask1,mxadjmsk)
      call adj_mem_ptr(endint,imask2,mxadjmsk)

      if (master) write(6, '(a,6X,a,i10)')
     .     '|','Integer memory needed by ADJ     = ', endint-startI
      startI = endint
c
c     --- THE LOCAL NONBOND NEEDS TO BE MAPPED LAST DUE TO NB LIST
c
      call loc_nb_mem(numatoms,
     $     startR,startI,endreal,endint,cutoffnb,mxlstmsk,
     $     enumtasks,maxnptrs,nucgmax,maximage,mempage,nimgmax,
     $     ltr_img,limgcrds,lfrction,lsavfrac,ldelcrd,ldfrac,
     $     lpforce,lscalars,
     $     inumatg,iindatg,iatmcell,iindoff,iimagptr,inghbptr,
     $     ibckptr,iucptr,inlogrid,inhigrid,inumimg,
     $     inummask,imaskptr,imask,iatmlist,ilist,iscratch,
     $     iiwa,iiwh,inumlist,ifailtsk,iimgind,inumvdw,inumhbnd,
     $     iipack)
      if (master) write(6, '(a,6X,a,i10)')
     .     '|','Integer memory used by local nonb= ', endint-startI
      if (master) write(6, '(a,6X,a,i10)')
     .     '|','Real memory used by local nonb   = ', endreal-startR
#ifdef MEM_ALLOC
c maybe calc to alloc max nb size here?
#else
      iserr = 0
      if ( endint .gt. maxint )then
         if (master) write(6, '(/5X,a,i12,a,i12,a)') 
     .      'INTEGER ARRAY OVERFLOW!  MAXINT IS ',maxint,
     .      ' NEED AT LEAST ', endint, ' INTEGER WORDS'
         iserr = 1
      endif
      if ( endreal .gt. maxreal )then
         if (master) write(6, '(/5X,a,i12,a,i12,a)') 
     .        'REAL ARRAY OVERFLOW!  MAXREAL IS ',maxreal,
     .      ' NEED AT LEAST ', endreal, ' REAL WORDS'
         iserr = 1
      endif
      if (iserr.ne.0) call mexit(6,1)
      maxnblst = maxpr
      if (master) write(6, '(/,a,4X,a,i10)') 
     .     '|','MAX NONBOND PAIRS = ', maxnblst
#endif
      return
      end


c-------------------------------------------------------------------
c
c     --- EW_STARTUP ---
c
c-------------------------------------------------------------------
c     ...called after load_ewald_info() and locmem() to fill in some 
c     arrays and perform other initial chores...
c
      subroutine ew_startup(numatoms,iblo,inb,X,IX,charge)
      implicit none
#     include "ew_pme_recip.h"
#     include "ew_erfc_spline.h"
#     include "ew_adj.h"
#     include "ew_localnb.h"
#     include "extra.h"
      integer numatoms,numbad,iblo(*),inb(*)
      _REAL_ X(*),charge(*)
      integer IX(*)
CCC  April 4,2014, modified by Jianfa Chen
cc      _REAL_ tim1,tim2
      REAL(KIND=4) :: tim1,tim2
cc   end modified by Jianfa
#ifdef MPI
# include "parallel.h"
#  include "ew_parallel.h"
#endif

      call second(tim1)
      call check_neutral(charge,numatoms,ischrgd)
      call init_profiles()
      call fill_erf_table(ERFTBDNS,mxerftab,
     $       X(lerf_arr),X(ltau))
      call load_adj_mask(iblo,inb,numatoms,
     $     mxadjmsk,IX(imask1),IX(imask2),numadjst)
      call load_list_mask(iblo,inb,numatoms,
     $     IX(inummask),mxlstmsk,IX(imaskptr),IX(imask))
#ifdef MPI
      if(i_do_recip)then
         call MPI_COMM_SIZE(recip_comm,numtasks,ierr)
         call MPI_COMM_RANK(recip_comm,mytaskid,ierr)
#endif
      call pmesh_kspace_setup(
     $    X(lbspmod1),X(lbspmod2),X(lbspmod3),X(lfftable),X(lffwork),
     $    nfft1,nfft2,nfft3,order,sizfftab,sizffwrk)
#ifdef MPI
         call MPI_COMM_SIZE(world_comm,numtasks,ierr)
         call MPI_COMM_RANK(world_comm,mytaskid,ierr)
      endif
#endif
      call second(tim2)
      if (master) then
        write(6, '(a,4x,a,f12.8)') 
     .     '|','Total Ewald setup time = ', tim2-tim1
        write(6, 1)
 1      format(t2,78('-'),/)
      endif
      return
      end

c-------------------------------------------------------------------
c
c     --- FILL_ERF_TABLE ---
c
c-------------------------------------------------------------------
      subroutine fill_erf_table(erftbdns,mxerftab,
     $       erf_arr,tau)

      implicit none
      integer mxerftab
      _REAL_ erftbdns,erf_arr(4,*),tau(*)

      _REAL_ del,x,pi,fac,erf
      integer i,ibcbeg,ibcend

      del = 1.d0 / erftbdns
      pi = 3.14159265358979323846d0
      fac = 2.d0 / sqrt(pi)

      erf_arr(2,1) = -fac
      x = (mxerftab-1)*del
      erf_arr(2,mxerftab) = -fac * exp(-x*x)
      do 100 i = 1,mxerftab
        x = del*(i-1)
#ifdef DPREC
        call derfcfun(x,erf)
#else
        call erfcfun(x,erf)
#endif
        tau(i) = x
        erf_arr(1,i) = erf
100   continue
      ibcbeg = 1
      ibcend = 1
      call cubspl ( tau, erf_arr, mxerftab, ibcbeg, ibcend )
      return
      end



c-------------------------------------------------------------------
c
c     --- FIND_EWALDCOF ---
c
c-------------------------------------------------------------------
      subroutine find_ewaldcof(cutoff,dtol,ewaldcof)
      implicit none
      _REAL_ cutoff,dtol,ewaldcof
      integer i,n
      _REAL_ term,x,xlo,xhi,y,erfc

c first get direct sum tolerance. How big must ewaldcof be to get
c terms outside the cutoff below tol

      x = 0.5d0
      i = 0
10    continue
          x = 2.d0 * x
          i = i + 1
          y = x * cutoff
#ifdef DPREC
          call derfcfun(y,erfc)
#else
          call erfcfun(y,erfc)
#endif
          term = erfc/cutoff
      if ( term .ge. dtol) goto 10
c binary search tolerance is 2 to the -50th
      n = i + 50
      xlo = 0.d0
      xhi = x
      do 20 i = 1,n
        x = (xlo+xhi)/2
        y = x * cutoff
#ifdef DPREC
        call derfcfun(y,erfc)
#else
        call erfcfun(y,erfc)
#endif
        term = erfc/cutoff
        if ( term .ge. dtol )then
           xlo = x
        else 
           xhi = x
        endif
20    continue
      ewaldcof = x

      return
      end


c-------------------------------------------------------------------
c
c     --- INIT_PROFILES ---
c
c-------------------------------------------------------------------
      subroutine init_profiles()
      implicit none
#     include "ew_time.h"
      integer i,ind
cc   April 4,2014,modified by Jianfa Chen
c      _REAL_ tim1
      REAL(KIND=4):: tim1
cc end modification by Jianfa
      call second(tim1)
      time1 = tim1
      do 10 i = 1,10
       listtime(i) = 0.d0
10    continue
      ind = 1
      call assign_ind(IMAP,ind)
      call assign_ind(ISETGRD,ind)
      call assign_ind(IGRDUC,ind)
      call assign_ind(IGRDIM,ind)
      call assign_ind(IBLDLST,ind)
      call assign_ind(ILSTTIM,ind)
      liststr(IMAP) =    'map, save_t     = '
      liststr(ISETGRD) = 'grid setup_t    = '
      liststr(IGRDUC) =  'grid ucell_t    = '
      liststr(IGRDIM) =  'grid image_t    = '
      liststr(IBLDLST) = 'build list_t    = '
      liststr(ILSTTIM) = 'total list_t    = '
      do 20 i = 1,20
        forctime(i) = 0.d0
20    continue
      ind = 1
      call assign_ind(ICLEAR,ind)
      call assign_ind(IADJMAP,ind)
      call assign_ind(ISELF,ind)
      call assign_ind(IBSPL,ind)
      call assign_ind(IFILLQ,ind)
      call assign_ind(ISCSUM,ind)
      call assign_ind(IGRSUM,ind)
      call assign_ind(IFFT,ind)
      call assign_ind(IDIRSUM,ind)
      call assign_ind(IADJST,ind)
      call assign_ind(IACCFRC,ind)
      call assign_ind(INBVIR,ind)
      call assign_ind(IFRCTIM,ind)
c
      forcestr(ICLEAR) =  'zero ene,force  = '
      forcestr(IADJMAP) = 'map,adjust      = '
      forcestr(ISELF) =   'self energy     = '
      forcestr(IBSPL) =   '1-dim b-spline  = '
      forcestr(IFILLQ) =  'grid charges    = '
      forcestr(ISCSUM) =  'scalar sum      = '
      forcestr(IGRSUM) =  'grad   sum      = '
      forcestr(IFFT) =    'FFT             = '
      forcestr(IDIRSUM) = 'direct force    = '
      forcestr(IADJST) =  'adjust masked   = '
      forcestr(IACCFRC) = 'accum  force    = '
      forcestr(INBVIR) =  'finish virial   = '
      forcestr(IFRCTIM) = 'total  force    = '

      return
      end




c-------------------------------------------------------------------
c
c     --- LOAD_ADJ_MASK ---
c
c-------------------------------------------------------------------
      subroutine load_adj_mask(iblo,inb,numatoms,
     $     maxmask,mask1,mask2,numadjst)
      implicit none
      integer iblo(*),inb(*),numatoms
      integer mask1(*),mask2(*),numadjst,maxmask

      integer ji,i,j,nx,k
#include "extra.h"
#ifdef MPI
# include "parallel.h"
#endif

c PASS 1 check sizes
      ji = 0
      numadjst = 0
      do 50 i = 1,numatoms-1
       nx = iblo(i)
       do 25 j = 1,nx
        k = inb(ji+j)
c use extra arrays mask1,mask2 to speed the nb_adjust routine
        if ( k .gt. i )then
            numadjst = numadjst + 1
        endif
25     continue
       ji = ji + nx
50    continue
      if (master) then
        write(6, '(5x,a,i10)') 'Total number of mask terms = ', numadjst
      endif
      if ( numadjst .gt. MAXMASK)
     $   call ewald_bomb('load_mask','MAXMASK not big enough!!',' ')

c PASS 2 fill mask array
      ji = 0
      numadjst = 0
      do 150 i = 1,numatoms-1
       nx = iblo(i)
       do 125 j = 1,nx
        k = inb(ji+j)
c use extra arrays mask1,mask2 to speed the nb_adjust routine
        if ( k .gt. i )then
            numadjst = numadjst + 1
            mask1(numadjst) = i
            mask2(numadjst) = k
        endif
125    continue
       ji = ji + nx
150   continue

      return
      end

c-------------------------------------------------------------------
c
c     --- LOAD_EWALD_INFO ---
c
c-------------------------------------------------------------------
c     ...routine which sets up some sizes and parameters necessary
c     for locmem() and must be called prior to calling locmem().
c     Currently this is called inside mdread().
c
      subroutine load_ewald_info(nf)
      implicit none
      integer nf
c
#     include "ew_cntrl.h"
#     include "ew_unitcell.h"
#     include "ew_pme_recip.h"
#     include "box.h"

      write(6, 9001) 
#ifdef FUJFFT
      write(6, '(5X,A)')
     .     'Using FUJITSU specific Fast Fourier Transform'
#endif
#ifdef PUBFFT
      write(6, '(A)')
     .     '|     Using PUBFFT Fast Fourier Transform'
#endif
#ifdef CRAYFFT
      write(6, '(A)')
     .     '|   Using the CRAY specific (CFFT3D) Fast Fourier Transform'
#endif
#ifdef SGIFFT
      write(6, '(A)')
     .     '|    Using the SGI specific (ZFFT3D) Fast Fourier Transform'
#endif
#ifdef MPI
      write(6, '(A)')
     .     '|    Using the T3D specific (FFT3D0) Fast Fourier Transform'
#endif
      call start_numtasks()
      cutoffnb = cut
      call read_ewald(nf,a,b,c,alpha,beta,gamma,
     $         nfft1,nfft2,nfft3,order,ischrgd,verbose,
     $         checkacc,dsum_tol)
      box(1) = a
      box(2) = b
      box(3) = c
c
      write(6, 9002) a, b, c
      write(6, 9003) alpha, beta, gamma
      write(6, 9004) nfft1, nfft2, nfft3
      write(6, 9006) cutoffnb, dsum_tol
      write(6, 9005) order
c
c     tec3: disable checkacc (exact) temporarilly until exact
c     ewald code is dropped in...
c
      if (checkacc .ne. 0) then
         write(6, '(5X,A)')
     .        'The exact Ewald option is current not implemented!'
         call mexit(6,1)
      endif
c
c     get some related quantities
      call get_ucell(a,b,c,alpha,beta,gamma,
     $     ucell,recip,dirlng,reclng,sphere,volume)
      call find_ewaldcof(cutoffnb,dsum_tol,ew_coeff)
      write(6, 9007) ew_coeff

 9001 format (/,3X,'EWALD SPECIFIC INPUT:',/)
 9002 format (5X,'Box X =',F9.3,3X,'Box Y =',F9.3,3X,'Box Z =',F9.3)
 9003 format (5X,'Alpha =',F9.3,3X,'Beta  =',F9.3,3X,'Gamma =',F9.3)
 9004 format (5X,'NFFT1 =',I5  ,7X,'NFFT2 =',I5  ,7X,'NFFT3 =',I5)
 9005 format (5X,'Interpolation order =',I5)
 9006 format (5X,'Cutoff=',F9.3,3X,'Tol   =',e9.3)
 9007 format (5X,'Ewald Coefficient =',F9.5,/)
      return
      end


c-------------------------------------------------------------------
c
c     --- LOAD_LIST_MASK ---
c
c-------------------------------------------------------------------
      subroutine load_list_mask(iblo,inb,numatoms,
     $     nummask,maxmask,maskptr,mask)
      implicit none
      integer iblo(*),inb(*),numatoms
      integer nummask(*),maxmask,maskptr(*),mask(*)

c double the mask to deal with our list generator
      integer ji,i,j,m,nx,off,k,tot
#include "extra.h"
#ifdef MPI
# include "parallel.h"
#endif
c PASS 1 get pointers, check size
      ji = 0
      do 10 i = 1,numatoms
        nummask(i) = 0
10    continue
      do 50 i = 1,numatoms - 1
       nx = iblo(i)
       do 25 j = 1,nx
        k = inb(ji+j)
        if ( k .gt. 0 )then
          nummask(k) = nummask(k) + 1
          nummask(i) = nummask(i) + 1
        endif
25     continue
       ji = ji + nx
50    continue
      tot = 0
      do 75 i = 1,numatoms
        tot = tot + nummask(i)
75    continue
      if (master) write(6, '(5x,a,i10)') 
     .     'Total number of mask terms = ', tot
      if ( tot .gt. MAXMASK)
     $   call ewald_bomb('load_mask','MAXMASK not big enough!!',' ')
      off = 0
      do 100 i = 1,numatoms
       maskptr(i) = off
       off = off + nummask(i)
100   continue

c PASS 2 fill mask array
      ji = 0
      do 110 i = 1,numatoms
        nummask(i) = 0
110   continue
      do 150 i = 1,numatoms-1
       nx = iblo(i)
       do 125 j = 1,nx
        k = inb(ji+j)
        if ( k .gt. 0 )then
          nummask(k) = nummask(k) + 1
          m = maskptr(k) + nummask(k)
          mask(m) = i
          nummask(i) = nummask(i) + 1
          m = maskptr(i) + nummask(i)
          mask(m) = k
        endif
125    continue
       ji = ji + nx
150   continue

      return
      end


c-------------------------------------------------------------------
c
c     --- LOC_NB_MEM ---
c
c-------------------------------------------------------------------
c     ...patterned after locmem in amber code...
c
      subroutine loc_nb_mem(numatoms,
     $      startreal,startint,endreal,endint,cutoffnb,maxmask,
     $      enumtasks,maxnptrs,nucgmax,maximage,mempage,nimgmax,
     $      ltr_img,limgcrds,lfrction,lsavfrac,ldelcrd,ldfrac,
     $      lpforce,lscalars,
     $      inumatg,iindatg,iatmcell,iindoff,iimagptr,inghbptr,
     $      ibckptr,iucptr,inlogrid,inhigrid,inumimg,
     $      inummask,imaskptr,imask,iatmlist,ilist,iscratch,
     $      iiwa,iiwh,inumlist,ifailtsk,iimgind,inumvdw,inumhbnd,
     $      iipack)
      implicit none
      _REAL_ cutoffnb
      integer numatoms,startreal,endreal,startint,endint,maxmask,
     $      enumtasks,maxnptrs,nucgmax,maximage,mempage,nimgmax
C REAL POINTERS
      integer ltr_img,limgcrds,lfrction,lsavfrac,ldelcrd,ldfrac,
     $      lpforce,lscalars
C INTEGER POINTERS
      integer inumatg,iindatg,iatmcell,iindoff,iimagptr,inghbptr,
     $      ibckptr,iucptr,inlogrid,inhigrid,inumimg,
     $      inummask,imaskptr,imask,iatmlist,ilist,iscratch,
     $      iiwa,iiwh,inumlist,ifailtsk,iimgind,inumvdw,inumhbnd,
     $      iipack

#ifdef MPI
# include "ew_parallel.h"
# include "parallel.h"
#endif

      integer mem_ptr

#ifdef MPI
      enumtasks = 1
#endif

c do real array offsets
      mem_ptr = startreal

      call adj_mem_ptr(mem_ptr,ltr_img,3*nimgmax)
      call adj_mem_ptr(mem_ptr,limgcrds,3*maximage)
      call adj_mem_ptr(mem_ptr,lfrction,3*numatoms)
      call adj_mem_ptr(mem_ptr,lsavfrac,3*numatoms)

c scratch from here on

#ifndef MPI
      call adj_mem_ptr(mem_ptr,ldelcrd,3*numatoms*enumtasks)
#endif
      call adj_mem_ptr(mem_ptr,ldfrac,3*numatoms)
      call adj_mem_ptr(mem_ptr,lpforce,3*numatoms*enumtasks)
c general multiprocessing scalar scratch
c used for direct energies
      call adj_mem_ptr(mem_ptr,lscalars,mempage*enumtasks)

      endreal =  mem_ptr
c do integer array  offsets
      mem_ptr = startint

      call adj_mem_ptr(mem_ptr,inumatg,nucgmax)
      call adj_mem_ptr(mem_ptr,iindatg,numatoms)
      call adj_mem_ptr(mem_ptr,iatmcell,numatoms)
      call adj_mem_ptr(mem_ptr,iindoff,nucgmax)
      call adj_mem_ptr(mem_ptr,iimagptr,nucgmax)
      call adj_mem_ptr(mem_ptr,inghbptr,maxnptrs*nucgmax)
      call adj_mem_ptr(mem_ptr,ibckptr,maximage)
      call adj_mem_ptr(mem_ptr,iucptr,nimgmax)
      call adj_mem_ptr(mem_ptr,inlogrid,nimgmax)
      call adj_mem_ptr(mem_ptr,inhigrid,nimgmax)
      call adj_mem_ptr(mem_ptr,inumimg,nimgmax)
      call adj_mem_ptr(mem_ptr,inummask,numatoms)
      call adj_mem_ptr(mem_ptr,imaskptr,numatoms)
      call adj_mem_ptr(mem_ptr,imask,maxmask)

#ifndef NEW_NUMAT
#  define NEW_NUMAT numatoms
#endif
      call adj_mem_ptr(mem_ptr,iatmlist,enumtasks*NEW_NUMAT)
      call adj_mem_ptr(mem_ptr,ilist,enumtasks*NEW_NUMAT)
      call adj_mem_ptr(mem_ptr,iscratch,enumtasks*numatoms)
      call adj_mem_ptr(mem_ptr,iiwa,enumtasks*NEW_NUMAT)
      call adj_mem_ptr(mem_ptr,iiwh,enumtasks*NEW_NUMAT)
#ifdef MPI
      enumtasks = numtasks
#endif
      call adj_mem_ptr(mem_ptr,inumlist,enumtasks)
      call adj_mem_ptr(mem_ptr,ifailtsk,enumtasks)
      call adj_mem_ptr(mem_ptr,iimgind,maximage)

c these are always last
      call adj_mem_ptr(mem_ptr,inumvdw,numatoms)
      call adj_mem_ptr(mem_ptr,inumhbnd,numatoms)
      iipack = mem_ptr

      endint = mem_ptr

      return
      end

c-------------------------------------------------------------------
c
c     --- LOCAL_NB_GET_SIZES ---
c
c-------------------------------------------------------------------
      subroutine local_nb_get_sizes(enumtasks,nghb,maxnptrs,numbad,
     $      cutoffnb,reclng,dirlng,nucgrd1,nucgrd2,nucgrd3,nghb1,
     $      nghb2,nghb3,nimgrd1,nimgrd2,nimgrd3,verbose,nucgmax,
     $      nimgmax,numatoms,maximage,maxmask)
      implicit none
      integer enumtasks,nghb,maxnptrs,numbad,nucgrd1,nucgrd2,nucgrd3,
     $      nghb1,nghb2,nghb3,nimgrd1,nimgrd2,nimgrd3,verbose,nucgmax,
     $      nimgmax,numatoms,maximage,maxmask
      _REAL_ cutoffnb,reclng(3),dirlng(3)

      _REAL_ density
      integer nucgmin
#ifdef MPI
# include "parallel.h"
#endif
#include "extra.h"

      maxnptrs = ((2*nghb+1)*(2*nghb+1)+1)/2
      if (master) write(6, '(a,4X,a,i10)') 
     .     '|','Max number of pointers           = ', maxnptrs
      maxmask = 2*numbad
      if (master) write(6, '(a,4X,a,i10)')
     .     '|','List build maxmask               = ', maxmask
      call setup_grids(nghb,cutoffnb,reclng,dirlng,
     $    nucgrd1,nucgrd2,nucgrd3,nghb1,nghb2,nghb3,nimgrd1,
     $    nimgrd2,nimgrd3,verbose)

c set upper bounds for ucell grid dimensions to allow for volume fluctuations
      nucgmax = int(nucgrd1*nucgrd2*nucgrd3 * 1.33d0)
      nucgmin = int(nucgrd1*nucgrd2*nucgrd3 / 1.33d0)
      nimgmax = int(nimgrd1*nimgrd2*nimgrd3 * 1.33d0)
c density is average number of atoms per subcell
      density = numatoms
      density = density / (nucgmin)
c the maximum image is density times nimgrd1*nimgrd2*nimgrd3
c times a fudge factor of 1.1
      maximage = int(density * nimgrd1*nimgrd2*nimgrd3 * 1.1d0)
      if (master) write(6, '(a,4X,a,i10,/)')
     .     '|','Maximage                         = ', maximage
      return
      end

c-------------------------------------------------------------------
c
c     --- PME_MEM ---
c
c-------------------------------------------------------------------
c     ...patterned after locmem in amber code...
c
      subroutine pme_mem(numatoms,startreal,endreal,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizscr,
     $     nfft1,nfft2,nfft3,order,
     $     lfftable,lbspmod1,lbspmod2,lbspmod3,
     $     lQarray,lffwork,ltheta1,ltheta2,ltheta3,ldtheta1,
     $     ldtheta2,ldtheta3,lfr1,lfr2,lfr3,lvscr)

      implicit none
      integer startreal,endreal
      integer nfft1,nfft2,nfft3,numatoms,order,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack,sizscr,
     $     lfftable,lbspmod1,lbspmod2,lbspmod3,
     $     lQarray,lffwork,ltheta1,ltheta2,ltheta3,ldtheta1,
     $     ldtheta2,ldtheta3,lfr1,lfr2,lfr3,lvscr

      integer mem_ptr
c     
      mem_ptr = startreal

c permanent or heap REAL storage
      call adj_mem_ptr(mem_ptr,lfftable,sizfftab)
      call adj_mem_ptr(mem_ptr,lbspmod1,nfft1)
      call adj_mem_ptr(mem_ptr,lbspmod2,nfft2)
      call adj_mem_ptr(mem_ptr,lbspmod3,nfft3)

c stack REAL storage
      call adj_mem_ptr(mem_ptr,lQarray,siz_Q)
      call adj_mem_ptr(mem_ptr,lffwork,sizffwrk)
      call adj_mem_ptr(mem_ptr,ltheta1,siztheta)
      call adj_mem_ptr(mem_ptr,ltheta2,siztheta)
      call adj_mem_ptr(mem_ptr,ltheta3,siztheta)
      call adj_mem_ptr(mem_ptr,ldtheta1,siztheta)
      call adj_mem_ptr(mem_ptr,ldtheta2,siztheta)
      call adj_mem_ptr(mem_ptr,ldtheta3,siztheta)
      call adj_mem_ptr(mem_ptr,lfr1,numatoms)
      call adj_mem_ptr(mem_ptr,lfr2,numatoms)
      call adj_mem_ptr(mem_ptr,lfr3,numatoms)
      call adj_mem_ptr(mem_ptr,lvscr,sizscr)

      endreal =  mem_ptr

      return
      end


c-------------------------------------------------------------------
c     --- PMESH_KSPACE_GET_SIZES ---
c-------------------------------------------------------------------
c     ...this routine computes the parameters needed for heap or 
c     stack allocation specified below...
c
      subroutine pmesh_kspace_get_sizes(
     $     nfft1,nfft2,nfft3,numatoms,order,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack,sizscr)
      implicit none
      integer nfft1,nfft2,nfft3,numatoms,order,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack,sizscr

c INPUT  
c      nfft1,nfft2,nfft3,numatoms,order
c      nfft1,nfft2,nfft3 are the dimensions of the charge grid array
c      numatoms is number of atoms
c      order is the order of B-spline interpolation

c OUTPUT
c      sizfftab,sizffwrk,siztheta,siz_Q
c      sizfftab is permanent 3d fft table storage
c      sizffwrk is temporary 3d fft work storage
c      siztheta is size of arrays theta1-3 dtheta1-3
c      sizscr is size of temporary arrays needed by vector code
c      sizheap is total size of permanent storage
c      sizstack is total size of temporary storage

      integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork
#include "extra.h"
#ifdef MPI
# include "parallel.h"
# include "ew_parallel.h"
#endif

#ifdef VECTOR
      sizscr = 4*numatoms
#else
      sizscr = 1
#endif
      call get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,
     $       sizfftab,sizffwrk)
      siztheta = numatoms*order
#ifndef MPI
      siz_Q = 2*nfftdim1*nfftdim2*nfftdim3
#else
      siz_Q = max( ntxyslab*(nxyslab(0)), ntxzslab*nxzslab(0))
#endif
      sizheap = nfft1+nfft2+nfft3+sizfftab
      sizstack = siz_Q+6*siztheta+sizffwrk+3*numatoms+sizscr
      if (master) then
        write(6, '(a,4X,a,i10)')
     .     '|','Total heap storage needed        = ', sizheap
        write(6, '(a,4X,a,i10)')
     .     '|','Total stack storage needed       = ', sizstack
      endif
      return
      end


c-------------------------------------------------------------------
c     --- PMESH_KSPACE_SETUP ---
c-------------------------------------------------------------------
      subroutine pmesh_kspace_setup(
     $    bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork,
     $    nfft1,nfft2,nfft3,order,sizfftab,sizffwrk)
      implicit none

c  see DO_PMESH_KSPACE for explanation of arguments

      integer nfft1,nfft2,nfft3,order,sizfftab,sizffwrk
      _REAL_ bsp_mod1(*),bsp_mod2(*),bsp_mod3(*)
      _REAL_ fftable(sizfftab),ffwork(sizffwrk)
   
      _REAL_ dummy
      integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw

      call get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
      call load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,
     $   nfft1,nfft2,nfft3,order)
      call fft_setup(dummy,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
      return
      end


c-------------------------------------------------------------------
c     --- PROF_EWTIME ---
c-------------------------------------------------------------------
      subroutine prof_ewtime(timtot)
      implicit none
#     include "ew_time.h"
      _REAL_ timtot

      _REAL_ p
      integer i

      write(6, '(/,2x,a)') ' Nonbond Ewald pairlist time components:'
      do 10 i = 1,ILSTTIM
         p = listtime(i) * 100.d0 / timtot
         write(6,15)liststr(i),listtime(i),p
         listtime(ILSTTIM) = listtime(ILSTTIM) + listtime(i)
10    continue
      write(6, '(/,2x,a)') ' Nonbond Ewald force time components:'
      do 20 i = 1,IFRCTIM
         p = forctime(i) * 100.d0 / timtot
c                    | in column 1 is for test/glsfdif
         write(6,15)forcestr(i),forctime(i),p
         forctime(IFRCTIM) = forctime(IFRCTIM) + forctime(i)
20    continue
15    format('|',5x,a,f12.2,f6.2)
      return
      end

c-------------------------------------------------------------------
c     --- READ_EWALD ---
c-------------------------------------------------------------------
      subroutine read_ewald(nf,a,b,c,alpha,beta,gamma,
     $         nfft1,nfft2,nfft3,order,ischrgd,verbose,
     $         checkacc,dsum_tol)
      implicit none
      integer nf
      _REAL_ a,b,c,alpha,beta,gamma,dsum_tol
      integer nfft1,nfft2,nfft3,order,
     $        ischrgd,verbose,checkacc

c read unit cell parameters
      read(nf,*, end=10)a,b,c,alpha,beta,gamma
c read the grid dimensions
      read(nf,*, end=11)nfft1,nfft2,nfft3,order,
     $          ischrgd,verbose,checkacc
c read the per-term error in the tail of the direct sum
      read(nf,*, end=12)dsum_tol
      return
   10 write(6,*) ' ** EOF reading Ewald unit cell parameters'
      call mexit(6,1)
   11 write(6,*) ' ** EOF reading Ewald grid dimensions'
      call mexit(6,1)
   12 write(6,*) ' ** EOF reading Ewald per-term error'
      call mexit(6,1)
      end




c-------------------------------------------------------------------
c     --- SET_TIMER ---
c-------------------------------------------------------------------
      subroutine set_timer()
      implicit none
#     include "ew_time.h"

cc   April 4,2014,modified by Jianfa Chen
cc      _REAL_ tim
      REAL(KIND=4)::tim
cc end modified by Jianfa
      call second(tim)
      time1 = tim
      return
      end







c-------------------------------------------------------------------
c     --- SETUP_GRIDS ---
c-------------------------------------------------------------------
c     ...this routine checks on the necessary resources for the unit 
c     cell and image cell grids used for short range particle pair 
c     calculations.  It is assumed that unit cell setup has already 
c     occurred.  This routine then checks the short range cutoff. 
c     The unit cell will be split into NUCGRD1 x NUCGRD2 x NUCGRD3 
c     geometrically similar subcells of size dirlng(1)/NUCGRD1  by  
c     dirlng(2)/NUCGRD2  by  dirlng(3)/NUCGRD3. 
c     The short range interactions will involve pairs in the subcell 
c     neighborhood of  +- NGHB1  by  +- NGHB2  by  +- NGHB3  subcells., 
c     about any given subcell.
c     The distances between  parallel faces of the unit cell are ,
c     respectively reclng(1), reclng(2) and reclng(3). 
c     Thus these subcell neighborhoods are guaranteed to contain all points 
c     within the minimum of (NGHB1/NUCGRD1)*reclng(1),
c     (NGHB2/NUCGRD2)*reclng(2),and (NGHB3/NUCGRD3)*1.0/reclng(3). 
c     This minimum is taken to be the short range cutoff...
c
      subroutine setup_grids(nghb,cutoffnb,reclng,dirlng,
     $     nucgrd1,nucgrd2,nucgrd3,nghb1,nghb2,nghb3,nimgrd1,
     $     nimgrd2,nimgrd3,verbose)
      implicit none
#include "extra.h"
#ifdef MPI
# include "ew_parallel.h"
# include "parallel.h"
#endif
      _REAL_ cutoffnb,reclng(3),dirlng(3)
      integer nghb,nucgrd1,nucgrd2,nucgrd3,nghb1,nghb2,nghb3,
     $     nimgrd1,nimgrd2,nimgrd3,verbose

      _REAL_ dc1,dc2,dc3,cut

      nghb1 = 3
      nghb2 = 3
      nghb3 = 3
      dc1 = cutoffnb / nghb1
      dc2 = cutoffnb / nghb2
      dc3 = cutoffnb / nghb3
      nucgrd1 = int(reclng(1) / dc1)
      nucgrd2 = int(reclng(2) / dc2)
      nucgrd3 = int(reclng(3) / dc3)
      nimgrd1 = nucgrd1 + 2 * nghb1
      nimgrd2 = nucgrd2 + 2 * nghb2
      nimgrd3 = nucgrd3 + 2 * nghb3

c check the short range cutoff
      dc1 = reclng(1)/nucgrd1
      dc2 = reclng(2)/nucgrd2
      dc3 = reclng(3)/nucgrd3
      cut = nghb1*dc1
      if ( nghb2*dc2 .lt. cut ) cut = nghb2*dc2
      if ( nghb3*dc3 .lt. cut ) cut = nghb3*dc3
      if (master .and. verbose .ge. 1) then
          write(6, '(5X,a,/,5X,i9,1X,i9,1X,i9)')
     .        'Number of grids per unit cell in each dimension:',
     .         nucgrd1, nucgrd2, nucgrd3
          write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)')
     .       'Unit cell edge lengths in each dimension:',
     .       dirlng(1), dirlng(2), dirlng(3)
          write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)')
     .       'Distance between parallel faces of unit cell:',
     .       reclng(1), reclng(2), reclng(3)
          write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)')
     .      'Distance between faces of short range grid subcells:',
     .       dc1, dc2, dc3
          write(6, '(5X,a,F9.3)')
     .       'Resulting cutoff from subcell neighborhoods is ', cut
      endif
      if ( cut .lt. cutoffnb )then
        call ewald_bomb('setup_grids',
     $    'Resulting cutoff is too small for your lower limit',' ')
      endif
      return
      end




c-------------------------------------------------------------------
c     --- START_NUMTASKS ---
      subroutine start_numtasks()
      implicit none
#     include "ew_numtasks.h"
#ifdef MPI
# include "parallel.h"
#endif
#ifdef SGI_MP
      integer mp_numthreads
      enumtasks = mp_numthreads()
#else
# ifdef MPI
#   ifdef PUBFFT
      enumtasks = 1
#else
      enumtasks = numtasks
#   endif
# else
      enumtasks = 1
# endif
#endif
      return
      end


c-------------------------------------------------------------------
c     --- UPDATE_TIME ---
c-------------------------------------------------------------------
      subroutine update_time(time)
      implicit none
      _REAL_ time
#     include "ew_time.h"


cc   April 4,2014,modified by Jianfa Chen
c      _REAL_ tim
      REAL(KIND=4):: tim
c     end modification by Jianfa
      call second(tim)
      time2 = tim
      time = time + time2-time1
      time1 = time2
      return
      end


c-------------------------------------------------------------------
c     --- ZERO_ARRAY ---
c-------------------------------------------------------------------
      subroutine zero_array(array,num)
      implicit none
      _REAL_ array(*)
      integer num

      integer i
#ifdef SGI_MP
C$DOACROSS LOCAL(i),SHARE(array,num)
#endif
      do 10 i = 1,num
        array(i) = 0.d0
10    continue
      return
      end
#ifdef MPI
c-------------------------------------------------------------------
c     --- STARTUP_GROUPS ---
c-------------------------------------------------------------------
c     --- STARTUP_GROUPS_mpi ---
      subroutine startup_groups(err)
      implicit none
      integer i,j,err
# include "mpif.h"
c # include "ew_numtasks.h"
# include "parallel.h"
# include "ew_parallel.h"

      call mpi_bcast(num_recip,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(NUM_RECIP.gt.numtasks)then
         print *, num_recip,'NUM_RECIP exceeds number of available PEs'
     $        ,numtasks
         num_recip = numtasks
      endif

      do i=1,num_recip
         ranks(i)=i-1
      enddo
      if(mytaskid.lt.num_recip)then
         i_do_recip = .TRUE.
      else
         i_do_recip = .FALSE.
      endif
      call MPI_GROUP_INCL(world_group,num_recip,ranks,recip_group
     $     ,err)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_COMM_CREATE(MPI_COMM_WORLD,recip_group,recip_comm
     $     ,err)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c         print *,mytaskid,"> numtasks ",numtasks
c         call flush(6)
c      do j=1,numtasks
c         if(mytaskid.eq.j-1)then
c            print *,"Recip Group Creation:  PE ",mytaskid,i_do_recip
c            print *,"num_recip ",num_recip
c     $           ,"ranks ",(ranks(i),i=1,numtasks)
c            print *,"rgroup,rcomm ",recip_group,recip_comm
c            call flush(6)
c         endif
c         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c         print *,mytaskid,"> numtasks ",numtasks,j
c         call flush(6)
c      enddo

c         print *,"rgroup,rcomm ",recip_group,recip_comm,mytaskid
c         call flush(6)
      if(num_recip.eq.numtasks)then
         i_do_direct = .TRUE.
         num_direct = numtasks
         do i=1,numtasks
            ranks(i)=i-1
         enddo
      else
         i_do_direct = .FALSE.
         num_direct = numtasks-num_recip
         do i=1,numtasks-num_recip
            ranks(i)=num_recip-1+i
            if(mytaskid.eq.num_recip-1+i)i_do_direct = .TRUE.
         enddo
      endif
      call MPI_GROUP_INCL(world_group,num_direct,ranks,direct_group
     $     ,err)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_COMM_CREATE(MPI_COMM_WORLD,direct_group,direct_comm
     $     ,err)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      

      do j=1,numtasks
      if(mytaskid.eq.j-1)then
c         print *,"Direct Group Creation: PE ",mytaskid,i_do_direct
c         print *,"num_direct ",num_direct
c     $        ,"ranks ",(ranks(i),i=1,numtasks)
c         print *,"dgroup,dcomm ",direct_group,direct_comm,mytaskid
c         call flush(6)
      endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      enddo

c         print *,'Leaving g start PE',mytaskid,numtasks,direct_comm
c         call flush(6)
      return
      end
#endif MPI
