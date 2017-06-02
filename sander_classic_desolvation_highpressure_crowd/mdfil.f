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
      subroutine mdfil
c     Author: George Seibel
*     implicit none
c
c     OUTPUT: (to common)
c
#include "files.h"
#   include "ew_numtasks.h"
#ifdef MPI
#   include "parallel.h"
#   include "ew_parallel.h"
#endif
c
c     INTERNAL:
c
      character*80 arg,nrecip
c        ... temp for each of the whitespace delimited command line words
      integer iarg
c        ... arg pointer, final number of arguments
c
c     --- default file names ---
c
      mdin   = 'mdin'
      mdout  = 'mdout'
      inpcrd = 'inpcrd'
      parm   = 'prmtop'
      restrt = 'restrt'
      refc   = 'refc'
      mdvel  = 'mdvel'
      mden   = 'mden'
      mdcrd  = 'mdcrd'
      mdinfo = 'mdinfo'
      vecs   = 'vecs'
      freqe   = 'dummy'
c Margaret add
      twhb = 'twhb'
      twvdw = 'twvdw'
      twchi = 'twchi'
      fenpc = 'fenpc'
      fencc = 'fencc'
c
c     --- default status of output: new
c
#ifdef MPI
      num_recip=numtasks
      num_direct=numtasks
#endif
#ifdef ARCHIVE
      owrite = 'U'
#else
      owrite = 'N'
#endif

#ifdef ARCHIVE
#  define NOGETARG
#endif
#ifndef NOGETARG
c
c     --- get com line arguments ---
c
# ifdef HITACHI_GETARG
      iarg = 1
# else
      iarg = 0
# endif
      indx = iargc()
# ifdef MPI
#  if TCGMSG
      iarg = iarg + 1
      indx = indx - 7
#  endif
# endif
      if (indx .eq. 0) goto 20
   10 continue
           iarg = iarg + 1
           call getarg(iarg,arg)
           if (arg .eq. '-O') then
                owrite = 'U'
#ifdef CSPP
           elseif (arg .eq. '-np') then
                iarg = iarg + 1
#endif
           elseif (arg .eq. '-i') then
                iarg = iarg + 1
                call getarg(iarg,mdin)
           elseif (arg .eq. '-o') then
                iarg = iarg + 1
                call getarg(iarg,mdout)
           elseif (arg .eq. '-p') then
                iarg = iarg + 1
                call getarg(iarg,parm)
           elseif (arg .eq. '-c') then
                iarg = iarg + 1
                call getarg(iarg,inpcrd)
           elseif (arg .eq. '-vecs') then
                iarg = iarg + 1
                call getarg(iarg,vecs)
           elseif (arg .eq. '-f') then
                iarg = iarg + 1
                call getarg(iarg,freqe)
           elseif (arg .eq. '-r') then
                iarg = iarg + 1
                call getarg(iarg,restrt)
           elseif (arg .eq. '-ref' .or. arg .eq.'-z') then
                iarg = iarg + 1
                call getarg(iarg,refc)
           elseif (arg .eq. '-e') then
                iarg = iarg + 1
                call getarg(iarg,mden)
           elseif (arg .eq. '-v') then
                iarg = iarg + 1
                call getarg(iarg,mdvel)
           elseif (arg .eq. '-x'.or.arg .eq.'-t') then
                iarg = iarg + 1
                call getarg(iarg,mdcrd)
           elseif (arg .eq. '-inf') then
                iarg = iarg + 1
                call getarg(iarg,mdinfo)
c Margaret add 
           elseif (arg .eq. '-b') then
                iarg = iarg + 1
                call getarg(iarg,twhb)
           elseif (arg .eq. '-d') then
                iarg = iarg + 1
                call getarg(iarg,twvdw)
           elseif (arg .eq. '-chiral') then
                iarg = iarg + 1
                call getarg(iarg,twchi)
           elseif (arg .eq. '-pc') then
                iarg = iarg + 1
                call getarg(iarg,fenpc)
           elseif (arg .eq. '-cc') then
                iarg = iarg + 1
                call getarg(iarg,fencc)
c end Margaret

           elseif (arg .eq. '-help') then
                write(0,9000)
                write(6,9000)
                call mexit(6,1)
#ifdef MPI
           elseif (arg .eq. '-p4pg') then
                iarg = iarg+1
           elseif (arg .eq. '-p4wd') then
                iarg = iarg+1
           elseif (arg .eq. '-np') then
                iarg = iarg+1
           elseif (arg .eq. '-mpedbg') then
              continue
           elseif (arg .eq. '-dbx') then
              continue
           elseif (arg .eq. '-gdb') then
              continue
#endif
#ifdef T3D
           elseif (arg .eq. '-freeze') then
                iarg = iarg+1
           elseif (arg .eq. '-debug') then
                iarg = iarg+1
           elseif (arg .eq. '-pool') then
                iarg = iarg+2
#endif
#ifdef MPI
           elseif ( arg .eq. '-nrecip') then
                iarg = iarg + 1
                call getarg(iarg,arg)
                read(arg,'(i5)',err=91)num_recip
                if(num_recip.eq.numtasks)then
                   num_direct=numtasks
                else
                   num_direct=numtasks-num_recip
                endif
                print *,"NUM_RECIP ",num_recip
     $               ,"  NUM_DIRECT ",num_direct
#endif
           else
                if (arg .eq. ' ') go to 20
                write(6,'(/,5x,a,a)') 'unknown flag: ',arg
                write(6,9000)
                call mexit(6, 1)
           endif
      if (iarg .lt. indx) go to 10
c
   20 continue
c
#endif
#ifdef MPI
      if(num_recip.eq.0)num_recip=numtasks
#endif
      return
 91   write(6,*)'....flag:  "-nrecip" expects an integer argument'
      write(6,*)'                   '
      call mexit(6, 1)
 9000 format(/,5x,
     +   'usage: sander  [-O] -i mdin -o mdout -p prmtop -c inpcrd ',
     +   '-r restrt',/19x,'[-ref refc -x mdcrd -v mdvel -e mden ',
     +   '-inf mdinfo]')
      end
