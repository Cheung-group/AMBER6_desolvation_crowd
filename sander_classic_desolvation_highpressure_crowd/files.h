      character*80 mdin, mdout, inpcrd, parm, restrt, 
     +             refc, mdvel, mden, mdcrd, mdinfo, nmr, mincor,
     +           vecs, freqe,redir,twhb,twvdw,twchi,fenpc,fencc
      character owrite
      common /files/ mdin, mdout, inpcrd, parm, restrt, 
     +                refc, mdvel, mden, mdcrd, mdinfo, nmr, mincor,
     +             vecs, freqe, owrite,twhb,twvdw,twchi,fenpc,fencc
      COMMON/RUNHED/ITITL(20),ITITL1(20)
      COMMON/HULP/NTPR,NTWR,NTWX,NTWXM,NTWV,NTWVM,NTWE,NTWEM,NTPP,
     *            IOUTFM,NTWPRT,NTWPR0
C      NMRRDR : Contains information about input/output file redirection
C               REDIR(6) and IREDIR(6) contain information regarding
C               LISTIN, LISTOUT, READNMR, NOESY, SHIFTS, DUMPAVE,
C               and PCSHIFT respectively. If IREDIR(I) > 0,
C               then that input/output has been redirected.
      COMMON/NMRRDR/REDIR(7),IREDIR(7)
