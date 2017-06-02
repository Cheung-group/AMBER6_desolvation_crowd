      parameter (NUMSTS=30)
c
c ... integer variables:
c
      COMMON/MDI/NRUN,NRP,NSM,NRAM,NSPM,IG,
     +   NTX,NTCX,NTXO,NTT,NTP,NTR,INIT,NTCM,NSCM,NTU,
     +      ISOLVP,NSOLUT,NTC,NTCC,NTF,NTID,NTN,NTNB,NSNB,NDFMIN,
     +      NSTLIM,NRC,NTRX,NPSCAL,IMIN,MAXCYC,NCYC,NTMIN,
     +      IREST,JFASTW,IWTNM,IOWTNM,IHWTNM(2),IBGWAT,IENWAT,IORWAT,
     +      IWATPR,NSOLW
c
c ... floats:
c
      COMMON/MDR/T,DT,TEMP0,TAUTP,TAUTS,PRES0,COMP,TAUP,TEMP,TEMPI,
     +      HEAT,TOL,DTEMP,TAUR,DX0,DXM,DRMS,TIMLIM,TIMTOT,TIMDEL,
     +      TIMRUN,TIMSTS(NUMSTS),TAUV,TAUV0,VZERO,VLIMIT,RBTARG(8),
     +      tmass,tmassinv
cc April 4,2014, modified by Jianfa Chen
       real(kind=4):: wtime0,wtime1,wtim0,wtim1,time0,time1,UTIM,tim1,tim2,tim
cc end modified by Jianfa
