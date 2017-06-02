      COMMON/RPARMS/RK(5000),REQ(5000),TK(5000),TEQ(5000),
     +     PK(3000),
     +     PN(3000),PHASE(3000),CN1(480000),CN2(480000),SOLTY(480000),
     +     GAMC(3000),GAMS(3000),FMN(3000),
     +     ASOL(480000),BSOL(480000),HBCUT(480000)
      common/iparms/IPN(3000)

#ifdef CHARMM
      common/p14/ cn114(15000),cn214(15000)
      common/ub/rkub(900),rub(900)
#endif
C
C NPHB is the number of h-bond parameters. NIMPRP is the number of
C improper torsional parameters (NPTRA-NIMPRP is the number of regular
C torsional parameters).
C
      COMMON/PRMLIM/NUMBND,NUMANG,NPTRA,NPHB,NIMPRP
