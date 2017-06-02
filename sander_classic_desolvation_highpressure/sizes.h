c
c  --- following four lines give "master" sizes for sander.  These
c      generally are what need to be changed when resizing the code.
c
c  Note that when MEM_ALLOC is set in the MACHINE file, only
c  MAXDUP (of the params in in this file) affects actual dimensioning.
c
c
c  --- "standard" (pretty-big) sizes:
c
cc March 24,2014, modified by Jianfa Chen
c      parameter (MAXINT=2000000)
c      parameter (MAXPR=4500000)
c      parameter (MAXREA=3700000)
c      parameter (MAXHOL=600000)
c      parameter (MAXDUP=8000)
      parameter (MAXINT=3000000)
      parameter (MAXPR=4500000)
      parameter (MAXREA=3700000)
      parameter (MAXHOL=600000)
      parameter (MAXDUP=150000)
c   end modified by Jianfa Chen
c
c  --- If the sizes above are too small, the program will complain.
c      Here are some typcial values:
c
c      Use       Parameter      Typical Values for Large Systems
c   floating pt    MAXREA       20 or 23 * Natom (23 for const pressure)
c   integers       MAXPR        cutoff-dependent
c   integers       MAXINT      5 - 10 * Natom??, data dependent
c   holerith       MAXHOL       25 * Natom
c   dihedral dup   MAXDUP      0  - 5000,  data dependent
c
c Integer needs are increased by the presence of large residues and
c use of the all-atom model instead of the united atom model.
c See the users manual for more information.
c
c  Note: if you are carrying out NMR-based refinement, you may also
c need to set some sizes in the file "nmr.h".
c
c  Packing scheme:
c
#ifdef CRAYFISH
c     --- Cray packing scheme ---
      parameter (IPACK=1)
#else
c    ---- generic packing scheme ---
      parameter (IPACK=0)
#endif
