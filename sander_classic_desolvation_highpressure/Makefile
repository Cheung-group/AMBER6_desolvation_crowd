#
#************************************************************************
#                             AMBER                                    **
#                                                                      **
#                Copyright (c) 1986, 1991, 1995 , 1995                 **
#             Regents of the University of California                  **
#                       All Rights Reserved.                           ** 
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
# Makefile for UCSF Amber 5   SANDER 
#
SHELL=/bin/sh
VPATH=
PREC=DPREC
include ../sfx.h

SRC= 	sander.f cshf.f noecalc.f caldis.f nonbon.f \
	calrate.f dinten.f drates.f indexn.f kmat.f pearsn.f \
	erfcc.f plane.f remarc.f nmrcal.f nmrred.f \
	restal.f getnat.f nmrnrg.f modwt.f disnrg.f angnrg.f \
	tornrg.f nmrprt.f nmrgrp.f nmrcms.f nmrcmf.f \
	impnum.f nmrrad.f nmrsht.f at2res.f chklin.f opnmrg.f \
	printe.f runmin.f ndvprt.f force.f rdparm.f \
	mdread.f locmem.f runmd.f newvel.f getcor.f \
	r6ave.f r6drv.f aveint.f degcnt.f decnvh.f corf.f \
	efield.f indip.f threeb.f tripl.f politr.f polder.f\
	fastwt.f solpnt.f echoin.f jnrg.f parallel.f \
	misc.f shake.f qiktip.f ene.f mdwrit.f minrit.f \
	resnba.f set.f setmm.f dynlib.f period.f mdfil.f nmlsrc.f \
	ew_force.f ew_setup.f ew_box.f ew_extra.f \
	ew_fft.f ew_direct.f ew_recip.f cosmo.f pcshift.f \
        wrtmny.f wrtint.f intclc.f redarc.f lgandb.f lgcall.f \
        lgclc.f lgread.f lgrst.f lgset.f lgwrt.f ovrlap.f lgcoor.f \
        lgnbad.f lgrfcl.f avefrc.f popbck.f

OBJ=  sander.o cshf.o noecalc.o caldis.o nonbon.o \
    calrate.o dinten.o drates.o indexn.o kmat.o pearsn.o \
    erfcc.o plane.o remarc.o nmrcal.o nmrred.o \
    restal.o getnat.o nmrnrg.o modwt.o disnrg.o angnrg.o \
    tornrg.o nmrprt.o nmrgrp.o nmrcms.o nmrcmf.o \
    impnum.o nmrrad.o nmrsht.o at2res.o chklin.o opnmrg.o \
    printe.o runmin.o ndvprt.o force.o rdparm.o \
    mdread.o locmem.o runmd.o newvel.o getcor.o \
    r6ave.o r6drv.o aveint.o degcnt.o decnvh.o corf.o \
    efield.o indip.o threeb.o tripl.o politr.o polder.o\
    fastwt.o solpnt.o echoin.o parallel.o jnrg.o \
    misc.o shake.o qiktip.o ene.o mdwrit.o minrit.o \
    resnba.o set.o setmm.o dynlib.o period.o mdfil.o nmlsrc.o \
    ew_force.o ew_setup.o ew_box.o ew_extra.o \
    ew_fft.o ew_direct.o ew_recip.o cosmo.o pcshift.o \
    wrtmny.o wrtint.o intclc.o redarc.o  lgandb.o lgcall.o \
    lgclc.o lgread.o lgrst.o lgset.o lgwrt.o ovrlap.o lgcoor.o \
    lgnbad.o lgrfcl.o avefrc.o popbck.o

LESOBJ=  sander.o cshf.o noecalc.o caldis.o nonbon.LES.o \
    calrate.o dinten.o drates.o indexn.o kmat.o pearsn.o \
    erfcc.o plane.o remarc.o nmrcal.o nmrred.o \
    restal.o getnat.o nmrnrg.o modwt.o disnrg.o angnrg.o \
    tornrg.o nmrprt.o nmrgrp.o nmrcms.o nmrcmf.o \
    impnum.o nmrrad.o nmrsht.o at2res.o chklin.o opnmrg.o \
    printe.o runmin.o ndvprt.o force.o rdparm.LES.o \
    mdread.LES.o locmem.o runmd.o newvel.o getcor.o \
    r6ave.LES.o r6drv.o aveint.o degcnt.o decnvh.o corf.o \
    efield.o indip.o threeb.o tripl.o politr.o polder.o\
    fastwt.o solpnt.o echoin.o parallel.LES.o jnrg.LES.o \
    misc.o shake.o qiktip.o ene.LES.o mdwrit.o minrit.o \
    resnba.o set.o setmm.o dynlib.o period.o mdfil.o nmlsrc.o \
    ew_force.o ew_setup.o ew_box.o ew_extra.o \
    ew_fft.o ew_direct.o ew_recip.o cosmo.o pcshift.o \
    wrtmny.o wrtint.o intclc.o redarc.o  lgandb.o lgcall.o \
    lgclc.o lgread.o lgrst.o lgset.o lgwrt.o ovrlap.o lgcoor.o \
    lgnbad.o lgrfcl.o avefrc.o popbck.o


# LIBSRC & LIBPSRC tell us the remaining lib files req'd for bigsource


LIBSRC=../lib/grdmax.f ../lib/rfree.f \
	../lib/rgroup.f ../lib/ddotp.f  ../lib/eis.f   \
	../lib/random.f ../lib/matinv.f ../lib/lsqfit.f \
    	../lib/util.f ../lib/amopen.f ../lib/namlst90.f ../lib/mexit.f

LIBOBJ=../lib/grdmax.o ../lib/rfree.o  \
	../lib/rgroup.o ../lib/ddotp.o  ../lib/eis.o  \
	../lib/random.o ../lib/matinv.o ../lib/lsqfit.o \
    	../lib/util.o ../lib/amopen.o ../lib/namlst90.o ../lib/mexit.o \
	../lib/malloc.o

all::    sander_classic$(SFX) sander.LES_classic$(SFX)

sander_classic$(SFX):		$(OBJ) libobj syslib lapack
		SYSLIB=`../sysdir lib` ; ../Compile LOAD -o sander_classic$(SFX) \
			$(OBJ) ../lapack/lapack.a   $(LIBOBJ) $$SYSLIB;
		
sander.LES_classic$(SFX):		$(LESOBJ) libobj syslib lapack
		SYSLIB=`../sysdir lib` ; ../Compile LOAD -o sander.LES_classic$(SFX) \
			$(LESOBJ) ../lapack/lapack.a $(LIBOBJ) $$SYSLIB;
		

bigsource: 	
	SYSDIR=`../sysdir dir`  ;  cd $$SYSDIR ; make sys.f
	SYSSRC=`../sysdir src` ;  ../Compile CPPONLY -DDPREC \
			-o sander_all.for  \
			$(SRC) $(LIBSRC) $$SYSSRC
	@echo "SANDER source in sander_all.for"

locsource:
	-/bin/rm *.for
	../Compile CPPONLY  $(SRC) 
	@echo "SANDER modules source now in xxx.for"

#----------HEADER DEPENDENCIES (if only 'make depend' worked w/ fortran..)

#
# box.h
#
dynlib.o efield.o ene.o force.o indip.o locmem.o mdread.o nonbon.o \
	parallel.o period.o polder.o politr.o rdparm.o resnba.o \
	runmd.o runmin.o sander.o sander_init.o setmm.o shake.o \
	threeb.o tripl.o ew_setup.o: box.h

#
# cosmo.h
#
cosmo.o: cosmo.h

#
# ew_adj.h
#
ew_force.o ew_setup.o parallel.o: ew_adj.h

#
# ew_cntrl.h
#
ew_box.o ew_direct.o ew_extra.o ew_force.o ew_setup.o \
	parallel.o: ew_cntrl.h

#
# ew_erfc_spline.h
#
ew_force.o ew_setup.o parallel.o: ew_erfc_spline.h

#
# ew_localnb.h
#
ew_direct.o ew_force.o ew_setup.o parallel.o: ew_localnb.h

#
# ew_nbe.h
#
ew_direct.o: ew_nbe.h

#
# ew_numtasks.h
#
ew_direct.o ew_force.o ew_setup.o parallel.o mdfil.o: ew_numtasks.h

#
# ew_parallel.h
#
ew_box.o ew_direct.o ew_fft.o ew_force.o ew_recip.o \
	ew_setup.o parallel.o: ew_parallel.h

#
# ew_pme_recip.h
#
ew_force.o ew_setup.o parallel.o: ew_pme_recip.h

#
# ew_time.h
#
ew_direct.o ew_force.o ew_setup.o: ew_time.h

#
# ew_unitcell.h
#
ew_box.o ew_direct.o ew_force.o ew_setup.o parallel.o:	ew_unitcell.h

#
# extra.h
#
dynlib.o ew_box.o ew_direct.o ew_direct.o:		extra.h

ew_force.o ew_setup.o ew_setup.o ew_setup.o:		extra.h

force.o forcemp.o mdread.o misc.o parallel.o:		extra.h

resnba.o runmd.o runmin.o sander.o set.o shake.o:	extra.h

#
# files.h
#
cshf.o dynlib.o getcor.o mdfil.o mdread.o mdwrit.o minrit.o \
	ndvprt.o nmrnrg.o nmrred.o noecalc.o parallel.o \
	runmd.o runmin.o sander.o sander_init.o \
	sander_min.o set.o pcshift.o: files.h

#
# fujfft.h
#
ew_fft.o ew_recip.o: fujfft.h

#
# iewald.h
#
dynlib.o ew_direct.o ew_force.o force.o getcor.o locmem.o \
	mdread.o mdwrit.o misc.o parallel.o rdparm.o resnba.o \
	runmd.o runmin.o sander.o: iewald.h

#
# les.h
#
ene.LES.o nonbon.LES.o parallel.LES.o r6ave.LES.o rdparm.LES.o \
   jnrg.LES.o:  les.h
#
# lgcom.h
#
lgandb.o lgcall.o lgclc.o lgread.o lgrst.o lgset.o lgwrt.o ovrlap.o \
	lgcoor.o lgnbad.o force.o mdread.o sander.o runmd.o lgrfcl.o:   lgcom.h

#
# avfrc.h
#
sander.o runmd.o force.o avefrc.o mdread.o:   avfrc.h

#
# md.h
#
corf.o cshf.o dynlib.o ene.o fastwt.o force.o locmem.o mdread.o \
	minrit.o misc.o noecalc.o parallel.o politr.o remarc.o \
	resnba.o rdparm.o runmd.o runmin.o sander.o sander_force.o \
	sander_init.o sander_min.o shake.o pcshift.o: md.h

#
# memory.h
#
force.o indexn.o locmem.o mdread.o noecalc.o parallel.o rdparm.o \
	remarc.o resnba.o runmd.o runmin.o sander.o \
	sander_force.o sander_init.o sander_min.o set.o setmm.o: memory.h

#
# mpif.h
#
dynlib.o ew_fft.o ew_force.o ew_recip.o fastwt.o force.o mdread.o \
    misc.o nmrcal.o nonbon.o parallel.o resnba.o runmd.o sander.o \
    sander_init.o set.o: ../Machines/mpi/mpif.h



#
# nmr.h
#
caldis.o calrate.o corf.o cshf.o dinten.o drates.o \
	force.o indexn.o kmat.o locmem.o mdread.o modwt.o \
	ndvprt.o newvel.o nmrprt.o noecalc.o parallel.o \
	rdparm.o remarc.o runmd.o runmin.o sander.o \
	sander_init.o pcshift.o: nmr.h

#
# other.h
#

mdread.o runmin.o runmd.o sander.o:   other.h

#
# ../lib/parallel.h
#
dynlib.o ene.o ew_box.o ew_direct.o ew_fft.o ew_force.o \
	ew_recip.o ew_setup.o fastwt.o force.o locmem.o \
	misc.o nonbon.o parallel.o qiktip.o \
	resnba.o runmd.o runmin.o sander.o \
	sander_init.o set.o shake.o: 		../lib/parallel.h


#
# parms.h
#
ene.o force.o mdread.o parallel.o rdparm.o sander.o sander_init.o \
	set.o: parms.h

#
# pcshift.h
#
pcshift.o: pcshift.h

#
# pol.h
#
force.o locmem.o mdread.o rdparm.o threeb.o tripl.o: pol.h

#
# sizes.h
#
force.o locmem.o resnba.o sander.o sander_force.o sander_init.o \
	sander_min.o set.o : sizes.h

#
# vector.h
#
period.o force.o qiktip.o indip.o polder.o nonbon.o fastwt.o cshf.o \
	remarc.o resnba.o shakemp.o pcshift.o:	../lib/vector.h

#
# shared memory file dependencies
#

force.o: force.f forcemp.f 

resnba.o resnbamp.o: resnba.f resnbamp.f

shake.o: shake.f shakemp.f 


#-----------LIBS

libobj:
		cd ../lib; ./Makelist $(LIBOBJ)
		
syslib:
		( SYSDIR=`../sysdir dir` ; echo sysdir is $$SYSDIR ; \
			cd $$SYSDIR ; make sys.a )

#-----------

lapack:
		cd ../lapack; make
 
#-----------

malloc.o:	malloc.c
		../Compile CC malloc.c

#-----------

rdparm.o:		rdparm.f
			../Compile L1 -P -D$(PREC) rdparm.f 

#----------


getcor.o:	getcor.f
		../Compile L1 -P -D$(PREC) -o $@  getcor.f


#-----------OBJ L3 : compile at highest optimization level

nonbon.o:	nonbon.f 
		../Compile L3 -P -D$(PREC) -o $@  nonbon.f 

nonbon.LES.o:	nonbon.f 
		../Compile L3 -P -DLES -D$(PREC) -o $@  nonbon.f 

caldis.o:	caldis.f
		../Compile L3 -P  -D$(PREC) $*.f

ene.o:		ene.f
		../Compile L3 -P  -D$(PREC) $*.f

ene.LES.o:		ene.f
		../Compile L3 -P  -DLES -D$(PREC) -o $@ ene.f

period.o:	period.f
		../Compile L3 -P  -D$(PREC) $*.f

calrate.o:	calrate.f
		../Compile L3 -P  -D$(PREC) $*.f

cshf.o:		cshf.f
		../Compile L3 -P -D$(PREC) cshf.f

dinten.o:	dinten.f
		../Compile L3 -P  -D$(PREC) $*.f

drates.o:	drates.f
		../Compile L3 -P  -D$(PREC) $*.f

fastwt.o:	fastwt.f
		../Compile L3 -o $@ -P -D$(PREC) $*.f

indexn.o:     	indexn.f
		../Compile L3 -P  -D$(PREC) $*.f

kmat.o:		kmat.f
		../Compile L3 -P  -D$(PREC) $*.f

remarc.o:	remarc.f
		../Compile L3 -DINVWGT  -D$(PREC) $*.f

qiktip.o:	qiktip.f
		../Compile L3 -P  -D$(PREC) qiktip.f

shake.o:	shake.f
		../Compile L3 -P  -D$(PREC) shake.f

efield.o:	efield.f
		../Compile L3 -P  -D$(PREC) efield.f

indip.o:	indip.f
		../Compile L3 -P  -D$(PREC) indip.f

pcshift.o:		pcshift.f
		../Compile L3 -P -D$(PREC) pcshift.f

politr.o:	politr.f
		../Compile L3 -P  -D$(PREC) politr.f

polder.o:	polder.f
		../Compile L3 -P  -D$(PREC) polder.f

tripl.o:	tripl.f
		../Compile L3 -P  -D$(PREC) tripl.f

threeb.o:	threeb.f
		../Compile L3 -P  -D$(PREC) threeb.f

resnba.o:	resnba.f
		../Compile L3 -P  -D$(PREC) resnba.f

ew_box.o:	ew_box.f 
		../Compile L3 -P  -D$(PREC) ew_box.f

ew_direct.o:	ew_direct.f 
		../Compile L3 -P  -D$(PREC) ew_direct.f

ew_extra.o:	ew_extra.f 
		../Compile L3 -P  -D$(PREC) ew_extra.f

ew_fft.o:	ew_fft.f 
		../Compile L3 -P  -D$(PREC) ew_fft.f

ew_force.o:	ew_force.f 
		../Compile L3 -P  -D$(PREC) ew_force.f

ew_recip.o:	ew_recip.f 
		../Compile L3 -P  -D$(PREC) ew_recip.f

ew_setup.o:	ew_setup.f 
		../Compile L3 -P  -D$(PREC) ew_setup.f

lgclc.o:	lgclc.f    
		../Compile L3 -P  -D$(PREC) lgclc.f   

lgnbad.o:	lgnbad.f    
		../Compile L3 -P  -D$(PREC) lgnbad.f   

lgrfcl.o:	lgrfcl.f    
		../Compile L3 -P  -D$(PREC) lgrfcl.f   

avefrc.o:	avefrc.f    
		../Compile L3 -P  -D$(PREC) avefrc.f   

popbck.o:	popbck.f    
		../Compile L3 -P  -D$(PREC) popbck.f   

#-----------ordinary compilations:

#
# default rule

.f.o:		$<
		../Compile L2 -P -D$(PREC) $<


#do not optimize nmrcal, due to bug on mips/SI compilers 
nmrcal.o:	nmrcal.f
		../Compile L0 -P -D$(PREC) nmrcal.f

#do not optimize nmrgrp, due to bug on HP compiler
nmrgrp.o:	nmrgrp.f
		../Compile L0 -P -D$(PREC) -o $@  nmrgrp.f

#locmem mdfil mdread don't need optimization
locmem.o:	locmem.f
		../Compile L0 -P -D$(PREC) locmem.f

mdfil.o:	mdfil.f
		../Compile L0 -P -D$(PREC) mdfil.f

mdread.o:	mdread.f
		../Compile L0 -P -D$(PREC) mdread.f

mdread.LES.o:	mdread.f
		../Compile L0 -P -DLES -D$(PREC) -o $@ mdread.f

jnrg.LES.o:  jnrg.f
		../Compile L2 -P -DLES -D$(PREC) -o $@ jnrg.f

parallel.LES.o:  parallel.f
		../Compile L2 -P -DLES -D$(PREC) -o $@ parallel.f

r6ave.LES.o:  r6ave.f
		../Compile L2 -P -DLES -D$(PREC) -o $@ r6ave.f

rdparm.LES.o:  rdparm.f
		../Compile L2 -P -DLES -D$(PREC) -o $@ rdparm.f

wrtmny.o:   wrtmny.f
		../Compile L2 -DDPREC -o $@  wrtmny.f

wrtint.o:   wrtint.f
		../Compile L2 -DDPREC -o $@  wrtint.f

intclc.o:   intclc.f
		../Compile L2 -DDPREC -o $@  intclc.f

redarc.o:   redarc.f
		../Compile L2 -DDPREC -o $@  redarc.f


#-----------

install:	all 
		/bin/mv sander_classic$(SFX) ../../exe
		/bin/mv sander.LES_classic$(SFX) ../../exe

clean:
		../Clean sander_classic$(SFX) sander.LES_classic$(SFX)

