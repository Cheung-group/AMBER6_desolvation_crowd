c margaret added. 9.18.01 Hydrogen Bond. real.
c input GO hydrogen bonds
c store atom number of theta,phi angles 
c H-B between JTHEI-JPHII, JTHEI<JPHII
      COMMON/HYDROB/NXI(92000),NXI1(92000),NXI2(92000),
     +     NXI3(92000)
c store the native Go phi,psi angles 
      COMMON/HBANG/THE0(10000),PHI0(10000)
c     check contact HB map
      COMMON/MAP/maphb(10000,10000)
c     output ener 
      COMMON/ENERTW/EHBA,EHBV,ECHI,XKCHI
c     input Amp 
	COMMON/DISPLACE/Amp
