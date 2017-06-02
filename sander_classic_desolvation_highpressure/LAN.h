c margaret added 12.15.00 for langevin dynamics
c generate random numbers.
      common/rand/JSEED(4)   
cc eta=6*pi*radius*wvisco
cc wvisco=0.1439 kcal*ps/mol*(AA**3) * 20.445  
cc 20.445(convert to internal unit)
c boltzman: 1.9872x10^(-3) kcal/moleK  
      COMMON/VISCO/eta,WVISCO,RRADIUS,BOLTZ
