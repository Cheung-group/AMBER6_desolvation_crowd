c margaret added. 10.3.00 double well potential. real.
c 9.15.00 r<=r0
c      COMMON/WID/kpower,npower,mpower
cc March 25,2014, modified by Jianfa
c      COMMON/DOU/DB0(200,200),DB1(200,200),CB1(200,200),
c     +     B2(200,200),h1(200,200),h2(200,200)
      COMMON/DOU/DB0(2000,2000),DB1(2000,2000),CB1(2000,2000),
     +     B2(2000,2000),h1(2000,2000),h2(2000,2000)
c margaret 9.16.00 widened douwell well potential. real.
c      COMMON/DIST/r0(200,200),r1(200,200),r2(200,200)
      COMMON/DIST/r0(2000,2000),r1(2000,2000),r2(2000,2000)
c
      COMMON/POWER/kpower,npower,mpower
c     check contact map
c      COMMON/MAP/mapgo(100,100)
      COMMON/GOMAP/mapgo(1000,1000)
     
cc 	end modified by Jianfa
