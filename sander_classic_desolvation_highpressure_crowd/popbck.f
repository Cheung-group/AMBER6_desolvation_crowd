      SUBROUTINE POPBCK(X,F,NATOM,SAFETY,IOUT)
C
C Subroutine POP BaCK
C
C This routine calculates the center of geometry of the system for
C a vacuum simulation. If the COG is within SAFETY Angstroms of
C + 9999.999 or SAFETY Angstroms of -999.999 in any dimension, 
C then the system is translated to relocate the COG at (0,0,0). 
C This is done so that
C the system doesn't translate so much that coordinates will
C overflow the f8.3 fields in the output coordinate files.
C
C Author: David Pearlman
C Date: 7/98

#ifdef DPREC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#endif

      DIMENSION X(3,*),COG(3)
C
      COG(1) = 0.0D0
      COG(2) = 0.0D0
      COG(3) = 0.0D0
      DO 10 I = 1,NATOM
        COG(1) = COG(1) + X(1,I)
        COG(2) = COG(2) + X(2,I)
        COG(3) = COG(3) + X(3,I)
   10 CONTINUE
      COG(1) = COG(1)/FLOAT(NATOM)
      COG(2) = COG(2)/FLOAT(NATOM)
      COG(3) = COG(3)/FLOAT(NATOM)

      IOK = 0
      DO 20 I = 1,3
         IF (COG(I).LT.-999.999D0+SAFETY) IOK = 1
         IF (COG(I).GT.9999.999D0-SAFETY) IOK = 1
   20 CONTINUE
      IF (IOK.EQ.1) THEN
         DO 30 J = 1,3
            DO 40 I = 1,NATOM
	    write(82,*)j,i,X(J,I),F(J,I)
               X(J,I) = X(J,I) - COG(J)
   40       CONTINUE
   30    CONTINUE
         WRITE(IOUT,9000)
 9000    FORMAT('COG too close to -999 or +9999; ',/,
     *          'COG translated back to [0,0,0]')
      END IF

      RETURN
      END
      


