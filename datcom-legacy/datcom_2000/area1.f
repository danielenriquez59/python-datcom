      SUBROUTINE AREA1 (X,Y,AREA,NSUM)
C
C***  CALCULATES INCREMEMTAL AREA OF VERTICAL TAIL SHADOWED BY MACH LINE
C
      DIMENSION X(4),Y(4)
      AREA2=0.
      A=SQRT((X(2)-X(1))**2+(Y(2)-Y(1))**2)
      B=SQRT((X(3)-X(2))**2+(Y(3)-Y(2))**2)
      C=SQRT((X(1)-X(3))**2+(Y(1)-Y(3))**2)
      S=(A+B+C)/2.
      AREA=SQRT(S*(S-A)*(S-B)*(S-C))
      IF(NSUM.NE.4.AND.NSUM.NE.6) GO TO 40
      A2=SQRT((X(4)-X(3))**2+(Y(4)-Y(3))**2)
      B2=SQRT((X(1)-X(4))**2+(Y(1)-Y(4))**2)
      S2=(A2+B2+C)/2.
      AREA2=SQRT(S2*(S2-A2)*(S2-B2)*(S2-C))
   40 AREA=AREA+AREA2
      RETURN
      END
