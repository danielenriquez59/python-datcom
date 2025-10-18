      SUBROUTINE AREA2(X,Y,AREA,INUM,AX,AY)
C
C***  CALCULATES INCREMENTAL AREA OF BODY SHADOWED BY MACH LINE
C
      DIMENSION X(4),Y(4)
C
C**** HERE FOR TRIANGLE OF TYPE 1 ****
C
      A=SQRT((X(2)-X(1))**2+(Y(2)-Y(1))**2)
      B=SQRT((X(3)-X(2))**2+(Y(3)-Y(2))**2)
      C=SQRT((X(1)-X(3))**2+(Y(1)-Y(3))**2)
      S=(A+B+C)/2.
      AREA=SQRT(S*(S-A)*(S-B)*(S-C))
      X1=(X(1)+X(2))/2.
      Y1=(Y(1)+Y(2))/2.
      X2=(X(1)+X(3))/2.
      Y2=(Y(1)+Y(3))/2.
      DNUM1=(Y(3)-Y1)/(X(3)-X1)
      DNUM2=(Y2-Y(2))/(X2-X(2))
      ANUM1=DNUM1*X1-Y1
      ANUM2=DNUM2*X(2)-Y(2)
      XBAR=(ANUM1-ANUM2)/(DNUM1-DNUM2)
      YBAR=(XBAR-X1)*DNUM1+Y1
      AX=AREA*XBAR
      AY=AREA*YBAR
      IF(INUM.EQ.3) GO TO 70
   30 A2=SQRT((X(4)-X(3))**2+(Y(4)-Y(3))**2)
C
C**** HERE FOR TRIANGLE OF TYP2 1 ****
C
      B2=SQRT((X(1)-X(4))**2+(Y(1)-Y(4))**2)
      C2=SQRT((X(3)-X(1))**2+(Y(3)-Y(1))**2)
      S2=(A2+B2+C2)/2.
      AREA22=SQRT(S2*(S2-A2)*(S2-B2)*(S2-C2))
      YBAR2=Y(3)/3.
      XBAR2=(YBAR2-Y(1))*(X(3)-(X(1)+X(4))/2.)/(Y(3)-Y(1))+(X(1)+X(4))/2
     1.
      AX2=AREA22*XBAR2
      AY2=AREA22*YBAR2
      AREA=AREA+AREA22
      AX=AX+AX2
      AY=AY+AY2
   70 CONTINUE
      RETURN
      END
