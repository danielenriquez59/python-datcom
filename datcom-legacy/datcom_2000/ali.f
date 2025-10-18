      FUNCTION ALI(Q,R,XST,XRT,XSLT)
C
C***  CPMPUTES VORTEX INTERFERENCE FACTORS
C
      DIMENSION AL(4)
      FF=R
      F=R
      H=Q
      DO 1030 MK=1,4
         AB=ALOG((H**2+(F-XST)**2)/(H**2+(F-XRT)**2))
C
C     ----MUST TEST FOR H=0.
C
         IF(H.NE.0.0) GO TO 1000
         BB=XST-XRT
         GO TO 1010
 1000    FST=(F-XST)/H
         FRT=(F-XRT)/H
         BB=(FST/ABS(FST))*H*ATAN(ABS(FST))
         BB=BB-(FRT/ABS(FRT))*H*ATAN(ABS(FRT))+(XST-XRT)
 1010    AL(MK)=((XST-XRT*XSLT)-F*(1.0-XSLT))/(2.*(XST-XRT))
         AL(MK)=AL(MK)*AB-((1.0-XSLT)/(XST-XRT))*BB
         F=-F
         IF(MK-2)1030,1020,1030
 1020    RT2=XRT*XRT
         H2=H*H
         F=(F*RT2)/((F**2)+H2)
         H=(H*RT2)/((FF**2)+H2)
 1030 CONTINUE
      ALI=2.0/(1.+XSLT)*(AL(1)-AL(2)-AL(3)+AL(4))
C
      RETURN
      END
