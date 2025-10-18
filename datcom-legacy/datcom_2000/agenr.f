      SUBROUTINE AGENR(BOAK,SB,A)
C
C***  GENERATES COEFFICIENTS FOR G/DELTA CALCULATIONS BY GDELTA
C
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      DIMENSION F1(4),F2(4),F3(4),F4(4),BOCO16(4),BOAK(4)
      DIMENSION A(16),AL(16),C1(16),C2(16),C3(4)
      DATA C1 /-.07612,.00001  ,.21677,.54120  , -.29289,-.21677,.00001
     1     ,.32443 , -.61732,-.54120,-.32443,.00001  ,
     2     -1.0,-.92388,-.70711,-.38267/
      DATA C2/1.92388,1.84776,1.63099,1.30656 , 1.70711,1.63099,1.41422,
     1        1.08979 , 1.38268,1.30655,1.08979,.76536 ,
     2         1.0,.92388,.70711,.38267/
      DATA F1 /5.2262,1.0360,0.,.11208/
      DATA F2 /1.91433,2.8284,.91418,0.0/
      DATA F3 /0.0,1.1944,2.1634,1.5772/
      DATA F4 /.14645,0.0,.85357,2.0/
      DATA C3 /.92388,.70711,.38268,.00001 /
C
      K=1
      KP=0
      TANSB=TAN(DEG*SB)
      DO 1010 I=1,16
         IF(KP.EQ.K)GO TO 1000
         KP=K
         BOCH=BOAK(K)
         BOCO16(K)=-BOCH/16.
         RCPLBC=1.0/BOCH
         BOC2=BOCH**2
         BOCTN=BOCH  *TANSB
         TBOCTN=2.0*BOCTN
         DENOM=1.+C3(K)*TBOCTN
 1000    TMP=(1.+C1(I)*BOCTN)**2
         AL(I)=RCPLBC/C1(I)*(SQRT(TMP+BOC2*(C1(I)**2))-1.)
     1         -RCPLBC/C2(I)*(SQRT(TMP+BOC2*(C2(I)**2))/DENOM-1.0)
     2         -2.*TANSB*SQRT((1.+C3(K)*BOCTN)**2+BOC2*(C3(K)**2))/DENOM
         IF(4*(I/4).EQ.I)K=K+1
 1010 CONTINUE
C
      K=1
      L=1
      DO 1020 I=1,4
         FST=-2.*F1(I)
         IF(I.EQ.1)FST=-FST
         A(L)=FST+BOCO16(I)*(2.6131*AL(K)+2.*(-.70711*AL(K+1)-.76537*
     1        AL(K+2)+.20711*AL(K+3)))
         L=L+1
         FST=-2.*F2(I)
         IF(I.EQ.2)FST=-FST
         A(L)=FST+BOCO16(I)*(-1.4142*AL(K)+2.*(1.8478*AL(K+1)-.50000*
     1   AL(K+2)-.76537*AL(K+3)))
         L=L+1
         FST=-2.*F3(I)
         IF(I.EQ.3)FST=-FST
         A(L)=FST+BOCO16(I)*(1.0824*AL(K)+2.*(-1.2071*AL(K+1)+1.8478
     1       *AL(K+2)-.70711*AL(K+3)))
         L=L+1
         FST=-2.*F4(I)
         IF(I.EQ.4)FST=-FST
         A(L)=FST+BOCO16(I)*(-.5*AL(K)+1.0824*AL(K+1)-1.4142*AL(K+2)
     1        +2.6131*AL(K+3))
         L=L+1
 1020    K=L
      RETURN
      END
