      SUBROUTINE CORD4M
C
C***** CALCULATE NACA FOUR DIGIT (MODIFIED) AIRFOIL COORDINATES
C
      COMMON /IWING/ PW, X(60)
      COMMON / IHT / PHT, XU(60),XL(60)
      COMMON / IVT / PVT, YUN(60),YLN(60)
      COMMON / IBW / PBW,L,I,J,K,II,JJ,KK,III,JJJ
      COMMON / IBH / PBH, THN(60),CAM(60)
      COMMON /IBWHV/ PBWHV, RHO,T,DELTAY,XOVC,TOVC,ZM,ZP
      COMMON /CONSNT/PI,DEG,UNUSED,RAD
      DIMENSION A(2,3),B(2),XN(2)
      DIMENSION AA(3,4),BA(3),XM(3)
      AI=I
      AJ=J
      AK=K
      AII=II
      AJJ=JJ
      AKK=KK
      AIII=III
      ZM=AI*.01
      ZP=AJ*.1
      ZT=AIII*.1
      T=AK*.1+AII*.01
      NN=2
      IF(ZT.GT..18.AND.ZT.LT..22) D1=T
      IF(ZT.GT..28.AND.ZT.LT..32) D1=1.17*T
      IF(ZT.GT..38.AND.ZT.LT..42) D1=1.575*T
      IF(ZT.GT..48.AND.ZT.LT..52) D1=2.325*T
      IF(ZT.GT..58.AND.ZT.LT..62) D1=3.5*T
      D0=.01*T
      A(1,1)=-2.*(1.-ZT)
      A(1,2)=-3.*((1.-ZT)**2)
      A(2,1)=(1.-ZT)**2
      A(2,2)=(1.-ZT)**3
      B(1)=D1
      B(2)=T/2.-D0-D1*(1.-ZT)
      LL=NN+1
      CALL SLEQ(A,XN,B,NN,LL)
      D2=XN(1)
      D3=XN(2)
      A0=SQRT(2.*1.1019*((T*AKK/6.)**2))
      RHO=.5*A0**2
      AA(1,1)=0.0
      AA(1,2)=2.0
      AA(1,3)=6.*ZT
      BA(1)=2.*D2+6.*D3*(1.-ZT)+A0/(4.*ZT**1.5)
      AA(2,1)=1.0
      AA(2,2)=2.*ZT
      AA(2,3)=3.*ZT**2
      BA(2)=-A0/(2.*ZT**.5)
      AA(3,1)=ZT
      AA(3,2)=ZT**2
      AA(3,3)=ZT**3
      BA(3)=-A0*ZT**.5+T/2.
      NN=3
      LL=NN+1
      CALL SLEQ(AA,XM,BA,NN,LL)
      A1=XM(1)
      A2=XM(2)
      A3=XM(3)
      DO 1000 M=1,L
        IF(X(M).EQ.ZT)YT=T/2.
        IF(X(M).EQ.ZP)YC=ZM
        IF(X(M).EQ.ZP)ALPHA=0.0
        IF(X(M).LT.ZP)YC=(ZM/ZP**2)*(2.*ZP*X(M)-X(M)**2)
        IF(X(M).LT.ZP)ALPHA=ATAN((2.*ZM/ZP**2)*(ZP-X(M)))
        IF(X(M).LT.ZT)YT=A0*X(M)**.5+A1*X(M)+A2*X(M)**2+A3*X(M)**3
        IF(X(M).GT.ZP)YC=(ZM/((1.-ZP)**2))*(1.-2.*ZP+2.*ZP*X(M)-X(M)**2)
        IF(X(M).GT.ZP)ALPHA=ATAN((2.*ZM/((1.-ZP)**2))*(ZP-X(M)))
        IF(X(M).GT.ZT)YT=D0+D1*(1.-X(M))+D2*(1.-X(M))**2+D3*(1.-X(M))**3
        XU(M)=X(M)-YT*SIN(ALPHA)
        YUN(M)=YC+YT*COS(ALPHA)
        XL(M)=X(M)+YT*SIN(ALPHA)
        YLN(M)=YC-YT*COS(ALPHA)
        CAM(M)=YC
        IF(CAM(M) .LT. 1.E-05) CAM(M)=0.0
        THN(M)=YT
 1000 CONTINUE
      THN(1)=0.0
      THN(L)=0.0
      CAM(1)=0.0
      CAM(L)=0.0
      XU(L)=1.
      YUN(L)=0.0
      XL(L)=1.
      YLN(L)=0.0
      XU(1)=0.0
      YLN(1)=0.0
      RETURN
      END
