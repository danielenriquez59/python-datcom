      SUBROUTINE COORD1
C
C***** CALCULATE NACA 1-SERIES AIRFOIL COORDINATES
C
      COMMON /IWING/ PW, X(60)
      COMMON / IHT / PHT, XU(60),XL(60)
      COMMON / IVT / PVT, YUN(60),YLN(60)
      COMMON / IBW / PBW,L,I,J,K,II,JJ,KK,III,JJJ
      COMMON / IBH / PBH, THN(60),CAM(60)
      COMMON /IBWH/  PBWH,ALPHAI,ALPHAO,AII
      COMMON /IBWHV/ PBWHV, RHO,T
      COMMON /CONSNT/PI,DEG,UNUSED,RAD
      DIMENSION A(2,3),B(2),XN(2)
      DIMENSION AA(3,4),BA(3),XM(3)
C
      AI=I
      AJ=J
      AK=K
      AII=II
      AJJ=JJ
      AKK=KK
      AIII=III
      AJJJ=JJJ
      ZT=AJ*.1-.1
      T=AJJ*.1+AKK*.01+AIII*.001+AJJJ*.0001
      NN=2
      IF(J.EQ.6)ZT=AJ*.1-.2
      D0=0.0
      IF(J.EQ.6)D1=2.157*T
      IF(J.EQ.6)SM=4.
      IF(J.EQ.8)D1=3.6833*T
      IF(J.EQ.8)SM=3.
      IF(J.EQ.9)D1=5.5283*T
      IF(J.EQ.9)SM=3.
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
      A0=SQRT(2.*1.1019*((T*SM/6.)**2))
      RHO=.5*A0**2
      AA(1,1)=0.0
      AA(1,2)=2.
      AA(1,3)=6.*ZT
      BA(1)=2.*D2+6.*D3*(1.-ZT)+A0/(4.*ZT**1.5)
      AA(2,1)=1.
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
      LL=L-1
      AII=AII*.1
      DO 1000 M=2,LL
        YC=-(AII/(4.*PI))*((1.-X(M))*ALOG(1.-X(M))+X(M)*ALOG(X(M)))
        ALPHA=ATAN((-AII/(4.*PI))*(ALOG(X(M))-ALOG(1.-X(M))))
        IF(X(M).EQ.ZT)YT=T/2.
        IF(X(M).LT.ZT)YT=A0*X(M)**.5+A1*X(M)+A2*X(M)**2+A3*X(M)**3
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
      XL(1)=0.0
      YLN(1)=0.0
      XL(1)=0.0
      YLN(1)=0.0
      ALPHAI=0.0
      ALPHAO=-RAD*AII/(2.*PI)
      RETURN
      END
