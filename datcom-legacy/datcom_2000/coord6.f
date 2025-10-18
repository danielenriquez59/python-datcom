      SUBROUTINE COORD6
C
C***** CALCULATE NACA 6-SERIES AIRFOIL COORDINATES
C
      COMMON /IWING/ PW, X(60)
      COMMON / IHT / PHT, XU(60),XL(60)
      COMMON / IVT / PVT, YUN(60),YLN(60)
      COMMON / IBW / PBW,L,I,J,K,II,JJ,KK,III,JJJ
      COMMON / IBH / PBH, THN(60),CAM(60)
      COMMON /IBWH/  PBWH,ALPHAI,ALPHAO,AJJ
      COMMON /IBWHV/ PBWHV, RHO,T
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
      AJJJ=JJJ
      T=AKK*.1+AIII*.01
      NN=2
      IF(J.EQ.3)GO TO 1000
      IF(J.EQ.4)GO TO 1010
      IF(J.EQ.5)GO TO 1020
      ZT=.45
      SM=-1.268
      R=SM*(T-.06)+.873
      D1=R*T
      GO TO 1030
 1000 ZT=.35
      SM=-.6116
      R=SM*(T-.06)+.46
      D1=R*T
      IF(II.GT.0)D1=T
      GO TO 1030
 1010 ZT=.40
      SM=-.6888
      R=SM*(T-.06)+.523
      D1=R*T
      IF(II.GT.0)D1=1.04*T
      GO TO 1030
 1020 ZT=.4
      SM=-.8833
      R=SM*(T-.06)+.65
      D1=R*T
      IF(II.GT.0)D1=1.17*T
 1030 D0=0.0
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
      RLE=.01*(68.682*T**2+.0182*T+.0014)
      RHO=RLE
      A0=SQRT(2.*RLE)
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
      ZA=AJJJ*.1
      IF(AJJJ.LT.1.)ZA=1.
      IF(ZA.EQ.1.)GO TO 1040
      C=1.-ZA
      G=(-1./C)*((ZA**2)*(.5*ALOG(ZA)-.25)+.25)
      H=(1./C)*((.5*C**2)*ALOG(C)-.25*C**2)+G
 1040 AJJ=AJJ*.1
      LL=L-1
      DO 1100 M=2,LL
         IF(ZA.EQ.1.)GO TO 1050
         S=1.-X(M)
         D=ZA-X(M)
      IF(D.EQ.0.0) D=1.0E-10
         YC=(AJJ/(2.*PI*(ZA+1.)))*((1./C)*((.5*D**2)*ALOG(ABS(D))-
     1    (.5*S**2)*ALOG(S)+.25*S**2-.25*D**2)-X(M)*ALOG(X(M))+G-X(M)*H)
         ALPHA=ATAN((AJJ/(2.*PI*(1.+ZA)))*((1./C)*(-D*ALOG(ABS(D))
     1         +S*ALOG(S))-ALOG(X(M))-1.-H))
         GO TO 1060
 1050    YC=-(AJJ/(4.*PI))*((1.-X(M))*ALOG(1.-X(M))+X(M)*ALOG(X(M)))
         ALPHA=ATAN((-AJJ/(4.*PI))*(ALOG(X(M))-ALOG(1.-X(M))))
 1060    IF(X(M).EQ.ZT)YT=T/2.
         IF(X(M).LT.ZT)YT=A0*X(M)**.5+A1*X(M)+A2*X(M)**2+A3*X(M)**3
        IF(X(M).GT.ZT)YT=D0+D1*(1.-X(M))+D2*(1.-X(M))**2+D3*(1.-X(M))**3
         XU(M)=X(M)-YT*SIN(ALPHA)
         YUN(M)=YC+YT*COS(ALPHA)
         XL(M)=X(M)+YT*SIN(ALPHA)
         YLN(M)=YC-YT*COS(ALPHA)
         IF(XU(M).GE..80.AND.II.GT.0)GO TO 1080
         NO=1
 1070    CONTINUE
         CAM(M)=YC
         IF(CAM(M) .LT. 1.E-05) CAM(M)=0.0
         THN(M)=YT
         GO TO 1100
 1080    IF(NO.NE.1)GO TO 1090
         SXU=XU(M)
         SXL=XL(M)
         SYU=YUN(M)
         SYL=YLN(M)
         SMU=-SYU/(1.-SXU)
         SML=-SYL/(1.-SXL)
         NO=2
         GO TO 1070
 1090    XU(M)=XU(M)-SXU
         XL(M)=XL(M)-SXL
         YUN(M)=SMU*XU(M)+SYU
         YLN(M)=SML*XL(M)+SYL
         XU(M)=XU(M)+SXU
         XL(M)=XL(M)+SXL
         GO TO 1070
 1100 CONTINUE
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
      XL(1)=0.0
      YUN(1)=0.0
      ALPHAI=0.0
      ALPHAO=-RAD*AJJ/(2.*PI)
      IF(ZA .EQ. 1.) RETURN
      ALPHAI=ALPHAO*H/(ZA+1.)
      ALPHAO=ALPHAO+ALPHAI
      RETURN
      END
