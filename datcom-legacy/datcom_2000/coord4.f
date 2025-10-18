      SUBROUTINE COORD4
C
C***** CALCULATE NACA FOUR DIGIT AIRFOIL COORDINATES
C
      COMMON /IWING/ PW, X(60)
      COMMON / IHT / PHT, XU(60),XL(60)
      COMMON / IVT / PVT, YUN(60),YLN(60)
      COMMON / IBW / PBW,L,I,J,K,II,JJ,KK,III,JJJ
      COMMON / IBH / PBH, THN(60),CAM(60)
      COMMON /IBWHV/ PBWHV, RHO,T,DELTAY,XOVC,TOVC,ZM,ZP
      AI=I
      AJ=J
      AK=K
      AII=II
      AJJ=JJ
      AKK=KK
      ZM=AI*.01
      ZP=AJ*.1
      T=AK*.1+AII*.01+AJJ*.001+AKK*.0001
      RHO=1.1019*T**2
      DO 1000 M=1,L
         YT=5.*T*(.2969*SQRT(X(M))-.126*X(M)-.3516*X(M)**2+.2843*X(M)**3
     1      -.1015*X(M)**4)
         IF(X(M).EQ.ZP)YC=ZM
         IF(X(M).EQ.ZP)ALPHA=0.0
         IF(X(M).LT.ZP)YC=(2.*ZP*X(M)-X(M)**2)*ZM/ZP**2
         IF(X(M).LT.ZP)ALPHA=ATAN((2.*ZM/(ZP**2))*(ZP-X(M)))
         IF(X(M).GT.ZP)
     1   YC=(ZM/((1.-ZP)**2))*(1.-2.*ZP+2.*ZP*X(M)-X(M)**2)
         IF(X(M).GT.ZP)ALPHA=ATAN((2.*ZM/((1.-ZP)**2))*((ZP-X(M))))
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
