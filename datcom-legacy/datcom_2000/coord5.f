      SUBROUTINE COORD5
C
C***** CALCULATE NACA FIVE DIGIT AIRFOIL COORDINATES
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
      AIII=III
      T=AII*.1+AJJ*.01+AKK*.001+AIII*.0001
      RHO=1.1019*T**2
      ZP=AJ*.1/2.
      A=6.*ZP-3.
      B=-2.+6.*ZP-3.*ZP**2
      G=B*B/4.+A*A*A/27.
      IF(G.LT.0.) GO TO 1000
      D=(-B/2.+G**.5)**(1./3.)
      E=(-B/2.-G**.5)**(1./3.)
      ZM=D+E+1.
      GO TO 1010
 1000 PHI=ARCCOS((-B/2.)/((-A**3/27.)**.5))
      ZM=1.+2.*((-A/3.)**.5)*COS(PHI/3.+4.18879)
 1010 XK=(6.*AI*.01)/(ZP**3-3.*ZM*ZP**2+ZM**2*(3.-ZM)*ZP)
      DO 1040 M=1,L
         IF(AK.NE.0.)GO TO 1020
         YT=5.*T*(.2969*SQRT(X(M))-.126*X(M)-.3516*X(M)**2+.2843*X(M)**3
     1      -.1015*X(M)**4)
         IF(X(M).LT.ZM)YC=(1./6.)*XK*(X(M)**3-3.*ZM*X(M)**2+ZM**2*(3.
     1      -ZM)*X(M))
         IF(X(M).LT.ZM)ALPHA=ATAN((1./6.)*XK*(3.*X(M)**2-6.*ZM*X(M)+
     1     ZM**2*(3.-ZM)))
         IF(X(M).EQ.ZP)YC=AI*.01
         IF(X(M).EQ.ZP)ALPHA=0.0
         IF(X(M).GT.ZM)YC=(1./6.)*XK*ZM**3*(1.-X(M))
         IF(X(M).GT.ZM)ALPHA=ATAN(-(1./6.)*XK*ZM**3)
         IF(AK.EQ.0.0)GO TO 1030
 1020    RK=(3.*((ZM-ZP)**2)-ZM**3)/((1.-ZM)**3)
         ZMX=AI*.01
         XK=(6.*ZMX)/((ZP-ZM)**3-RK* (1.-ZM)**3*ZP-(ZM**3)*ZP+ZM**3)
         YT=5.*T*(.2969*SQRT(X(M))-.126*X(M)-.3516*X(M)**2+.2843*X(M)**3
     1      -.1015*X(M)**4)
         IF(X(M).LT.ZM)YC=(1./6.)*XK*((X(M)-ZM)**3-RK*X(M)*(1.-ZM)**3
     1      -ZM**3*X(M)+ZM**3)
         IF(X(M).LT.ZM)ALPHA=ATAN((1./6.)*XK*(3.*(X(M)-ZM)**2-RK*(1.
     1      -ZM)**3-ZM**3))
         IF(X(M).EQ.ZP)YC=AI*.01
         IF(X(M).EQ.ZP)ALPHA=0.0
         IF(X(M).GT.ZM)YC=(1./6.)*XK*(RK*(X(M)-ZM)**3-RK*X(M)*(1.-ZM)**3
     1      -X(M)*ZM**3+ZM**3)
         IF(X(M).GT.ZM)ALPHA=ATAN((1./6.)*XK*(3.*RK*(X(M)-ZM)**2-RK*
     1      (1.-ZM)**3-ZM**3))
 1030    XU(M)=X(M)-YT*SIN(ALPHA)
         YUN(M)=YC+YT*COS(ALPHA)
         XL(M)=X(M)+YT*SIN(ALPHA)
         YLN(M)=YC-YT*COS(ALPHA)
         CAM(M)=YC
         IF(CAM(M) .LT. 1.E-05) CAM(M)=0.0
         THN(M)=YT
 1040 CONTINUE
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
