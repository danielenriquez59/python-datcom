      SUBROUTINE BDAREA(ABORT)
C
C***  EXECUTIVE FOR BODY PARTS SHADOWED BY MACH LINES
C
      LOGICAL ABORT
      DIMENSION A(3),X(4),Y(4),XL(42),YL(42),XLI(42),YLI(42),XT(42),
     1YT(42),XTI(42),YTI(42)
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG
      COMMON /FLGTCD/ FLC(73)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /HTI/    HTIN(154)
      COMMON /HTDATA/ AHT(195)
      COMMON /BODYI/  XNX,XB(20),SB(20),PB(20),RB(20)
      COMMON /SYNTSS/ SYNA(10)
C
      NX=XNX
      ILI=0
      ILIE=0
      ITI=0
      ITIE=0
      SAXB=0
      SAYB=0
      AREA=0.
      AREAB=0.
 1000 AMUU=ATAN(1./SQRT(FLC(I+2)**2-1.))
      DO 1050 K=1,2
         IF(K.NE.1) GO TO 1020
 1010    A(1)=SYNA(6)+(HTIN(4)-HTIN(3))*AHT(62)*COS(SYNA(8)/RAD)
         A(2)=SYNA(7)-(A(1)-SYNA(6))*SIN(SYNA(8)/RAD)/COS(SYNA(8)/RAD)
         A(3)=SYNA(8)
         GO TO 1030
 1020    A(1)=A(1)+AHT(10)*COS(SYNA(8)/RAD)
         A(2)=A(2)-AHT(10)*SIN(SYNA(8)/RAD)
         GO TO 1040
 1030    CALL PTINT2(XB,RB,A,AMUU,XL,YL,XLI,YLI,NX,ILUI,ILLI,ILUIE,
     1               ILLIE,ABORT)
         IF(ABORT) RETURN
         GO TO 1050
 1040    CALL PTINT2(XB,RB,A,AMUU,XT,YT,XTI,YTI,NX,ITUI,ITLI,ITUIE,
     1              ITLIE,ABORT)
         IF(ABORT) RETURN
 1050 CONTINUE
      DO 1210 K=1,2
         SAXB=0.
         SAYB=0.
         AREA=0.
         AREAB=0.
C
C     CHECK FOR UPPER OR LOWER AREA
C
         IF(K.NE.1) GO TO 1060
         ILI=ILUI
         ILIE=ILUIE
         ITI=ITUI
         ITIE=ITUIE
         GO TO 1070
 1060    ILI=ILLI
         ILIE=ILLIE
         ITI=ITLI
         ITIE=ITLIE
         NX=NX+21
 1070    IF(ITIE.NE.0) GO TO 1080
C
C**** NEITHER MACH LINE INTERSECTS BODY PROFILE EXTENSION ****
C
         KIND=1
         GO TO 1100
 1080    IF(ILIE.NE.0) GO TO 1090
C
C**** ONLY THE TRAILING EDGE MACH LINE INTERSECTS BODY PROFILE EXTENSION
C
         KIND=2
         GO TO 1100
C
C**** BOTH MACH LINES INTERSECT BODY PROFILE EXTENSION ****
C
 1090    KIND=3
 1100    GO TO (1110,1150,1180) KIND
 1110    IF(ILI.NE.ITI) GO TO 1120
         X(1)=0.
         X(2)=XLI(ILI)
         X(3)=XTI(ITI)+AHT(10)
         X(4)=AHT(10)
         Y(1)=0.
         Y(2)=YLI(ILI)
         Y(3)=YTI(ITI)
         Y(4)=0.
         CALL AREA2(X,Y,SAREA,2,AX,AY)
         AREA=SAREA
         AREAB=AREA
         SAXB=AX
         SAYB=AY
         GO TO 1190
 1120    II=1
         X(1)=0.
         X(2)=XLI(ILI)
         X(3)=XL(ILI+II)
         X(4)=AHT(10)
         Y(1)=0.
         Y(2)=YLI(ILI)
         Y(3)=YL(ILI+II)
         Y(4)=0.
         CALL AREA2(X,Y,SAREA,2,AX,AY)
         AREA=AREA+SAREA
         SAXB=SAXB+AX
         SAYB=SAYB+AY
 1130    II=II+1
         IF(XL(ILI+II).GT.(XTI(ITI)+AHT(10))) GO TO 1140
         X(1)=XL(ILI+II-1)
         X(2)=XL(ILI+II)
         X(3)=AHT(10)
         Y(1)=YL(ILI+II-1)
         Y(2)=YL(ILI+II)
         Y(3)=0.
         CALL AREA2(X,Y,SAREA,3,AX,AY)
         AREA=AREA+SAREA
         SAXB=SAXB+AX
         SAYB=SAYB+AY
         GO TO 1130
 1140    X(1)=XL(ILI+II-1)
         X(2)=XTI(ITI)+AHT(10)
         X(3)=AHT(10)
         Y(1)=YL(ILI+II-1)
         Y(2)=YTI(ITI)
         Y(3)=0.
         CALL AREA2(X,Y,SAREA,3,AX,AY)
         AREA=AREA+SAREA
         AREAB=AREA
         SAXB=SAXB+AX
         SAYB=SAYB+AY
         GO TO 1190
 1150    II=1
         X(1)=0.
         X(2)=XLI(ILI)
         X(3)=XL(ILI+II)
         X(4)=AHT(10)
         Y(1)=0.
         Y(2)=YLI(ILI)
         Y(3)=YL(ILI+II)
         Y(4)=0.
         CALL AREA2(X,Y,SAREA,2,AX,AY)
         AREA=AREA+SAREA
         SAXB=SAXB+AX
         SAYB=SAYB+AY
 1160    IF((ILI+II).EQ.NX) GO TO 1170
         X(1)=XL(ILI+II)
         X(2)=XL(ILI+II+1)
         X(3)=AHT(10)
         Y(1)=YL(ILI+II)
         Y(2)=YL(ILI+II+1)
         Y(3)=0.
         CALL AREA2(X,Y,SAREA,3,AX,AY)
         AREA=AREA+SAREA
         SAXB=SAXB+AX
         SAYB=SAYB+AY
         II=II+1
         GO TO 1160
 1170    X(1)=XL(ILI+II)
         X(2)=XTI(ITIE)+AHT(10)
         X(3)=AHT(10)
         Y(1)=YL(ILI+II)
         Y(2)=YTI(ITIE)
         Y(3)=0.
         CALL AREA2(X,Y,SAREA,3,AX,AY)
         AREA=AREA+SAREA
         SAXB=SAXB+AX
         SAYB=SAYB+AY
         X(3)=XTI(ITI)+AHT(10)
         Y(3)=YTI(ITI)
         CALL AREA2(X,Y,EXTARE,3,AX,AY)
         AREAB=AREA-EXTARE
         SAXB=SAXB-AX
         SAYB=SAYB-AY
         GO TO 1190
 1180    X(1)=0.
         X(2)=XLI(ILIE)
         X(3)=XTI(ITIE)+AHT(10)
         Y(1)=0.
         Y(2)=YLI(ILIE)
         Y(3)=YTI(ITIE)
         X(4)=AHT(10)
         Y(4)=0.
         CALL AREA2(X,Y,SAREA,2,AX,AY)
         AREA=AREA+SAREA
         X(2)=XLI(ILI)
         X(3)=XTI(ITI)+AHT(10)
         Y(2)=YLI(ILI)
         Y(3)=YTI(ITI)
         CALL AREA2(X,Y,SAREA,2,AX,AY)
         AREAB=AREAB+SAREA
         SAXB=SAXB+AX
         SAYB=SAYB+AY
 1190    IF(K.EQ.2) GO TO 1200
         AREAU=AREA
         AREAUB=AREAB
         SAXUB=SAXB
         SAYUB=SAYB
         GO TO 1210
 1200    AREAL=AREA
         AREALB=AREAB
         SAXLB=SAXB
         SAYLB=SAYB
 1210 CONTINUE
      HTIN(114+I)=AREAUB+AREALB
      HTIN(134+I)=AREAU+AREAL
      AXTB=SAXUB+SAXLB
      AYTB=SAYUB+SAYLB
      XCB=AXTB/HTIN(114+I)
      YCB=AYTB/HTIN(114+I)
      XBARB=XCB*COS(-SYNA(8)/RAD)-YCB*SIN(-SYNA(8)/RAD)+SYNA(6)
     1      +(HTIN(4)-HTIN(3))*AHT(62)*COS(SYNA(8)/RAD)
      HTIN(94+I)=XBARB-SYNA(1)
      RETURN
      END
