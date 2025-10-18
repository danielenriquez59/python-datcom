      SUBROUTINE BODOPT
C
C***  COMPUTES ASYMMETRICAL BODY AERODYNAMICS
C
      COMMON /WINGD/  A(195), B(48)
      COMMON /BDATA/  BD(762)
      COMMON /IBODY/  PBODY, BODY(400)
      COMMON /BODYI/  XNX,X(20),S(20),P(20),R(20),ZU(20),ZL(20),
     1                BNOSE,BTAIL,BLN,BLA,DS
      COMMON /WINGI/  CHRDTP,SSPNOP,SSPNE,SSPN,CHRDBP,CHRDR,SAVSI,SAVSO,
     1                CHSTAT,SWAFP,TWISTA,SSPNDD,DHDADI,DHDADO,
     2                TYPE,
     3                TOVC,DELTAY,XOVC,CLI,ALPHAI,CLALPA(20),
     4                CLMAX(20),CMO,LERI,LERO,CAMBER,TOVCO,XOVCO,CMOT,
     5                CLMAXL,CLAMO,TCEFF,KSHARP,SLOPE(6)
      COMMON /OVERLY/ NLOG,NMACH,M,NALPHA,IG
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4                VFPL,VFSC,CTAB
      COMMON /FLGTCD/ FLC(73)
      COMMON /OPTION/ SREF,LREF,ROUGFC
      COMMON /SYNTSS/ XCG,XW,ZW,ALIW,ZCG,XH,ZH,ALIH,XV,VERTUP,HINAX,
     1                XVF,SCALE,ZV,ZVF,YV,YF,PHIV,PHIF
      COMMON /CONSNT/ CONST(4)
C
      EQUIVALENCE (PI,CONST(1)), (RAD,CONST(4))
      EQUIVALENCE (UNUSED,CONST(3))
      LOGICAL PBODY
      LOGICAL VERTUP
      LOGICAL FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1        HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2        TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3        HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4        VFPL,VFSC,CTAB
      REAL LREF
C
      DIMENSION ROUTID(2),Q15127(3),Q95MX(3),QKPARB(3),QKVARB(3)
     1,QALFR(3)
      DIMENSION ALPH(20),CMP(20)
      DIMENSION D(20)
      DIMENSION XE(20),XOLE(20),ZPOLE(20),SE(20),DSEDX(20),YE(20),PE(20)
      DIMENSION CY95(6)
      DIMENSION XBA1(9),YBA1(9),CBA1(6)
      DIMENSION XBA2(9),YBA2(9),CBA2(6)
      DIMENSION X1BA3(5),X2BA3(12),YBA3(60)
      DIMENSION TX95(20),TY95(20)
      DIMENSION XE95(20),YE95(20),PE95(20),SE95(20),DYEDX(20)
      DIMENSION RP(20),RS(20)
      DIMENSION DELA(20)
      DIMENSION X27M(4),X27I(4),C27(6)
C
      DATA ROUTID /4HBODO,4HPT  /,   Q15127 /4H4.1.,4H5.1-,4H27  /,
     1 Q95MX/4HFIND,4H .95,4HYMAX/,  QKPARB /4H4.2.,4H1.2-,4H36A /,
     2 QKVARB/4H4.2.,4H1.2-,4H36B /, QALFR/4H4.2.,4H1.2-,4H37  /
      DATA X27I /1.57780,1.67221,1.98509,2.28874/
      DATA X27M /0.0,1.0,2.0,3.0/
      DATA I27  /0/
      DATA RP   /20*1./,  RS/20*1./
C
C---------     FIGURE  BA1   ( KP VS ASPECT RATIO )         ------------
C
      DATA XBA1 /0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0/
      DATA YBA1 /0.0,.72,1.33,1.83,2.225,2.575,2.88,3.125,3.37/
C
C---------     FIGURE BA2    ( KV VS ASPECT RATIO )         ------------
C
      DATA XBA2 /0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0/
      DATA YBA2 /3.14,3.14,3.15,3.19,3.21,3.626,3.31,3.9,3.47/
C
C---------  FIGURE 4.2.1.2-37
C
      DATA X1BA3 /1.,2.,3.,4.,100./
      DATA X2BA3 /0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,100./
      DATA YBA3  /34.85,31.20,27.55,23.90,20.25,16.60,12.95,11.25,
     1     10.30,9.75,9.40,0.,  30.65,25.90,21.15,16.40,11.65,9.10,
     2     7.60,6.70,6.10,5.70,5.50,0.,  26.70,20.00,13.30,9.10,
     3     7.50,4.90,4.00,3.50,3.10,2.90,2.75,0.,  22.70,16.00,9.30,
     4     7.20,5.00,3.95,3.00,2.40,2.05,1.85,1.70,0.,  12*0./
C
      IY95=0
      IBA1=0
      IBA2=0
C
C     ------CALCULATE CD0 -------------
C
      NX=XNX+.5
      BD(1)=X(NX)
      BD(57)=S(NX)
      CALL GETMAX(X,S,NX,BD(2),BD(3),IMAX)
      BD(56)=BD(3)
C
C*****ADJUST BASE AREA IF BODY IS BOAT TAILED ********
C
      IF(BD(57).LE.0.3*BD(56)) BD(57)=0.30*BD(56)
      BD(55)=12.0*BD(1)/ROUGFC
C
C     ------FIGURE 4.1.5.1-27-----------
C
      CALL TBFUNX(B(1),CEPT,DYDX,4,X27M,X27I,C27,I27,MI,NG,
     1           0,0,Q15127,3,ROUTID)
      BD(91)=BD(55)**1.0482*10.0**CEPT
      BD(90)=FLC(M+42)*BD(1)
      RUSE=BD(90)
      IF(BD(91).LT.BD(90))RUSE=BD(91)
C
C     ------FIGURE 4.1.5.1-26-----------
C
      CALL FIG26(RUSE,B(1),BD(92))
C
C     ------CALCULATE SS=PIN -----------
C
      CALL EQSPCE(X,R,P,S,NX,XE,YE,PE,SE,20,DSEDX)
      CALL TRAPZ (PE,20,XE,PIN,1)
      BD(93)=PIN
      BD(86)=SQRT(BD(57)*4./PI)
      BD(85)=SQRT(BD(56)*4./PI)
      BD(75)=BD(1)/BD(85)
      BD(59)=BD(92)*(1.+60./(BD(75)**3)+.0025*BD(75))*BD(93)/BD(56)
      BD(60)=.029*((BD(86)/BD(85))**3/SQRT(BD(59))*BD(56)/SREF)
      BD(59)=BD(59)*BD(56)/SREF
      BD(61)=BD(59)+BD(60)
      BD(10)=BD(61)
C
C---------     COMPUTE SLOPE OF LINE CONNECTING END MIDPOINTS    -------
C
      N=1
      Z0= 0.5*(ZU(N)+ZL(N))
      N=NX
      ZNX= 0.5*(ZU(N)+ZL(N))
      DZDX= (ZNX-Z0)/(X(N)-X(1))
 1020 TEMPD=ATAN(DZDX)*RAD
C
      DO 1050 N=1,NX
C
C---------     COMPUTE BODY DEPTH ARRAY, D                       -------
C
         D(N)= ZU(N)-ZL(N)
         IF(D(N).EQ.0. .AND. N.NE.1)D(N)=UNUSED
C
C---------     COMPUTE LOCUS OF BODY MIDPOINTS (Z-PRIME ARRAY ) --------
C
 1030    BD(N+740)=ZL(N)+ 0.5*D(N) - (X(N)*DZDX + Z0 )
C
C---------     COMPUTE Z-PRIME OVER REFERENCE LENGTH (ZPOL-ARRAY ) -----
C
 1040    BD(N+720)=BD(N+740)/X(NX)
C
C---------     COMPUTE (X)N OVER REFERENCE LENGTH   ( XOL-ARRAY ) ------
C
 1050 BD(N+700)=X(N)/X(NX)
C
C---------     COMPUTE (X)CG/L,  XOVLCG
C
      XOVLCG= XCG/X(NX)
C
C---------     EQUAL SPACE X,XOL,ZPOL, AND S                ------------
C
      CALL EQSPCE (X,BD(701),BD(721),S,NX,XE,XOLE,ZPOLE,SE,
     1             20,DSEDX)
C
C---------     COMPUTE NEW ARRAYS OF VALUES TO BE INTEGRATED     -------
C
      DO 1080 L=1,20
         IF(L.GT.1) GO TO 1060
         ALPH(L)=0.0
         CMP(L)=0.0
         GO TO 1080
 1060    CONTINUE
         IF(L.LT.20) GO TO 1070
         ALPH(L)=0.0
         CMP(L)=0.0
         GO TO 1080
 1070    CONTINUE
         ALPH(L) = ZPOLE(L)/((1.0-XOLE(L))*SQRT(XOLE(L)-XOLE(L)**2))
         CMP(L)= ZPOLE(L)*((1.-2.*XOLE(L))/SQRT(XOLE(L)-XOLE(L)**2))
 1080 CONTINUE
C
C---------     COMPUTE ALPHA-ZERO                      -----------------
C
      CALL TRAPZ(ALPH,19,XOLE,ALZERO,1)
 1090 ALPHZR=-ALZERO/PI
 1100 ALPHZD=ALPHZR*RAD
      BD(81)=ALPHZD
C
C---------     COMPUTE TOTAL BODY PLANFORM AREA, SP              -------
C
 1110 CALL EQSPCE(X,R,BD(701),S,NX,XE,YE,XOLE,SE,20,
     1            DSEDX)
      CALL TRAPZ(YE,20,XE,SP,1)
 1120 SP = 2.0 * SP
C
C---------     COMPUTE CM-ZERO                              ------------
C
      CALL TRAPZ(CMP,20,XOLE,CMZP,1)
 1130 CM0=2.*CMZP*(SP/SREF)*(X(NX)/LREF)
      BD(62)=CM0
C
C---------     COMPUTE ASPECT RATIO FOR BODY ALONE, ARB     ------------
C              FIND Y-MAX OF PLANFORM VIEW
C
      CALL GETMAX(X,R,NX,XMAX,YMAX,IMAX)
 1140 ARB= (4.0*YMAX*YMAX)/SP
C
C     -------    LOCATE BODY STATION OF .96 YMAX  ------
C
 1150 Y95=0.96*YMAX
      IM=IMAX-1
      X95=X(IM)+(X(IMAX)-X(IM))*(Y95-R(IM))/(R(IMAX)-R(IM))
C
C---------     FIGURE BA-1 ( KP VS ARB )                    ------------
C
      CALL TBFUNX(ARB,AKP,DKPDAR,9,XBA1,YBA1,CBA1,IBA1,MI,NG,
     1            0,2,QKPARB,3,ROUTID)
C
C---------     FIGURE BA-2 (KV VS ARB )                     ------------
C
 1160 CALL TBFUNX(ARB,AKV,DKVDAR,9,XBA2,YBA2,CBA2,IBA2,MI,NG,
     1             0,2,QKVARB,3,ROUTID)
C
C---------     FIND SMAX, XSMAX, AND NMAX AT DSDX= ZERO          -------
C
      SMAX=S(NX)
      DEFF=SQRT(4.0*SMAX/PI)
      IF(DEFF.EQ.0.)GO TO 1190
C
C---------     COMPUTE FINENESS RATIO, FR                   ------------
C
 1170 FR= X(NX)/DEFF
C
C---------     FIND DEPTH AND WIDTH AT SMAX                 ------------
C
      YSMAX=R(NX)
      DSMAX=D(NX)
 1180 RATIO= 2.0*YSMAX/DSMAX
      GO TO 1200
 1190 RATIO=3.0
      FR=10.0
 1200 CONTINUE
C
C     -------    COMPUTE PLANFORM AREA FOR BODY LENGTH X95 AT .96YMAX --
C     -------    FIND INDEX OF LOCATION OF .96XMAX
C
      DO 1210  I= 1,NX
         TX95(I)= X(I)
         TY95(I)= R(I)
         IF(X95.GT.X(I)) GO TO 1210
         N95= I
         TX95(I)= X95
         TY95(I)= Y95
         GO TO 1220
 1210 CONTINUE
 1220 CONTINUE
C
C---------     EQUAL SPACE X AND Y ARRAYS FROM X(1) TO X95       -------
C
      CALL EQSPCE(TX95,TY95,RP,RS,N95,XE95,YE95,PE95,SE95,20,DYEDX)
C
C---------     COMPUTE PLANFORM AREA FOR BODY LENGTH X95         -------
C
 1230 CALL TRAPZ(YE95,20,XE95,A95,1)
 1240 A95MAX=2.0*A95
C
C---------     FIGURE BA-3 ( ALPHV VS FR FOR SEVERAL A/B )       -------
C
      CALL TLINEX(X1BA3,X2BA3,YBA3,5,12,RATIO,FR,ALPHV,
     1            2,2,2,2,QALFR,3,ROUTID)
C
C---------     COMPUTE EXPONENT FOR (X/L) CNP CALCULATIONS       -------
C
 1250 EXPON= (2.0*Y95*X95)/A95MAX-1.00
C
C---------     COMPUTE ( X/L )CNP,  XOLCNP                       -------
C
      XOLCNP= 2.0*EXPON/(2.0*EXPON+1.0)
C
C---------     COMPUTE   DELTA(X)/L , DELXOL                     -------
C
      DELXOL= XOVLCG-XOLCNP*(X95  /X(NX))
C
C              COMPUTE AREA ARRAY, DELA
C
      CALL TRAPZ(YE,20,XE,DELA,0)
      DO 1260 J= 1,19
 1260 DELA(J)= DELA(J+1)-DELA(J)
      SUM=0.0
C
C---------     COMPUTE PLANFORM AREA MOMENT OF BODY, SUM         -------
C
      DO 1270 NUM=1,19
 1270 SUM = SUM + 2.0*DELA(NUM) * (XE(NUM)+(XE(NUM+1)-XE(NUM))/2.0)
C
C---------     COMPUTE (X/L).CA.,  XOVLCA                        -------
C
      XOVLCA= (SUM/SP)/X(NX)
C
C********           ALPHA LOOP
C
      BD(94)=FLC(2)
      NALPHA=BD(94)+.5
      DO 1300 J=1,NALPHA
C
C     ---------     REVISED ANGLE-OF-ATTACK  (29 MARCH 71 )
C
 1280    ANGD= FLC(J+22) - ALPHZD
         ANGR= ANGD/RAD
         SA= SIN(ANGR)
         CA=COS(ANGR)
C
C---------     COMPUTE (CNP)J,   BD(680+J)                       -------
C
         BD(J+680)=AKP*SA*CA
C
C---------     COMPUTE (CNV)J,   BD(660+J)                       -------
C
         BD(J+660)=AKV*(SIN((FLC(J+22)-ALPHV)/RAD))**2
         IF(FLC(J+22).LE.ALPHV) BD(J+660)=0.0
C
C---------     COMPUTE (CNTOT)J   BODY(J+40)                     -------
C
 1290    BODY(J+60)=(BD(J+680)+BD(J+660))*(SP/SREF)
C
C-----------------COMPUTE (CD)J -----------------------
C
         BODY(J)=BODY(J+60)*SA + BD(10)
C
C---------     COMPUTE  (CM)CNP ,   BD(234+J)                    -------
C
         BD(J+234)= BD(J+680)*(DELXOL       )*(SP/SREF)*(X(NX)/LREF)
C
C---------     COMPUTE (CM) CNV     BD(114+J)                    -------
C
         BD(J+114)= BD(J+660)*(XOVLCG-XOVLCA)*(SP/SREF)*(X(NX)/LREF)
C
C---------     COMPUTE (CM)J,  BODY(40+J)                        -------
C
         BODY(J+40)= CM0 + BD(J+234) + BD(J+114)
C
C     ------CALCULATE CL AND CA --------
C
         BODY(J+20)=(BODY(J+60)-BODY(J)*SIN(FLC(J+22)/RAD))/COS
     1              (FLC(J+22)/RAD)
         BODY(J+80)=BODY(J)*COS(FLC(J+22)/RAD)-BODY(J+20)*SIN(FLC(J+22)/
     1              RAD)
         B(J+22)=FLC(J+22)+ALIW
         BD(J+214)=BODY(J+60)*SA
 1300 CONTINUE
      BODY(121)=AKP*DELXOL*(SP/SREF)*(X(NX)/LREF)/RAD
      BODY(101)=(AKP*SP)/(RAD*SREF)
      RETURN
      END
