      SUBROUTINE MAXCL
C
C  CALCULATE CLMAX AT FLIGHT CONDITIONS
C
      COMMON /IBODY/ PB, NACA(80), BF(232), CBAR
      COMMON /IBWH/   PBWH,AI,ALO,CLI,ASEP,CMCO4,CLA0,CLA(20),CLMAX0,
     1                CLMAX(20)
      COMMON /IBWHV/  PBWHV, RHO,TMAX,DELTAY,XOVC,TOVC,COVC,CXVC,CAMBER
      COMMON /FLGTCD/ FLC(93)
C
      DIMENSION PARM58(26),EVAL58(52),PARM59(30),EVAL59(160),
     1          PARM60(20),EVAL60(30),PARM61(16),EVAL61(32),
     2          PARM62(8),EVAL62(8),VAR(4),LEN(4),ROUT(2),C(6)
      LOGICAL CAMBER
C
      DATA PARM58 /0.,1.,1.1,1.25,1.50,2.00,2.25,2.50,3.00,3.50,4.00,
     1    4.50,5.00,.3,.35,.4,.45,9*0./
      DATA EVAL58 /2*.8,.81,.85,.98,1.2,1.31,1.43,1.58,1.59,1.55,1.47,
     1    1.42,2*.8,.81,.85,.98,1.2,1.31,1.43,1.51,1.52,1.48,1.41,1.41,
     2    2*.8,.81,.85,.98,1.2,1.31,1.39,1.44,1.45,1.43,1.35,1.35,
     3    2*.8,.81,.85,.98,1.2,1.29,1.33,1.35,1.35,1.35,1.35,1.35/
      DATA PARM59 /0.,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,0.,.02,.04,.06,
     1    6*0.,.15,.30,.40,.50,6*0./
      DATA EVAL59 /10*0.,
     1 0.,.12,.23,.32,.23,.12,.07,.07,.06,.03,
     2 0.,.27,.36,.40,.30,.16,3*.1,.06,
     3 0.,.27,.36,.40,.30,.16,3*.1,.06,
     4 10*0.,
     5 2*0.,.17,.22,.25,.15,.08,.06,.08,.02,
     6 0.,.15,.27,.40,.33,.18,.08,.06,.07,.01,
     7 0.,.53,.64,.56,.38,.22,.10,2*.06,0.,
     8 10*0.,
     9 2*0.,.05,.19,.23,.14,.07,2*.06,.02,
     A 0.,.16,.27,.37,.34,.19,.10,.07,.09,.05,
     B 0.,.39,.50,.50,.40,.25,.13,.10,.13,.09,
     C 10*0.,
     D 2*0.,.04,.10,.15,.14,.09,.07,.10,.08,
     E 0.,.10,.19,.27,.30,.22,.15,.13,.18,.12,
     F 0.,.26,.36,.43,.45,.30,.21,.20,.23,.20 /
      DATA PARM60 /0.,.25,1.,1.75,2.,2.75,3.,3.5,4.,4.5,.35,.4,.45,7*0./
      DATA EVAL60 /0.,.16,.19,.20,.18,.08,.05,-.01,-.03,-.04,
     1             0.,.16,.19,.20,.18,.08,.05,.03,.04,.02,
     2             0.,.16,.19,.20,.18,.08,.08,.10,.12,.08 /
      DATA PARM61 /0.,1.,1.5,2.,2.5,3.,3.5,4.,3.E06,6.E06,9.E06,2.5E07,
     1    4*0./
      DATA EVAL61 /0.,-.09,2*-.13,-.09,-.11,-.20,-.24,
     1             2*0.,-.08,-.07,2*-.03,-.07,-.10,
     2             8*0.,0.,.12,.14,.07,.01,.11,.18,.20 /
      DATA PARM62 /0.,1.,1.5,2.,2.5,3.,3.5,4./
      DATA EVAL62 /2*0.,-.19,-.37,-.46,3*-.5/
      DATA ROUT /4HMAXC, 4HL   /
C
C  CLMAX BASE
C
      IF(ABS(XOVC-.3) .LE. 1.E-5) XOVC = .3
      VAR(1)=DELTAY
      VAR(2)=XOVC
      LEN(1)=13
      LEN(2)=4
      CALL INTERX(2,PARM58,VAR,LEN,EVAL58,CLBASE,13,52,0,0,0,0,0,0,0,0,
     1            4HBASE,1,ROUT)
C
C  DELTA CLMAX
C
      DEL1=0.0
      DEL2=0.0
      IF(.NOT.CAMBER) GO TO 1000
      VAR(2)=COVC
      VAR(3)=CXVC
      LEN(1)=10
      LEN(2)=4
      LEN(3)=4
      CALL INTERX(3,PARM59,VAR,LEN,EVAL59,DEL1,10,200,0,0,0,0,0,0,0,0,
     1            4HDEL1,1,ROUT)
      IF(XOVC .EQ. 0.30) GO TO 1000
      VAR(2)=XOVC
      LEN(1)=10
      LEN(2)=3
      CALL INTERX(2,PARM60,VAR,LEN,EVAL60,DEL2,10,30,0,0,0,0,0,0,0,0,
     1            4HDEL2,1,ROUT)
 1000 NMACH=FLC(1)+0.5
      IN=0
      DO 1010 I=1,NMACH
         RN=FLC(I+42)*CBAR
         VAR(2)=RN
         IF(RN .LT. 1.0) VAR(2)=9.0E6
         LEN(1)=8
         LEN(2)=4
         CALL INTERX(2,PARM61,VAR,LEN,EVAL61,DEL3,8,32,0,0,0,0,0,0,0,0,
     1               4HDEL3,1,ROUT)
         CALL TBFUNX(DELTAY,DEL4,DYDX,8,PARM62,EVAL62,C,IN,MI,NG,0,0,
     1               4HDEL4,1,ROUT)
         CLMAX0=CLBASE+DEL1+DEL2+DEL3
         CLMAX(I)=CLMAX0
 1010 CONTINUE
      RETURN
      END
