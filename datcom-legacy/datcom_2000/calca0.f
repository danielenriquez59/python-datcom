      SUBROUTINE CALCA0(A,B,AIN)
C
C***  COMPUTES LIFTING SURFACE ALPHA ZERO LIFT
C
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      DIMENSION ROUTID(2)
      DIMENSION A(195),B(49),AIN(64)
      DIMENSION TR(3),AR(4),SAQC(22),DA0OT(22,10),DA(220)
      DIMENSION DAP(132),DAPP(88)
      EQUIVALENCE (DAP(1),DA(1)),(DAPP(1),DA(133))
      DIMENSION NPT(7),LOCX(7),LOCY(7),CC(6,7),IN(7),MI(7)
      DIMENSION TOC(7) , CX(12),CY(33)
      LOGICAL CAMBER,NOINTP,NOINT2 , NOINT3
      EQUIVALENCE (CAMBER,CAFAKE)
      EQUIVALENCE (DA0OT(1,1) , DA(1))
      DIMENSION Q41315(3)
      DATA Q41315 /4H4.1.,4H3.1-,4H5   /
      DATA ROUTID /4HCALC,4HA0  /
      DATA TR     /0.,.5,1./ , AR/1.5,3.5,6.0,10./
      DATA SAQC   /-45.,-40.,-35.,-30.,-25.,-20.,-15.,-10.,-5.,0.,
     1    5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60./
      DATA DAP/-.399  ,-.3995 ,-.400  ,-.400  ,-.4005 ,-.4005 ,-.4005 ,
     1 -.4005 ,-.4005 ,-.400  ,-.400  ,-.400  ,-.3995 ,-.3990 ,-.3985 ,
     2 -.398  ,-.397  ,-.3955 ,-.3940 ,-.385  ,-.380  ,-.372  ,
     3         -.384  ,-.385  ,-.3855 ,-.386  ,-.386  ,-.3865 ,-.386  ,
     4 -.3855 ,-.3845 ,-.384  ,-.3825 ,-.3815 ,-.380  ,-.378  ,-.375  ,
     5 -.372  ,-.369  ,-.364  ,-.358  ,-.350  ,-.3435 ,-.335  ,
     6         -.375  ,-.375  ,-.375  ,-.3745 ,-.374  ,-.373  ,-.372  ,
     7 -.371  ,-.370  ,-.3685 ,-.367  ,-.365  ,-.362  ,-.359  ,-.355  ,
     8 -.3515 ,-.347  ,-.342  ,-.336  ,-.331  ,-.325  ,-.318  ,
     9         -.417  ,-.4155 ,-.414  ,-.413  ,-.412  ,-.411  ,-.4105 ,
     A -.410  ,-.4095 ,-.409  ,-.408  ,-.4075 ,-.407  ,-.4065 ,-.406  ,
     B -.405  ,-.404  ,-.4025 ,-.401  ,-.396  ,-.393  ,-.387  ,
     C         -.430  ,-.427  ,-.424  ,-.422  ,-.420  ,-.4175 ,-.4155 ,
     D -.414  ,-.412  ,-.410  ,-.4085 ,-.407  ,-.405  ,-.403  ,-.401  ,
     E -.399  ,-.396  ,-.393  ,-.390  ,-.385  ,-.381  ,-.375  ,
     F         -.437  ,-.434  ,-.431  ,-.428  ,-.4245 ,-.422  ,-.419  ,
     G -.417  ,-.414  ,-.412  ,-.409  ,-.407  ,-.404  ,-.402  ,-.399  ,
     H -.3965 ,-.394  ,-.391  ,-.3885 ,-.3875 ,-.386  ,-.380/
      DATA DAPP/-.419 ,-.417  ,-.416  ,-.414  ,-.413  ,-.413  ,-.412  ,
     1 -.411  ,-.4105 ,-.410  ,-.4095 ,-.409  ,-.408  ,-.4075 ,-.407  ,
     2 -.406  ,-.405  ,-.4035 ,-.402  ,-.400  ,-.395  ,-.390  ,
     3         -.4405 ,-.4365 ,-.433  ,-.430  ,-.427  ,-.425  ,-.422  ,
     4 -.4205 ,-.419  ,-.4175 ,-.416  ,-.415  ,-.413  ,-.412  ,-.410  ,
     5 -.409  ,-.407  ,-.406  ,-.406  ,-.4055 ,-.405  ,-.405  ,
     6         -.456  ,-.451  ,-.447  ,-.442  ,-.439  ,-.436  ,-.433  ,
     7 -.431  ,-.4285 ,-.426  ,-.424  ,-.422  ,-.420  ,-.419  ,-.417  ,
     8 -.416  ,-.415  ,-.414  ,-.415  ,-.416  ,-.417  ,-.418  ,
     9         -.469  ,-.465  ,-.460  ,-.456  ,-.452  ,-.449  ,-.445  ,
     A -.442  ,-.439  ,-.437  ,-.434  ,-.432  ,-.429  ,-.428  ,-.426  ,
     B -.425  ,-.424  ,-.423  ,-.423  ,-.423  ,-.424  ,-.425/
      DATA IN   /7*0/
      DATA TOC  /16.,14.,12.,10.,9.,8.,7./
      DATA NPT  /7,6,4,4,3,4,5/
      DATA LOCX /1,3,5,6,8,8,8/
      DATA LOCY /1,8,14,18,22,25,29/
      DATA CX   /.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95/
      DATA CY /1.,.95,.85,.65,.35,-.15,-1. , 1.,.95,.75,.45,-.06,-1.25 ,
     1 1.,.9,.6,-.5 , 1.,1.,.85,.1 , 1.,.8,0. , 1.,.95,.65,0. ,
     2 1.,1.,.95,.85,.55/
C
      CAFAKE = AIN(64)
      IF(AIN(10).EQ.UNUSED) GO TO 1000
      A(134)=AIN(10)
      GO TO 1010
 1000 A(134)=AIN(20)-AIN(19)/A(131)
C
C     ----TEST FOR TWIST
C
 1010 IF(ABS(AIN(11)).LT.0.5) GO TO 1160
C
C     ----HAVE TWIST-CALCULATE CORRECTED ALPHA0
C
      NOINTP = .FALSE.
      NOINT2 = .FALSE.
      NOINT3  = .FALSE.
      DO 1030 I=1,3
         II = I
         TEMP = A(27) - TR(I)
         IF(ABS(TEMP).LT.2.E-2) GO TO 1040
         IF(TEMP.LT.0.) GO TO 1050
         TEMPP = TEMP
 1030 CONTINUE
 1040 NOINTP = .TRUE.
 1050 ITR = 3*(II-1)
      IF(II.EQ.1) NOINTP =.TRUE.
      IF(.NOT.NOINTP) TRR = TEMPP/(TEMPP-TEMP)
      DO 1060 I=1,4
         II = I
         TEMP2 = A(7) - AR(I)
         IF(ABS(TEMP2).LT.2.E-2) GO TO 1070
         IF(TEMP2.LT.0.) GO TO 1080
         TEMPP2 = TEMP2
 1060 CONTINUE
 1070 NOINT2 = .TRUE.
 1080 IAR = II
      IF(II.EQ.1) NOINT2 =.TRUE.
      IF(.NOT.NOINT2) ARR = TEMPP2/(TEMPP2-TEMP2)
      DO 1090 I=1,22
         II = I
         TEMP3 = A(40) - SAQC(I)
         IF(ABS(TEMP3).LT.2.E-2) GO TO 1100
         IF(TEMP3.LT.0.0) GO TO 1110
         TEMPP3 = TEMP3
 1090 CONTINUE
 1100 NOINT3= .TRUE.
 1110 IF(II.EQ.1) NOINT3 = .TRUE.
      IF(.NOT.NOINT3) SAR = TEMPP3/(TEMPP3-TEMP3)
      IP1=1
      IF(IAR.EQ.4.AND.A(7).LE.0.5)IAR=3
 1120 IDX=ITR+IAR
      D2=DA0OT(II,IDX)
      IF(.NOT.NOINT3)D2=DA0OT(II-1,IDX)+SAR*(D2-DA0OT(II-1,IDX))
      IF(NOINT2) GO TO 1130
      IDX=IDX-1
      D1=DA0OT(II,IDX)
      IF(.NOT.NOINT3)D1=DA0OT(II-1,IDX)+SAR*(D1-DA0OT(II-1,IDX))
      D2=D1+ARR*(D2-D1)
 1130 IF(NOINTP) GO TO 1150
      IF(IP1.EQ.2) GO TO 1140
      IP1=IP1+1
      ITR=ITR-3
      DSV=D2
      IF(IAR.EQ.4) IAR=3
      GO TO 1120
 1140 D2=D2+TRR*(DSV-D2)
 1150 A(135)=D2
      A(134) = AIN(11) * A(135) + A(134)
C
C     ----TEST FOR CAMBER
C
 1160 CONTINUE
      IF(.NOT.CAMBER) GO TO 1230
C
C     ----HAVE CAMBER-CALCULATE ALPHA0 FOR EACH MACH NUMBER
C
      NOINTP = .FALSE.
      DO 1170 I=1,7
         II = I
         TEMP = 100.*AIN(16) - TOC(I)
         IF(ABS(TEMP).LT.2.E-2) GO TO 1180
         IF(TEMP.GT.0.0) GO TO 1190
         TEMPP = TEMP
 1170 CONTINUE
 1180 NOINTP = .TRUE.
 1190 IF(II.EQ.1.OR.TEMP.LT.0.) NOINTP = .TRUE.
      IF(.NOT.NOINTP)FACT=TEMPP/(TEMPP-TEMP)
      XAG = A(43) * B(1)
      I=1
 1200 NP = NPT(II)
      IX = LOCX(II)
      IY = LOCY(II)
      CALL TBFUNX(XAG,YAG1,DYDX,NP,CX(IX),CY(IY),CC(1,II),IN(II),MI(II),
     1            NG,-1,0,Q41315,3,ROUTID)
      IF(I.EQ.2) GO TO 1210
      IF(NOINTP) GO TO 1220
      I=I+1
      II=II-1
      YSV = YAG1
      GO TO 1200
 1210 YAG1=YAG1+(YSV-YAG1)*FACT
      II=II+1
 1220 B(49) = YAG1 * A(134)
      A(137) = YAG1
      RETURN
 1230 CONTINUE
      B(49) = A(134)
      RETURN
      END
