      SUBROUTINE WTLIFT(A,B,AIN,AOUT)
C
C *** CALCULATES WING OR TAIL LIFT CHARACTERISTICS
C
      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF
      COMMON /OVERLY/ NLOG,NMACH,IM,NALPHA,IG
      COMMON /CONSNT/ PI, DEG, UNUSED, RAD
      DIMENSION ROUTID(2),Q3417B(3),Q13242(3),Q3418A(3),Q3418B(3)
      DIMENSION Q13415(3),Q3414A(3),Q3414B(3)
      DIMENSION AMACH(5) , DYA(6) , SALE4(4), DCLTB(5,6,4) , DCAR(120)
      DIMENSION A(195),B(49),AOUT(101),AIN(101)
      DIMENSION DACLMX(13,4),DACLL(52)
      DIMENSION SALE(13) , DELTAY(7) ,CLOCL(13,7),CLL(91) , DY(4)
      DIMENSION TR(11),C2(11),C(6)
      DIMENSION C1ABCS(10),ACLMX(10),C17(6)
      DIMENSION DMN(3), C1TAB(20), DACL(60)
      DIMENSION ACLE(10),DACLO(100) , C1TABO(10)
      DIMENSION  BA(12),CLOVCL(12) , CBA(6)
      DIMENSION WTYPE(4)
C
      EQUIVALENCE (DCLTB(1,1,1),DCAR(1))
      EQUIVALENCE (CLOCL(1,1),CLL(1)) ,(DACLMX(1,1),DACLL(1))
C
      LOGICAL ITYPTT
      LOGICAL LOARAT
C
      DATA ROUTID /4HWTLI,4HFT  /
      DATA Q3417B /4H4.1.,4H3.4-,4H24B /,
     1     Q13242 /4H4.1.,4H3.2-,4H52  /,
     2     Q3418A /4H4.1.,4H3.4-,4H25A /,
     3     Q3418B /4H4.1.,4H3.4-,4H25B /,
     4     Q13415 /4H4.1.,4H3.4-,4H22  /,
     5     Q3414A /4H4.1.,4H3.4-,4H21A /,
     6     Q3414B /4H4.1.,4H3.4-,4H21B /
C
C     ----4.1.3.4-22
C
      DATA AMACH/.2,.3,.4,.5,.6/
      DATA DYA  /2.,2.25,2.5,3.,4.,4.5/
      DATA SALE4/0.,20.,40.,60./
      DATA DCAR/0.,2*-.02,2*0., 0.,-.13,-.19,-.2,-.19 , 0.,-.185,-.32,
     1-.4,-.445 , 0.,-.21,-.36,-.45,-.5 , 0.,-.24,-.42,-.545,-.64 ,
     2-0.,-.24,-.455,-.605,-.72 , 0.,-.04,-.045,-.02,0. , 0.,-.095,-.165
     3,-.22,-.25 , 0.,-.105,-.19,-.26,-.29 , 0.,-.125,-.225,-.3,-.345 ,
     4 0.,-.15,-.265,-.37,-.45 , 0.,-.15,-.29,-.41,-.51 , 0.,-.04,-.07,
     5-.095,-.1 , 0.,-.043,-.080,-.108,-.123 , 0.,-.045,-.085,-.12,-.145
     6 , 0.,-.05,-.095,-.14,-.19 , 0.,-.07,-.125,-.19,-.26 , 0.,-.07,
     7-.14,-.215,-.3,4*0.,-.02,4*0.,-.022,4*0.,-.03,3*0.,-.02,-.07,2*0.,
     8-.02,-.055,-.085,0.,-.04,-.085,-.14,-.2/
C
C     ----4.1.3.4-21 A
C
      DATA SALE/0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60./
      DATA DELTAY/1.4,1.6,1.8,2.0,2.2,2.4,2.5/
      DATA CLL /.9,.915,.93,.95,.975,1.,1.03,1.07,1.1,1.15,1.19,1.25,1.3
     1 , .9,.91,.92,.94,.96,.98,1.,1.02,1.05,1.08,1.11,1.15,1.19 ,   .9,
     2 .91,.92,.93,.94,.95,.96,.965,.975,.99,1.,1.01,1.03 , 6*.9, 2*.89,
     3.885,2*.88,.875,.87 , .9,.895,.89,.88,.87,.86,.85,.835,.82,.8,.78,
     4.755,.73 , .9,.89,.88,.87,.85,.83,.81,.79,.76,.725,.695,.65,.59 ,
     5 .9,.89,.87,.86,.84,.82,.795,.77,.73,.7,.65,.59,.52 /
C
C     ----4.1.3.4-21 B
C
      DATA DY/1.2,2.,3.,4./
      DATA DACLL/1.75,1.9,2.2,2.7,3.4,4.15,5.1,6.1,7.3,8.7,10.15,11.75,
     113.3 , .1,.5,1.05,1.65,2.3,3.1,3.9,4.7,5.7,6.7,7.7,8.75,9.8,
     21.2,1.4,1.7,2.,2.4,2.85,3.35,3.7,4.25,4.7,5.3,5.9,6.65,
     32.2,2.1,2*2.,2.1,2.15,2.3,2.4,2.55,2.7,2.9,3.05,3.3/
C
C     ----4.1.3.4-24(B)
C
      DATA IN/0/
      DATA TR/0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1./
      DATA C2/0.,.21,.5,.9,1.08,1.05,1.,.94,.9,.86,.85/
C
C     ----4.1.3.4-25(A)
C
      DATA C1ABCS/0.,.4,.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6/
      DATA ACLMX/3*35.,32.,28.,25.,23.2,22.,21.5,21./
C
C     ----4.1.3.4-25(B)
C
      DATA I17/0/
      DATA DMN/ .2,.4,.6/
      DATA C1TAB/4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,10.5,11.,
     111.5,12.,12.5,13.,13.5,14./
      DATA DACL /0.,.5,.9,1.4,1.9,2.5,3.3,4.,4.6,5.6,6.4,7.3,8.2,9.2,10.
     1,11.,12.,13.,14.,15. , 0.,.2,.4,.7,1.2,1.7,2.4,3.,3.7,4.6,5.2,6.,
     26.9,7.8,8.6,9.5,10.4,11.4,12.3,13.5 , 0.,0.,.1,.2,.5,.7,1.,1.3,1.6
     3,2.,2.5,3.,3.6,4.3,4.9,5.5,6.2,7.,7.6,8.5/
      DATA ACLE /0.,2.,3.,4.,5.,6.,7.,8.,9.,30./
      DATA DACLO /                                    10.,8.5,6.9,5.5,4.
     1,2.6,1.5,.7,.1,0. , 8.7,7.3,5.3,4.2,2.6,1.4,0.5,-.2,-.5,0. , 7.5,
     25.9,4.2,2.5,1.2,0.,-.7,-1.1,-.8,0. , 5.5,3.4,1.6,0.,-1.3,-2.1,-2.5
     3,-2.,-.8,0. , 3.,.7,-1.4,-3.3,-4.3,-4.3,-3.1,-2.,-.8,0. , .3,-2.5,
     4-4.7,-5.8,-5.3,-4.3,-3.1,-2.,-.8,0. ,-2.2,-5.,-6.7,-6.3,-5.3,-4.3,
     5-3.1,-2.,-0.8,0. , -3.3,-6.6,-7.2,-6.3,-5.3,-4.3,-3.1,-2.,-.8,0.,
     6  -4.2,-7.,-7.2,-6.3,-5.3,-4.3,-3.1,-2.,-.8,0. , -8.5,-7.9,-7.2,
     7 -6.3,-5.3,-4.3,-3.1,-2.,-.8,0. /
      DATA C1TABO /.0,.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5/
C
C     ----4.1.3.2-52
C
      DATA BA/1.15,1.4,2.,2.2,3.,3.6,4.,5.,6.,7.,8.,9./
      DATA CLOVCL /1.25,1.2,1.12,1.1,1.04,1.,.985,.96,.95,3*.94/
      DATA IN4132 /0/
      DATA WTYPE /4HSTRA ,4HDOUB ,4HCRAN ,4HCURV /
C
C     ----INITIALIZATION ENTRY
C
      ITYPTT = .TRUE.
      LOARAT = .FALSE.
      IF(A(7).LT.A(125)) LOARAT = .TRUE.
      TYPE = AIN(15)
C
C     ----COMPUTE CLA0 + CLAI FOR OTHER THAN STRAIGHT TAPERED
C
      IF(TYPE.EQ.WTYPE(1)) GO TO 1000
C
C     ----CLA0
C
      TEMP1 = 2.*PI*A(168)*DEG
      TEMP2 =(A(168)*DEG*2.*PI/A(131))**2
      TEMP3 = A(98)**2
      TEMP4=B(2)**2
      A(172)=TEMP1/(2.+SQRT(TEMP2*(1.+TEMP3/TEMP4)+4.))
C
C     ----CLAI
C
      TEMP1 = 2.*PI*A(5)*DEG
      TEMP2 = (A(5)*DEG*2.*PI/A(131))**2
      TEMP3 = A(74)**2
      A(171)=TEMP1/(2.+SQRT(TEMP2*(1.+TEMP3/TEMP4)+4.))
 1000 CONTINUE
C
C     ----CURVED WING TEST
C
      IF(TYPE.EQ.WTYPE(4)) GO TO 1070
      ITYPTT = .FALSE.
C
C     ----4.1.3.4-24(B)
C
      CALL TBFUNX(A(27),A(159),DCDT,11,TR,C2,C,IN,MI,NG,
     1            0,0,Q3417B,3,ROUTID)
      A(160) = (A(159)+1.)*A(38)*A(7)
C
C     ----COMPUTE BETA*S
C
      TEMP1 = 2.*PI*A(7)*DEG
      TEMP2=A(7)*DEG*2.*PI/A(131)
      TEMP2 = TEMP2**2
      TEMP3 = A(50)**2
      CALL TLINEX(DELTAY, SALE,CLL,7,13,AIN(17),A(34),A(145),
     1           -1,0,-1,2,Q3414A,3,ROUTID)
      CALL TLINEX(DY,SALE,DACLL,4,13,AIN(17),A(34),A(144),
     1           -1,0,-1,2,Q3414B,3,ROUTID)
 1010 IF(ITYPTT) GO TO 1070
      TEMP4 = B(2)**2
      AOUT(101)=(TEMP1/(2.+SQRT(TEMP2*(1.+TEMP3/TEMP4)+4.)))*A(3)/
     1          SREF
      IF(TYPE.NE.WTYPE(3)) GO TO 1020
C
C     ----HERE FOR CRANKED HIGH ASPECT RATIO WING-FINISH CLA COMPUTATION
C     ----FIGURE 4.1.3.2-52
C
      BAARG = A(7)* B(2)
      CALL TBFUNX(BAARG,CLRAT,DYDX,12,BA,CLOVCL,CBA ,IN4132,M4132,N4132,
     1            2,0,Q13242,3,ROUTID)
      AOUT(101) = CLRAT * AOUT(101)
 1020 CONTINUE
      IF(LOARAT) GO TO 1030
      CALL TLIN3X(DYA,AMACH,SALE4,DCAR,6,5,4,AIN(17),B(1),A( 34),TEMP5
     1            ,0,-1,0,2,2,2,Q13415,3,ROUTID)
      A(146)=A(145)*A(132)
      B(44)=(A(145)*A(132)+TEMP5)*A(3)/SREF
      B(43) = B(44)/AOUT(101)+B(49)+A(144)
      GO TO 1060
C
C     ----HERE FOR LOW ASPECT RATIO
C
 1030 C1P1AC = A(7)*A(124)/B(2)
      CALL CLMXBS(C1P1AC,B(44),A,B,AIN)
C
C     ----4.1.3.4-25(A)
C
      CALL TBFUNX(C1P1AC,ACLBAS,DYDX,10,C1ABCS,ACLMX,C17,I17,MI,NG,
     1            0,0,Q3418A,3,ROUTID)
      IF(A(160).LE.4.5) GO TO 1040
      CALL TLINEX(DMN,C1TAB,DACL,3,20,B(1),A(160),DACMAX,
     1            -1,0,2,0,Q3418B,3,ROUTID)
      GO TO 1050
 1040 TMP=A(7)*A(37)*(1.+4.*A(27)**2)
      CALL TLINEX(ACLE,C1TABO,DACLO,10,10,TMP,A(160),DACMAX,
     1            0,0,2,0,Q3418B,3,ROUTID)
 1050 B(43)=ACLBAS+DACMAX
 1060 CONTINUE
      RETURN
 1070 WRITE(6,1080) TYPE,LOARAT
 1080 FORMAT(33H NO CLALPHA COMPUTATION- WING IS  ,A4,18H LOW ASPECT RAT
     1IO= ,L1)
      RETURN
      END
