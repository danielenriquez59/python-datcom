      SUBROUTINE TRSONI(NZ)
C
C*** COMPUTES TRANSONIC WING LIFT SLOPE,CLMAX,ALPHA CLMAX,AND CD0
C***      BODY LIFT AND MOMENT SLOPES,DRAG AT ANGLE OF ATTACK
C  *** TRANSONIC LATERAL STABILITY CALCULATED IN SUBLAT ***
C
      DIMENSION ROUTID(2)
      DIMENSION Q3243A(3),Q3243B(3),Q3244B(3),Q3244A(3),Q3244C(3),Q34171
     1  (4),Q34172(4),Q3418A(3),Q34181(3),Q34182(3),Q3419B(3),Q3419C(3),
     2  Q3419A(3),Q15127(4),Q5129R(3),Q5129L(3),Q23124(3),Q23126(3)
      REAL K,MACH,MFB,MT(5)
      REAL MFB0
      DIMENSION T429L(15),D429L(56),T429R(14),D429R(35),T424(16),
     1          D424(39)
      DIMENSION WTYPE(4),C17(6) ,TR(11),DR(11)
      DIMENSION T426(16),D426(55)
      DIMENSION X(11),Y(11),C(6)
      DIMENSION T418A(10),D418A(10), T18B1(3),D18B1(20),C18B1(60),
     1          T18B2(10), D18B2(10) ,C18B2(100)
      DIMENSION T419A(22),D419A(77) ,T419B(13),D419B(13),T419C(26),
     1          D419C(130)
      DIMENSION T43A(15),D43A(56),T43B(13),D43B(40),T44A(13),D44A(42),
     1          T44B(13),D44B(42),T44C(6),D44C(6)
      DIMENSION XMTD(15),CDW1(15),CDW2(15),VAR(4),LGH(4)
      DIMENSION CLAMT(5),ALPHA(20) ,CDJB(20)
      DIMENSION X27M(4),X27I(4),C27(6)
C
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /OVERLY/ NLOG,NMACH,M,NALPHA,IG,NF
      COMMON /OPTION/ SR,CBARR,RUFF,BLREF
      COMMON /FLGTCD/ FLC(93)
      COMMON /BDATA/  BD(762)
      COMMON /WINGD/  A(195),B(49)
      COMMON /WINGI/  WINGIN(77)
      COMMON /SUPBOD/ SBD(227)
      COMMON /SBETA/  STB(135),TRA(108)
      COMMON /IBODY/  PB, BODY(400)
      COMMON /IWING/  PW, WING(400)
      COMMON /IBW/    PBW, BW(380)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4                VFPL,VFSC,CTAB
C
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4         VFPL,VFSC,CTAB
      EQUIVALENCE (TOC,WINGIN(16)),(ARSTAR,A(7)),(SWEPE2,A(71)),
     1            (COSLE,A(61)),(TAPR,A(27)),(CBAR,A(16)),(RNFS,A(129)),
     2          (XOC,WINGIN(18)),(SRSTAR,A(3)),(COSC4,A(67)),(SS,BD(93))
     3,           (DB,SBD(6)),(DMAX,SBD(120)),(XLB,SBD(2)),(SB,BD(57)),
     4            (TANC2,A(74))
      EQUIVALENCE (CLA14,TRA(1)),(CLAMO,WINGIN(69)),(K,TRA(3))
     1  ,(MFB0,TRA(5)),(MFB,TRA(6)),(AOC,TRA(7)),(CFBCT,TRA(8)),(BETAFB,
     2  TRA(9)),(CLAFBT,TRA(10)),(MACH,TRA(4)),(CLAFB,TRA(12)),(CLAA,TRA
     3  (13)),(BOC,TRA(14)),(CLAB,TRA(15)),(MT(1),TRA(16)),(CLAMT(1),TRA
     4(21)),(CLA,WING(101)),(C1,TRA(27)),(ARATIO,TRA(28)),(BU4,TRA(29)),
     5   (CLMAX6,TRA(30)),(ACLBAS,TRA(31)),(DACMA6,TRA(32)),(C3,TRA(33))
      EQUIVALENCE (CLAB14,SBD(18)),(CMAB14,SBD(110)),(CD014,SBD(124)),
     1     (CMAB6,BODY(121)),(CFB6,BD(92)),(CD0B6,BD(61)),(CDBB6,BD(60))
      EQUIVALENCE (DALCM,TRA(34)),(DCLMAX,TRA(35)),(ALCLM6,TRA(36)),
     1        (ALCLMT,TRA(37)),(CLMAXT,TRA(38)),(RLCOFF,TRA(39)),
     2        (RNN,TRA(40)),(RL,TRA(41)),(CF,TRA(42)),(CDW2(1),TRA(43)),
     3        (CDW,TRA(67)),(CDF,TRA(68)),(CD0W,WING(1)),(CLAW6,TRA(70))
     4        ,(CLABM6,BODY(101)),(CLABD,BODY(101)),(CMAB,BODY(121)),
     5        (CDPB,TRA(80)),(CDBFIG,TRA(81)),(CDBB,TRA(76)),
     6        (CDWB,TRA(77)),(CD0B,TRA(78)),(CDJB(1),BODY(1)),
     7        (CD0WB,BW(1)),(CDFB,TRA(79)),(DCLA14,TRA(82))
      EQUIVALENCE (TANLE,A(62))
      LOGICAL LOWAR
      DATA ROUTID /4HTRSO,4HNI  /
      DATA Q3243A /4H4.1.,4H3.2-,4H53A /,Q3243B/4H4.1.,4H3.2-,4H53B /,
     1Q3244B/4H4.1.,4H3.2-,4H54B /,Q3244A/4H4.1.,4H3.2-,4H54A /,Q3244C/
     24H4.1.,4H3.2-,4H54C /,QCLA/4HCLA /,Q34171/4H4.1.,4H3.4-,4H24B , 4H
     3C1  /,Q34172/4H4.1.,4H3.4-,4H24B ,4HC2  /,Q3418A/4H4.1.,4H3.4-, 4H
     425A /,Q34181/4H4.1.,4H3.4-,4H25B1/,Q34182/4H4.1.,4H3.4-,4H25B2/,
     5Q3419B/4H4.1.,4H3.4-,4H26B /,Q3419C/4H4.1.,4H3.4-,4H26C /,Q3419A
     6/4H4.1.,4H3.4-,4H26A /,Q15127/4H4.1.,4H5.1-,4H27  ,4HCEPT/,Q5129R
     7/4H4.1.,4H5.1-,4H29R /,Q5129L/4H4.1.,4H5.1-,4H29L /,QCDW/4HCDW /,
     8Q23124/4H4.2.,4H3.1-,4H24  /,Q23126/4H4.2.,4H3.1-,4H26  /
C
C              FIGURE 4.1.5.1-27
C
      DATA X27I  /1.57780,1.67221,1.98509,2.28874/
      DATA X27M  /0.0,1.0,2.0,3.0/
      DATA I27   /0/
      DATA WTYPE /4HSTRA,4HDOUB,4HCRAN,4HCURV/
      DATA I17,IN/0,0/
      DATA NP/11/
C
C              FIGURE 4.1.3.4-24B
C
      DATA TR/0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1./
      DATA DR/0.,.21,.5,.9,1.08,1.05,1.,.94,.9,.86,.85/
      DATA Y/0.,.225,.47,.496,.43,.32,.21,.125,.075,.0475,.0/
      DATA X/0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1./
C
C     ----4.1.3.4-25(A)
C
      DATA T418A /0.,.4,.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6/
      DATA D418A/3*35.,32.,28.,25.,23.2,22.,21.5,21./
C
C     ----4.1.3.4-25(B)
C
      DATA T18B1/.2,.4,.6/
      DATA D18B1/4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,10.5,11.,
     1           11.5,12.,12.5,13.,13.5,14./
      DATA C18B1/0.,.5,.9,1.4,1.9,2.5,3.3,4.,4.6,5.6,6.4,7.3,8.2,9.2,10.
     1,11.,12.,13.,14.,15. , 0.,.2,.4,.7,1.2,1.7,2.4,3.,3.7,4.6,5.2,6.,
     26.9,7.8,8.6,9.5,10.4,11.4,12.3,13.5 , 0.,0.,.1,.2,.5,.7,1.,1.3,1.6
     3,2.,2.5,3.,3.6,4.3,4.9,5.5,6.2,7.,7.6,8.5/
      DATA T18B2/0.,2.,3.,4.,5.,6.,7.,8.,9.,30./
      DATA C18B2 /                                    10.,8.5,6.9,5.5,4.
     1,2.6,1.5,.7,.1,0. , 8.7,7.3,5.3,4.2,2.6,1.4,0.5,-.2,-.5,0. , 7.5,
     25.9,4.2,2.5,1.2,0.,-.7,-1.1,-.8,0. , 5.5,3.4,1.6,0.,-1.3,-2.1,-2.5
     3,-2.,-.8,0. , 3.,.7,-1.4,-3.3,-4.3,-4.3,-3.1,-2.,-.8,0. , .3,-2.5,
     4-4.7,-5.8,-5.3,-4.3,-3.1,-2.,-.8,0. ,-2.2,-5.,-6.7,-6.3,-5.3,-4.3,
     5-3.1,-2.,-0.8,0. , -3.3,-6.6,-7.2,-6.3,-5.3,-4.3,-3.1,-2.,-.8,0.,
     6  -4.2,-7.,-7.2,-6.3,-5.3,-4.3,-3.1,-2.,-.8,0. , -8.5,-7.9,-7.2,
     7 -6.3,-5.3,-4.3,-3.1,-2.,-.8,0. /
      DATA D18B2  /.0,.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5/
C
C              FIGURE 4.1.3.4-26B
C
      DATA T419B
     1 / 0., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 16. /
      DATA D419B
     1 / 1., 1., .93, .8, .55, .4, .3, .24, .2, .16, .14, .12, .12 /
C
C              FIGURE 4.1.3.2-53A
C
      DATA T43A/
     1 0.,     2.,     4.,     6.,     8.,     10.,    12.,    14. ,
     2 0., 1., 2., 3., 4., 6., 8./
      DATA D43A
     1/8*1.0,
     2 6*1.0,.98,.9,
     3 5*1.0,.983,.89,.855,
     4 4*1.0,.939,.888,.857,.835,
     5 2*1.0,.958,.904,.869,.842,.82,.802,
     6 1.0,    .952,   .908,   .872,   .842,   .816,   .80,    .787,
     7 1.0,    .952,   .908,   .872,   .842,   .816,   .80,    .787 /
C
C              FIGURE 4.1.3.2-53B
C
      DATA T43B
     1  / 0.0, 0.17453, 0.34907 , 0.52360, 0.69813,
     2         1.04720, 1.22173 , 1.57080    ,
     3 .799,.85,.9,.95,1.000/
      DATA D43B/
     1 .8,     .807,   .82,    .841,   .867,   .918,   .945,   1.0,
     2 .85,    .855,   .865,   .88,    .9,     .938,   .958,   1.0,
     3 .90,    .903,   .91,    .92,    .932,   .958,   .974,   1.0,
     4 .95,    .952,   .955,   .96,    .967,   .98,    .987,   1.0,
     5    1.0,  1.0,  1.0,  1.0,  1.0,  1.0,1.0,  1.0  /
C
C              FIGURE 4.1.3.2-54A
C
      DATA T44A/
     1 0.,     2.,     4.,     6.,     8.,     12.,    14. ,
     2 1.,     2.,     3.,     4.,     6.,     8. /
      DATA D44A/
     1 1.09,   1.088,  1.056,  1.008,  .962,   .872,   .828,
     2 .98,    1.088,  1.1,    1.05,   .992,   .882,   .83,
     3 .85,    1.,     1.075,  1.06,   1.005,  .89,    .835,
     4 .74,    .9,     1.015,  1.07,   1.025,  .895,   .837,
     54*1.15,1.07,.915,.915,
     65*1.12,.94,.94/
C
C              FIGURE 4.1.3.2-54B
C
      DATA T44B/
     1 0.,     4.,     6.,     8.,     10.,    12.,    14. ,
     2 1.,     2.,     3.,     4.,     6.,     8. /
      DATA D44B/
     1 -.14,   -.14,   -.14,   -.14,   -.14,   0.,     .44,
     2 -.07,   -.07,   -.07,   -.07,   0.,     .31,    .64,
     3 .05,    .05,    .05,    .06,    .23,    .49,    .75,
     4 .06,    .06,    .06,    .15,    .35,    .58,    .79,
     5 .08,    .08,    .14,    .29,    .48,    .66,    .84,
     6 .09,    .15,    .24,    .39,    .555,   .72,    .89 /
C
C              FIGURE 4.1.3.2-54C
C
      DATA T44C
     1/0.,     2.,     4.,     6.,     8.,     16. /
      DATA D44C
     1  / -.04, .01, .075, .13, .15, .15 /
C
C              FIGURE 4.1.3.4-26A
C
      DATA T419A
     1 /          0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 3.0, 4*0.,
     2   .6, .7, .8, .9, .95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.4 /
      DATA D419A
     1 /              0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,
     2               -.050,-.050,-.050,-.030,-.015,-.010, 0.0 ,
     3               -.115,-.100,-.083,-.047,-.018, 0.0 , .025,
     4               -.145,-.110,-.080,-.020, .022, .050, .097,
     5               -.120,-.080,-.030, .050, .100, .140, .204,
     6               -.060, 0.0 , .060, .150, .210, .250, .325,
     7                .090, .130, .180, .257, .310, .330, .360,
     8                .113, .165, .210, .280, .320, .345, .360,
     9                .108, .160, .205, .284, .320, .344, .360,
     A                .096, .145, .190, .273, .320, .343, .355,
     B                .020, .050, .100, .200, .300, .330, .340/
C
C              FIGURE 4.1.3.4-26C
C
      DATA T419C
     1 / 0.0, 0.4, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 3.0, 5.0  ,3*0.,
     2   .6,.7,.75,.8,.85,.9,.95,1.0,1.05,1.1,1.15,1.2,1.4  /
      DATA D419C
     1 /   0.0, 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0,  0.0 , 0.0 ,
     2   -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,
     3   -1.25,-1.22,-1.20,-1.20,-1.20,-1.10,-1.10,-1.00,-1.00,-0.85,
     4   -2.05,-1.75,-1.50,-1.40,-1.30,-1.20,-1.10,-1.00,-1.00,-0.95,
     5   -2.20,-1.90,-1.60,-1.50,-1.40,-1.20,-1.10,-1.00,-0.90,-0.63,
     6   -1.85,-1.60,-1.40,-1.30,-1.20,-1.00,-0.80,-0.70,-0.60,-0.40,
     7   -1.80,-1.38,-1.00,-0.90,-0.70,-0.50,-0.30,-0.10, 0.0 , 0.10,
     8   -0.60,-0.20, 0.20, 0.30, 0.50, 0.80, 1.00, 1.20, 1.40, 2.00,
     9    1.90, 2.40, 3.00, 3.20, 3.50, 3.90, 4.30, 4.60, 5.00, 6.20,
     A    2.40, 4.05, 5.50, 6.20, 7.00, 8.00, 8.90, 10.0, 10.9, 13.9,
     B    1.60, 4.04, 6.90, 8.00, 9.10, 10.6, 11.9, 13.0, 14.3, 18.2,
     C   -1.20, 3.90, 7.60, 8.90, 10.2, 12.1, 13.5, 14.7, 15.9, 18.8,
     D   -11.0,-0.60, 8.10, 10.5, 12.0, 14.4, 16.2, 17.5, 18.8, 21.3 /
C
C              FIGURE 4.1.5.1-29 LEFT SIDE (SUBSONIC)
C
      DATA T429L/
     10.0,0.4,0.6,0.8,1.0,1.2,1.4,1.475,
     20.,.5,1.0,1.5,2.0,3.0,4.0/
      DATA D429L/
     A8*0.0,
     10.90,0.71,0.55,0.39,0.21,0.10,0.00,0.0,
     22.34,2.03,1.70,1.05,0.51,0.21,0.01,0.0,
     33.01,2.92,2.68,1.99,1.00,0.32,0.01,0.0,
     43.30,3.27,3.12,2.50,1.26,0.42,0.03,0.0,
     53.49,3.43,3.38,2.80,1.52,0.51,0.06,0.0,
     63.61,3.58,3.52,3.14,1.83,0.52,0.09,0.0/
C
C              FIGURE 4.1.5.1-29 RIGHT SIDE (SUPERSONIC)
C
      DATA T429R/
     10.0,.4,.8,1.2,1.4,0.0,0.0,
     20.,.5,1.0,1.5,2.0,3.0,4.0/
      DATA D429R/
     A5*0.0,
     10.90,1.08,1.18,1.28,1.31,
     22.34,2.53,2.68,2.75,2.78,
     33.01,3.08,3.10,3.12,3.12,
     43.30,3.32,3.32,3.30,3.28,
     53.49,3.49,3.48,3.42,3.39,
     63.61,3.61,3.60,3.56,3.51/
C
C              FIGURE 4.2.3.1-24
C
      DATA T424/
     1.75,.85,.90,.95,1.0,1.015,1.030,1.050,1.075,1.125,1.15,1.2,1.3,
     2.011,.067,.122/
      DATA D424/
     1.011,.014,.018,.031,.069,.040,.020,6*0.0,
     2.067,.071,.078,.092,.140,.150,.120,.086, .067, .050,.047,.045,.05,
     3.122,.130,.139,.155,.204,.220,.236,.210, .168, .110, .101,.100,.1/
      DATA T426/
     16.0,7.0,8.0,10.,12.,14.,16.,18.,20.,24.,26.,
     21.0,1.025,1.05,1.1,1.2/
      DATA D426/
     1.118,.084,.065,.041,.029,.020,.016,.012,.010,.007,.006,
     2.147,.110,.088,.061,.043,.033,.026,.022,.019,.015,.013,
     3.166,.120,.097,.068,.050,.039,.031,.026,.022,.015,.013,
     4.186,.138,.110,.076,.057,.044,.034,.028,.022,.015,.013,
     5.200,.152,.120,.080,.059,.046,.036,.028,.022,.015,.013/
C
C  *** WING ALONE ***
C      WING LIFT CURVE SLOPE
C
      MACH=FLC(NZ+2)
      IF(WINGIN(15).NE.WTYPE(1))GO TO 1100
      CALL TRANWG(CLA14, DCLA14)
      K=CLAMO*RAD/(2.0*PI)
      BETA6=0.80
C
C              FIGURE 4.1.3.2-53A FORCE BREAK MACH NUMBER=MFB ( 0 SWEEP)
C
      VAR(1)=TOC*100.
      VAR(2)=ARSTAR
      LGH(1)=8
      LGH(2)=7
      CALL INTERX(2,T43A,VAR,LGH,D43A,MFB0,8,56,
     1            2,2,0,0,2,-1,0,0,Q3243A,3,ROUTID)
C
C              FIGURE 4.1.3.2-53B MFB0 CORRECTED FOR SWEEP
C
      VAR(1)=SWEPE2
      VAR(2)=MFB0
      LGH(1)=8
      LGH(2)=5
      CALL INTERX(2,T43B,VAR,LGH,D43B,MFB,8,40,
     1            0,2,0,0,2,0,0,0,Q3243B,3,ROUTID)
      IF(NF .LT. 0) RETURN
C
C              FIGURE 4.1.3.2-54B (A/C)
C
      VAR(1)=TOC*100.
      VAR(2)=ARSTAR
      LGH(1)=7
      LGH(2)=6
      CALL INTERX(2,T44B,VAR,LGH,D44B,AOC,7,42,
     1            0,2,0,0,1,2,0,0,Q3244B,3,ROUTID)
C
C              FIGURE 4.1.3.2-54A CORRECTION TO LIFT CURVE SLOPE
C
      CALL INTERX(2,T44A,VAR,LGH,D44A,CFBCT,7,42,
     1            2,1,0,0,2,2,0,0,Q3244A,3,ROUTID)
C
      BETAFB = SQRT(1.0-MFB**2)
      IF(MFB .GT. 0.98) BETAFB = 0.0
      BB   = BETAFB
      ARG1 = 2.0*PI*ARSTAR/RAD
      ARG2=(ARSTAR/K)**2
      ARG3=BB**2+TANC2**2
      ARG4 = 2.0+SQRT(ARG2*ARG3+4.0)
      CLAFBT = ARG1/ARG4
C
      BB   = 0.8
      ARG2 = (BB*ARSTAR/K)**2
      ARG3 = 0.0
      IF(ARG2 .GT. 0.0) ARG3 = 1.0+(TANC2/BB)**2
      ARG4 = 2.0+SQRT(ARG2*ARG3+4.0)
      CLAW6 = ARG1*SRSTAR/(ARG4*SR)
C
      XM = 0.75
      IF(XM .GT. (MFB-0.1)) XM = MFB-0.1
      IF(XM .GT. MACH) XM=MACH
      BB   = SQRT(1.0-XM**2)
      ARG2=(ARSTAR/K)**2
      ARG3=BB**2+TANC2**2
      ARG4 = 2.0+SQRT(ARG2*ARG3+4.0)
      CLAW7 = ARG1*SRSTAR/(ARG4*SR)
      DCLA7 = XM*ARG1*ARG2*SRSTAR/SR
      DCLA7 = DCLA7/(ARG4**2*SQRT(ARG2*ARG3+4.0))
      CLAFB=CLAFBT*CFBCT *SRSTAR/SR
      CLAA=(1.-AOC)*CLAFB
C
C              FIGURE 4.1.3.2-54C (B/C)
C
      VAR(1)=TOC*100.
      LGH(1)=6
      CALL INTERX(1,T44C,VAR,LGH,D44C,BOC,6,6,
     1            2,0,0,0,1,0,0,0,Q3244C,3,ROUTID)
      CLAB=(1.-BOC)*CLAFB
      MT(1)=XM
      MT(2)=MFB
      MT(3)=MFB+.07
      MT(4)=MFB+.14
      MT(5)=1.40
      CLAMT(1)=CLAW7
      CLAMT(2)=CLAFB
      CLAMT(3)=CLAA
      CLAMT(4)=CLAB
      CLAMT(5)=CLA14
      CALL TRANF(5, MT, CLAMT, DCLA7, DCLA14, MACH, CLA)
C
C  *** CLASSIFICATION OF ASPECT RATIO ***
C
      CALL TBFUNX(A(27), C1   ,DYDX,NP,X,Y,C,IN,MI,NG,
     1            0,0,Q34171,4,ROUTID)
      ARG1=(C1+1.)*COSLE
      ARATIO=A(128)/ARG1
C
C              FIGURE 4.1.3.4-24B (C1) AND (C2)
C
      VAR(1)=TAPR
      LGH(1)=11
      CALL INTERX(1,TR,VAR,LGH,DR,C2,11,11,
     1            0,0,0,0,0,0,0,0,Q34172,4,ROUTID)
      LOWAR=.FALSE.
      IF(ARSTAR.LT.ARATIO)LOWAR=.TRUE.
C
C  *** IF LOW ASPECT RATIO CALCULATE CLMAX AND (ALPHA)CLMAX ***
C
      IF(.NOT.LOWAR)GO TO 1070
      BU4=(C1+1.)*ARSTAR*COSLE/BETA6
      ARG=(C2+1.)*ARSTAR*TANLE
      A(160)=ARG
      CALL CLMXB1(BU4,CLMAX6,A,WINGIN)
C
C  *** CALCULATE (ALPHA)CLMAX AT MACH=0.60 ***
C              FIGURE 4.1.3.4-25A
C
      CALL TBFUNX(BU4,ACLBAS,DYDX,10,T418A,D418A,C17,I17,MI,NG,
     1            0,0,Q3418A,3,ROUTID)
      IF(ARG.LE.4.5)GO TO 1050
C
C              FIGURE 4.1.3.4-25B
C
      CALL TLINEX(T18B1,D18B1,C18B1,3,20,MT(1),ARG,DACMA6,
     1            0,0,2,1,Q34181,3,ROUTID)
      GO TO 1060
 1050 ARG1=ARSTAR*COSLE*(1.+4.*TAPR**2)
      CALL TLINEX(T18B2,D18B2,C18B2,10,10,ARG1,ARG,DACMA6,
     1            0,2,2,0,Q34182,3,ROUTID)
 1060 ALCLM6=ACLBAS+DACMA6
C
C  *** ADJUST PARAMETERS TO TRANSONIC MACH NUMBER ***
C              FIGURE 4.1.3.4-26B
C
      VAR(1)=ARG
      LGH(1)=13
      CALL INTERX(1,T419B,VAR,LGH,D419B,C3,13,13,
     1            0,0,0,0,0,0,0,0,Q3419B,3,ROUTID)
C
C              FIGURE 4.1.3.4-26C  DELTA (ALPHA)CLMAX
C
      VAR(1)=BU4*BETA6
      VAR(2)=MACH
      LGH(1)=10
      LGH(2)=13
      CALL INTERX(2,T419C,VAR,LGH,D419C,DALCM ,13,130,
     1            2,0,0,0,2,2,0,0,Q3419C,3,ROUTID)
C
C              FIGURE 4.1.3.4-26A DELTA CLMAX
C
      VAR(1)=C3*VAR(1)
      LGH(1)=7
      LGH(2)=11
      CALL INTERX(2,T419A,VAR,LGH,D419A,DCLMAX,11,77,
     1            2,0,0,0,2,2,0,0,Q3419A,3,ROUTID)
      ALCLMT=ALCLM6+DALCM
      CLMAXT=CLMAX6+DCLMAX*(SRSTAR/SR)
 1070 CONTINUE
C
C  *** ZERO ANGLE OF ATTACK DRAG,CD0 ***
C              FIGURE 4.1.5.1-27 (EQUATION FOR RLCOFF)
C
      ARG=12.*CBAR/RUFF
      CALL TBFUNX(0.6,CEPT,DYDX,4,X27M,X27I,C27,I27,MI,NG,
     1            0,0,Q15127,4,ROUTID)
      RLCOFF=ARG**1.0482*10.0**CEPT
      RNN=CBAR*RNFS
      IF(RLCOFF.LT.RNN)RNN=RLCOFF
      CALL FIG26(RNN,0.60,CF)
      RL=1.2
      IF(XOC.LT.0.30)RL=2.0
      CDF=CF*(1.+RL*TOC)*2.*SRSTAR/SR
      ARG=SQRT(COSC4)
      XMTD(1)=0.60/ARG
      XMTD(2)=0.65/ARG
      XMTD(3)=0.70/ARG
      XMTD(4)=0.75/ARG
      XMTD(5)=0.80/ARG
      XMTD(6)=0.85/ARG
      XMTD(7)=0.90/ARG
      XMTD(8)=0.925/ARG
      XMTD(9)=0.950/ARG
      XMTD(10)=0.975/ARG
      XMTD(11)=1.000/ARG
      XMTD(12)=1.050/ARG
      XMTD(13)=1.100/ARG
      XMTD(14)=1.200/ARG
      XMTD(15)=1.400/ARG
      ARG2=TOC**.3333
      VAR(2)=ARSTAR*ARG2
      ARG3=COSC4**2.5
      DO 1090 I=1,15
         ARG4=ABS((XMTD(I)*ARG)**2-1.)
         VAR(1)=SQRT(ARG4)/ARG2
C
C              FIGURE 4.1.5.1-29
C
         IF(XMTD(I)*ARG.LE.1.)GO TO 1080
         LGH(1)=5
         LGH(2)=7
         CALL INTERX(2,T429R,VAR,LGH,D429R,CDW1(I),7,35,
     1              0,0,0,0,0,2,0,0,Q5129R,3,ROUTID)
         GO TO 1090
 1080    LGH(1)=8
         LGH(2)=7
         CALL INTERX(2,T429L,VAR,LGH,D429L,CDW1(I),8,56,
     1               0,0,0,0,-1,2,0,0,Q5129L,3,ROUTID)
         CDW2(I)=CDW1(I)*TOC**1.666*ARG3*SRSTAR/SR
 1090 CONTINUE
      CALL TRANF(15, XMTD, CDW2, 0.0, 0.0, MACH, CDW)
      CD0W=CDW+CDF
C
C  *** BODY ALONE ***
C
 1100 IF(.NOT.BO)RETURN
C
C  ***BODY LIFT AND MOMENT SLOPES ***
C
      ARG1=(MACH-.60)/.80
      CLABD=CLABM6+(CLAB14-CLABM6)*ARG1
      CMAB=CMAB6+(CMAB14-CMAB6)*ARG1
C
C  *** BODY DRAG ***
C
      CDFB=CFB6*SS/SR
      CDPB=CD0B6-CDBB6-CDFB
      IF(MACH.GE.1.0.AND.MACH.LE.1.2)CDPB=CDPB*(1.-(MACH-1.0)/.2)
      IF(MACH.GT.1.2)CDPB=0.0
      IF(DB.LT.0.3*DMAX)DB=0.3*DMAX
      ARG=(DB/DMAX)**2
C
C              FIGURE 4.2.3.1-24 (CDBT)
C
      VAR(1)=MACH
      VAR(2)=CDBB6/ARG*4.0*SR/(PI*DMAX**2)
      LGH(1)=13
      LGH(2)=3
      CALL INTERX(2,T424,VAR,LGH,D424,CDBFIG,13,39,
     1            0,2,0,0,0,2,0,0,Q23124,3,ROUTID)
      CDBB=CDBFIG*ARG*PI*DMAX**2/(4.0*SR)
      IF(MACH.LT.1.0)GO TO 1110
      VAR(1)=XLB/DMAX
      VAR(2)=MACH
      LGH(1)= 11
      LGH(2)=5
      CALL INTERX(2,T426,VAR,LGH,D426,CDWB,11,55,
     1            0,0,0,0,2,2,0,0,Q23126,3,ROUTID)
      CDWB=CDWB*PI*DMAX**2/(4.*SR)
      GO TO 1120
 1110 CDWB=0.0
 1120 CD0B=CDFB+CDPB+CDBB+CDWB
      IF(MACH.GT.1.2)GO TO 1130
      GO TO 1140
 1130 CD0B=CD0B+(CD014-CD0B)*(MACH-1.2)/0.2
C
C  *** BODY DRAG AT ANGLE OF ATTACK ***
C
 1140 DO 1150 J=1,NALPHA
         CDJB(J)=CD0B+(FLC(J+22)/RAD)**2*SB/SR
 1150 CONTINUE
      CD0WB=CD0B+CD0W
      TRA(73)=CD0WB
      RETURN
      END
