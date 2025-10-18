      SUBROUTINE SUPLNG
C
C***  CALCULATES SUPERSONIC WING CL, CLA AND CMA
C
      DIMENSION ROUTID(2),Q3250A(3),Q3250B(3),Q3246A(3),Q3246G(3),
     1 Q3359A(3),Q3359B(3),Q3361A(3),Q3360A(3),Q246AF(4),
     2 Q13252(3),Q13253(3),Q13251(3),Q222AF(4),Q15127(3),Q15258(3)
      REAL MACH,KSHARP
      REAL LHS,LERBW,LERI,LERO
      DIMENSION T13246(52),D13246(1104),VAA(4),LGA(4),VAB(4),LGB(4),
     1 A1350(17), B1350(19),DA50(72),DB50(88),G13246(9),DG3246(9),
     2A13359(22),D13359(99),ALSCHD(20),BR3359(21),DR3359(110),BL3359(20)
     3,DL3359(50),ALPHAJ(20),CNAAA(20),T61A(9),
     4SLOPE(9)
      DIMENSION DS5258(5),DR5258(5),T15258 (5)
      DIMENSION R13252(12),DRND52(12),S13252(12),DSHP52(12),T13251(12),
     1D13251(12)
      DIMENSION WTYPE(4)
      DIMENSION CAF(20)
      DIMENSION CN(20),CL(20),DUMY4(184),DR360A(156),SUPAF(216),SUBAF(21
     16)
      DIMENSION DUMYC(108),DUMYD(108)
      DIMENSION TL360A(24),DL360A(72),TR360A(25),T13253(24),D13253(140),
     1T422AF(18),DUMY7(65),DUMY8(91),DUMYA(108),DUMYB(108)
      DIMENSION T455(10),D455(10),CD(20),CDL(20)
      DIMENSION DUMY1(184),DUMY2(184),DUMY3(184),DUMY5(184),DUMY6(184)
      DIMENSION EPSLON(20),CM(20)
      DIMENSION X27M(4),X27I(4),C27(6)
      COMMON /OVERLY/ NLOG,NMACH,MIDX,NALPHA
      COMMON /FLGTCD/ FLC(95)
      COMMON /SYNTSS/ SYNA(19)
      COMMON /OPTION/ SW,CBARR,RUFF,BLREF
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /SUPWH/  SLG(141)
      COMMON /WINGD/  A(195)
      COMMON /WINGI/  HTIN(100)
C
C********* HTIN CORRESPONDS TO WINGIN ELSEWHERE ***********************
C
      COMMON /IWING/ PWING,HT(400)
C
C********* HT CORRESPONDS TO WING ELSEWHERE ***************************
C
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1    HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,TRANSN,HYPERS
      LOGICAL       FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1   HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,TRANSN,HYPERS
      EQUIVALENCE (CLADEG,HT(101)),(CMADEG,HT(121)),(CAF(1),HT(81))
      EQUIVALENCE (RKL,SLG(133)),(S2,SLG(74)),(A2,SLG(73))
      EQUIVALENCE  (CNAAAP,SLG(75))
      EQUIVALENCE(CM,HT(41))
      EQUIVALENCE (SRSTAR,A(3)),(P,SLG(82)),(DRAGC,SLG(81)),
     1  (CDO,SLG(80)),(CDL(1),SLG(53)),(CD(1),HT(1))
      EQUIVALENCE (SBW,SLG(119)),(CNTBW,SLG(77)),(XACCRI,SLG(76))
      EQUIVALENCE (TAPEXP,A(27)),(SO,A(167))
      EQUIVALENCE  (CBP,A(166)),              (BSTARB,A(23))
      EQUIVALENCE (TANLE,A(62)),(BETA,SLG(1)),(SWEPLE,A(58)),(AR,A(7  ))
     1,(TRATIP,A(25)),(BOVERT,SLG(2)),(CNNNT,SLG(3)),(COSLE,A(61)),
     2(DELTAY,HTIN(17)),(BCNA,SLG(4)),(DELTYT,SLG(8)),(DELTDT,SLG(9)),
     3(E,SLG(11)),(CC,SLG(12)),(ALSCHD(1),FLC(23)),(ALPHAI,SYNA(4)),
     4(ALPHAJ(1),SLG(33)),(CNAAA(1),SLG(13)),(SINLE,A(60))
     5,(DETACH,SLG(94)),(CNAAST,SLG(116)),(DETALP,SLG(117))
      EQUIVALENCE (CNTHRY,SLG(5)),(CNA1,SLG(7)),(TLE192,SLG(10)),
     1(RMACH,SLG(93)),(CN(1),HT(61)),(CL(1),HT(21)),(DETANG,SLG(115)),
     2(XACCRO,SLG(78))
      EQUIVALENCE (SPAN,HTIN(4)),(SPANO,HTIN(2)),(CRBW,SLG(118))
     1,(TANTEO,A(104)),(CB,HTIN(5)),(CT,HTIN(1)),(TAPBW,SLG(121))
     2,(TANLEI,A(62)),(TANTEI,A(80)),(CN1,SLG(127)),(CN2,SLG(128)),
     3(ARBW,SLG(120)),(CLEBW,SLG(122)),(CRGLV,SLG(123)),(SGLV,SLG(124))
     4,(ARGLV,SLG(125)),(BE,SLG(126)),(CNAE,SLG(129)),(CNAGLV,SLG(130))
     5,(CNABW,SLG(131)),(CLEGLV,SLG(132)),(TANLEO,A(86)),(DETXCG,A(173))
      EQUIVALENCE (ARSTAR,A(7)),(XACCR,SLG(134)),(DCMCL,SLG(135)),
     1(CRSTAR,A(10)),(CMA,SLG(136)),(ARIP,A(5))
     2 ,(TAPIP,A(26)),(CNCNTI,SLG(137)),(TAPOP,A(169)),(COSLEO,A(85)),
     3(AROP,A(168)),(SI,A(1)),          (CNCNTO,SLG(138)),(CNATI,SLG(139
     4)),(CNATO,SLG(140)),(CNAI,SLG(92)),(CNAO,SLG(91))
      EQUIVALENCE (D13246(1),DUMY1(1)),(D13246(185),DUMY2(1)),(D13246(36
     19),DUMY3(1)),(D13246(553),DUMY4(1)),(D13246(737),DUMY5(1)) ,
     2(D13246(921),DUMY6(1))
      EQUIVALENCE (SUBAF(1),DUMYA(1)),(SUBAF(109),DUMYB(1))
      EQUIVALENCE (SUPAF(1),DUMYC(1)),(SUPAF(109),DUMYD(1))
      EQUIVALENCE (DR360A(1),DUMY7(1)),(DR360A(66),DUMY8(1))
      EQUIVALENCE (SWTEI,A(76)),(SWTEO,A(100))
      EQUIVALENCE (LERO,HTIN(63)),(LERI,HTIN(62)),(KSHARP,HTIN(71))
      EQUIVALENCE (SPANS,HTIN(3)),(XW,SYNA(2)),(XCG,SYNA(1))
      EQUIVALENCE (SWEPLO,A(82))
      EQUIVALENCE (SLG(141),RKT)
      LOGICAL CNPRME
      LOGICAL SUPLE,ROUND,DETACH, GLOVE,PITCH  ,TRANS
      DATA ROUTID/4HSUPL,4HNG  /,Q3250A/4H4.1.,4H3.2-,4H60A /,
     1                           Q3250B/4H4.1.,4H3.2-,4H60B /,
     2 Q3246A/4H4.1.,4H3.2-,4H56A /,Q3359A/4H4.1.,4H3.3-,4H59A /,
     3 Q3246G/4H4.1.,4H3.2-,4H56G /,Q3359B/4H4.1.,4H3.3-,4H59B /,
     4 Q3361A/4H4.1.,4H3.3-,4H61A /,
     5 Q3360A/4H4.1.,4H3.3-,4H60A /,Q246AF/4H4.1.,4H3.2-,4H56A-,4HF   /,
     6 Q13252/4H4.1.,4H3.2-,4H62  /,Q13253/4H4.1.,4H3.2-,4H63  /,
     7 Q15127/4H4.1.,4H5.1-,4H27  /,Q13251/4H4.1.,4H3.2-,4H61  /,
     8 Q15258/4H4.1.,4H5.2-,4H58  /,Q222AF/4H4.1.,4H4.2-,4H26A-,4HF   /
C
C                  FIGURE 4.1.5.1-27
C
      DATA X27M /0.0,1.0,2.0,3.0/
      DATA X27I /1.57780,1.67221,1.98509,2.28874/
      DATA I27  /0/
      DATA WTYPE/4HSTRA,4HDOUB,4HCRAN,4HCURV/
C
C                FIGURE 4.1.3.2-65
C
      DATA T455
     1/0.,.2,.4,.8,1.,1.25,1.6667,2.5,5.,10./
      DATA D455
     1/1.,1.27,1.58,2.28,2.70,2.56,2.49,2.44,2.405,2.4/
C
C              FIGURE 4.1.3.2-60A (CNA/CNA)T FOR SUBSONIC L.E.
C
      DATA A1350/
     1 0.,     .3,     .4,     .5,     .6,     .7,     .8,     .9,
     2 1.0001 ,
     3 0.,     .41,    .82,    1.24,   2.12,   3.18,   6.95,   16.1 /
      DATA DA50/
     1 1.,     1.,     1.,     1.,     1.,     1.,     1.,     1.,
     2         1.,
     3 1.05,   1.05,   1.05,   1.05,   .985,   .945,   .915,   .9,
     4 .9,
     5 1.04,   1.04,   1.04,   1.04,   .965,   .908,   .87,    .85,
     6 .84,
     7 1.12,   1.12,   1.12,   1.015,  .94,    .88,    .838,   .81,
     8 .796,
     9 1.11,   1.11,   1.11,   1.,     .903,   .84,    .795,   .765,
     A .75,
     B 1.08,   1.08,   1.08,   .954,   .865,   .8,     .75,    .72,
     C .7,
     D 1.2,    1.2,    1.043,  .907,   .817,   .75,    .707,   .675,
     E .66,
     F 1.14,   1.14,   .975,    .857,  .772,   .717,   .675,   .65,
     G .632 /
C
C              FIGURE 4.1.3.2-60B (CNA/CNA)T FOR SUPERSONIC L.E.
C
      DATA B1350/
     1 0.,     .1,     .2,      .3,    .4,     .5,     .6,     .7,
     2 .8,     .9,     1.0001 ,
     3 0.,     4.,     8.,      12.,   20.,    30.,    50.,    70./
      DATA DB50/
     1 1.,     1.,     1.,      1.,    1.,     1.,     1.,     1.,
     2 1.,     1.,     1.,
     3 1.15,   1.15,   1.15,    1.15,  1.095,  1.04,   .99,    .96,
     4 .935,   .915,   .9,
     5 1.12,   1.12,   1.12,    1.12,  1.05,   .985,   .94,    .905,
     6 .88,    .857,   .84,
     7 1.15,   1.15,   1.15,    1.08,  1.005,  .945,   .902,   .87,
     8 .842,   .82,    .796,
     9 1.22,   1.14,   1.05,    .98,   .93,    .89,    .853,   .823,
     A .795,   .77,    .75,
     B 1.13,   1.05,   .98,    .925,   .88,    .845,   .81,    .782,
     C .752,   .73,    .7,
     D 1.02,   .942,   .895,   .855,   .82,    .79,    .76,    .735,
     E .71,    .685,   .66,
     F 1.,     .92,    .87,    .825,   .79,    .755,   .728,   .7,
     G .678,   .655,   .632 /
C
C              FIGURE 4.1.3.2-56G (CNA/A) FOR A*BETA LESS THAN 1.0
C
      DATA G13246/0.,.2,.4,.45,.5,.6,.8,.9,1.0/
      DATA DG3246/1.61,1.58,1.55,1.57,1.62,1.75,1.94,2.0,2.0/
C
C              FIGURE 4.1.3.2-61 (KL)
C
      DATA T13251/.015,.03,.05,.075,.11,.16,.23,.45,.90,1.2,1.7,2.0/
      DATA D13251/.963,.94,.92,.90,.88,.86,.84,.80,.76,.74,.72,.709/
C
C              FIGURE 4.1.5.2-58
C
      DATA T15258/0.0,.28,.40,.50, 1.87/
      DATA DS5258/.54,.54,.559,.60,1.95/
      DATA DR5258/.54,.54,.593,.68,2.0/
C
C              FIGURE 4.1.3.2-56A-F (WING SUPERSONIC CNA)
C
      DATA T13246
     1  / 0.,.1,.2,.3,.33,.4,.5,.6,.7,.8,.9,1.0,1.111,1.25,1.429,1.667,
     2    2.0,2.5,2.941,4.167,7.143,14.286, 30. ,
     3    .25,.5,1.,2.,3.,4.,5.,6. , 15*0.,
     4    0.,.2,.25,.3333,.5,1.0/
C
C     FIG 4.1.3.2-56A
C
      DATA DUMY1
     1  /  .39, .39, .39, .39, .39, .39, .39, .39, .39, .39, .39, .41,
     2     .44, .5 , .58, .65, .80,1.0 ,1.15,1.55,2.40,3.92,4.0 ,
     3     .77, .77, .78, .78, .79, .79, .80, .80, .80, .81, .82, .84,
     4     .92,1.0 ,1.18,1.32,1.60,1.9 ,2.18,2.80,3.84,3.92,4.0 ,
     5    1.55,1.56,1.57,1.57,1.59,1.59,1.60,1.63,1.66,1.68,1.70,1.75,
     6    1.88,2.09,2.30,2.6 ,2.94,3.35,3.70,3.80,3.98,3.98,4.0 ,
     7    3.15,3.15,3.15,3.15,3.15,3.15,3.17,3.19,3.23,3.27,3.33,3.4 ,
     8    3.46,3.54,3.60,3.7 ,3.75,3.8 ,3.86,3.91,4.0 ,4.0 ,4.0 ,
     9    4.71,4.74,4.83,5.09,5.25,5.05,4.80,4.55,4.30,4.09,3.90,3.72,
     A    3.78,3.8 ,3.83,3.88,3.90,3.94,3.92,3.96,4.0 ,4.0 ,4.0 ,
     B    6.29,6.2 ,5.99,5.72,5.61,5.42,5.18,4.9 ,4.64,4.42,4.2 ,4.0 ,
     C    4.0 ,4.0 ,4.0 ,4.0 ,4.0 ,4.0 ,4.0 ,4.0 ,4.0 ,4.0 ,4.0 ,
     D    6.3 ,6.34,6.39,6.08,5.97,5.78,5.5 ,5.2 ,4.95,4.7 ,4.48,4.25,
     E    4.2 ,4.18,4.14,4.11,4.07,4.05,4.04,4.02,4.  ,4.  ,4.  ,
     F    6.32,6.4 ,6.4 ,6.34,6.32,6.13,5.82,5.51,5.25,4.99,4.73,4.5 ,
     G    4.4 ,4.32,4.27,4.2 ,4.13,4.11,4.07,4.04,4.  ,4.  ,4.  /
C
C     FIG 4.1.3.2-56B
C
      DATA DUMY2
     1  /  .41, .41, .4 , .4 , .4 , .4 , .4 , .4 , .4 , .4 , .4 , .41,
     2     .41, .48, .56, .67, .84,1.08,1.29,1.83,2.80,3.8 ,4.0 ,
     3     .8 , .8 , .8 , .8 , .8 , .8 , .8 , .8 , .8 , .8 , .83, .85,
     4     .98,1.09,1.27,1.48,1.80,2.25,2.53,3.22,3.79,3.92,4.  ,
     5    1.57,1.59,1.6 ,1.6 ,1.6 ,1.61,1.63,1.65,1.69,1.73,1.77,1.8 ,
     6    2.  ,2.26,2.57,2.9 ,3.27,3.53,3.65,3.8 ,3.91,3.97,4.  ,
     7    3.17,3.17,3.22,3.4 ,3.5 ,3.6 ,3.7 ,3.67,3.58,3.49,3.37,3.23,
     8    3.38,3.52,3.62,3.71,3.80,3.88,3.9 ,3.95,3.97,3.99,4.0 ,
     9    4.72,5.  ,5.1 ,5.02,4.97,4.84,4.64,4.45,4.24,4.05,3.89,3.70,
     A    3.84,3.91,3.97,3.99,4.0 ,4.  ,4.  ,4.  ,3.99,4.  ,4.  ,
     B    5.57,5.6 ,5.62,5.61,5.6 ,5.42,5.15,4.91,4.7 ,4.48,4.25,4.08,
     C    4.13,4.18,4.16,4.11,4.09,4.05,4.02,4.01,4.  ,4.  ,4.  ,
     D    5.73,5.77,5.79,5.79,5.77,5.75,5.6 ,5.33,5.08,4.83,4.59,4.4 ,
     E    4.42,4.4 ,4.32,4.23,4.17,4.08,4.05,4.03,4.01,4.  ,4.  ,
     F    5.83,5.88,5.91,5.93,5.93,5.91,5.88,5.68,5.4 ,5.14,4.9 ,4.7 ,
     G    4.68,4.59,4.46,4.33,4.23,4.13,4.1 ,4.06,4.02,4.01,4.  /
C
C     FIG 4.1.3.2-56C
C
      DATA DUMY3
     1  /  .4 , .4 , .4 , .4 , .4 , .4 , .4 , .4 , .4 , .4 , .4 , .42,
     2     .43, .5 , .6 , .71, .9 ,1.12,1.35,1.97,2.95,3.80,4.  ,
     3     .8 , .8 , .8 , .8 , .8 , .8 , .8 , .81, .82, .83, .83, .88,
     4     .95,1.06,1.23,1.46,1.83,2.3 ,2.65,3.32,3.76,3.92,4.  ,
     5    1.59,1.59,1.59,1.59,1.59,1.6 ,1.62,1.65,1.68,1.72,1.78,1.82,
     6    2.  ,2.29,2.65,3.  ,3.27,3.5 ,3.6 ,3.78,3.90,3.98,4.  ,
     7    3.14,3.2 ,3.41,3.62,3.68,3.75,3.7 ,3.62,3.51,3.42,3.3 ,3.18,
     8    3.34,3.5 ,3.6 ,3.7 ,3.8 ,3.85,3.88,3.93,3.98,4.0 ,4.  ,
     9    4.7 ,4.85,5.02,4.98,4.91,4.8 ,4.6 ,4.4 ,4.2 ,4.  ,3.82,3.68,
     A    3.8 ,3.9 ,3.98,4.  ,4.  ,4.  ,4.  ,4.  ,4.  ,4.  ,4.  ,
     B    5.4 ,5.44,5.45,5.45,5.42,5.35,5.12,4.9 ,4.69,4.48,4.22,4.02,
     C    4.12,4.18,4.17,4.15,4.1 ,4.05,4.04,4.01,4.  ,4.  ,4.  ,
     D    5.61,5.63,5.64,5.62,5.61,5.6 ,5.54,5.3 ,5.02,4.8 ,4.6 ,4.39,
     E    4.42,4.44,4.38,4.26,4.16,4.1 ,4.07,4.03,4.01,4.  ,4.  ,
     F    5.72,5.75,5.76,5.78,5.77,5.76,5.72,5.69,5.4 ,5.15,4.9 ,4.68,
     G    4.68,4.6 ,4.49,4.33,4.22,4.15,4.10,4.06,4.02,4.01,4.  /
C
C      FIG 4.1.3.2-56D
C
      DATA DUMY4
     1  /  .41, .41, .41, .41, .41, .41, .41, .41, .41, .41, .41, .41,
     2     .43, .50, .6 , .71, .89,1.11,1.35,2.00,3.00,3.72,4.00,
     3     .82, .82, .81, .8 , .8 , .81, .82, .83, .84, .85, .87, .89,
     4     .98,1.1 ,1.3 ,1.5 ,1.82,2.23,2.68,3.3 ,3.7 ,3.91,4.  ,
     5    1.6 ,1.59,1.58,1.58,1.58,1.59,1.60,1.62,1.66,1.73,1.81,1.92,
     6    2.2 ,2.45,2.7 ,2.98,3.22,3.45,3.59,3.75,3.7 ,3.91,4.  ,
     7    3.13,3.18,3.32,3.6 ,3.64,3.72,3.70,3.64,3.54,3.45,3.3 ,3.14,
     8    3.32,3.48,3.6 ,3.7 ,3.79,3.87,3.9 ,3.95,3.99,4.  ,4.  ,
     9    4.71,4.76,4.8 ,4.83,4.84,4.7 ,4.53,4.34,4.18,4.  ,3.82,3.63,
     A    3.78,3.9 ,3.98,4.  ,4.  ,4.  ,4.  ,4.  ,4.  ,4.  ,4.  ,
     B    5.2 ,5.22,5.22,5.2 ,5.2 ,5.15,5.1 ,4.88,4.63,4.42,4.21,4.0 ,
     C    4.15,4.21,4.21,4.17,4.14,4.07,4.05,4.03,4.01,4.  ,4.  ,
     D    5.45,5.47,5.47,5.45,5.45,5.41,5.36,5.3 ,5.02,4.81,4.58,4.36,
     E    4.45,4.47,4.41,4.30,4.19,4.12,4.08,4.02,4.  ,4.  ,4.  ,
     F    5.58,5.59,5.59,5.59,5.59,5.58,5.57,5.53,5.41,5.14,4.91,4.65,
     G    4.7 ,4.65,4.53,4.4 ,4.27,4.16,4.12,4.05,4.02,4.01,4.  /
C
C     FIG 4.1.3.2-56E
C
      DATA DUMY5
     1 /   .4 , .4 , .4 , .4 , .4 , .4 , .4 , .4 , .4 , .4 , .4 , .4 ,
     2     .44, .52, .63, .78, .98,1.21,1.4 ,2.07,3.07,3.6 ,4.  ,
     3     .8 , .8 , .79, .79, .79, .79, .8 , .81, .83, .85, .88, .9 ,
     4    1.01,1.17,1.38,1.62,2.00,2.50,2.74,3.22,3.66,3.87,4.  ,
     5    1.58,1.59,1.62,1.71,1.76,1.77,1.8 ,1.85,1.88,1.91,1.92,1.93,
     6    2.13,2.4 ,2.6 ,2.88,3.12,3.39,3.51,3.71,3.87,3.99,4.  ,
     7    3.15,3.21,3.38,3.61,3.7 ,3.67,3.59,3.5 ,3.4 ,3.29,3.17,3.03,
     8    3.18,3.37,3.51,3.67,3.78,3.85,3.88,3.93,3.98,4.  ,4.  ,
     9    4.42,4.42,4.41,4.40,4.39,4.38,4.35,4.25,4.08,3.9 ,3.72,3.57,
     A    3.70,3.86,3.93,3.99,4.  ,4.  ,4.  ,4.  ,4.  ,4.  ,4.  ,
     B    4.88,4.85,4.81,4.80,4.79,4.77,4.73,4.68,4.53,4.35,4.15,3.99,
     C    4.08,4.20,4.23,4.19,4.1 ,4.07,4.05,4.01,4.  ,4.  ,4.  ,
     D    5.08,5.09,5.08,5.08,5.07,5.06,5.03,4.99,4.92,4.74,4.53,4.3 ,
     E    4.43,4.48,4.43,4.32,4.20,4.12,4.08,4.03,4.01,4.  ,4.  ,
     F    5.19,5.21,5.22,5.23,5.23,5.23,5.23,5.21,5.19,5.09,4.88,4.6 ,
     G    4.72,4.69,4.58,4.41,4.28,4.18,4.12,4.07,4.02,4.01,4.  /
C
C     FIG 4.1.3.2-56F
C
      DATA DUMY6
     1 /   .4 , .4 , .4 , .4 , .4 , .41, .44, .49, .51, .54, .59, .61,
     2     .69, .75, .81, .9 ,1.01,1.26,1.53,2.1 ,2.9 ,3.42,4.  ,
     3     .81, .81, .81, .83, .84, .86, .88, .91, .99,1.02,1.1 ,1.19,
     4    1.3 ,1.42,1.58,1.75,2.  ,2.34,2.6 ,3.  ,3.46,3.71,4.   ,
     5    1.6 ,1.58,1.59,1.6 ,1.6 ,1.61,1.67,1.7 ,1.79,1.86,1.92,2.  ,
     6    2.13,2.32,2.51,2.75,2.97,3.2 ,3.35,3.58,3.8 ,3.9 ,4.  ,
     7    3.13,3.1 ,3.08,3.05,3.04,3.03,3.01,2.99,2.95,2.93,2.91,2.89,
     8    3.02,3.2 ,3.37,3.52,3.67,3.77,3.82,3.88,3.92,3.97,4.  ,
     9    3.79,3.84,3.86,3.8 ,3.79,3.72,3.63,3.58,3.5 ,3.47,3.43,3.41,
     A    3.57,3.75,3.89,3.98,4.  ,3.98,3.94,3.93,3.92,3.99,4.  ,
     B    4.12,4.2 ,4.19,4.09,4.06,4.03,3.99,3.95,3.92,3.9 ,3.88,3.87,
     C    4.01,4.1 ,4.21,4.21,4.12,4.05,4.  ,3.96,3.99,4.  ,4.  ,
     D    4.39,4.45,4.44,4.39,4.35,4.31,4.27,4.25,4.25,4.25,4.25,4.25,
     E    4.37,4.50,4.49,4.39,4.23,4.1 ,4.03,3.99,4.  ,4.  ,4.  ,
     F    4.6 ,4.68,4.65,4.59,4.57,4.54,4.5 ,4.45,4.48,4.53,4.59,4.61,
     G    4.71,4.75,4.69,4.5 ,4.31,4.18,4.08,4.01,4.  ,4.  ,4.  /
C
C              FIGURE 4.1.3.3-59A
C
      DATA A13359
     1 / 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0 , 2*0.,
     2   .7,.8,.9,1.0,1.1,1.25,1.5,2.0,2.5,3.0,3.25  /
      DATA D13359
     1 / 0.,3.50,4.95,5.85,6.40,6.83,7.25,7.57,7.85,
     2   0.,3.00,4.19,4.95,5.58,6.04,6.49,6.80,7.00,
     3   0.,2.34,3.40,4.00,4.49,4.91,5.28,5.63,5.92,
     4   0.,1.86,2.78,3.20,3.53,3.90,4.25,4.60,4.85,
     5   0.,1.42,2.20,2.68,3.05,3.34,3.62,3.87,4.06,
     6   0.,1.10,1.62,1.94,2.19,2.40,2.64,2.80,2.95,
     7   0.,0.92,1.35,1.65,1.84,2.00,2.17,2.28,2.35,
     8   0.,0.73,1.10,1.35,1.54,1.67,1.76,1.83,1.85,
     9   0.,0.62,0.94,1.00,1.10,1.12,1.12,1.10,1.09,
     A   0.,0.54,0.72,0.68,0.42,0.0 ,0.0 ,0.0, 0.0 ,
     B   0.,0.45,0.63,0.51,0.23,0.0 ,0.0 , 0.0 ,0.0 /
C
C              FIGURE 4.1.3.3-59B (LEFT SIDE)
C
      DATA BL3359
     1 / 0.5,0.6,0.7,0.8,1.0,5*0,
     2   0.6,0.8,1.0,1.2,1.6,2.0,2.4,2.6,2.8,3.0/
      DATA DL3359
     1 / 2.03,1.95,1.90,1.86,1.78,
     2   1.94,1.87,1.82,1.77,1.70,
     3   1.82,1.75,1.69,1.64,1.57,
     4   1.71,1.63,1.58,1.54,1.46,
     5   1.51,1.44,1.37,1.32,1.24,
     6   1.32,1.24,1.16,1.10,0.99,
     7   0.98,0.89,0.82,0.75,0.65,
     8   0.58,0.52,0.47,0.42,0.37,
     9   0.18,0.15,0.12,0.13,0.17,
     A   0.00,0.00,0.00,0.00,0.00/
C
C              FIGURE 4.1.3.3-59B (RIGHT SIDE)
C
      DATA BR3359
     1 / 0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0 ,
     2   0.6,0.8,1.0,1.2,1.6,2.0,2.4,2.6,2.8,3.0/
      DATA DR3359
     1 / 7*1.69,1.7,1.72,1.75,1.77,
     2   1.69,1.66,1.64,1.62,1.60,1.60,1.61,1.63,1.65,1.67,1.69,
     3   1.69,1.64,1.59,1.56,1.52,1.51,1.51,1.51,1.52,1.54,1.57,
     4   1.69,1.61,1.55,1.49,1.44,1.41,1.39,1.40,1.41,1.43,1.46,
     5   1.69,1.56,1.45,1.36,1.28,1.22,1.18,1.17,1.18,1.20,1.23,
     6   1.69,1.52,1.34,1.20,1.08,0.99,0.94,0.91,0.91,0.94,0.98,
     7   1.69,1.40,1.14,0.94,0.80,0.70,0.65,0.62,0.61,0.63,0.65,
     8   1.69,1.32,1.08,0.89,0.74,0.64,0.54,0.47,0.42,0.38,0.37,
     9   1.69,1.32,1.04,0.84,0.67,0.54,0.45,0.34,0.27,0.20,0.17,
     A   1.69,1.31,1.02,0.79,0.61,0.46,0.33,0.21,0.12,0.03,0.0/
C
C                FIGURE 4.1.3.3-60A (LEFT SIDE)
C
      DATA TL360A
     1 /  0, 0.2, 0.4, 0.6, 0.8, 1.0   , 6*0.,
     2   1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.26, 2.4, 2.6, 2.8, 3.0 /
      DATA DL360A
     1 / 1.60, 1.35, 1.12, 0.94, 0.79, 0.68,
     2   1.60, 1.35, 1.12, 0.94, 0.79, 0.68,
     3   1.60, 1.35, 1.12, 0.94, 0.79, 0.68,
     4   1.60, 1.35, 1.12, 0.94, 0.79, 0.68,
     5   1.60, 1.35, 1.12, 0.94, 0.79, 0.68,
     6   1.60, 1.35, 1.12, 0.94, 0.79, 0.68,
     7   1.60, 1.35, 1.12, 0.94, 0.79, 0.68,
     8   1.60, 1.35, 1.12, 0.94, 0.79, 0.68,
     9   1.48, 1.25, 1.03, 0.85, 0.70, 0.60,
     A   0.84, 0.63, 0.50, 0.42, 0.37, 0.37,
     B   0.34, 0.25, 0.17, 0.13, 0.11, 0.16,
     C   0.0,  0.0,  0.0,  0.0,  0.0,  0.0    /
C
C                FIGURE 4.1.3.3-60A (RIGHT SIDE)
C
      DATA TR360A
     1 / 0.0, 0.2, 0.220, 0.24, 0.270, 0.3, 0.32, 0.375, 0.4, 0.471,
     2   0.6, 0.8, 1.0  ,
     3   1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.26, 2.4, 2.6, 2.8, 3.0 /
      DATA DUMY7
     1 / 1.68, 1.46, 1.45, 1.40, 1.31, 1.23, 1.16, 1.05, 1.00, 0.88,
     2                                             0.72, 0.64, 0.68,
     3   1.68, 1.41, 1.40, 1.40, 1.31, 1.23, 1.16, 1.05, 1.00, 0.88,
     4                                             0.72, 0.64, 0.68,
     5   1.68, 1.37, 1.35, 1.33, 1.31, 1.23, 1.16, 1.05, 1.00, 0.88,
     6                                             0.72, 0.64, 0.68,
     7   1.68, 1.30, 1.26, 1.24, 1.21, 1.18, 1.16, 1.05, 1.00, 0.88,
     8                                             0.72, 0.64, 0.68,
     9   1.68, 1.24, 1.20, 1.16, 1.14, 1.11, 1.10, 1.05, 1.00, 0.88,
     A                                             0.72, 0.64, 0.68/
      DATA DUMY8/
     1   1.68, 1.16, 1.11, 1.10, 1.07, 1.04, 1.02, 0.96, 0.94, 0.88,
     2                                             0.72, 0.64, 0.68,
     3   1.68, 1.11, 1.07, 1.03, 1.00, 0.97, 0.95, 0.89, 0.87, 0.81,
     4                                             0.72, 0.64, 0.68,
     5   1.68, 1.10, 1.03, 1.02, 0.98, 0.95, 0.93, 0.87, 0.85, 0.79,
     6                                             0.72, 0.64, 0.68,
     7   1.68, 1.07, 1.00, 0.98, 0.95, 0.90, 0.87, 0.81, 0.79, 0.72,
     8                                             0.63, 0.56, 0.60,
     9   1.68, 1.04, 0.98, 0.94, 0.90, 0.86, 0.83, 0.76, 0.74, 0.66,
     A                                             0.54, 0.42, 0.37,
     B   1.68, 1.03, 0.97, 0.92, 0.86, 0.77, 0.76, 0.69, 0.65, 0.57,
     C                                             0.42, 0.25, 0.16,
     D   1.68, 1.02, 0.96, 0.90, 0.84, 0.775,0.75, 0.67, 0.63, 0.52,
     E                                             0.34, 0.10, 0.0  /
C
C              FIGURE 4.1.3.2-62 LEADING-EDGE-EFFECT FACTORS FOR NORMAL-
C                                FORCE-CURVE SLOPE AT SUPERSONIC SPEEDS
C              FIGURE 4.1.3.2-62 (ROUND LEADING EDGE)
C
      DATA R13252
     1 / 0., .5, .94, 1., 1.05, 1.5, 2., 2.5, 3., 3.5, 4., 20./
      DATA DRND52
     1 / 1.043, .97, .893, .89, .894, .987, 1.005, 1.02, 1.025, 1.03,
     2 1.034, 1.034 /
C
C              FIGURE 4.1.3.2-62 (SHARP LEADING EDGE)
C
      DATA S13252
     1 / 0., .5, .94, 1., 1.05, 1.5, 2., 2.5, 3., 3.5, 4., 20./
      DATA DSHP52
     1 / 1.08, 1.013, .94, .942, .948, 1.005, 1.03, 1.04, 1.048, 1.05,
     2 1.055, 1.055 /
C
C              FIGURE 4.1.3.2-63 WING SUPERSONIC NORMAL-FORCE-CURVE
C              SLOPE, TAPER RATIO ZERO
C
      DATA T13253
     1  / .0,.2,.4,.6,.8,1.,1.2,1.4,1.6,1.8,2.,2.4,2.8,3.2 ,
     2     -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8 /
      DATA D13253
     1  / 1.57 ,1.69 ,1.752,1.78 ,1.774,1.7  ,1.472,1.288,1.15 ,1.034,
     2     .938, .795, .69 , .605,
     3    1.57 ,1.708,1.76 ,1.778,1.727,1.56 ,1.35 ,1.18 ,1.047, .947,
     4     .86 , .72 , .617, .543,
     5    1.57 ,1.708,1.77 ,1.737,1.575,1.418,1.233,1.07 , .94 , .842,
     6     .768, .646, .558, .49 ,
     7    1.57 ,1.708,1.763,1.618,1.442,1.278,1.104, .955, .845, .758,
     8     .68 , .57 , .49 , .428,
     9    1.57 ,1.708,1.548,1.41 ,1.275,1.133, .98 , .84 , .738, .66 ,
     A     .594, .495, .425, .37 ,
     B    1.57 ,1.49 ,1.367,1.245,1.12 , .995, .845, .725, .628, .555,
     C     .5  , .423, .363, .312,
     D    1.254,1.273,1.169,1.062, .958, .848, .7  , .59 , .518, .465,
     E     .415, .345, .29 , .253,
     F     .94 , .943, .948, .875, .788, .695, .55 , .465, .404, .36 ,
     G     .325, .268, .227, .196,
     H     .62 , .623, .63 , .648, .597, .537, .403, .335, .287, .25 ,
     I     .22 , .18 , .15 , .13 , .31 , .315, .327, .347, .391, .35 ,
     J     .23 , .177, .147, .129, .112, .09 , .073, .062 /
C
C              FIGURE 4.1.3.3-61A
C
      DATA T61A/1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0/
      DATA SLOPE/59.375,41.667,33.0,27.778,23.611,16.667,12.5,8.75,6.0/
C
C           ****** FIGURE 4.1.4.2-26(A),(B),(C),(D),(E),(F) ******
C              ****SUPERSONIC****
C
      DATA T422AF
     1/0.,.2,.4,.6,.8,1.,
     21.,2.,3.,4.,5.,6.,
     30.,.2,.25,.33,.5,1./
      DATA DUMYA/
     1  .165,   .210,   .250,   .290,   .310,   .345,   .335,   .365,
     2  .390,   .415,   .445,   .470,   .500,   .540,   .560,   .560,
     3  .560,   .560,   .670,   .670,   .670,   .670,   .670,   .670,
     4  .830,   .775,   .775,   .775,   .775,   .775,   .990,   .930,
     5  .895,   .895,   .895,   .895,
     6  .200,   .230,   .280,   .305,   .335,   .360,   .400,   .445,
     7  .485,   .500,   .520,   .530,   .600,   .630,   .650,   .660,
     8  .665,   .665,   .795,   .800,   .800,   .805,   .810,   .815,
     9  .970,   .965,   .955,   .955,   .955,   .955,   1.150,  1.135,
     A  1.120,  1.100,  1.100,  1.105,
     B  .230,   .275,   .300,   .330,   .350,   .370,   .415,   .470,
     C  .500,   .530,   .545,   .550,   .630,   .670,   .680,   .685,
     D  .690,   .690,   .830,   .835,   .835,   .840,   .845,   .850,
     E  1.030,  1.015,  1.005,  1.000,  1.005,  1.010,  1.250,  1.225,
     F  1.200,  1.170,  1.165,  1.160/
      DATA DUMYB/
     1  .220,   .280,   .315,   .345,   .375,   .390,   .440,   .500,
     2  .535,   .560,   .570,   .580,   .670,   .700,   .720,   .725,
     3  .740,   .740,   .880,   .885,   .895,   .900,   .900,   .900,
     4  1.070,  1.070,  1.075,  1.075,  1.080,  1.080,  1.270,  1.260,
     5  1.260,  1.255,  1.255,  1.255,
     6  .250,   .300,   .330,   .380,   .415,   .445,   .500,   .560,
     7  .600,   .620,   .635,   .640,   .750,   .780,   .800,   .820,
     8  .820,   .825,   .980,   .990,   1.000,  1.020,  1.020,  1.020,
     9  1.190,  1.200,  1.200,  1.210,  1.220,  1.225,  1.380,  1.390,
     A  1.400,  1.410,  1.420,  1.420,
     B  .340,   .380,   .410,   .460,   .500,   .540,   .680,   .700,
     C  .730,   .770,   .790,   .840,   .950,   .990,   1.010,  1.050,
     D  1.080,  1.120,  1.200,  1.240,  1.290,  1.330,  1.370,  1.420,
     E  1.440,  1.500,  1.550,  1.610,  1.670,  1.720,  1.680,  1.760,
     F  1.820,  1.890,  1.950,  2.020/
      DATA DUMYC/
     1  .415,   .410,   .400,   .385,   .370,   .345,   .500,   .500,
     2  .495,   .485,   .480,   .470,   .585,   .580,   .580,   .575,
     3  .570,   .560,   .670,   .670,   .670,   .670,   .670,   .670,
     4  .750,   .750,   .755,   .760,   .765,   .775,   .830,   .840,
     5  .845,   .855,   .870,   .895,
     6  .460,   .455,   .445,   .420,   .390,   .360,   .575,   .575,
     7  .570,   .560,   .545,   .530,   .695,   .695,   .690,   .685,
     8  .680,   .665,   .800,   .805,   .805,   .810,   .815,   .815,
     9  .920,   .930,   .935,   .945,   .970,   .955,   1.040,  1.045,
     A  1.050,  1.075,  1.110,  1.105,
     B  .475,   .465,   .450,   .430,   .400,   .370,   .600,   .600,
     C  .595,   .585,   .575,   .550,   .725,   .730,   .730,   .725,
     D  .715,   .690,   .850,   .850,   .855,   .865,   .870,   .850,
     E  .970,   .975,   .980,   1.000,  1.020,  1.010,  1.110,  1.110,
     F  1.110,  1.130,  1.180,  1.160/
      DATA DUMYD/
     1  .500,   .490,   .470,   .450,   .425,   .390,   .640,   .635,
     2  .630,   .620,   .600,   .580,   .770,   .775,   .780,   .775,
     3  .765,   .740,   .920,   .915,   .920,   .930,   .935,   .900,
     4  1.050,  1.055,  1.060,  1.080,  1.105,  1.080,  1.195,  1.200,
     5  1.205,  1.225,  1.265,  1.255,
     6  .550,   .535,   .525,   .500,   .475,   .445,   .720,   .715,
     7  .710,   .690,   .670,   .640,   .890,   .890,   .890,   .885,
     8  .870,   .825,   1.060,  1.050,  1.050,  1.060,  1.065,  1.020,
     9  1.215,  1.215,  1.220,  1.245,  1.270,  1.225,  1.380,  1.380,
     A  1.395,  1.420,  1.470,  1.420,
     B  .760,   .730,   .700,   .650,   .600,   .540,   1.000,  1.000,
     C  .970,   .930,   .890,   .840,   1.240,  1.230,  1.230,  1.220,
     D  1.190,  1.120,  1.500,  1.480,  1.480,  1.490,  1.470,  1.420,
     E  1.750,  1.720,  1.730,  1.760,  1.780,  1.720,  2.000,  1.970,
     F  1.980,  2.020,  2.070,  2.020/
C
      SQP3=12.4707672
      CNPRME=.FALSE.
      MACH=FLC(MIDX+2)
      IF(TANLEO.EQ.0.0)TANLEO=.00001
      IF(TANLEI.EQ.0.0)TANLEI=.00001
      BETA=SQRT(MACH**2-1.)
      DO 1000 J=1,NALPHA
      ALPHAJ(J)= (ALSCHD(J)+ALPHAI)/RAD
 1000 CONTINUE
      IF(HTIN(15).NE.WTYPE(1))GO TO 1310
C
C       **STRAIGHT TAPERED SUPERSONIC NORMAL FORCE SLOPE**
C
      BOVERT=BETA/TANLE
      CA=COSLE
      IF( BOVERT.GT.1.0)GO TO 1010
      SUPLE=.FALSE.
      GO TO 1020
 1010 SUPLE=.TRUE.
 1020 CONTINUE
      PITCH=.FALSE.
C
C                   FIGURE 4.1.3.2-60A (CNA/CNA)T FOR SUBSONIC L.E.
C
 1030 DELTYT=DELTAY/CA
      IF(SUPLE)GO TO 1040
      VAB(1)=BOVERT
      VAB(2)=DELTYT
      LGB(1)=9
      LGB(2)=8
      CALL INTERX(2,A1350,VAB,LGB,DA50,CNCNT,9,72,
     1            0,0,0,0,0,2,0,0,Q3250A,3,ROUTID)
      IF(CNCNT.GT.1.)CNCNT=1.
      IF(PITCH)GO TO 1320
      CNNNT=CNCNT
      GO TO 1050
C
C                   FIGURE 4.1.3.2-60B (CNA/CNA)T FOR SUPERSONIC L.E.
C
 1040 CONTINUE
      ARG=DELTAY/(5.85*CA)
      DELTDT=ATAN(ARG)*RAD
      VAB(1)=1./BOVERT
      VAB(2)=DELTDT
      LGB(1)=11
      LGB(2)=8
      CALL INTERX(2,B1350,VAB,LGB,DB50,CNCNT,11,88,
     1            0,0,0,0,0,2,0,0,Q3250B,3,ROUTID)
      IF(CNCNT.GT.1.)CNCNT=1.
      IF(PITCH)GO TO 1320
      CNNNT=CNCNT
 1050 IF(TRATIP.EQ.1.0.AND.SWEPLE.EQ.0.0) GO TO 1060
C
C                   NON-RECTANGULAR WING ANALYSIS
C                   FIGURES 4.1.3.2-56A THROUGH F
C
      VAA(1)=BOVERT
      VAA(2)=AR*TANLE
      VAA(3)=TAPEXP
      LGA(1)=23
      LGA(2)=8
      LGA(3)=6
      CALL INTERX(3,T13246,VAA,LGA,D13246,BCNA,23,1104,
     1            0,2,0,0,0,2,0,0,Q3246A,3,ROUTID)
      IF(SUPLE)GO TO 1080
      CNTHRY=BCNA/TANLE
      GO TO 1090
 1060 CONTINUE
C
C                        RECTANGULAR WING ANALYSIS
C
      IF(AR*BETA.GT.1.0)GO TO 1070
C
C                        FIGURE 4.1.3.2-56G(CNA/A)
C
      VAB(1)=AR*BETA
      LGB(1)=9
      CALL INTERX(1,G13246,VAB,LGB,DG3246,CNAA,9,9,
     1            0,0,0,0,0,0,0,0,Q3246G,3,ROUTID)
      CNTHRY=CNAA*AR
      GO TO 1090
 1070 BCNA= 4.-2.*(1./(AR*BETA))
 1080 CNTHRY=BCNA/BETA
 1090 CNA=CNTHRY*CNCNT
      CNA1=CNA*SRSTAR/SW
      TLE192=TANLE/1.92
      IF(TLE192.LE.1.)GO TO 1100
      ARG=DELTAY/(5.85*COSLE)
      DELTDT=ATAN(ARG)*RAD
C
C                   FIGURE 4.1.3.3-59A (VALUE OF CC)
C
      VAB(1)=DELTDT
      VAB(2)=CNA*TLE192
      LGB(1)=9
      LGB(2)=11
      CALL INTERX(2,A13359,VAB,LGB,D13359,CC,11,99,
     1            0,2,0,0,2,2,0,0,Q3359A,3,ROUTID)
      E=CNA*(TLE192+CC*(TLE192-1.))
      GO TO 1110
 1100 E=CNA
 1110 CONTINUE
      IF(SUPLE)GO TO 1160
C
C     WING LIFT IN THE NONLINEAR ANGLE OF ATTACK RANGE SUBSONIC L.E.
C                         STRAIGHT TAPERED WING
C           COMPUTATION OF LIFT COEFICIENT VS ANGLE OF ATTACK
C
      ARG=E*BOVERT
      DO 1150 J=1,NALPHA
         IF(ALPHAJ(J).NE.0.0)GO TO 1120
         CNAAA(J)=0.0
         GO TO 1140
 1120    ARG1=CNCNT/(BETA*ABS(TAN(ALPHAJ(J))))
         IF(ARG1.GT.1.)GO TO 1130
C
C              FIGURE 4.1.3.3-59B (CN(AA)) RIGHT SIDE
C
         VAB(1)=ARG1
         VAB(2)=ARG
         LGB(1)=11
         LGB(2)=10
         CALL INTERX(2,BR3359,VAB,LGB,DR3359,CNAAA(J),11,110,
     1               2,2,0,0,0,0,0,0,Q3359B,3,ROUTID)
         GO TO 1140
 1130    CONTINUE
C
C              FIGURE 4.1.3.3-59B(CN(AA)) LEFT SIDE
C
         VAB(1)=1./ARG1
         VAB(2)=ARG
         LGB(1)=5
         LGB(2)=10
         CALL INTERX(2,BL3359,VAB,LGB,DL3359,CNAAA(J),10,50,
     1               2,2,0,0,0,0,0,0,Q3359B,3,ROUTID)
 1140    CONTINUE
         CN(J)=CNA* SIN(2.*ALPHAJ(J))/2.+CNAAA(J)*SIN(ALPHAJ(J))*ABS
     1         (SIN(ALPHAJ(J)))
         CN(J)=CN(J)*SRSTAR/SW
         CL(J)=CN(J)*COS(ALPHAJ(J))
 1150 CONTINUE
      GO TO 1420
C
C     WING LIFT IN THE NONLINEAR ANGLE OF ATTACK RANGE SUPERSONIC L.E.
C                         STRAIGHT TAPERED WING
 1160 CONTINUE
C
C ****IF SHOCK IS ATTACHED AT ZERO ANGLE OF ATTACK SET DETACH=FALSE
C
      TRANS=.FALSE.
      RMACH=MACH*SQRT(1.-SINLE **2)
      IF(RMACH.LT.1.)GO TO 1250
      CALL ANGDET(RMACH,DELMAX)
      IF(DELTDT/RAD.GT.DELMAX)GO TO 1250
      DETACH=.FALSE.
C
C         FIND ANGLE OF ATTACH AT WHICH SHOCK WILL START TO DETACH
C
      ANGLE=1.
      ANGMAX=ABS(ALPHAJ(1))
      IF(ANGMAX.LT.ABS(ALPHAJ(NALPHA)))ANGMAX=ABS(ALPHAJ(NALPHA))
 1170 DEG=ANGLE/RAD
      RRMACH=MACH* (1.-(SINLE *COS(DEG))**2)**.5
      IF(RRMACH.LT.1.)GO TO 1180
      CALL ANGDET(RRMACH,DELMAX)
      IF((DELTDT/RAD+TAN(DEG)/COSLE).LT.DELMAX)GO TO 1180
      DETANG=DEG
      GO TO 1200
 1180 ANGLE=ANGLE+1.
      IF(ANGLE.GT.ANGMAX)GO TO 1190
      GO TO 1170
 1190 DETANG=90./RAD
      GO TO 1220
 1200 CONTINUE
C
C              FIGURE 4.1.3.3-60B
C
      BTANA = BETA*TAN(DETANG)
      CALL FIG60B(BETA, BTANA, CNAAST)
C
C              FIGURE 4.1.3.3-61A
C
      VAB(1)=BETA
      LGB(1)=9
      CALL INTERX(1,T61A,VAB,LGB,SLOPE,DELCN,9,9,
     1            0,0,0,0,2,0,0,0,Q3361A,3,ROUTID)
      DETALP=DELCN*CNAAST/RAD
      ALPHAP=DETANG+DETALP
      CNPRME = .TRUE.
      J=1
      ARG2=BETA*ABS(TAN(ALPHAP))
      GO TO 1280
 1210 CNAAAP=CNAAA(J)
      CNPRME=.FALSE.
 1220 J=0
 1230 J=J+1
      IF(ABS(ALPHAJ(J)).LT.DETANG)GO TO 1240
      TRANS=.TRUE.
      RATIO=(ABS(ALPHAJ(J))-DETANG)/(ALPHAP-DETANG)
      ARG2=BETA*ABS(TAN(ALPHAJ(J)))
      IF(ABS(ALPHAJ(J)).LT.ALPHAP)GO TO 1300
      TRANS=.FALSE.
      GO TO 1270
 1240 CONTINUE
C
C              FIGURE 4.1.3.3-60B
C
      BTANA = BETA*ABS(TAN(ALPHAJ(J)))
      CALL FIG60B(BETA, BTANA, CNAAA(J))
      TRANS=.FALSE.
      GO TO 1300
 1250 CONTINUE
      DETACH=.TRUE.
C
C              SUPERSONIC WING LIFT VS ANGLE OF ATTACK WITH DETACHED
C                   SHOCK AT ZERO ANGLE OF ATTACK
C
      J=0
 1260 J=J+1
      ARG2=BETA*ABS(TAN(ALPHAJ(J)))
 1270 IF(ALPHAJ(J).NE.0.0)GO TO 1280
      CNAAA(J)=0.0
      GO TO 1300
 1280 ARG1=CNCNT/ARG2
C
C                   FIGURE 4.1.3.3-60A (CNAAA) RIGHT SIDE
C
      VAB(2)=E
      IF(ARG1.GT.1.) GO TO 1290
      VAB(1)=ARG1
      LGB(1)=13
      LGB(2)=12
      CALL INTERX(2,TR360A,VAB,LGB,DR360A,CNAAA(J),13,156,
     1            2,0,0,0,0,0,0,0,Q3360A,3,ROUTID)
      GO TO 1300
 1290 CONTINUE
C
C                   FIGURE 4.1.3.3-60A (CNAAA) LEFT SIDE
C
      VAB(1)=1./ARG1
      LGB(1)=6
      LGB(2)=12
      CALL INTERX(2,TL360A,VAB,LGB,DL360A,CNAAA(J),12,72,
     1            0,0,0,0,0,0,0,0,Q3360A,3,ROUTID)
 1300 CONTINUE
      IF(CNPRME) GO TO 1210
      IF(TRANS) CNAAA(J)=CNAAST-RATIO*(CNAAST-CNAAAP)
      CN(J)=CNA*SIN(2.*ALPHAJ(J))/2.+CNAAA(J)*SIN(ALPHAJ(J))*ABS(SIN(ALP
     1       HAJ(J)))
      CN(J)=CN(J)* SRSTAR/SW
      CL(J)=CN(J)*COS(ALPHAJ(J))
      IF(J.GE.NALPHA)GO TO 1420
      IF(DETACH)GO TO 1260
      GO TO 1230
C
C    ***SUPERSONIC NORMAL FORCE SLOPE,NON-STRAIGHT TAPERED WINGS***
C                    BASIC WING COMPONENT
C
 1310 CONTINUE
      SPANIN=SPANS-SPANO
      CRBW=CB+SPANIN*(TANLEO-TANTEO)
      SBW=(CRBW+CT)*SPANS
      ARBW=4.*SPANS**2/SBW
      TAPBW=CT/CRBW
      BOVERT=BETA/TANLEO
      SUPLE=.FALSE.
      IF(BOVERT.GT.1.)SUPLE=.TRUE.
C
C                   FIGURE 4.1.3.2-56A-F (CNABW)
C
      VAB(1)=BOVERT
      VAB(2)=ARBW*TANLEO
      VAB(3)=TAPBW
      PITCH=.FALSE.
      I=0
      GLOVE=.FALSE.
      TA=TANLEO
 1320 LGB(1)=23
      LGB(2)=8
      LGB(3)=6
      IF(I.NE.0) VAB(2)=ARP*TA
      CALL INTERX(3,T13246,VAB,LGB,D13246,BCNA,23,1104,
     1            0,2,0,0,0,2,0,0,Q246AF,4,ROUTID)
      CNTHRY=BCNA/TA
      IF(SUPLE)CNTHRY=BCNA/BETA
      IF(PITCH.AND.I.EQ.1)GO TO 1520
      IF(PITCH.AND.I.EQ.2)GO TO 1530
C
C                   LEADING-EDGE-EFFECT PARAMETERS
C              FIGURE 4.1.3.2-62 (CLE)BW
C
      ROUND=.FALSE.
      IF(KSHARP.EQ.UNUSED)ROUND=.TRUE.
      VAB(1)=BOVERT
      LGB(1)=12
C
C              SHARP LEADING EDGE
C
      IF(ROUND)GO TO 1330
      CALL INTERX(1,S13252,VAB,LGB,DSHP52,CLE  ,12,12,
     1            0,0,0,0,0,0,0,0,Q13252,3,ROUTID)
      GO TO 1340
C
C              ROUND LEADING EDGE
C
 1330 CALL INTERX(1,R13252,VAB,LGB,DRND52,CLE  ,12,12,
     1            0,0,0,0,0,0,0,0,Q13252,3,ROUTID)
 1340 IF(GLOVE)GO TO 1370
      CNTBW=CNTHRY
      CLEBW=CLE
      IF(TANLEO.GE.TANLEI)GO TO 1350
      CNABW=CNTBW*SBW/SW *CLEBW
      GO TO 1360
 1350 CONTINUE
      S2=SPANIN**2*TANLEO
      A2=4.*SPANIN**2/S2
      VAB(2)=A2*TANLEO
      VAB(3)=0.0
      CALL INTERX(3,T13246,VAB,LGB,D13246,BCNA,23,1104,
     1            0,2,0,0,0,2,0,0,Q246AF,4,ROUTID)
      CNT2=BCNA/BETA
      IF(VAB(1).GT.1.)CNT2  =BCNA/TANLEO
      CNABW=(CNTHRY*SBW/SW-CNT2  *  S2/SW )*CLEBW
 1360 CONTINUE
C
C                        GLOVE COMPONENT
C
      IF(ABS(SWEPLE-SWEPLO).LT.4.)GO TO 1380
      CRGLV=TANLEI*SPANIN
      SGLV=CRGLV*SPANIN
      ARGLV=4.*SPANIN**2/SGLV
      VAB(1)=BETA/TANLEI
      VAB(2)=ARGLV*TANLEI
      VAB(3)=0.0
      GLOVE=.TRUE.
      BOVERT=VAB(1)
      SUPLE=.FALSE.
      TA=TANLEI
      IF(BOVERT.GT.1.)SUPLE=.TRUE.
      GO TO 1320
 1370 CONTINUE
      CLEGLV=CLE
      CNAGLV=CNTHRY*CLEGLV*SGLV/SW
C
C              FIGURE 4.1.3.2-63 (CN1=CNA/A)
C
      GO TO 1390
 1380 CNAGLV=0.0
C
C                        EXTENSION COMPONENT
C
 1390 IF(ABS(SWTEI-SWTEO).LT.4.)GO TO 1400
      BE= 2.*SPANIN
      VAB(1)=BETA/TANLEI
      VAB(2)=TANTEI/TANLEI
      LGB(1)=14
      LGB(2)=10
      CALL INTERX(2,T13253,VAB,LGB,D13253,CN1,14,140,
     1            2,2,0,0,2,2,0,0,Q13253,3,ROUTID)
      VAB(2)=TANTEO/TANLEI
      CALL INTERX(2,T13253,VAB,LGB,D13253,CN2,14,140,
     1            2,2,0,0,2,2,0,0,Q13251,3,ROUTID)
      CNAE=(CN1-CN2)*BE**2/SW
      GO TO 1410
 1400 CNAE=0.0
C
C**SUPERSONIC NORMAL FORCE SLOPE,NON-STRAIGHT TAPERED WINGS,TOTAL WING
C
 1410 CONTINUE
C
C              FIGURE 4.1.3.2-61
C
      VAB(1)=CNAGLV/(CLEGLV*BETA)*SW/A(4)
      LGB(1)= 12
      CALL INTERX(1,T13251,VAB,LGB,D13251,RKL,12,12,
     1            2,0,0,0,2,0,0,0,Q13251,3,ROUTID)
      CNA1=RKL*(CNABW+CNAGLV+CNAE)
 1420 CONTINUE
C
C      ***PITCHING MOMENT SLOPES AT LOW ANGLE OF ATTACK,STRAIGHT TAPERED
C
      DXW=(SPAN-SPANS)*TANLEI*COS(ALPHAI/RAD)
      DETXCG=XCG-(XW+DXW)
      IF(HTIN(15).NE.WTYPE(1))GO TO 1500
C
C              FIGURE 4.1.4.2-26A-F (XAC/CRSTAR)
C
      TA=TANLE
      ARP=ARSTAR
      TAPR=TAPEXP
      PITCH=.FALSE.
      VAB(2)=ARP*TA
      VAB(3)=TAPR
 1430 LGB(1)=6
      LGB(2)=6
      LGB(3)=6
      IF(BOVERT.GT.1.)GO TO 1440
      VAB(1)=BOVERT
      CALL INTERX(3,T422AF,VAB,LGB,SUBAF,XACCR,6,216,
     1            0,2,0,0,0,2,0,0,Q222AF,4,ROUTID)
      GO TO 1450
 1440 VAB(1)=1./BOVERT
      CALL INTERX(3,T422AF,VAB,LGB,SUPAF,XACCR,6,216,
     1 0,2,0,0,0,2,0,0,Q222AF,4,ROUTID)
 1450 CONTINUE
      IF(BOVERT.GE.0.) GO TO 1455
      VAB(1)=1./BOVERT
      CALL FWDXAC(VAB(2),VAB(3),VAB(1),MACH,ROUTID,XACCR)
 1455 CONTINUE
      IF(SWEPLE .LT. 1.E-10 .AND. HTIN(MIDX+71) .NE. UNUSED)
     1      XACCR=HTIN(MIDX+71)*A(16)/CRSTAR
      IF(PITCH)GO TO 1460
      GO TO 1480
 1460 IF(I.EQ.1)XACCRI=XACCR
      IF(I.EQ.2)GO TO 1470
      GO TO 1030
 1470 ARG1=XACCR*CBP/CRSTAR
      ARG2=BSTARB/(CRSTAR*2.)*TANLEO
      ARG3=BSTARB*TANLEI/CRSTAR
      XACCRO=ARG1-ARG2+ARG3
      GO TO 1030
 1480 DCMCL=(DETXCG/CRSTAR-XACCR)*CRSTAR/CBARR
      CMA=DCMCL*CNA1
      DO 1490 J=1,NALPHA
         IF(ALPHAJ(J).GT.0.175)GO TO 1490
         CM(J)=CMA*ALPHAJ(J)
 1490 CONTINUE
      GO TO 1540
C
C   ***SUPERSONIC WING PITCHING MOMENT SLOPE,NON-STRAIGHT TAPERED WINGS*
C
 1500 CONTINUE
      ARP=ARIP
      I=1
      CA=COSLE
      TA=TANLE
      TAPR=TAPIP
 1510 BOVERT=BETA/TA
      VAB(1)=BOVERT
      VAB(2)=ARP*TA
      VAB(3)=TAPR
      SUPLE=.FALSE.
      IF(BOVERT.GT.1.)SUPLE=.TRUE.
      PITCH=.TRUE.
      GO TO 1430
 1520 CNCNTI=CNCNT
      CNATI=CNTHRY
      CNAI=CNCNTI*CNATI
      CA=COSLEO
      TA=TANLEO
      TAPR=TAPOP
      I=2
      ARP=AROP
      GO TO 1510
 1530 CNCNTO=CNCNT
      CNATO=CNTHRY
      CNAO=CNCNTO*CNATO
      ARG1=CNAI*SI*XACCRI + CNAO*SO*XACCRO
      ARG2=CNAI*SI+CNAO*SO
      XACCR=ARG1/ARG2
      GO TO 1480
 1540 CONTINUE
      ARG=(1.+P)/(PI*ARSTAR*P)
C
C  ***SUPERSONIC WING DRAG,CN AND CA VS ANGLE OF ATTACK***
C
      DO 1550 J=1,NALPHA
         CDL(J)=DRAGC*ARG*SW/SRSTAR*CL(J)**2
         CD(J)=CDO + CDL(J)
         COSA=COS(FLC(J+22)/RAD)
         SINA=SIN(FLC(J+22)/RAD)
         IF(CL(J).NE.UNUSED) CN(J)=CL(J)*COSA+CD(J)*SINA
         IF(CL(J).NE.UNUSED) CAF(J)=CD(J)*COSA-CL(J)*SINA
 1550 CONTINUE
      IF(CMA.NE.UNUSED) CMADEG=CMA/RAD
      IF(CNA1.NE.UNUSED) CLADEG=CNA1/RAD
      IF(TANLEO.EQ.0.00001)TANLEO=0.0
      IF(TANLEI.EQ.0.00001)TANLEI=0.0
      RETURN
      END
