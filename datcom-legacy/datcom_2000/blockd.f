      BLOCK DATA BLOCKD
C
C     SET PROGRAM CONSTANTS
C
      COMMON /IBODY/   PBODY, LEN2(18),LDM2(18),FLTCON(81),DUM1(283)
      COMMON /BODYI/          LEN3(4),LDM3(4),OPTINS(20),DUM2(101)
      COMMON /IHT/     PHT,   LEN4(15),LDM4(15),BODY(44),DUM3(306)
      COMMON /IVT/     PVT,   LEN5(21),LDM5(21),PLNF(103),DUM4(235)
      COMMON /IBW/     PBW,   LEN6(36),LDM6(36),SCHR(170),DUM5(138)
      COMMON /IBH/     PBH,   LEN7(19),LDM7(19),SYNTHS(60),DUM6(282)
      COMMON /IBV/     PBV,   LENC(14),LDMC(14),PRPOWR(78),DUM7(274)
      COMMON /IBWH/    PBWH,  LEND(15),LDMD(15),JETPWR(89),DUM8(261)
      COMMON /IBWV/    PBWV,  LENE(21),LDME(21),LARWB(89),DUM9(249)
      COMMON /IPOWER/  PPOWR, LENF(2),LDMF(2),GRNDEF(8),DUM10(188)
      COMMON /IDWASH/  PDWAS, LENG(8),LDMG(8),TVTPAN(23),DUM11(21)
      COMMON /IBWHV/   PBWHV, LENH(32),LDMH(32),EXPR(131),DUM12(185)
      COMMON /IWING/   PWING, LENI(27),LDMI(27),SYMFLP(139),DUM13(207)
      COMMON /VTI/            LENJ(14),LDMJ(14),ASYFLP(78),DUM14(210)
      COMMON /WINGI/          LENK(7),LDMK(7),HYPEFF(32),DUM15(55)
      COMMON /HTI/            LENL(12),LDML(12),TRNJET(38),DUM16(92)
      COMMON /IVF/     PVF,   LENM(22),LDMM(22),CONTAB(76),DUM17(260)
      COMMON /CONSNT/         CONST(5)
C
      INTEGER FLTCON, OPTINS, BODY, PLNF, SCHR, PRPOWR, SYNTHS,
     1        JETPWR, LARWB, GRNDEF, TVTPAN, EXPR, SYMFLP,
     2        ASYFLP, HYPEFF, TRNJET, CONTAB
C
      DIMENSION SCH1(158),SCH2(12)
      EQUIVALENCE (SCH1(1),SCHR(1)),(SCH2(1),SCHR(159))
C
C
      DATA CONST / 3.141592654,0.01745329,1.E-30,57.2957795,4H$   /
C
      DATA LEN2 / 5,4,6,6,5,6,4,6,6,2,3,4,4,2,5,4,4,5 /
      DATA LDM2 / 1,20,1,20,20,-1,20,3*1,3*20,4*1,20 /
      DATA FLTCON /   4HN   ,4HM   ,4HA   ,4HC   ,4HH   ,4HM   ,
     1  4HA   ,4HC   ,4HH   ,4HN   ,4HA   ,4HL   ,4HP   ,4HH   ,
     2  4HA   ,4HA   ,4HL   ,4HS   ,4HC   ,4HH   ,4HD   ,4HR   ,
     3  4HN   ,4HN   ,4HU   ,4HB   ,4HH   ,4HY   ,4HP   ,4HE   ,
     4  4HR   ,4HS   ,4HP   ,4HI   ,4HN   ,4HF   ,4HS   ,4HT   ,
     5  4HM   ,4HA   ,4HC   ,4HH   ,4HT   ,4HS   ,4HM   ,4HA   ,
     6  4HC   ,4HH   ,4HT   ,4HR   ,4HA   ,4HL   ,4HT   ,4HT   ,
     7  4HI   ,4HN   ,4HF   ,4HV   ,4HI   ,4HN   ,4HF   ,4HW   ,
     8  4HT   ,4HG   ,4HA   ,4HM   ,4HM   ,4HA   ,4HN   ,4HA   ,
     9  4HL   ,4HT   ,4HL   ,4HO   ,4HO   ,4HP   ,4HA   ,4HL   ,
     A  4HP   ,4HH   ,4HA   /
C
      DATA LEN3 / 4,5,6,5 /
      DATA LDM3 / 4*1 /
      DATA OPTINS /   4HS   ,4HR   ,4HE   ,4HF   ,4HC   ,4HB   ,
     1  4HA   ,4HR   ,4HR   ,4HR   ,4HO   ,4HU   ,4HG   ,4HF   ,
     2  4HC   ,4HB   ,4HL   ,4HR   ,4HE   ,4HF   /
C
      DATA LEN4 / 2,4*1,2,2,5,5,3,3,2,5,6,5 /
      DATA LDM4 / 1,6*20,8*1 /
      DATA BODY /     4HN   ,4HX   ,4HX   ,4HS   ,4HP   ,4HR   ,
     1  4HZ   ,4HU   ,4HZ   ,4HL   ,4HB   ,4HN   ,4HO   ,4HS   ,
     2  4HE   ,4HB   ,4HT   ,4HA   ,4HI   ,4HL   ,4HB   ,4HL   ,
     3  4HN   ,4HB   ,4HL   ,4HA   ,4HD   ,4HS   ,4HI   ,4HT   ,
     4  4HY   ,4HP   ,4HE   ,4HM   ,4HE   ,4HT   ,4HH   ,4HO   ,
     5  4HD   ,4HE   ,4HL   ,4HL   ,4HI   ,4HP   /
C
      DATA LEN5 / 6,5,6,6,4,5,6,3*5,6,4,3*6,3,3*4,3,4 /
      DATA LDM5 / 15*1,6*20 /
      DATA PLNF /     4HC   ,4HH   ,4HR   ,4HD   ,4HB   ,4HP   ,
     1  4HC   ,4HH   ,4HR   ,4HD   ,4HR   ,4HC   ,4HH   ,4HR   ,
     2  4HD   ,4HT   ,4HP   ,4HC   ,4HH   ,4HS   ,4HT   ,4HA   ,
     3  4HT   ,4HS   ,4HS   ,4HP   ,4HN   ,4HS   ,4HS   ,4HP   ,
     4  4HN   ,4HE   ,4HS   ,4HS   ,4HP   ,4HN   ,4HO   ,4HP   ,
     5  4HS   ,4HA   ,4HV   ,4HS   ,4HI   ,4HS   ,4HA   ,4HV   ,
     6  4HS   ,4HO   ,4HS   ,4HW   ,4HA   ,4HF   ,4HP   ,4HT   ,
     7  4HW   ,4HI   ,4HS   ,4HT   ,4HA   ,4HT   ,4HY   ,4HP   ,
     8  4HE   ,4HS   ,4HS   ,4HP   ,4HN   ,4HD   ,4HD   ,4HD   ,
     9  4HH   ,4HD   ,4HA   ,4HD   ,4HI   ,4HD   ,4HH   ,4HD   ,
     A  4HA   ,4HD   ,4HO   ,4HS   ,4HH   ,4HB   ,4HS   ,4HE   ,
     B  4HX   ,4HT   ,4HR   ,4HL   ,4HP   ,4HH   ,4HS   ,4HV   ,
     C  4HW   ,4HB   ,4HS   ,4HV   ,4HB   ,4HS   ,4HV   ,4HH   ,
     D  4HB   /
C
      DATA LEN6 / 4,6,4,3,6,6,5,6,3,5,4,4,4,5,3,4,5,6,6,5,5,5,3,
     1            5,6,4,5,6,6,4,5,3,3,4,6,6 /
      DATA LDM6 / 5*1,20,20,-1,11*1,6,1,1,20,3*1,5*50,5*1 /
      DATA SCH1 /     4HT   ,4HO   ,4HV   ,4HC   ,4HD   ,4HE   ,
     1  4HL   ,4HT   ,4HA   ,4HY   ,4HX   ,4HO   ,4HV   ,4HC   ,
     2  4HC   ,4HL   ,4HI   ,4HA   ,4HL   ,4HP   ,4HH   ,4HA   ,
     3  4HI   ,4HC   ,4HL   ,4HA   ,4HL   ,4HP   ,4HA   ,4HC   ,
     4  4HL   ,4HM   ,4HA   ,4HX   ,4HC   ,4HA   ,4HM   ,4HB   ,
     5  4HE   ,4HR   ,4HC   ,4HM   ,4H0   ,4HX   ,4HO   ,4HV   ,
     6  4HC   ,4HO   ,4HC   ,4HM   ,4H0   ,4HT   ,4HL   ,4HE   ,
     7  4HR   ,4HI   ,4HL   ,4HE   ,4HR   ,4HO   ,4HT   ,4HO   ,
     8  4HV   ,4HC   ,4HO   ,4HC   ,4HM   ,4HO   ,4HC   ,4HM   ,
     9  4HO   ,4HT   ,4HT   ,4HC   ,4HE   ,4HF   ,4HF   ,4HK   ,
     A  4HS   ,4HH   ,4HA   ,4HR   ,4HP   ,4HC   ,4HL   ,4HM   ,
     B  4HA   ,4HX   ,4HL   ,4HS   ,4HL   ,4HO   ,4HP   ,4HE   ,
     C  4HC   ,4HL   ,4HA   ,4HM   ,4HO   ,4HC   ,4HL   ,4HA   ,
     D  4HM   ,4H0   ,4HX   ,4HA   ,4HC   ,4HD   ,4HW   ,4HA   ,
     E  4HS   ,4HH   ,4HT   ,4HY   ,4HP   ,4HE   ,4HI   ,4HN   ,
     F  4HN   ,4HP   ,4HT   ,4HS   ,4HX   ,4HC   ,4HO   ,4HR   ,
     G  4HD   ,4HY   ,4HU   ,4HP   ,4HP   ,4HE   ,4HR   ,4HY   ,
     H  4HL   ,4HO   ,4HW   ,4HE   ,4HR   ,4HM   ,4HE   ,4HA   ,
     I  4HN   ,4HT   ,4HH   ,4HI   ,4HC   ,4HK   ,4HY   ,4HC   ,
     J  4HM   ,4HC   ,4HL   ,4HD   ,4HA   ,4HR   ,4HC   ,4HL   /
      DATA SCH2 /     4HA   ,4HL   ,4HP   ,4HH   ,4HA   ,4HO   ,
     1  4HA   ,4HL   ,4HP   ,4HH   ,4HA   ,4H0   /
C
      DATA LEN7 / 3,2,2,4,3,2,2,4,2,6,5,3,5,2,3,2,2,4,4 /
      DATA LDM7 / 9*1,-1,9*1 /
      DATA SYNTHS /   4HX   ,4HC   ,4HG   ,4HX   ,4HW   ,4HZ   ,
     1  4HW   ,4HA   ,4HL   ,4HI   ,4HW   ,4HZ   ,4HC   ,4HG   ,
     2  4HX   ,4HH   ,4HZ   ,4HH   ,4HA   ,4HL   ,4HI   ,4HH   ,
     3  4HX   ,4HV   ,4HV   ,4HE   ,4HR   ,4HT   ,4HU   ,4HP   ,
     4  4HH   ,4HI   ,4HN   ,4HA   ,4HX   ,4HX   ,4HV   ,4HF   ,
     5  4HS   ,4HC   ,4HA   ,4HL   ,4HE   ,4HZ   ,4HV   ,4HZ   ,
     6  4HV   ,4HF   ,4HY   ,4HV   ,4HY   ,4HF   ,4HP   ,4HH   ,
     7  4HI   ,4HV   ,4HP   ,4HH   ,4HI   ,4HF   /
C
      DATA LENC / 12*6,4,2 /
      DATA LDMC / 12*1,-1,1 /
      DATA PRPOWR /   4HA   ,4HI   ,4HE   ,4HT   ,4HL   ,4HP   ,
     1  4HN   ,4HE   ,4HN   ,4HG   ,4HS   ,4HP   ,4HT   ,4HH   ,
     2  4HS   ,4HT   ,4HC   ,4HP   ,4HP   ,4HH   ,4HA   ,4HL   ,
     3  4HO   ,4HC   ,4HP   ,4HH   ,4HV   ,4HL   ,4HO   ,4HC   ,
     4  4HP   ,4HR   ,4HP   ,4HR   ,4HA   ,4HD   ,4HE   ,4HN   ,
     5  4HG   ,4HF   ,4HC   ,4HT   ,4HB   ,4HW   ,4HA   ,4HP   ,
     6  4HR   ,4H3   ,4HB   ,4HW   ,4HA   ,4HP   ,4HR   ,4H6   ,
     7  4HB   ,4HW   ,4HA   ,4HP   ,4HR   ,4H9   ,4HN   ,4HO   ,
     8  4HP   ,4HB   ,4HP   ,4HE   ,4HB   ,4HA   ,4HP   ,4HR   ,
     9  4H7   ,4H5   ,4HC   ,4HR   ,4HO   ,4HT   ,4HY   ,4HP   /
C
      DATA LEND / 14*6,5 /
      DATA LDMD / 15*1 /
      DATA JETPWR /   4HA   ,4HI   ,4HE   ,4HT   ,4HL   ,4HJ   ,
     1  4HN   ,4HE   ,4HN   ,4HG   ,4HS   ,4HJ   ,4HT   ,4HH   ,
     2  4HS   ,4HT   ,4HC   ,4HJ   ,4HJ   ,4HI   ,4HA   ,4HL   ,
     3  4HO   ,4HC   ,4HJ   ,4HE   ,4HV   ,4HL   ,4HO   ,4HC   ,
     4  4HJ   ,4HE   ,4HA   ,4HL   ,4HO   ,4HC   ,4HJ   ,4HI   ,
     5  4HN   ,4HL   ,4HT   ,4HA   ,4HJ   ,4HE   ,4HA   ,4HN   ,
     6  4HG   ,4HL   ,4HJ   ,4HE   ,4HV   ,4HE   ,4HL   ,4HO   ,
     7  4HA   ,4HM   ,4HB   ,4HT   ,4HM   ,4HP   ,4HJ   ,4HE   ,
     8  4HS   ,4HT   ,4HM   ,4HP   ,4HJ   ,4HE   ,4HL   ,4HL   ,
     9  4HO   ,4HC   ,4HJ   ,4HE   ,4HT   ,4HO   ,4HT   ,4HP   ,
     A  4HA   ,4HM   ,4HB   ,4HS   ,4HT   ,4HP   ,4HJ   ,4HE   ,
     B  4HR   ,4HA   ,4HD   /
C
      DATA LENE / 2,4,6,6,2,6,6,1,4,6,5,2,2,3,3,6,6,3,5,6,5 /
      DATA LDME / 13*1,-1,1,1,-1,4*1 /
      DATA LARWB  /   4HZ   ,4HB   ,4HS   ,4HR   ,4HE   ,4HF   ,
     1  4HD   ,4HE   ,4HL   ,4HT   ,4HE   ,4HP   ,4HS   ,4HF   ,
     2  4HR   ,4HO   ,4HN   ,4HT   ,4HA   ,4HR   ,4HR   ,4H3   ,
     3  4HL   ,4HE   ,4HO   ,4HB   ,4HD   ,4HE   ,4HL   ,4HT   ,
     4  4HA   ,4HL   ,4HL   ,4HS   ,4HW   ,4HE   ,4HT   ,4HP   ,
     5  4HE   ,4HR   ,4HB   ,4HA   ,4HS   ,4HS   ,4HB   ,4HA   ,
     6  4HS   ,4HE   ,4HH   ,4HB   ,4HB   ,4HB   ,4HB   ,4HL   ,
     7  4HF   ,4HX   ,4HC   ,4HG   ,4HT   ,4HH   ,4HE   ,4HT   ,
     8  4HA   ,4HD   ,4HR   ,4HO   ,4HU   ,4HN   ,4HD   ,4HN   ,
     9  4HS   ,4HB   ,4HS   ,4HS   ,4HB   ,4HS   ,4HL   ,4HB   ,
     A  4HX   ,4HC   ,4HE   ,4HN   ,4HS   ,4HB   ,4HX   ,4HC   ,
     B  4HE   ,4HN   ,4HW   /
C
      DATA LENF / 3,5 /
      DATA LDMF / 1,20 /
      DATA GRNDEF /   4HN   ,4HG   ,4HH   ,4HG   ,4HR   ,4HD   ,
     1  4HH   ,4HT   /
C
      DATA LENG / 3,2,3,2,2,6,3,2 /
      DATA LDMG / 8*1 /
      DATA TVTPAN /   4HB   ,4HV   ,4HP   ,4HB   ,4HV   ,4HB   ,
     1  4HD   ,4HV   ,4HB   ,4HH   ,4HS   ,4HV   ,4HV   ,4HP   ,
     2  4HH   ,4HI   ,4HT   ,4HE   ,4HV   ,4HL   ,4HP   ,4HZ   ,
     3  4HP   /
C
      DATA LENH / 3*3,4,4,3*3,4,4,3*3,5*4,5,5,6,6,5,3,5*5,4,5,4 /
      DATA LDMH / 23*20,9*1 /
      DATA EXPR /     4HC   ,4HD   ,4HB   ,4HC   ,4HL   ,4HB   ,
     1  4HC   ,4HM   ,4HB   ,4HC   ,4HL   ,4HA   ,4HB   ,4HC   ,
     2  4HM   ,4HA   ,4HB   ,4HC   ,4HD   ,4HW   ,4HC   ,4HL   ,
     3  4HW   ,4HC   ,4HM   ,4HW   ,4HC   ,4HL   ,4HA   ,4HW   ,
     4  4HC   ,4HM   ,4HA   ,4HW   ,4HC   ,4HD   ,4HH   ,4HC   ,
     5  4HL   ,4HH   ,4HC   ,4HM   ,4HH   ,4HC   ,4HL   ,4HA   ,
     6  4HH   ,4HC   ,4HM   ,4HA   ,4HH   ,4HC   ,4HD   ,4HW   ,
     7  4HB   ,4HC   ,4HL   ,4HW   ,4HB   ,4HC   ,4HM   ,4HW   ,
     8  4HB   ,4HC   ,4HL   ,4HA   ,4HW   ,4HB   ,4HC   ,4HM   ,
     9  4HA   ,4HW   ,4HB   ,4HQ   ,4HO   ,4HQ   ,4HI   ,4HN   ,
     A  4HF   ,4HE   ,4HP   ,4HS   ,4HL   ,4HO   ,4HN   ,4HD   ,
     B  4HE   ,4HO   ,4HD   ,4HA   ,4HC   ,4HD   ,4HV   ,4HA   ,
     C  4HL   ,4HP   ,4HO   ,4HW   ,4HA   ,4HL   ,4HP   ,4HL   ,
     D  4HW   ,4HA   ,4HL   ,4HP   ,4HO   ,4HH   ,4HA   ,4HL   ,
     E  4HP   ,4HL   ,4HH   ,4HA   ,4HC   ,4HL   ,4HM   ,4HW   ,
     F  4HC   ,4HL   ,4HM   ,4HW   ,4HA   ,4HC   ,4HL   ,4HM   ,
     G  4HH   ,4HC   ,4HL   ,4HM   ,4HH   /
C
      DATA LENI / 6*6,3*5,4,2,2,4,5,6,6,4,4,3,8*6 /
      DATA LDMI / 13*1,5*10,1,10,1,4*10,2*1 /
      DATA SYMFLP /   4HC   ,4HH   ,4HR   ,4HD   ,4HF   ,4HI   ,
     1  4HC   ,4HH   ,4HR   ,4HD   ,4HF   ,4HO   ,4HS   ,4HP   ,
     2  4HA   ,4HN   ,4HF   ,4HI   ,4HS   ,4HP   ,4HA   ,4HN   ,
     3  4HF   ,4HO   ,4HN   ,4HD   ,4HE   ,4HL   ,4HT   ,4HA   ,
     4  4HP   ,4HH   ,4HE   ,4HT   ,4HE   ,4HP   ,4HP   ,4HH   ,
     5  4HE   ,4HT   ,4HE   ,4HF   ,4HT   ,4HY   ,4HP   ,4HE   ,
     6  4HN   ,4HT   ,4HY   ,4HP   ,4HE   ,4HS   ,4HC   ,4HH   ,
     7  4HA   ,4HC   ,4HB   ,4HT   ,4HC   ,4HS   ,4HC   ,4HH   ,
     8  4HD   ,4HD   ,4HE   ,4HL   ,4HT   ,4HA   ,4HC   ,4HP   ,
     9  4HR   ,4HM   ,4HE   ,4HI   ,4HC   ,4HP   ,4HR   ,4HM   ,
     A  4HE   ,4HO   ,4HS   ,4HC   ,4HL   ,4HD   ,4HS   ,4HC   ,
     B  4HM   ,4HD   ,4HC   ,4HM   ,4HU   ,4HD   ,4HE   ,4HL   ,
     C  4HJ   ,4HE   ,4HT   ,4HJ   ,4HE   ,4HT   ,4HF   ,4HL   ,
     D  4HP   ,4HE   ,4HF   ,4HF   ,4HJ   ,4HE   ,4HT   ,4HC   ,
     E  4HA   ,4HP   ,4HI   ,4HN   ,4HB   ,4HC   ,4HA   ,4HP   ,
     F  4HO   ,4HU   ,4HT   ,4HD   ,4HO   ,4HB   ,4HD   ,4HE   ,
     G  4HF   ,4HD   ,4HO   ,4HB   ,4HC   ,4HI   ,4HN   ,4HD   ,
     H  4HO   ,4HB   ,4HC   ,4HO   ,4HT   /
C
      DATA LENJ / 4*6,4,4,5,6*6,5 /
      DATA LDMJ / 6*10,8*1 /
      DATA ASYFLP /   4HD   ,4HE   ,4HL   ,4HT   ,4HA   ,4HL   ,
     1  4HD   ,4HE   ,4HL   ,4HT   ,4HA   ,4HR   ,4HD   ,4HE   ,
     2  4HL   ,4HT   ,4HA   ,4HD   ,4HD   ,4HE   ,4HL   ,4HT   ,
     3  4HA   ,4HS   ,4HX   ,4HS   ,4HO   ,4HC   ,4HH   ,4HS   ,
     4  4HO   ,4HC   ,4HS   ,4HT   ,4HY   ,4HP   ,4HE   ,4HX   ,
     5  4HS   ,4HP   ,4HR   ,4HM   ,4HE   ,4HN   ,4HD   ,4HE   ,
     6  4HL   ,4HT   ,4HA   ,4HC   ,4HH   ,4HR   ,4HD   ,4HF   ,
     7  4HI   ,4HC   ,4HH   ,4HR   ,4HD   ,4HF   ,4HO   ,4HS   ,
     8  4HP   ,4HA   ,4HN   ,4HF   ,4HI   ,4HS   ,4HP   ,4HA   ,
     9  4HN   ,4HF   ,4HO   ,4HP   ,4HH   ,4HE   ,4HT   ,4HE   /
C
      DATA LENK / 5,3,5,2,5,6,6 /
      DATA LDMK / 4*1,-1,1,10 /
      DATA HYPEFF /   4HA   ,4HL   ,4HI   ,4HT   ,4HD   ,4HX   ,
     1  4HH   ,4HL   ,4HT   ,4HW   ,4HO   ,4HT   ,4HI   ,4HC   ,
     2  4HF   ,4HL   ,4HA   ,4HM   ,4HN   ,4HR   ,4HH   ,4HN   ,
     3  4HD   ,4HL   ,4HT   ,4HA   ,4HH   ,4HD   ,4HE   ,4HL   ,
     4  4HT   ,4HA   /
C
      DATA LENL / 4,2,5,6,2,2,3,4,3,2,2,3 /
      DATA LDML / 3*10,-10,8*1 /
      DATA TRNJET /   4HT   ,4HI   ,4HM   ,4HE   ,4HF   ,4HC   ,
     1  4HA   ,4HL   ,4HP   ,4HH   ,4HA   ,4HL   ,4HA   ,4HM   ,
     2  4HN   ,4HR   ,4HJ   ,4HN   ,4HT   ,4HM   ,4HE   ,4HI   ,
     3  4HS   ,4HP   ,4HS   ,4HP   ,4HA   ,4HN   ,4HP   ,4HH   ,
     4  4HE   ,4HG   ,4HP   ,4HC   ,4HC   ,4HL   ,4HF   ,4HP   /
C
      DATA LENM / 5*5,4*4,7*2,2*5,2,2,3,4 /
      DATA LDMM / 22*1 /
      DATA CONTAB /   4HT   ,4HT   ,4HY   ,4HP   ,4HE   ,4HC   ,
     1  4HF   ,4HI   ,4HT   ,4HC   ,4HC   ,4HF   ,4HO   ,4HT   ,
     2  4HC   ,4HC   ,4HF   ,4HI   ,4HT   ,4HT   ,4HC   ,4HF   ,
     3  4HO   ,4HT   ,4HT   ,4HB   ,4HI   ,4HT   ,4HC   ,4HB   ,
     4  4HO   ,4HT   ,4HC   ,4HB   ,4HI   ,4HT   ,4HT   ,4HB   ,
     5  4HO   ,4HT   ,4HT   ,4HB   ,4H1   ,4HB   ,4H2   ,4HB   ,
     6  4H3   ,4HB   ,4H4   ,4HD   ,4H1   ,4HD   ,4H2   ,4HD   ,
     7  4H3   ,4HG   ,4HC   ,4HM   ,4HA   ,4HX   ,4H    ,4H    ,
     8  4H    ,4H    ,4H    ,4HK   ,4HS   ,4HR   ,4HL   ,4HB   ,
     9  4HG   ,4HR   ,4HD   ,4HE   ,4HL   ,4HR   /
C
      END
