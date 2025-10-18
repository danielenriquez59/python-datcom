      SUBROUTINE CCARD(KOL)
C
C***  CHECK FOR LEGAL CONTROL CARDS
C
      DIMENSION KOL(80), PART(4), IDIM(3), IDIMT(8), NMLIST(8)
      DIMENSION BUILD(5), NACA(4), DERD(9), DERR(9), KASE(6)
      DIMENSION SAVE(4), NEXT(9), IDMP(5), IDMC(4), TRIM(4)
      DIMENSION DAMP(4), NACAT(4), PLOT(4)
      DIMENSION KEYDP(224),  NDP(75)
      DIMENSION KEYARY(136), NDPARY(48)
      DIMENSION KEYINP(47),  NDPINP(13)
      DIMENSION KEYIOM(41),  NDPIOM(14)
      EQUIVALENCE (KEYDP( 1 ),KEYARY(1)), (NDP(1 ),NDPARY(1))
      EQUIVALENCE (KEYDP(137),KEYINP(1)), (NDP(49),NDPINP(1))
      EQUIVALENCE (KEYDP(184),KEYIOM(1)), (NDP(62),NDPIOM(1))
      INTEGER BUILD, DERD, DERR, SAVE, TRIM, DAMP, PART, PLOT
      LOGICAL NMTEST
      INTEGER EXTRAP
      DIMENSION EXTRAP(12)
C
C**   VALID CONTROL CARDS
C
      DATA EXTRAP /4HP   ,4HR   ,4HI   ,4HN   ,4HT   ,4H    ,
     1             4HE   ,4HX   ,4HT   ,4HR   ,4HA   ,4HP   /
      DATA BUILD  /4HB   ,4HU   ,4HI   ,4HL   ,4HD   /
      DATA NACA   /4HN   ,4HA   ,4HC   ,4HA   /
      DATA NACAT  /4HW   ,4HH   ,4HV   ,4HF   /
      DATA DERD   /4HD   ,4HE   ,4HR   ,4HI   ,4HV   ,4H    ,
     1             4HD   ,4HE   ,4HG   /
      DATA DERR   /4HD   ,4HE   ,4HR   ,4HI   ,4HV   ,4H    ,
     1             4HR   ,4HA   ,4HD   /
      DATA PART   /4HP   ,4HA   ,4HR   ,4HT   /
      DATA IDIM   /4HD   ,4HI   ,4HM    /
      DATA IDIMT  /4HF   ,4HT   ,4HI   ,4HN   ,4HM   ,4H    ,
     1             4HC   ,4HM   /
      DATA NMLIST /4HN   ,4HA   ,4HM   ,4HE   ,4HL   ,
     1             4HI   ,4HS   ,4HT   /
      DATA KASE   /4HC   ,4HA   ,4HS   ,4HE   ,4HI   ,4HD   /
      DATA SAVE   /4HS   ,4HA   ,4HV   ,4HE    /
      DATA NEXT   /4HN   ,4HE   ,4HX   ,4HT   ,4H    ,4HC   ,4HA   ,
     1             4HS   ,4HE   /
      DATA IDMP   /4HD   ,4HU   ,4HM   ,4HP   ,4H    /
      DATA IDMC   /4HC   ,4HA   ,4HS   ,4HE   /
      DATA TRIM   /4HT   ,4HR   ,4HI   ,4HM   /
      DATA DAMP   /4HD   ,4HA   ,4HM   ,4HP   /
      DATA PLOT   /4HP   ,4HL   ,4HO   ,4HT   /
C
      DATA IBLNK  /4H    /,KOMMA /4H,   /
      DATA KEYARY  /
     1       4HB   ,4HA   ,4HB   ,4HD   ,4HA   ,4HV   ,4HF   ,4HF   ,
     2       4HA   ,4HC   ,4HT   ,4HW   ,4HB   ,4HT   ,4HB   ,4HH   ,
     3       4HT   ,4HA   ,4HV   ,4HT   ,4HA   ,4HH   ,4HT   ,4HC   ,
     4       4HD   ,4HW   ,4HB   ,4HC   ,4HH   ,4HT   ,4HD   ,4HH   ,
     5       4HT   ,4HD   ,4HY   ,4HN   ,4HH   ,4HD   ,4HV   ,4HT   ,
     6       4HG   ,4HR   ,4HL   ,4HB   ,4HP   ,4HW   ,4HS   ,4HT   ,
     7       4HB   ,4HS   ,4HB   ,4HD   ,4HS   ,4HL   ,4HG   ,4HS   ,
     8       4HW   ,4HB   ,4HS   ,4HT   ,4HP   ,4HD   ,4HW   ,4HA   ,
     9       4HS   ,4HT   ,4HG   ,4HS   ,4HL   ,4HA   ,4HT   ,4HR   ,
     A       4HA   ,4HX   ,4HP   ,4HD   ,4HD   ,4HV   ,4HF   ,4HF   ,
     B       4HL   ,4HP   ,4HF   ,4HH   ,4HG   ,4HF   ,4HC   ,4HM   ,
     C       4HT   ,4HC   ,4HD   ,4HF   ,4HL   ,4HA   ,4HT   ,4HR   ,
     D       4HM   ,4HS   ,4HP   ,4HR   ,4HT   ,4HR   ,4HN   ,4HT   ,
     E       4HR   ,4HM   ,4H2   ,4HH   ,4HY   ,4HP   ,4HD   ,4HY   ,
     F       4HN   ,4HJ   ,4HE   ,4HT   ,4HH   ,4HB   ,4HS   ,4HH   ,
     G       4HB   ,4HT   ,4HR   ,4HA   ,4HH   ,4HS   ,4HT   ,4HB   ,
     H       4HH   ,4HS   ,4HE   ,4HC   ,4HS   ,4HL   ,4HA   ,4HH   /
      DATA NDPARY / 1,1,2,3,4,4*3,2*1,2,2*3,4,3,3*2,19*3,4,3*3,2,3,
     1              4,4,3,4 /
C
      DATA KEYINP  /
     1       4HI   ,4HN   ,4HP   ,4HT   ,4HF   ,4HL   ,4HC   ,4HO   ,
     2       4HP   ,4HT   ,4HN   ,4HS   ,4HY   ,4HN   ,4HA   ,4HB   ,
     3       4HD   ,4HI   ,4HN   ,4HW   ,4HG   ,4HI   ,4HN   ,4HV   ,
     4       4HT   ,4HI   ,4HN   ,4HT   ,4HV   ,4HT   ,4HV   ,4HF   ,
     5       4HI   ,4HN   ,4HH   ,4HT   ,4HI   ,4HN   ,4HP   ,4HW   ,
     6       4HI   ,4HN   ,4HL   ,4HB   ,4HI   ,4HN   ,4HF   /
      DATA NDPINP / 4,3,5*4,3,4*4,1 /
C
      DATA KEYIOM  /
     1       4HI   ,4HO   ,4HM   ,4HB   ,4HO   ,4HD   ,4HY   ,4HW   ,
     2       4HI   ,4HN   ,4HG   ,4HH   ,4HT   ,4HV   ,4HT   ,4HV   ,
     3       4HF   ,4HB   ,4HW   ,4HB   ,4HH   ,4HB   ,4HV   ,4HB   ,
     4       4HW   ,4HH   ,4HB   ,4HW   ,4HV   ,4HB   ,4HW   ,4HH   ,
     5       4HV   ,4HP   ,4HO   ,4HW   ,4HR   ,4HD   ,4HW   ,4HS   ,
     6       4HH   /
      DATA NDPIOM / 3,4,4,6*2,3,3,4,4,4 /
C
      NA = 0
      NB = 0
C
C***  TEST FOR - NEXT CASE
C
      DO 1000 I=1,9
         IF(KOL(I) .NE. NEXT(I)) GO TO 1010
 1000 CONTINUE
      GO TO 1380
C
C***  TEST FOR - DUMP
C
 1010 DO 1020 I=1,5
         IF(KOL(I) .NE. IDMP(I)) GO TO 1100
 1020 CONTINUE
C
C***  TEST FOR - DUMP CASE
C
      J=1
      DO 1030 L=6,9
         IF(KOL(L) .NE. IDMC(J)) GO TO 1040
 1030 J=J+1
      GO TO 1380
C
C***  TEST FOR ARRAY NAMES
C
 1040 I=5
 1050 I=I+1
      IF(I.GT.80) GO TO 1380
      IF(KOL(I) .EQ. IBLNK .OR. KOL(I) .EQ. KOMMA)GO TO 1050
      NS=I
      NCHAR=1
 1060 I=I+1
      IF(I.GT.80)GO TO 1380
      IF(KOL(I) .EQ. IBLNK .OR. KOL(I) .EQ. KOMMA) GO TO 1070
      NCHAR=NCHAR+1
      GO TO 1060
 1070 IF(NCHAR .LT. 7) GO TO 1080
      NA = NA + 1
      GO TO 1050
 1080 LOC = 1
      DO 1090 J=1,75
         IF(J .GT. 1) LOC = LOC + NDP(J-1)
         IF(NDP(J) .NE. NCHAR) GO TO 1090
         IF(NMTEST(KOL(NS),KEYDP(LOC),NDP(J))) GO TO 1050
 1090 CONTINUE
      NA = NA + 1
      GO TO 1050
C
C***  TEST FOR - SAVE
C
 1100 DO 1110 I=1,4
         IF(KOL(I) .NE. SAVE(I)) GO TO 1120
 1110 CONTINUE
      GO TO 1380
C
C***  TEST FOR - CASEID
C
 1120 DO 1130 I=1,6
         IF(KOL(I) .NE. KASE(I)) GO TO 1140
 1130 CONTINUE
      GO TO 1380
C
C***  TEST FOR - TRIM
C
 1140 DO 1150 I=1,4
         IF(KOL(I) .NE. TRIM(I)) GO TO 1160
 1150 CONTINUE
      GO TO 1380
C
C***  TEST FOR - PRINT EXTRAP
C
 1160 DO 1170 I=1,12
         IF(KOL(I) .NE. EXTRAP(I)) GO TO 1180
 1170 CONTINUE
      GO TO 1380
C
C***  TEST FOR - DAMP
C
 1180 DO 1190 I=1,4
         IF(KOL(I) .NE. DAMP(I)) GO TO 1200
 1190 CONTINUE
      GO TO 1380
C
C***  TEST FOR - BUILD
C
 1200 DO 1210 I=1,5
         IF(KOL(I) .NE. BUILD(I)) GO TO 1220
 1210 CONTINUE
      GO TO 1380
C
C***  TEST FOR - NACA
C
 1220 DO 1230 I=1,4
         IF(KOL(I) .NE. NACA(I)) GO TO 1250
 1230 CONTINUE
      DO 1240 I=1,4
         IF(KOL(6) .EQ. NACAT(I)) GO TO 1380
 1240 CONTINUE
      NB = 1
      GO TO 1380
C
C***  TEST FOR - DERIV DEG
C
 1250 DO 1260 I=1,9
         IF(KOL(I) .NE. DERD(I)) GO TO 1270
 1260 CONTINUE
      GO TO 1380
C
C***  TEST FOR - DERIV RAD
C
 1270 DO 1280 I=1,9
         IF(KOL(I) .NE. DERR(I)) GO TO 1290
 1280 CONTINUE
      GO TO 1380
C
C***  TEST FOR - PART
C
 1290 DO 1300 I=1,4
         IF(KOL(I) .NE. PART(I)) GO TO 1310
 1300 CONTINUE
      GO TO 1380
C
C***  TEST FOR DIMENSION CARD
C
 1310 DO 1320 I=1,3
         IF(KOL(I) .NE. IDIM(I)) GO TO 1330
 1320 CONTINUE
      IF((KOL(5) .EQ. IDIMT(1) .AND. KOL(6) .EQ. IDIMT(2)) .OR.
     1   (KOL(5) .EQ. IDIMT(3) .AND. KOL(6) .EQ. IDIMT(4)) .OR.
     2   (KOL(5) .EQ. IDIMT(5) .AND. KOL(6) .EQ. IDIMT(6)) .OR.
     3   (KOL(5) .EQ. IDIMT(7) .AND. KOL(6) .EQ. IDIMT(8)))
     4   GO TO 1380
C
C***  TEST FOR NAMELIST CARD
C
 1330 DO 1340 I=1,8
         IF(KOL(I) .NE. NMLIST(I)) GO TO 1350
 1340 CONTINUE
      GO TO 1380
C
C***  TEST FOR PLOT CARD
C
 1350 DO 1360 I=1,4
         IF(KOL(I) .NE. PLOT(I)) GO TO 1370
 1360 CONTINUE
      GO TO 1380
C
C***   CARD NOT FOUND
C
 1370 WRITE(6,1400) (KOL(I),I=1,80)
      GO TO 1390
C
C***  CARD FOUND
C
 1380 NC = NA+NB
      IF(NA .GT. 0) WRITE(6,1410) (KOL(I),I=1,80), NA
      IF(NB .GT. 0) WRITE(6,1420) (KOL(I),I=1,80)
      IF(NC .EQ. 0) WRITE(6,1430) (KOL(I),I=1,80)
 1390 WRITE(11,1440) (KOL(I),I=1,80)
C
      RETURN
C
C***  FORMAT
C
 1400 FORMAT(1X,80A1,3X,23H** ILLEGAL CONTROL CARD)
 1410 FORMAT(1X,80A1,3X,11H** ERROR **,I4,
     1                  22H INCORRECT ARRAY NAMES)
 1420 FORMAT(1X,80A1,3X,29H** ERROR ** INCORRECT LIFTING,/84X,
     1                    33H SURFACE DESIGNATION ON NACA CARD)
 1430 FORMAT(1X,80A1)
 1440 FORMAT(80A1)
      END
