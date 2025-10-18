      SUBROUTINE WBCD
C
C *** EXEC FOR WING-BODY AND HORIZONTAL TAIL-BODY CD
C
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /FLGTCD/ FLC(95), TR
      COMMON /BDATA/  BD(762)
      COMMON /WINGD/  A(195)
      COMMON /WHWB/   FACT(182), WB(39), HB(39)
      COMMON /IHT/    PHT,HT(380)
      COMMON /HTI/    HTIN(131)
      COMMON /HTDATA/ AHT(195)
      COMMON /WINGI/  WINGIN(100)
      COMMON /IBH/    PBH,BH(380)
      COMMON /SYNTSS/ XCG, XW, ZW, ALIW, ZCG, XH, ZH, ALIH, XV,
     1                VERTUP, HINAX, XVF, SCALE, ZV, ZVF, YV, YF,
     2                PHIV, PHIF
      COMMON /IBW/    PBW, BW(380)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST
      REAL LN, LA
      LOGICAL NA
      DIMENSION CDL(20)
      DATA STRA / 4HSTRA /
C
C***  WING-BODY DRAG
C
      IF(WINGIN(15) .NE. STRA .OR. .NOT. WGPL) GO TO 1020
      CD0 = WB(17)
      RN = FLC(I+42)*A(122)
      DB = 2.0*(WINGIN(4)-WINGIN(3))
      LN = (XW + 0.5*DB*A(38))/DB
      LA = (XW + WINGIN(6) + 0.5*DB*A(56))/DB
      IF(BD(2) .LT. (LN*DB)) LN = LN*DB/BD(85)
      IF(BD(2) .GT. (LA*DB)) LA = LA*DB/BD(85)
      LA = BD(1)/DB-LA
      CALL WBCDL(A(120), A(38), WINGIN(16), LN, LA, A(118),
     1          WINGIN(62), WINGIN(11), WINGIN(93), WINGIN(94), RN,
     1          TR, FLC(I+2), FLC(23), CDL, NA)
      IF(NA) GO TO 1020
        CALL EXSUBT
        DO 1010 J=1,NALPHA
           BW(J) = CD0 + CDL(J)
           CA = COS(FLC(J+22)/RAD)
           SA = SIN(FLC(J+22)/RAD)
           BW(J+60) = BW(J+20)*CA + BW(J)*SA
           BW(J+80) = BW(J)*CA - BW(J+20)*SA
           IF(CDL(J) .NE. UNUSED) GO TO 1010
           BW(J)    = -UNUSED
           BW(J+60) = -UNUSED
           BW(J+80) = -UNUSED
 1010   CONTINUE
 1020 CONTINUE
C
C***  H.T.-BODY DRAG
C
      IF(HTIN(15) .NE. STRA .OR. .NOT. HTPL) GO TO 1040
      CD0 = HB(17)
      RN = FLC(I+42)*AHT(122)
      DB = 2.0*(HTIN(4)-HTIN(3))
      LN = (XH + 0.5*DB*AHT(38))/DB
      LA = (XH + HTIN(6) + 0.5*DB*AHT(56))/DB
      IF(BD(2) .LT. (LN*DB)) LN = LN*DB/BD(85)
      IF(BD(2) .GT. (LA*DB)) LA = LA*DB/BD(85)
      LA = BD(1)/DB-LA
      CALL WBCDL(AHT(120), AHT(38), HTIN(16), LN, LA, AHT(118),
     1          HTIN(62), HTIN(11), HTIN(93), HTIN(94), RN,
     1          TR, FLC(I+2), FLC(23), CDL, NA)
      IF(NA) GO TO 1040
        CALL EXSUBT
        DO 1030 J=1,NALPHA
           BH(J) = CD0 + CDL(J)
           CA = COS(FLC(J+22)/RAD)
           SA = SIN(FLC(J+22)/RAD)
           BH(J+60) = BH(J+20)*CA + BH(J)*SA
           BH(J+80) = BH(J)*CA - BH(J+20)*SA
           IF(CDL(J) .NE. UNUSED) GO TO 1030
           BH(J)    = -UNUSED
           BH(J+60) = -UNUSED
           BH(J+80) = -UNUSED
 1030   CONTINUE
 1040 CONTINUE
      RETURN
      END
