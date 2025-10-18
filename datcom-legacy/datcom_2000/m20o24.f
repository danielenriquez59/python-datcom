      SUBROUTINE M20O24
C
C***  EXEC FOR OVERLAY 20, SUPERSONIC WING-BODY AERO
C
      COMMON /IBODY/   PBODY,  BODY(400)
      COMMON /IWING/   PWING,  WING(400)
      COMMON /IHT/     PHT,    HT(380)
      COMMON /IVT/     PVT,    VT(380)
      COMMON /IVF/     PVF,    VF(380)
      COMMON /IBW/     PBW,    BW(380)
      COMMON /IBH/     PBH,    BH(380)
      COMMON /IBV/     PBV,    BV(380)
      COMMON /IBWV/    PBWV,   BWV(380)
C
      COMMON /FLGTCD/ FLC(160)
      COMMON /WBHCAL/ STP(156)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG,NF,LF,L,NOVLY
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
C
      DIMENSION CR(6), DR(6), ROUTID(2)
      DATA ROUTID / 4HM20O, 4H24   /
C
      NOVLY=20
      IF(WGPL .AND. BO)CALL SUPWB
      IF(HTPL .AND. BO) CALL SUPHB
      DO 1000 J=1,NALPHA
 1000 BWV(J) = BW(J)
      IF(VTPL) CALL VRTCDO(I)
      IF(VFPL) CALL VFCDO(I)
      CALL SUPCM0
C
C     VT DATA
C
      IF(VT(1) .EQ. UNUSED) VT(1) = STP(1)
      VT(21) = 0.0
      VT(41) = 0.0
      VT(61) = 0.0
      VT(81) = VT(1)
      VT(101) = 0.0
      VT(121) = 0.0
      IF(VF(1) .EQ. UNUSED) VF(1) = STP(156)
      VF(21) = 0.0
      VF(41) = 0.0
      VF(61) = 0.0
      VF(81) = VF(1)
      VF(101) = 0.0
      VF(121) = 0.0
      DO 1010 J=1,NALPHA
         IF(J .EQ. 1) GO TO 1010
         VT(J) = -UNUSED
         VT(J+20) = -UNUSED
         VT(J+40) = -UNUSED
         VT(J+60) = -UNUSED
         VT(J+80) = -UNUSED
         VT(J+100) = -UNUSED
         VT(J+120) = -UNUSED
         VF(J) = -UNUSED
         VF(J+20) = -UNUSED
         VF(J+40) = -UNUSED
         VF(J+60) = -UNUSED
         VF(J+80) = -UNUSED
         VF(J+100) = -UNUSED
         VF(J+120) = -UNUSED
 1010 CONTINUE
      IF( .NOT. BO) GO TO 1060
      IF(.NOT. WGPL) GO TO 1040
C
C     B-W  DATA
C
      IN = 0
      DO 1020 J=1,NALPHA
         CA = COS(FLC(J+22)/RAD)
         SA = SIN(FLC(J+22)/RAD)
         BW(J+60) = BW(J+20)*CA + BW(J)*SA
         BW(J+80) = BW(J)*CA - BW(J+20)*SA
         IF(J .EQ. 1) GO TO 1020
        CALL TBFUNX(FLC(J+22),X,BW(J+100),NALPHA,FLC(23),BW(21),
     1               CR,IN,MI,NG,0,0,4HCLA ,1,ROUTID)
         BW(J+120) = -UNUSED
 1020 CONTINUE
C
C     B-W-V  DATA
C
      IN = 0
      DO 1030 J=1,NALPHA
         CA = COS(FLC(J+22)/RAD)
         SA = SIN(FLC(J+22)/RAD)
         BWV(J+20) = BW(J+20)
         BWV(J+40) = BW(J+40)
         BWV(J+60) = BWV(J+20)*CA + BWV(J)*SA
         BWV(J+80) = BWV(J)*CA - BWV(J+20)*SA
         BWV(J+100)= BW(J+100)
         BWV(J+120)= BW(J+120)
 1030 CONTINUE
 1040 IF(.NOT. HTPL) GO TO 1060
C
C     B-H  DATA
C
      IN = 0
      DO 1050 J=1,NALPHA
         CA = COS(FLC(J+22)/RAD)
         SA = SIN(FLC(J+22)/RAD)
         BH(J+60) = BH(J+20)*CA + BH(J)*SA
         BH(J+80) = BH(J)*CA - BH(J+20)*SA
         IF(J .EQ. 1) GO TO 1050
         CALL TBFUNX(FLC(J+22),X,BH(J+100),NALPHA,FLC(23),BH(21),
     1               CR,IN,MI,NG,0,0,4HCLA ,1,ROUTID)
         BH(J+120) = -UNUSED
 1050 CONTINUE
 1060 IF(.NOT. (BO .AND. VTPL)) GO TO 1080
C
C     B-V DATA
C
      DO 1070 J=1,NALPHA
         DO 1070 K=1,7
            KK=(K-1)*20
            IF(K .EQ. 1) BV(J+KK)=BODY(J+KK)+VT(1) + VF(1)
            IF(K .NE. 1) BV(J+KK)=BODY(J+KK)
 1070 CONTINUE
 1080 CONTINUE
      RETURN
      END
