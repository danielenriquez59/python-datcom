      SUBROUTINE WBTRA
C
C *** CALCULATES TRANSONIC WING-BODY CDL
C
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG
      COMMON /EXPER/  KLIST, NLIST(100), NNAMES, IMACH, MDATA,
     1                KBODY, KWING, KHT, KVT, KWB,
     2                KDWASH(3), ALPO, ALPL
      COMMON /CONSNT/ CONST(4)
      COMMON /IDWASH/ PDWASH,DWASH(60)
      COMMON /IWING/  PWING,WING(400)
      COMMON /IBODY/  PBODY,BODY(400)
      COMMON /IVT/    PVT,VT(380)
      COMMON /IVF/    PVF,VF(380)
      COMMON /IHT/    PHT,HT(380)
      COMMON /IBH/    PBH,BH(380)
      COMMON /IBW/    PBW,BW(380)
      COMMON /IBV/    PBV,BV(380)
      COMMON /IBWH/   PBWH,BWH(380)
      COMMON /IBWV/   PBWV,BWV(380)
      COMMON /IBWHV/  PBWHV,BWHV(380)
      DIMENSION KEY(3)
      LOGICAL KDWASH
      EQUIVALENCE (KEY(1),KQOQIN), (KEY(2),KEPSLN), (KEY(3),KDEODA)
      EQUIVALENCE(UNUSED,CONST(3))
C
      CALL EXSUBT
      DO 1000 J=1,3
         KEY(J) = 0
         IF(KDWASH(J)) KEY(J) = 1
 1000 CONTINUE
      CALL TRAWBT(KDEODA,KQOQIN)
C
C     SET IOM FOR DOWNWASH RESULTS AND HT, VT, B-W, B-W-H, B-W-V,
C                      AND B-W-H-V CL-ALPHA
C
      DWASH(21) = 0.0
      BWH(101) = BWHV(101)
      VT(21)   = 0.0
      VT(41)   = 0.0
      VT(61)   = 0.0
      VT(101)  = 0.0
      VT(121)  = 0.0
      VF(21)   = 0.0
      VF(41)   = 0.0
      VF(61)   = 0.0
      VF(101)  = 0.0
      VF(121)  = 0.0
      BWV(101) = BW(101)
      BWV(121) = BW(121)
      BV(101) = BODY(101)
      BV(121) = BODY(121)
      DO 1010 J=2,NALPHA
         DWASH(J)    = -UNUSED
         DWASH(J+20) = -UNUSED
         DWASH(J+40) = -UNUSED
 1010 CONTINUE
      DO 1040 II=1,9
         INDX = 20*(II-1)
         DO 1030 J=2,NALPHA
            WING(INDX+J) = -UNUSED
            HT(INDX+J)   = -UNUSED
            VT(INDX+J)   = -UNUSED
            VF(INDX+J)   = -UNUSED
            BV(INDX+J)   = -UNUSED
            BWH(INDX+J)  = -UNUSED
            IF(II .NE. 3) BWV(INDX+J)  = -UNUSED
            BWHV(INDX+J) = -UNUSED
            IF(II .EQ. 1) GO TO 1020
             BW(INDX+J) = -UNUSED
             BH(INDX+J) = -UNUSED
 1020       CONTINUE
            IF(II .NE. 1) GO TO 1030
            IF(BW(INDX+J) .EQ. UNUSED) BW(INDX+J) = -UNUSED
            IF(BH(INDX+J) .EQ. UNUSED) BH(INDX+J) = -UNUSED
 1030    CONTINUE
 1040 CONTINUE
      CALL EXSUBT
      RETURN
      END
