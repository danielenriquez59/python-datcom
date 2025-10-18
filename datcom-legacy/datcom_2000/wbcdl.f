      SUBROUTINE WBCDL(AR,TANLE,TOVC,LN,LA,TAPR,LER,TWIST,YCM,CLD,RN,
     1      TR,MACH,ALPHA,CDL,NA)
C
C *** CALCULATES THE WING-BODY AND  HORIZONTAL TAIL-BODY CDL
C
      COMMON /OVERLY/ NLOG, NMACH, I, NALPHA
      COMMON /CONSNT/ PI, DEG, UNUSED, RAD
C
      EQUIVALENCE (B(1),B0), (B(2),B1), (B(3),B2), (B(4),B3),
     1            (B(5),B4), (B(6),B5), (B(7),B6), (B(8),B7),
     2            (B(9),B8), (B(10),B9), (B(11),B10), (B(12),B11),
     3            (B(13),B12), (B(14),B13), (B(15),B14),
     4            (B(16),B15)
      DIMENSION ALPHA(20), CDL(20), B(16)
      REAL LA, LN, LER, MACH
      LOGICAL NA, NDM
C
      NA = .FALSE.
      IF(AR .LT. 1.6 .OR. AR .GT. 6.0)            GO TO 1020
      IF(TANLE .LT. 0.0 .OR. TANLE .GT. 2.74748)  GO TO 1020
      IF(TOVC .LT. 0.025 .OR. TOVC .GT. 0.100)    GO TO 1020
      IF(LN .LT. 2.2 .OR. LN .GT. 8.4)            GO TO 1020
      IF(LA .LT. 0.3 .OR. LA .GT. 5.6)            GO TO 1020
      IF(TAPR .LT. 0.0 .OR. TAPR .GT. 1.0)        GO TO 1020
      IF(LER .LT. 0.0 .OR. LER .GT. 0.015)        GO TO 1020
      IF(TWIST .LT. -9.4 .OR. TWIST .GT. UNUSED)  GO TO 1020
      IF(YCM .LT. 0.0 .OR. YCM .GT. 0.0263)       GO TO 1020
      IF(CLD .LT. 0.0 .OR. CLD .GT. 0.45)         GO TO 1020
        IF(RN .LT. 8.0E5) RN = 8.0E5
        IF(RN .GT. 8.0E6) RN = 8.0E6
        DO 1010 J=1,NALPHA
           CALL TABLES(B, MACH, ALPHA(J), NDM)
           CDL(J) = B0 + B1/AR + B2*AR + B3*SQRT(TANLE) + B4*TOVC
     1            + B5*LN + B6*LA + B7*TAPR + B8*TAPR**2 + B9*TAPR**3
     2            + B10*TR + B11*LER - B12*TWIST/RAD + B13*YCM
     3            + B14*CLD + B15*RN/1.E6
           IF(NDM) CDL(J) = UNUSED
           IF(NDM .AND. J .EQ. 1) NA = .TRUE.
 1010   CONTINUE
      GO TO 1030
 1020 CONTINUE
        NA = .TRUE.
 1030 CONTINUE
      RETURN
      END
