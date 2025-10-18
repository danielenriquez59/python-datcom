      SUBROUTINE WBCM0(AR,TANLE,TOVC,LN,LA,TAPR,LER,TWIST,YCM,CLD,RN,
     1                 TR,WL,VT,HD,DB,MACH,CM0)
C
C *** CALCULATES CM0 FOR WING-BODY USING REGRESSION METHOD
C
      COMMON /CONSNT/ PI, DEG, UNUSED, RAD
C
      EQUIVALENCE (C(1),C0), (C(2),C1), (C(3),C2), (C(4),C3),
     1            (C(5),C4), (C(6),C5), (C(7),C6), (C(8),C7),
     2            (C(9),C8), (C(10),C9), (C(11),C10), (C(12),C11),
     3            (C(13),C12), (C(14),C13), (C(15),C14),
     4            (C(16),C15), (C(17),C16), (C(18),C17), (C(19),C18)
      DIMENSION C(19)
      REAL LA, LN, LER, MACH
C
      IF(MACH .GT. 2.5)                           GO TO 1010
      IF(AR .LT. 1.6 .OR. AR .GT. 6.0)            GO TO 1010
      IF(TANLE .LT. 0.0 .OR. TANLE .GT. 2.74748)  GO TO 1010
      IF(TOVC .LT. 0.025 .OR. TOVC .GT. 0.100)    GO TO 1010
      IF(LN .LT. 2.2 .OR. LN .GT. 8.4)            GO TO 1010
      IF(LA .LT. 0.3 .OR. LA .GT. 5.6)            GO TO 1010
      IF(TAPR .LT. 0.0 .OR. TAPR .GT. 1.0)        GO TO 1010
      IF(LER .LT. 0.0 .OR. LER .GT. 0.015)        GO TO 1010
      IF(TWIST .LT. -9.4 .OR. TWIST .GT. UNUSED)  GO TO 1010
      IF(YCM .LT. 0.0 .OR. YCM .GT. 0.0263)       GO TO 1010
      IF(CLD .LT. 0.0 .OR. CLD .GT. 0.45)         GO TO 1010
        IF(RN .LT. 8.0E5) RN = 8.0E5
        IF(RN .GT. 8.0E6) RN = 8.0E6
        CALL TABLEC(C,MACH)
        CM0 = C0 + C1/AR + C2*AR + C3*TANLE + C4*TOVC
     1        + C5*LN + C6*LA + C7*TAPR + C8*TAPR**2 + C9*TR
     2        + C10*LER + C11*TWIST + C12*YCM + C13*CLD
     3        + C14*WL + C15*VT + C16*HD + C17*DB + C18*RN/1.0E6
 1010 CONTINUE
      RETURN
      END
