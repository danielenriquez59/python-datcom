      SUBROUTINE RVALUE(KOL,L,NDML,NF,BLANK,COMMA,NUMBER)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      DIMENSION KOL(80), NUMBER(15)
      INTEGER BLANK, COMMA, FIRST, POWR
      LOGICAL END, SIGN, EXP, STAR, DEC, FLAG
C
 1010 CONTINUE
      END  = .FALSE.
      SIGN = .FALSE.
      EXP  = .FALSE.
      STAR = .FALSE.
      DEC  = .FALSE.
      MULT = 0
      POWR = 0
 1020 CONTINUE
      FIRST = L
      IF(L .GE. 81) GO TO 1140
      IF(KOL(L) .NE. BLANK) GO TO 1030
        L = L+1
        GO TO 1020
 1030 CONTINUE
      DO 1040 J=1,15
      IF(KOL(L) .EQ. NUMBER(J)) GO TO 1050
 1040 CONTINUE
      J = 16
 1050 CONTINUE
      IF(J .GT. 10) GO TO 1060
        IF(.NOT. STAR) MULT = 10*MULT+J-1
        IF(EXP)        POWR = 10*POWR+J-1
        GO TO 1120
 1060 CONTINUE
      IF(J .GT. 12) GO TO 1070
        IF(SIGN) NF = NF+1
        FLAG = (L .EQ. FIRST) .OR. (KOL(L-1) .EQ. NUMBER(15)) .OR.
     1         (KOL(L-1) .EQ. NUMBER(14))
        IF(.NOT. FLAG) NF = NF+1
        SIGN = .TRUE.
        GO TO 1120
 1070 CONTINUE
      IF(J .GT. 13) GO TO 1080
        IF(DEC) NF = NF+1
        IF(EXP) NF = NF+1
        DEC = .TRUE.
        GO TO 1120
 1080 CONTINUE
      IF(J .GT. 14) GO TO 1090
        IF(STAR) NF = NF+1
        IF(SIGN) NF = NF+1
        IF(EXP)  NF = NF+1
        IF(DEC)  NF = NF+1
        IF(MULT .EQ. 0) NF = NF+1
        IF(MULT .EQ. 0) MULT = 1
        STAR = .TRUE.
        GO TO 1120
 1090 CONTINUE
      IF(J .GT. 15 .OR. L .EQ. FIRST) GO TO 1100
        IF(.NOT. DEC) NF = NF+1
        IF(EXP)       NF = NF+1
        EXP  = .TRUE.
        SIGN = .FALSE.
        GO TO 1120
 1100 CONTINUE
      IF(KOL(L) .NE. BLANK) GO TO 1110
        NF = NF+1
        GO TO 1120
 1110 CONTINUE
      IF(KOL(L) .EQ. KAND .AND. L .EQ. FIRST) GO TO 1140
      IF(KOL(L) .NE. COMMA .AND. KOL(L) .NE. KAND) GO TO 1140
        END = .TRUE.
C
C***  NUMERIC CHARACTER FOUND
C
 1120 CONTINUE
      IF(L .EQ. FIRST .AND. EXP) GO TO 1140
      IF(END) GO TO 1130
        L = L+1
        GO TO 1030
 1130 CONTINUE
        NDML = NDML+1
        IF(STAR) NDML = NDML-1+MULT
        IF(.NOT. DEC .AND. .NOT. EXP) NF = NF+1
        IF(EXP .AND. POWR .EQ. 0) NF = NF+1
        IF(KOL(L) .EQ. KAND) GO TO 1140
        L = L+1
        GO TO 1010
C
C***  NON-NUMERIC CHARACTER FOUND OR END OF CARD
C
 1140 CONTINUE
      IF(NDML .EQ. 0) NF = NF+1
      RETURN
      END
