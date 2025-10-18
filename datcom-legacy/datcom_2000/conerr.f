      SUBROUTINE CONERR
C
C***  EXECUTIVE FOR INPUT ERROR DIAGNOSTICS
C
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      DIMENSION KOL(80), NUM(10)
      DIMENSION NLNAME(138), LOC(24), LEN(24)
      LOGICAL NMTEST
C
      DATA NNAME / 24 /, IBLANK / 4H     /
      DATA NUM /4H0   ,4H1   ,4H2   ,4H3   ,4H4   ,4H5   ,4H6   ,
     1          4H7   ,4H8   ,4H9   /
      DATA LOC /1,4,10,16,20,26,32,38,44,50,56,62,68,74,79,85,91,97,
     1          103,109,115,121,127,133/
      DATA LEN /3,2*6,4,9*6,5,2*6,4,7*6/
      DATA NLNAME /
     1 4HE   ,4HN   ,4HD   ,4HF   ,4HL   ,4HT   ,4HC   ,4HO   ,
     2 4HN   ,4HO   ,4HP   ,4HT   ,4HI   ,4HN   ,4HS   ,4HB   ,4HO   ,
     3 4HD   ,4HY   ,4HW   ,4HG   ,4HP   ,4HL   ,4HN   ,4HF   ,4HW   ,
     4 4HG   ,4HS   ,4HC   ,4HH   ,4HR   ,4HS   ,4HY   ,4HN   ,4HT   ,
     5 4HH   ,4HS   ,4HH   ,4HT   ,4HP   ,4HL   ,4HN   ,4HF   ,4HH   ,
     6 4HT   ,4HS   ,4HC   ,4HH   ,4HR   ,4HV   ,4HT   ,4HP   ,4HL   ,
     7 4HN   ,4HF   ,4HV   ,4HT   ,4HS   ,4HC   ,4HH   ,4HR   ,4HP   ,
     8 4HR   ,4HO   ,4HP   ,4HW   ,4HR   ,4HJ   ,4HE   ,4HT   ,4HP   ,
     9 4HW   ,4HR   ,4HL   ,4HA   ,4HR   ,4HW   ,4HB   ,4HG   ,4HR   ,
     A 4HN   ,4HD   ,4HE   ,4HF   ,4HT   ,4HV   ,4HT   ,4HP   ,
     B 4HA   ,4HN   ,4HE   ,4HX   ,4HP   ,4HR   ,4H    ,4H    ,4HS   ,
     C 4HY   ,4HM   ,4HF   ,4HL   ,4HP   ,4HA   ,4HS   ,4HY   ,4HF   ,
     D 4HL   ,4HP   ,4HH   ,4HY   ,4HP   ,4HE   ,4HF   ,4HF   ,4HT   ,
     E 4HR   ,4HN   ,4HJ   ,4HE   ,4HT   ,4HC   ,4HO   ,4HN   ,4HT   ,
     F 4HA   ,4HB   ,4HV   ,4HF   ,4HS   ,4HC   ,4HH   ,4HR   ,4HV   ,
     G 4HF   ,4HP   ,4HL   ,4HN   ,4HF   /
C
C***  READ ONE CARD AND CHECK TYPE -
C***           K = 1 FIRST CARD OF NAMELIST
C***           K = 2 CONTINUATION CARD IN NAMELIST
C***           K = 3 CONTROL CARD
C
      WRITE(6,1120)
C
C***  THE FOLLOWING TWO CARDS ARE FOR FORTRAN-IV
C
C1000 READ(5,1020) (KOL(I),I=1,80)
C     IF(EOF(5)) 1110,1030
C
C***  THE FOLLOWING CARD IS FOR FORTRAN-V
C
 1000 READ(5,1020,END=1110) (KOL(I),I=1,80)
 1020 FORMAT(80A1)
 1030 K = 2
      IF(KOL(2) .EQ. KAND)   K = 1
      IF(KOL(1) .NE. IBLANK) K = 3
C
C***  CHECK NAMELIST NAME IF K = 1
C
      L = 2
      IF(K .NE. 1) GO TO 1100
        IERR = 0
        DO 1080 NAME = 1,NNAME
           I = LOC(NAME)
           L = LEN(NAME)
           IF(.NOT. NMTEST(KOL(3),NLNAME(I),L)) GO TO 1080
           IF(NAME .NE. 17) GO TO 1070
           L=L+2
           M=L+1
           MNUM = 0
 1040      DO 1050 J=1,10
              IF(KOL(M) .EQ. NUM(J)) GO TO 1060
 1050      CONTINUE
           IERR = 1
           GO TO 1080
 1060      M=M+1
           MNUM = 10*MNUM + J
           IF(M .EQ. (L+2)) GO TO 1040
           IF(MNUM .LT. 1 .OR. MNUM .GT. 20) IERR=1
 1070      CONTINUE
           IF(KOL(L+3) .EQ. IBLANK) GO TO 1090
 1080   CONTINUE
        IERR = 1
 1090   L = L+3
 1100 CONTINUE
C
C***  CHECK NAMELIST VARIABLES IF K = 1 OR 2
C***  CHECK CONTROL CARD IF K = 3
C
      IF(K .LE. 2) CALL NMLIST(KOL,L,NAME,K,IERR)
      IF(K .EQ. 3) CALL CCARD (KOL)
      GO TO 1000
 1110 CONTINUE
      REWIND 11
 1120 FORMAT(1H1,25X,29HCONERR - INPUT ERROR CHECKING/
     1    1H0,39H ERROR CODES - N* DENOTES THE NUMBER OF,
     2        25H OCCURENCES OF EACH ERROR/
     3    1H0,26H A - UNKNOWN VARIABLE NAME/
     4    1H0,47H B - MISSING EQUAL SIGN FOLLOWING VARIABLE NAME/
     5    1H0,44H C - NON-ARRAY VARIABLE HAS AN ARRAY ELEMENT,
     6        18H DESIGNATION - (N)/
     7    1H0,52H D - NON-ARRAY VARIABLE HAS MULTIPLE VALUES ASSIGNED/
     8    1H0,43H E - ASSIGNED VALUES EXCEED ARRAY DIMENSION/
     9    1H0,17H F - SYNTAX ERROR//
     A    1H0,30H******************************,
     B        20H  INPUT DATA CARDS  ,
     C        30H******************************/)
      RETURN
      END
