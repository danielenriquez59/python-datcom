      SUBROUTINE XNAM12(IOP)
C
      COMMON /POWER/  AD(50)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      LOGICAL EOF
C
      DIMENSION LEND(15),LDMD(15),JETPWR(89),LOCD(15)
      DIMENSION NLNAME(6)
C
C***  NAMELIST JETPWR
C
      DATA NLNAME / 4HJ   ,4HE   ,4HT   ,4HP   ,4HW   ,4HR   /
      DATA LEND / 14*6,5 /
      DATA LDMD / 15*1 /
      DATA LOCD / 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27 /
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
C**   IF IOP EQUAL ZERO READ NAMELIST JETPWR
C**   IF IOP EQUAL ONE WRITE NAMELIST JETPWR
C
      IF(IOP .EQ. 0)
     1  CALL NAMER(KAND,9,NLNAME,6,JETPWR,89,LEND,15,LDMD,AD,50,
     2             LOCD,EOF)
      IF(IOP .EQ. 1)
     1  CALL NAMEW(KAND,6,NLNAME,6,JETPWR,89,LEND,15,LDMD,AD,50,LOCD)
C
      RETURN
      END
