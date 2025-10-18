      SUBROUTINE XNAM11(IOP)
C
      COMMON /POWER/  AC(50)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      INTEGER PROPWR
      LOGICAL EOF
C
      DIMENSION LENC(14),LDMC(14),PROPWR(78),LOCC(14)
      DIMENSION NLNAME(6)
C
      DATA NLNAME / 4HP   ,4HR   ,4HO   ,4HP   ,4HW   ,4HR   /
      DATA LENC / 12*6,4,2 /
      DATA LDMC / 12*1,-1,1 /
      DATA LOCC / 1,2,3,4,5,6,7,8,9,10,11,12,29,28 /
      DATA PROPWR /   4HA   ,4HI   ,4HE   ,4HT   ,4HL   ,4HP   ,
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
C**   IF IOP EQUAL ZERO READ NAMELISR PROPWR
C**   IF IOP EQUAL ONE WRITE NAMELIST PROPWR
C
      IF(IOP .EQ. 0)
     1 CALL NAMER(KAND,9,NLNAME,6,PROPWR,78,LENC,14,LDMC,AC,50,LOCC,EOF)
      IF(IOP .EQ. 1)
     1 CALL NAMEW(KAND,6,NLNAME,6,PROPWR,78,LENC,14,LDMC,AC,50,LOCC)
C
      RETURN
      END
