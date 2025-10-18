      SUBROUTINE XNAM17(IOP)
C
      COMMON /FLAPIN/ AJ(137)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      INTEGER ASYFLP
      LOGICAL EOF
C
      DIMENSION LENJ(14),LDMJ(14),ASYFLP(78),LOCJ(14)
      DIMENSION NLNAME(6)
C
      DATA NLNAME / 4HA   ,4HS   ,4HY   ,4HF   ,4HL   ,4HP   /
      DATA LENJ / 4*6,4,4,5,6*6,5 /
      DATA LDMJ / 6*9,8*1 /
      DATA LOCJ / 19,29,1,39,49,60,18,59,16,12,13,14,15,11 /
      DATA ASYFLP /   4HD   ,4HE   ,4HL   ,4HT   ,4HA   ,4HL   ,
     1  4HD   ,4HE   ,4HL   ,4HT   ,4HA   ,4HR   ,4HD   ,4HE   ,
     2  4HL   ,4HT   ,4HA   ,4HD   ,4HD   ,4HE   ,4HL   ,4HT   ,
     3  4HA   ,4HS   ,4HX   ,4HS   ,4HO   ,4HC   ,4HH   ,4HS   ,
     4  4HO   ,4HC   ,4HS   ,4HT   ,4HY   ,4HP   ,4HE   ,4HX   ,
     5  4HS   ,4HP   ,4HR   ,4HM   ,4HE   ,4HN   ,4HD   ,4HE   ,
     6  4HL   ,4HT   ,4HA   ,4HC   ,4HH   ,4HR   ,4HD   ,4HF   ,
     7  4HI   ,4HC   ,4HH   ,4HR   ,4HD   ,4HF   ,4HO   ,4HS   ,
     8  4HP   ,4HA   ,4HN   ,4HF   ,4HI   ,4HS   ,4HP   ,4HA   ,
     9  4HN   ,4HF   ,4HO   ,4HP   ,4HH   ,4HE   ,4HT   ,4HE   /
C
C**   IF IOP EQUAL ZERO READ NAMELIST ASYFLP
C**   IF IOP EQUAL ONE WRITE NAMELIST ASYFLP
C
      IF(IOP .EQ. 0)
     1 CALL NAMER(KAND,9,NLNAME,6,ASYFLP,78,LENJ,14,LDMJ,AJ,137,
     2            LOCJ,EOF)
      IF(IOP .EQ. 1)
     1 CALL NAMEW(KAND,6,NLNAME,6,ASYFLP,78,LENJ,14,LDMJ,AJ,137,LOCJ)
C
      RETURN
      END
