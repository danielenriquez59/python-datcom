      SUBROUTINE XNAM22(IOP)
C
      COMMON /FLAPIN/ AM(137)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      INTEGER CONTAB
      LOGICAL EOF
C
      DIMENSION LENM(21),LDMM(21),CONTAB(71),LOCM(21)
      DIMENSION NLNAME(6)
C
      DATA NLNAME / 4HC   ,4HO   ,4HN   ,4HT   ,4HA   ,4HB   /
      DATA LENM / 5*5,4*4,7*2,5,2,2,3,4 /
      DATA LDMM / 21*1 /
      DATA LOCM / 117,118,119,122,123,120,121,124,125,126,127,
     1            128,129,130,131,132,133,134,135,136,137 /
      DATA CONTAB /   4HT   ,4HT   ,4HY   ,4HP   ,4HE   ,4HC   ,
     1  4HF   ,4HI   ,4HT   ,4HC   ,4HC   ,4HF   ,4HO   ,4HT   ,
     2  4HC   ,4HC   ,4HF   ,4HI   ,4HT   ,4HT   ,4HC   ,4HF   ,
     3  4HO   ,4HT   ,4HT   ,4HB   ,4HI   ,4HT   ,4HC   ,4HB   ,
     4  4HO   ,4HT   ,4HC   ,4HB   ,4HI   ,4HT   ,4HT   ,4HB   ,
     5  4HO   ,4HT   ,4HT   ,4HB   ,4H1   ,4HB   ,4H2   ,4HB   ,
     6  4H3   ,4HB   ,4H4   ,4HD   ,4H1   ,4HD   ,4H2   ,4HD   ,
     7  4H3   ,4HG   ,4HC   ,4HM   ,4HA   ,4HX   ,4HK   ,4HS   ,
     8  4HR   ,4HL   ,4HB   ,4HG   ,4HR   ,4HD   ,4HE   ,4HL   ,
     9  4HR   /
C
C**   IF IOP EQUAL ZERO READ NAMELIST CONTAB
C**   IF IOP EQUAL ONE WRITE NAMELIST CONTAB
C
      IF(IOP .EQ. 0)
     1 CALL NAMER(KAND,9,NLNAME,6,CONTAB,71,LENM,21,LDMM,AM,137,
     2            LOCM,EOF)
      IF(IOP .EQ. 1)
     1 CALL NAMEW(KAND,6,NLNAME,6,CONTAB,71,LENM,21,LDMM,AM,137,LOCM)
C
      RETURN
      END
