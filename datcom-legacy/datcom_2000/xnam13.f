      SUBROUTINE XNAM13(IOP)
C
      COMMON /POWER/  AE(50)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      COMMON /OPTION/ SR,CBARR,ROUGFC,BLREF
C
      LOGICAL EOF
C
      DIMENSION LENE(21),LDME(21),LARWB(89),LOCE(21)
      DIMENSION NLNAME(5)
C
      DATA NLNAME / 4HL   ,4HA   ,4HR   ,4HW   ,4HB   /
      DATA LENE / 2,4,6,6,2,6,6,1,4,6,5,2,2,3,3,6,6,3,5,6,5 /
      DATA LDME / 13*1,-1,1,1,-1,4*1 /
      DATA LOCE / 30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,
     1            45,46,47,48,49,50 /
      DATA LARWB  /   4HZ   ,4HB   ,4HS   ,4HR   ,4HE   ,4HF   ,
     1  4HD   ,4HE   ,4HL   ,4HT   ,4HE   ,4HP   ,4HS   ,4HF   ,
     2  4HR   ,4HO   ,4HN   ,4HT   ,4HA   ,4HR   ,4HR   ,4H3   ,
     3  4HL   ,4HE   ,4HO   ,4HB   ,4HD   ,4HE   ,4HL   ,4HT   ,
     4  4HA   ,4HL   ,4HL   ,4HS   ,4HW   ,4HE   ,4HT   ,4HP   ,
     5  4HE   ,4HR   ,4HB   ,4HA   ,4HS   ,4HS   ,4HB   ,4HA   ,
     6  4HS   ,4HE   ,4HH   ,4HB   ,4HB   ,4HB   ,4HB   ,4HL   ,
     7  4HF   ,4HX   ,4HC   ,4HG   ,4HT   ,4HH   ,4HE   ,4HT   ,
     8  4HA   ,4HD   ,4HR   ,4HO   ,4HU   ,4HN   ,4HD   ,4HN   ,
     9  4HS   ,4HB   ,4HS   ,4HS   ,4HB   ,4HS   ,4HL   ,4HB   ,
     A  4HX   ,4HC   ,4HE   ,4HN   ,4HS   ,4HB   ,4HX   ,4HC   ,
     B  4HE   ,4HN   ,4HW   /
C
C**  IF IOP EQUAL ZERO READ NAMELIST LARWB
C**  IF IOP EQUAL ONE WRITE NAMELIST LARWB
C
      IF(IOP .EQ. 0)
     1  CALL NAMER(KAND,9,NLNAME,5,LARWB,89,LENE,21,LDME,AE,50,LOCE,EOF)
      IF(IOP .EQ. 1)
     1  CALL NAMEW(KAND,6,NLNAME,5,LARWB,89,LENE,21,LDME,AE,50,LOCE)
C
      SR=AE(31)
      CBARR=AE(37)
      RETURN
      END
